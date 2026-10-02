test_that("the core is linked and reports its compatibility target", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  target <- ancombc2_compat_target()
  expect_type(target, "character")
  expect_match(target, "ANCOMBC 2\\.15\\.2", fixed = FALSE)
  # The exact commit is the claim; a target string without it is not a target.
  expect_match(target, "dc4febdf59badb3a8dfe0c767ef2186323c2199a", fixed = TRUE)
})

# A small, deterministic two-group problem. Fixed rather than random: a test that
# changes its own input cannot fail reproducibly.
make_counts <- function(n_taxa = 40, n_samp = 20) {
  set.seed(42)
  grp <- rep(0:1, length.out = n_samp)
  m <- matrix(0, n_taxa, n_samp,
              dimnames = list(paste0("t", seq_len(n_taxa)),
                              paste0("s", seq_len(n_samp))))
  for (i in seq_len(n_taxa)) {
    base <- 20 + 60 * ((i * 37) %% 11) / 11
    signal <- if (i %% 5 == 0) 1.2 else 0
    mu <- (base + signal * grp)^2
    v <- ((i * 13 + seq_len(n_samp) * 7) %% 97) / 97
    m[i, ] <- pmax(0, round(sqrt(mu) * v * 4))
  }
  m
}

make_meta <- function(n_samp = 20) {
  data.frame(
    grp = factor(rep(c("a", "b"), length.out = n_samp)),
    x1 = seq_len(n_samp) / n_samp,
    row.names = paste0("s", seq_len(n_samp))
  )
}

test_that("a clean two-group problem runs end to end", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  res <- ancombc2(counts, meta, ~ grp + x1, group = "grp")
  expect_s3_class(res, "ancombc2_result")
  expect_true(is.data.frame(res$res))
  # `grp` has two levels, so `model.matrix` drops the first and names the single
  # dummy after the *remaining* level.
  expect_identical(res$fix_eff, c("(Intercept)", "grpb", "x1"))
  expect_identical(nrow(res$res), length(res$taxa))
  expect_identical(colnames(res$res)[1:7],
                   c("taxon", "lfc_(Intercept)", "lfc_grpb", "lfc_x1",
                     "se_(Intercept)", "se_grpb", "se_x1"))
  expect_true(all(c("p_grpb", "q_grpb", "diff_grpb") %in% colnames(res$res)))
  expect_true(is.numeric(unlist(res$delta_em)))
  expect_length(res$delta_em, 3)
})

test_that("the group is required exactly where the reference requires it", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  expect_error(ancombc2(counts, meta, ~ grp + x1, global = TRUE),
               "Group variable is required")
  expect_error(ancombc2(counts, meta, ~ grp + x1, pairwise = TRUE),
               "Group variable is required")
  expect_error(ancombc2(counts, meta, ~ grp + x1, struc_zero = TRUE),
               "Group variable is required")
  # ... and not otherwise.
  expect_silent(ancombc2(counts, meta, ~ grp + x1))
})

test_that("unsupported arguments are named, not dropped", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  expect_warning(ancombc2(counts, meta, ~ grp + x1, rand_formula = ~ 1),
                 "rand_formula")
  expect_warning(ancombc2(counts, meta, ~ grp + x1, dunnet = TRUE), "dunnet")
  expect_warning(ancombc2(counts, meta, ~ grp + x1, trend = TRUE), "trend")
})

test_that("a single-level group is a typed error, and two levels deactivate", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  meta$grp <- factor(rep("only", ncol(counts)))
  expect_error(ancombc2(counts, meta, ~ grp, group = "grp"),
               ">= 2 categories", fixed = TRUE)

  # Two levels with `global = TRUE`: the reference *deactivates* the test with a
  # warning rather than erroring, and the result must not carry one.
  meta2 <- make_meta(ncol(counts))
  res <- suppressWarnings(ancombc2(counts, meta2, ~ grp + x1, group = "grp",
                                    global = TRUE))
  expect_null(res$global)
  expect_true(any(grepl("< 3 categories", res$warnings)))
})

test_that("a size-one group is rejected", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  lev <- as.character(meta$grp)
  lev[1] <- "c"                       # a third level with exactly one sample
  meta$grp <- factor(lev)
  expect_error(ancombc2(counts, meta, ~ grp + x1, group = "grp"),
               ">= 2")
})

test_that("the p-value adjustment method reaches the core", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  col <- function(m, field) {
    r <- ancombc2(counts, meta, ~ grp + x1, group = "grp", p_adj_method = m)
    r$res[[field]]
  }
  # `none` is the least adjusted, so its q *is* its p.
  expect_equal(col("none", "q_grpb"), col("none", "p_grpb"))
  # Every other method is at least as conservative, so its q dominates `none`'s
  # pointwise. Comparing whole columns rather than one element, because a single
  # p-value of 1 stays 1 under every adjustment.
  for (m in c("holm", "BH", "bonferroni")) {
    q <- col(m, "q_grpb")
    expect_true(all(q >= col("none", "q_grpb") - 1e-12),
                info = paste("q under", m, "should dominate q under none"))
  }
  # And the methods are not all the same, so the argument reaches the core.
  expect_false(isTRUE(all.equal(col("holm", "q_grpb"), col("BH", "q_grpb"))))
  expect_error(ancombc2(counts, meta, ~ grp + x1, p_adj_method = "nonesuch"))
})

test_that("the two compatibility modes are both reachable", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  a <- ancombc2(counts, meta, ~ grp + x1, group = "grp", compat = "ancombc2-2.15")
  b <- ancombc2(counts, meta, ~ grp + x1, group = "grp", compat = "strict")
  # They differ wherever a taxon is unobserved, which the 0.1 substitution turns
  # into a variance contribution -- and with zeros in the table that is nearly
  # everywhere. `StrictSpec` is never larger, because it drops terms rather than
  # adding a constant.
  expect_false(isTRUE(all.equal(a$res$se_grpb, b$res$se_grpb)))
  expect_true(all(b$res$se_grpb <= a$res$se_grpb + 1e-12))
})

test_that("the sensitivity analysis adds its columns", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  res <- ancombc2(counts, meta, ~ grp + x1, group = "grp",
                  pseudo_sens = TRUE, conservative = TRUE)
  expect_true(all(c("passed_ss_grpb", "diff_robust_grpb") %in% colnames(res$res)))
  expect_length(res$sensitivity$pseudo, 4)
  # The reference's conservative grid.
  expect_equal(res$sensitivity$pseudo, list(0, 0.1, 0.5, 1))

  res2 <- ancombc2(counts, meta, ~ grp + x1, group = "grp",
                   pseudo_sens = TRUE, conservative = FALSE)
  expect_length(res2$sensitivity$pseudo, 50)
})

test_that("an all-zero taxon is filtered rather than fitted", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  # `rbind` rather than `counts["allzero", ] <- 0`, which is a subscript error on
  # a matrix that has no such row.
  counts <- rbind(counts, allzero = rep(0, ncol(counts)))
  meta <- make_meta(ncol(counts))
  res <- ancombc2(counts, meta, ~ grp + x1, group = "grp")
  # Prevalence, not the structural-zero screen, is what removes it: a taxon with
  # no counts has prevalence 0, which is below `prv_cut`, so it never reaches
  # either table. The two-taxon-set behaviour needs a taxon that *is* observed
  # somewhere and zero everywhere else, which is the structural-zero case.
  expect_false("allzero" %in% res$taxa)
  expect_false("allzero" %in% res$taxa_bias)
})

test_that("a run is reproducible and thread-count independent", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  a <- ancombc2(counts, meta, ~ grp + x1, group = "grp")
  b <- ancombc2(counts, meta, ~ grp + x1, group = "grp")
  expect_identical(a$res, b$res)
})

test_that("permuting taxa permutes the results", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  p <- c(3L, 1L, 2L, seq_len(nrow(counts))[-c(1L, 2L, 3L)])
  a <- ancombc2(counts, meta, ~ grp + x1, group = "grp")
  b <- ancombc2(counts[p, , drop = FALSE], meta, ~ grp + x1, group = "grp")
  # Run b's row i is run a's row p[i], so a's *permuted* column is b's.
  expect_equal(a$res$p_grpb[p], b$res$p_grpb)
  expect_identical(a$taxa[p], b$taxa)
})

test_that("an over-parameterised design is rejected rather than fitted", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts(n_samp = 4)
  meta <- make_meta(4)
  # Four samples and four columns leaves zero residual degrees of freedom, which
  # is R's `lm_smoke$df.residual == 0` guard. A *smaller* table is caught earlier,
  # by the group-size check, because each group would hold a single sample.
  # `x2` must not be collinear with `x1`, or the *rank* check fires first and the
  # test would pass for the wrong reason.
  meta$x2 <- c(4, 1, 3, 2) / 4
  expect_error(ancombc2(counts, meta, ~ grp + x1 + x2, group = "grp"),
               "degrees of freedom|residual")
  counts2 <- make_counts(n_samp = 2)
  meta2 <- make_meta(2)
  expect_error(ancombc2(counts2, meta2, ~ grp + x1, group = "grp"),
               ">= 2")
})

test_that("the print and summary methods work", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  res <- ancombc2(counts, meta, ~ grp + x1, group = "grp")
  out <- capture.output(print(res))
  expect_true(any(grepl("ancombc2-rs result", out)))
  s <- summary(res)
  expect_s3_class(s, "summary.ancombc2_result")
  expect_equal(dim(s$beta), c(length(res$taxa), length(res$fix_eff)))
  expect_true(any(grepl("significant", capture.output(print(s)))))
})

test_that("a malformed request is an R error, not a crash", {
  skip_if_not(ancombc2_available(), "the Rust core is not linked")
  counts <- make_counts()
  meta <- make_meta(ncol(counts))
  expect_error(ancombc2(counts, as.list(meta), ~ grp + x1), "data.frame")
  expect_error(ancombc2(counts, meta), "required")
  expect_error(ancombc2(counts, meta, ~ grp + nosuchcolumn), "not found")
  # A missing covariate leaves an incomplete design row, which the core reports as
  # an unidentifiable covariate rather than silently dropping the sample.
  meta_na <- meta
  meta_na$x1[3] <- NA
  expect_error(ancombc2(counts, meta_na, ~ grp + x1, group = "grp"))
})
