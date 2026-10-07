# S07: the output transport's acceptance check.
#
# The claim under test is narrow and stated exactly: *the core's payloads survive
# a round trip out of R and back into R bit for bit, and the original's result can
# be rebuilt from them to the byte.* Nothing here claims the bridge computed
# anything -- it cannot yet -- so the test does not pretend otherwise.
#
# Run against the installed replacement library, so what is exercised is the
# compiled `init.c` and the compiled Rust bridge. `make r-output-selftest` does that
# and also passes the fixture path; under a bare `R CMD check` the fixture is not
# found and this exits 1 naming the path, rather than passing with nothing tested.

suppressWarnings(suppressMessages({
  library(ANCOMBC)
}))

fails <- 0L
checks <- 0L
ok <- function(label, cond) {
  checks <<- checks + 1L
  if (isTRUE(cond)) {
    cat(sprintf("ok    %s\n", label))
  } else {
    fails <<- fails + 1L
    cat(sprintf("FAIL  %s\n", label))
  }
}
bad <- function(label, cond, why) {
  checks <<- checks + 1L
  if (isTRUE(cond)) {
    cat(sprintf("ok    %s\n", label))
  } else {
    fails <<- fails + 1L
    cat(sprintf("FAIL  %s -- %s\n", label, why))
  }
}

# The exact comparison the whole project is judged by. Not `all.equal`, not a
# tolerance: `serialize()` of each, byte for byte.
same_bytes <- function(x, y) {
  identical(serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3),
            serialize(y, NULL, ascii = FALSE, xdr = TRUE, version = 3))
}

payload_path <- Sys.getenv("ANCOMBC_PAYLOAD_RDS", "validation/exact/payloads/tiny-defaults.rds")
if (!file.exists(payload_path)) {
  cat(sprintf("FAIL  the payload fixture %s does not exist\n", payload_path))
  cat("       These payloads are captured from the pinned original and are not\n",
      "       reproducible without it. Regenerate with:\n",
      "         Rscript scripts/capture_core_payloads.R\n")
  quit(status = 1L)
}
ok("the payload fixture exists", TRUE)
fixture <- readRDS(payload_path)
int <- fixture$internals
want <- fixture$expected

ok("the bridge is present", ANCOMBC:::rb_available())
bad("the fixture names the pinned original it came from",
    grepl("ANCOMBC 2.15.2 @ dc4febdf", fixture$provenance$oracle, fixed = TRUE),
    fixture$provenance$oracle)
ok("the fixture carries the warning that it is test-only",
   grepl("never be read by a candidate run", fixture$provenance$warning))

# ---------------------------------------------------------------------------
# 1. The plan
# ---------------------------------------------------------------------------
pay <- ANCOMBC:::rb_emit_payload(int)

bad("the plan has the eleven documented entries",
    length(pay$plan_names) == 11L,
    paste("got", length(pay$plan_names)))
bad("the plan is in the documented order",
    identical(pay$plan_names,
              c("theta_hat", "beta_hat", "var_hat", "dof", "vcov_hat", "delta_em",
                "delta_wls", "var_delta", "s02", "y_bias_crt", "O2")),
    paste(paste(pay$plan_names, collapse = ",")))
bad("`dof` comes back integer because the original's was",
    identical(pay$dof_is_int, 1L) && is.integer(pay$payload$dof),
    paste("dof_is_int =", pay$dof_is_int, "typeof =", typeof(pay$payload$dof)))
bad("`vcov_hat` comes back as a list, not a matrix",
    is.list(pay$payload$vcov_hat) && length(pay$payload$vcov_hat) == nrow(int$O2),
    paste("is.list =", is.list(pay$payload$vcov_hat)))

# ---------------------------------------------------------------------------
# 2. The payloads themselves
# ---------------------------------------------------------------------------
pl <- pay$payload
bad("`beta_hat` is byte-identical", same_bytes(pl$beta_hat, int$beta_hat), "differs")
bad("`var_hat` is byte-identical", same_bytes(pl$var_hat, int$var_hat), "differs")
bad("`dof` is byte-identical, integer included",
    same_bytes(pl$dof, int$dof), "differs")
bad("`theta_hat` is byte-identical", same_bytes(pl$theta_hat, int$theta_hat), "differs")
bad("`delta_em` is byte-identical", same_bytes(pl$delta_em, int$delta_em), "differs")
bad("`delta_wls` is byte-identical", same_bytes(pl$delta_wls, int$delta_wls), "differs")
bad("`var_delta` is byte-identical", same_bytes(pl$var_delta, int$var_delta), "differs")
bad("`s02` is byte-identical", same_bytes(pl$s02, int$s02), "differs")
bad("`O2` is byte-identical", same_bytes(pl$O2, int$O2), "differs")
bad("`y_bias_crt` is byte-identical as a matrix",
    same_bytes(as.matrix(pl$y_bias_crt), as.matrix(int$y_bias_crt)), "differs")
bad("every covariance matrix is byte-identical",
    identical(lapply(pl$vcov_hat, unname), lapply(int$vcov_hat, unname)),
    "differs")
bad("every covariance matrix keeps its dim", all(vapply(pl$vcov_hat, is.matrix, logical(1))),
    "a covariance came back without a dim")

# A transposed covariance is still a symmetric covariance, so symmetry cannot be
# the check. `p == 3` makes the transpose visible: compare the raw wire block
# against an explicitly row-major rebuild.
vc_wire <- ANCOMBC:::rb_flatten_vcov(int$vcov_hat)
bad("the wire block is row-major within each taxon",
    identical(vc_wire[1:9], as.vector(t(int$vcov_hat[[1L]]))),
    "the first taxon block is not the row-major rebuild")
bad("row-major in, column-major out, for a non-symmetric block",
    identical(pl$vcov_hat[[1L]][1, 2], int$vcov_hat[[1L]][1, 2]) &&
      identical(pl$vcov_hat[[1L]][2, 1], int$vcov_hat[[1L]][2, 1]) &&
      identical(pl$vcov_hat[[1L]][1, 3], int$vcov_hat[[1L]][1, 3]),
    "the (1,2)/(2,1)/(1,3) entries moved")

# The awkward cases: NA, NaN, Inf, -0.0. These are the bits a lossy transport drops
# first, and they are the bits `serialize()` compares.
special <- list(
  beta_hat = matrix(c(NA_real_, NaN, Inf, -Inf, -0.0, 0.0), nrow = 2),
  var_hat = matrix(c(NaN, NA_real_, Inf, -0.0, 0.0, -Inf), nrow = 2),
  dof = matrix(c(1L, 2L, 3L, 4L, 5L, 6L), nrow = 2),
  # 3x3, because `x` has three columns and every covariance is p x p.
  vcov_hat = list(matrix(c(NA_real_, 1, 2, 3, 4, NaN, 6, 7, 8), 3),
                  matrix(c(-0.0, 0, Inf, 12, -Inf, 14, 15, 16, 17), 3)),
  y_bias_crt = matrix(c(NA_real_, -0.0, 0, Inf, -Inf, NaN), nrow = 2),
  theta_hat = c(NA_real_, -0.0, NaN),
  # Length p, which `x` makes three.
  delta_em = c(NA_real_, -0.0, 0.0),
  delta_wls = c(NaN, 0.0, Inf),
  var_delta = c(Inf, -Inf, -0.0),
  s02 = c(-0.0, 0.0, NaN),
  O2 = matrix(c(0, -0.0, NA_real_, NaN, Inf, -Inf), nrow = 2),
  x = matrix(c(1, 0, 0, 1, 1, 1), nrow = 2),
  fix_eff = c("(Intercept)", "x1", "x2")
)
special$tax_name <- c("T1", "T2")
special$y_bias_crt <- data.frame(special$y_bias_crt)
sp <- ANCOMBC:::rb_emit_payload(special)
bad("NA, NaN, Inf, -Inf and -0.0 all survive `beta_hat`",
    identical(sp$payload$beta_hat[1:6], special$beta_hat[1:6]),
    "differs")
# `-0.0 == 0.0`, so equality cannot see the sign. Division can: `1/-0.0` is `-Inf`
# and `1/0.0` is `Inf`. That is the observable difference, and it is what
# `serialize()` compares.
bad("-0.0 keeps its sign bit through `beta_hat`",
    identical(1 / sp$payload$beta_hat[5], -Inf) && identical(1 / sp$payload$beta_hat[6], Inf),
    "the sign of zero was lost or flipped")
# `is.na()` is TRUE for NaN, so it cannot tell the two apart on its own.
# `is.nan()` can, and that is the only thing that distinguishes them.
bad("NA and NaN stay distinguishable in `var_hat`",
    is.nan(sp$payload$var_hat[1]) && is.na(sp$payload$var_hat[2]) &&
      !is.nan(sp$payload$var_hat[2]),
    "NA and NaN were conflated")
# The fixture's first covariance is `matrix(c(NA,1,2,3,4,NaN,6,7,8), 3)`, which in
# column-major order puts NA at (1,1) and NaN at (3,2). Naming those two positions
# explicitly matters: the transport is where a row-major/column-major mix-up would
# move them, and a symmetric covariance would not notice.
bad("NA survives in `vcov_hat` as NA, not NaN",
    is.na(sp$payload$vcov_hat[[1]][1, 1]) && !is.nan(sp$payload$vcov_hat[[1]][1, 1]) &&
      is.nan(sp$payload$vcov_hat[[1]][3, 2]),
    "the NA/NaN distinction was lost in a covariance")
bad("the two covariance matrices stay separate",
    is.na(sp$payload$vcov_hat[[1]][1, 1]) && !is.na(sp$payload$vcov_hat[[2]][1, 1]),
    "the taxa blocks were merged")
# `unname` and a bare dim: the transport *adds* dimnames the fixture did not have,
# and the test is about the type and the values, not about the absence of names.
bare <- function(x) {
  dimnames(x) <- NULL
  attr(x, "dimnames") <- NULL
  names(x) <- NULL
  x
}
bad("integer `dof` comes back integer",
    is.integer(sp$payload$dof) && identical(bare(sp$payload$dof), bare(special$dof)),
    paste("typeof =", typeof(sp$payload$dof)))

# The reverse asymmetry: `n_tax` is `nrow(O2)`, so `taxa_bias` may be longer.
asym <- special
asym$taxa <- c(1L, 3L)
asym$taxa_bias <- c(1L, 2L, 3L)
a <- ANCOMBC:::rb_emit_payload(asym)
ok("a reported set shorter than the bias set is carried, not rejected",
   identical(a$taxa, c(1L, 3L)) && identical(a$taxa_bias, c(1L, 2L, 3L)))
bad("a bias set *shorter* than the reported set is refused",
    inherits(try(ANCOMBC:::rb_emit_payload(
      modifyList(asym, list(taxa = c(1L, 2L, 3L, 4L), taxa_bias = c(1L, 2L)))),
      silent = TRUE), "try-error"),
    "a taxa_bias shorter than taxa was accepted: O1 is a superset of O2")

# ---------------------------------------------------------------------------
# 3. The assembly, against the original's own result
# ---------------------------------------------------------------------------
got <- ANCOMBC:::ancombc2_assemble_core(
  pay,
  p_adj_method = int$p_adj_method, alpha = int$alpha, tax_name = int$tax_name,
  global = int$global, pairwise = int$pairwise, dunnet = int$dunnet,
  trend = int$trend
)

bad("the assembled result has the original's ten names",
    identical(names(got), names(want)),
    paste(paste(names(got), collapse = ",")))
bad("`feature_table` is byte-identical to the original's",
    same_bytes(got$feature_table, want$feature_table), "differs")
bad("`bias_correct_log_table` is byte-identical to the original's",
    same_bytes(got$bias_correct_log_table, want$bias_correct_log_table), "differs")
bad("`samp_frac` is byte-identical", same_bytes(got$samp_frac, want$samp_frac), "differs")
bad("`delta_em` is byte-identical", same_bytes(got$delta_em, want$delta_em), "differs")
bad("`delta_wls` is byte-identical", same_bytes(got$delta_wls, want$delta_wls), "differs")
bad("`res` is byte-identical", same_bytes(got$res, want$res), "differs")
bad("`res_global` is the original's NULL", identical(got$res_global, want$res_global), "differs")
bad("`res_pair` is the original's NULL", identical(got$res_pair, want$res_pair), "differs")
bad("`res_dunn` is the original's NULL", identical(got$res_dunn, want$res_dunn), "differs")
bad("`res_trend` is the original's NULL", identical(got$res_trend, want$res_trend), "differs")
bad("the whole core result is byte-identical", same_bytes(got, want), "differs")

# `res` is where a rounding or naming slip shows up, so it gets its own report of
# where it first went wrong rather than a bare "differs".
if (!same_bytes(got$res, want$res)) {
  for (nm in names(want$res)) {
    if (!same_bytes(got$res[[nm]], want$res[[nm]])) {
      cat(sprintf("      res$%s differs: got %s want %s\n", nm,
                  paste(utils::head(format(got$res[[nm]]), 3), collapse = " "),
                  paste(utils::head(format(want$res[[nm]]), 3), collapse = " ")))
    }
  }
}

# ---------------------------------------------------------------------------
# 4. The claim this file does not make
# ---------------------------------------------------------------------------
bad("an untransported test is a named error, not a silent NULL",
    inherits(try(ANCOMBC:::ancombc2_assemble_core(pay, global = TRUE), silent = TRUE),
             "try-error"),
    "global = TRUE was accepted silently")

cat(sprintf("\n%d/%d checks passed\n", checks - fails, checks))
if (fails > 0L) quit(status = 1L)
