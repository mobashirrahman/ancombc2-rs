#!/usr/bin/env Rscript
# Generate the structural-zero and edge-case matrix (PLAN.md section 5.5).
#
#   Rscript --vanilla scripts/generate_edge_cases.R
#
# Each case is a small fixture plus a golden captured from the pinned oracle:
# which taxa were retained, the structural-zero table, and the reported
# coefficients. `crates/ancombc2-core/tests/edge_cases.rs` then runs the core on
# the same fixture and asserts both the case's own predicate and general parity.
#
# Why a golden rather than a hand-written expectation
# ---------------------------------------------------
# The plan's table is written as "Expected", and for most rows the only
# authoritative answer is what the reference does. `prv_cut` is
# `prevalence >= prv_cut` in `ancombc_prep.R`, so a taxon exactly at the cut is
# *kept*; whether a group of size 1 is dropped or errors is likewise a property
# of the reference, not a preference. Capturing the oracle's answer and asserting
# against it is what makes these cases parity rather than opinion.
#
# Two rows name parameters this oracle does not have
# ---------------------------------------------------
# The plan lists `keep_zero` and `perc_thres`. ANCOMBC 2.15.2 at
# dc4febdf has neither; the nearest real arguments are `struc_zero` (which
# controls whether the flagged taxa are removed from the primary fit) and
# `prv_cut` (the prevalence filter). Those two rows are implemented as
# `struc_zero_off` and `prv_cut_boundary`, and the substitution is recorded in
# `docs/reference_behavior.md` rather than papered over.

suppressWarnings(suppressMessages({
  library(utils)
}))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  args[[i + 1L]]
}
root <- arg("--out", "validation/edge")
oracle_dir <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC")
for (f in c("oracle.R", "stubs.R", "harness.R")) {
  if (!file.exists(file.path("reference/R", f))) {
    stop("cannot find reference/R/", f, "; run from the repository root")
  }
  source(file.path("reference/R", f))
}
H <- load_harness(oracle_dir)

# ---- writers ----------------------------------------------------------------

write_table <- function(m, path, ids) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  # R's `write.table(row.names = TRUE)` layout: an unnamed first header field.
  writeLines(paste(colnames(m), collapse = "\t"), con)
  for (i in seq_len(nrow(m))) {
    writeLines(paste(c(ids[i], format(m[i, ], trim = TRUE, scientific = FALSE,
                                           digits = 15)), collapse = "\t"), con)
  }
}

write_analysis <- function(a, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("key\tvalue", con)
  for (nm in names(a)) {
    v <- a[[nm]]
    if (is.null(v) || length(v) == 0L) next
    if (is.logical(v)) v <- if (v) "true" else "false"
    writeLines(paste(nm, v, sep = "\t"), con)
  }
}

# JSON, by hand: `jsonlite` is not a dependency of this repository and the row
# schema here is small.
jnum <- function(x) {
  if (length(x) == 0L) return("null")
  x <- as.numeric(x)
  if (is.na(x[1]) || is.infinite(x[1])) return("null")
  sprintf("%.17g", x)
}
jvec <- function(x) {
  if (length(x) == 0L) return("[]")
  paste0("[", paste(vapply(x, jnum, ""), collapse = ","), "]")
}
jstrs <- function(x) {
  x <- as.character(x)
  if (length(x) == 0L) return("[]")
  one <- function(s) paste0("\"", gsub("\"", "\\\\\"", s), "\"")
  paste0("[", paste(vapply(x, one, ""), collapse = ","), "]")
}
# A logical vector as JSON booleans, element by element rather than through
# `ifelse`: `ifelse` on a factor returns a *factor*, whose `as.character` gives
# the level labels -- so a column of `NA` came out as the literal text `NA`,
# which is not JSON and made the whole golden unparseable.
jbool <- function(x) {
  if (length(x) == 0L) return("[]")
  one <- function(v) {
    v <- as.logical(v)
    if (length(v) == 0L || is.na(v)) "null" else if (v) "true" else "false"
  }
  paste0("[", paste(vapply(as.list(x), one, ""), collapse = ","), "]")
}
# The zero_ind table as an object keyed by column name.
#
# The `taxon` column holds names, not flags, so it is written separately: passing
# it through `jbool` turned every label into `null` and lost the only thing that
# identifies which row a flag belongs to.
jmat <- function(m) {
  cols <- colnames(m)
  flag_cols <- which(vapply(m, function(c) is.logical(c) || is.numeric(c), TRUE))
  body <- vapply(flag_cols, function(j) {
    paste0(jstr(cols[j]), ":", jbool(m[, j]))
  }, "")
  paste0("{", paste(body, collapse = ","), "}")
}
jstr <- function(s) paste0("\"", gsub("\"", "\\\\\"", as.character(s)[1]), "\"")

# ---- the oracle's answer ----------------------------------------------------

# Run the oracle and capture everything the Rust side must reproduce. A case the
# oracle *refuses* is a first-class result: `error` is recorded and the Rust
# side must fail the same way.
capture <- function(counts, meta, cfg) {
  t0 <- Sys.time()
  res <- try(suppressWarnings(suppressMessages(
    H$ref_run(
      counts, meta, as.character(cfg$formula),
      p_adj_method = if (is.null(cfg$p_adj_method)) "BH" else cfg$p_adj_method,
      pseudo = if (is.null(cfg$pseudo)) 0.5 else cfg$pseudo,
      pseudo_sens = identical(cfg$pseudo_sens, "true"),
      conservative = if (is.null(cfg$conservative)) TRUE else identical(cfg$conservative, "true"),
      prv_cut = if (is.null(cfg$prv_cut)) 0 else cfg$prv_cut,
      lib_cut = if (is.null(cfg$lib_cut)) 0 else cfg$lib_cut,
      s0_perc = if (is.null(cfg$s0_perc)) 0.05 else cfg$s0_perc,
      group = if (is.null(cfg$group)) NULL else cfg$group,
      struc_zero = identical(cfg$struc_zero, "true"),
      neg_lb = identical(cfg$neg_lb, "true"),
      alpha = if (is.null(cfg$alpha)) 0.05 else cfg$alpha,
      global = identical(cfg$global, "true"),
      pairwise = identical(cfg$pairwise, "true")
    )
  )), silent = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "try-error")) {
    return(list(
      error = gsub("[\"\r\n\t]", " ", trimws(as.character(res))),
      elapsed = elapsed
    ))
  }
  taxa <- as.character(res$taxa_retained)
  coef_names <- colnames(res$beta)
  # The group contrast, named. See scripts/realdata_r.R for why a positional
  # fallback is wrong on a design with a covariate.
  col <- NA_integer_
  if (!is.null(cfg$group) && !is.null(cfg$group_levels)) {
    lv <- trimws(strsplit(as.character(cfg$group_levels), ",", fixed = TRUE)[[1]])
    for (i in seq_along(lv)[-1]) {
      cand <- match(paste0(cfg$group, lv[[i]]), coef_names)
      if (!is.na(cand)) { col <- cand; break }
    }
  }
  if (is.na(col)) col <- which(coef_names != "(Intercept)")[[1]]
  zi <- NULL
  if (!is.null(res$zero_ind)) zi <- as.data.frame(res$zero_ind)
  list(
    error = NULL,
    elapsed = elapsed,
    taxa_retained = taxa,
    fix_eff = coef_names,
    coefficient = coef_names[col],
    zero_ind = zi,
    beta = as.numeric(res$beta[, col]),
    se = as.numeric(res$se[, col]),
    p = as.numeric(res$p[, col]),
    q = as.numeric(res$q[, col]),
    diff_abn = as.logical(res$diff_abn[, col])
  )
}

emit_case <- function(c) {
  dir <- file.path(root, c$name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write_table(c$counts, file.path(dir, "counts.tsv"), rownames(c$counts))
  write_analysis(c$analysis, file.path(dir, "analysis.tsv"))
  m <- c$meta
  # The metadata's row names are the sample names, and the label column comes
  # from the *original* frame's row names: subsetting a data frame to one row
  # replaces explicit row names with automatic ones, so `rownames(m[i, ])` is
  # empty and the line comes out as values only.
  con <- file(file.path(dir, "meta.tsv"), open = "wt")
  writeLines(paste(colnames(m), collapse = "\t"), con)
  labels <- rownames(m)
  for (i in seq_len(nrow(m))) {
    row <- m[i, , drop = FALSE]
    values <- vapply(row, as.character, "")
    writeLines(paste(c(labels[i], values[colnames(m)]), collapse = "\t"), con)
  }
  close(con)

  g <- capture(c$counts, m, c$analysis)
  out <- c(
    sprintf('"name":%s', jstr(c$name)),
    sprintf('"plan_row":%s', jstr(c$plan_row)),
    sprintf('"about":%s', jstr(c$about)),
    sprintf('"n_taxa":%d,"n_samp":%d', nrow(c$counts), ncol(c$counts)),
    sprintf('"groups":%s', jstrs(c$groups))
    # Deliberately no wall-clock. The CI job regenerates these and diffs them
    # against the committed copies, and a measured duration differs every run --
    # so recording it would make the drift check fail on a healthy oracle. How
    # long the reference takes is a benchmark question, and
    # `benchmarks/results/results.jsonl` is where that belongs.
  )
  if (!is.null(g$error)) {
    out <- c(out, sprintf('"error":%s', jstr(g$error)), '"taxa_retained":[],"beta":[],"se":[],"p":[],"q":[],"diff_abn":[],"zero_ind":null')
  } else {
    out <- c(
      out, '"error":null',
      sprintf('"taxa_retained":%s', jstrs(g$taxa_retained)),
      sprintf('"fix_eff":%s', jstrs(g$fix_eff)),
      sprintf('"coefficient":%s', jstr(g$coefficient)),
      sprintf('"zero_ind":%s', if (is.null(g$zero_ind)) "null" else jmat(g$zero_ind)),
      sprintf('"beta":%s', jvec(g$beta)),
      sprintf('"se":%s', jvec(g$se)),
      sprintf('"p":%s', jvec(g$p)),
      sprintf('"q":%s', jvec(g$q)),
      sprintf('"diff_abn":%s', jbool(g$diff_abn))
    )
  }
  writeLines(paste0("{", paste(out, collapse = ","), "}"),
             file.path(dir, "golden.json"))
  status <- if (is.null(g$error)) {
    sprintf("%d taxa, %d flagged", length(g$taxa_retained),
            if (is.null(g$zero_ind)) 0L else sum(as.matrix(g$zero_ind[, -1, drop = FALSE])))
  } else {
    "the oracle refuses it"
  }
  cat(sprintf("  %-24s %3d x %-4d %s\n", c$name, nrow(c$counts), ncol(c$counts), status))
}

# ---- the cases --------------------------------------------------------------

# A small deterministic count table. `set.seed` is fixed so the cases are
# reproducible, and the values are integers so the file is a count table and not
# a question about rounding.
base_counts <- function(n_t, n_s, seed) {
  set.seed(seed)
  m <- matrix(as.numeric(rpois(n_t * n_s, 40)), n_t, n_s)
  rownames(m) <- paste0("T", seq_len(n_t))
  colnames(m) <- paste0("S", seq_len(n_s))
  m
}

# `sizes` gives each group's sample count, so a case can ask for an unbalanced
# split. `rep(levels, length.out = n)` -- the obvious form -- *alternates* the
# levels, so a "group of size 1" case built that way silently produced a 4/3
# split and tested nothing. The size is stated, and the counts are verified.
two_group_meta <- function(sizes = c(5, 5), levels = c("g1", "g2")) {
  if (length(sizes) != length(levels)) stop("one size per level")
  grp <- rep(levels, times = sizes)
  m <- data.frame(group = factor(grp, levels = levels),
                  row.names = paste0("S", seq_along(grp)))
  stopifnot(nrow(m) == sum(sizes))
  m
}

cat("generating the structural-zero and edge-case matrix\n")
dir.create(root, recursive = TRUE, showWarnings = FALSE)

cases <- list()

# 1. Completely absent in one group.
{
  cnt <- base_counts(6, 10, 101)
  cnt[1, 1:5] <- 0
  cases[[length(cases) + 1L]] <- list(
    name = "absent_in_group_a",
    plan_row = "Completely absent in group A",
    about = paste("taxon T1 is zero in every sample of g1 and present in g2, so",
                  "its group-1 prevalence is 0 and it must be flagged for that",
                  "group only"),
    counts = cnt,
    meta = two_group_meta(c(5, 5)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0, lib_cut = 0, struc_zero = "true",
                    neg_lb = "true", alpha = 0.05)
  )
}

# 2. Rare in one group: a single observation.
{
  cnt <- base_counts(6, 10, 102)
  cnt[2, 2:5] <- 0
  cnt[2, 1] <- 12
  cases[[length(cases) + 1L]] <- list(
    name = "rare_in_group_a",
    plan_row = "Rare in group A (< 1 obs)",
    about = paste("taxon T2 has exactly one observation in g1, so its g1",
                  "prevalence is 1/5. The reference's `neg_lb` widens the net",
                  "beyond a strict zero, so it should still be flagged -- which is",
                  "the point of the case: a `< 1 obs` taxon and a",
                  "structurally-zero one are not the same test"),
    counts = cnt,
    meta = two_group_meta(c(5, 5)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0, lib_cut = 0, struc_zero = "true",
                    neg_lb = "true", alpha = 0.05)
  )
}

# 3. Present in every group.
{
  cnt <- base_counts(6, 10, 103)
  cases[[length(cases) + 1L]] <- list(
    name = "present_in_every_group",
    plan_row = "Present in every group",
    about = paste("every taxon is observed in every sample of every group, so no",
                  "taxon may be flagged. This is the negative control for the two",
                  "cases above: a rule that flags too eagerly passes those and",
                  "fails this"),
    counts = cnt,
    meta = two_group_meta(c(5, 5)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0, lib_cut = 0, struc_zero = "true",
                    neg_lb = "true", alpha = 0.05)
  )
}

# 4. NA in the counts, on both sides of the zero/present boundary.
{
  cnt <- base_counts(6, 10, 104)
  # NA is *missing*, not zero. The two must be told apart: the reference's
  # sandwich replaces an NA term with 0.1, which is a documented quirk, and it
  # treats a zero as a structural absence. A case that conflates them would pass
  # while testing nothing.
  cnt[3, 1:3] <- NA
  cnt[4, 1:3] <- 0
  cases[[length(cases) + 1L]] <- list(
    name = "na_counts",
    plan_row = "NA in counts",
    about = paste("T3 is NA across g1 and T4 is zero across g1, so the two",
                  "differ only in how they are missing. The reference's",
                  "sandwich replaces an NA term with 0.1 rather than skipping it",
                  "(docs/reference_behavior.md section 2), and this case pins",
                  "that behaviour as a contract rather than a caveat"),
    counts = cnt,
    meta = two_group_meta(c(5, 5)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0, lib_cut = 0, struc_zero = "true",
                    neg_lb = "true", alpha = 0.05)
  )
}

# 5. A group of size 1.
{
  cnt <- base_counts(5, 7, 105)
  cases[[length(cases) + 1L]] <- list(
    name = "group_of_size_one",
    plan_row = "Group of size 1",
    about = paste("g2 has exactly one sample, which is the plan's case. The",
                  "reference validates group sizes before anything else, so this",
                  "is expected to be a refusal rather than a result; the golden",
                  "records whichever it is, and the Rust side must match it",
                  "either way. Written as sizes c(6, 1) rather than by alternating",
                  "the levels, which would have produced 4/3 and tested nothing"),
    counts = cnt,
    meta = two_group_meta(c(6, 1)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0, lib_cut = 0, struc_zero = "true",
                    neg_lb = "true", alpha = 0.05)
  )
}

# 6. `keep_zero = TRUE`, which this oracle spells `struc_zero = FALSE`.
#
# The flagged taxa stay in the primary fit instead of being dropped, and
# `zero_ind` is still reported. Two different questions in one case, and the
# second is the one a "just skip the removal" implementation would get wrong.
{
  cnt <- base_counts(6, 10, 106)
  cnt[1, 1:5] <- 0
  cases[[length(cases) + 1L]] <- list(
    name = "struc_zero_off",
    plan_row = "keep_zero = TRUE (as struc_zero = FALSE)",
    about = paste("the same fixture as `absent_in_group_a` with the structural-",
                  "zero removal switched off, so T1 is retained *and* still",
                  "reported in the zero_ind table. ANCOMBC 2.15.2 has no",
                  "`keep_zero` argument; `struc_zero` is the nearest one and this",
                  "is the documented substitution"),
    counts = cnt,
    meta = two_group_meta(c(5, 5)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0, lib_cut = 0, struc_zero = "false",
                    neg_lb = "true", alpha = 0.05)
  )
}

# 7. Prevalence exactly at the filter boundary, which this oracle spells
#    `prv_cut`.
{
  cnt <- base_counts(5, 10, 107)
  # T1 is observed in exactly 3 of 10 samples. With prv_cut = 0.3 that is
  # `prevalence >= prv_cut` exactly, and `ancombc_prep.R` uses `>=`, so T1 is
  # *kept*. A `<` would drop it. One taxon at the boundary is the whole test.
  cnt[1, 1:7] <- 0
  cnt[2, 1:7] <- 0   # 3/10 observed: exactly at the cut
  cnt[3, 1:8] <- 0   # 2/10 observed: just below
  cases[[length(cases) + 1L]] <- list(
    name = "prv_cut_boundary",
    plan_row = "perc_thres exactly at boundary (as prv_cut)",
    about = paste("T2 is observed in exactly 3 of 10 samples and prv_cut is 0.3,",
                  "so its prevalence equals the cut. `ancombc_prep.R` filters on",
                  "`prevalence >= prv_cut`, so T2 is kept and T3 (2/10) is not.",
                  "ANCOMBC 2.15.2 has no `perc_thres`; `prv_cut` is the nearest",
                  "one and this is the documented substitution"),
    counts = cnt,
    meta = two_group_meta(c(5, 5)),
    groups = c("g1", "g2"),
    analysis = list(formula = "group", group = "group",
                    group_levels = "g1,g2", p_adj_method = "BH", pseudo = 0.5,
                    prv_cut = 0.3, lib_cut = 0, struc_zero = "false",
                    neg_lb = "false", alpha = 0.05)
  )
}

for (c in cases) emit_case(c)

# An index, so the Rust test and a reader can enumerate the cases without
# globbing, and so the plan's row names are recorded once.
idx <- vapply(cases, function(c) {
  sprintf('  {"name":%s,"plan_row":%s,"about":%s,"dir":%s}',
          jstr(c$name), jstr(c$plan_row), jstr(c$about), jstr(c$name))
}, "")
writeLines(c("{", " \"matrix\": \"PLAN.md section 5.5, structural-zero and edge cases\",",
             " \"note\": \"Two of the plan's rows name parameters ANCOMBC 2.15.2 does",
             " not have (keep_zero, perc_thres); the substitutions used are recorded",
             " per case and in docs/reference_behavior.md.\",",
             " \"cases\": [", paste(idx, collapse = ",\n"), "  ]", "}"),
           file.path(root, "index.json"))
cat(sprintf("wrote %d case(s) to %s\n", length(cases), root))
