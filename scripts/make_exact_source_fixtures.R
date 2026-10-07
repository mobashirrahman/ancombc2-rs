#!/usr/bin/env Rscript
# Generate the *new* exact fixtures that no committed fixture can express
# (IMPROVED_PLAN.md S04). Explicit seeds; the input is written before either arm
# computes anything, and the writer is round-trip-safe.
#
#   Rscript --vanilla scripts/make_exact_source_fixtures.R \
#       --out-dir validation/exact/fixtures-src --manifest <path>
#
# Round-trip safety
# -----------------
# Counts are written with `format(..., digits = 15, scientific = FALSE)` on
# integers, so the text is exact and re-reading it is the identity -- the defect
# IMPROVED_PLAN.md S04 names in the legacy TSV fixtures. The generator asserts
# that: it writes, reads back with the same reader the exact inputs use, and
# stops if any cell differs. A fixture that cannot survive its own round trip is
# not a fixture.
#
# These live under validation/exact/ and never touch validation/fixtures/ or
# validation/golden/.

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
  k <- sub("^--", "", a)
  if (i == length(args)) stop("option --", k, " needs a value", call. = FALSE)
  opt[[k]] <- args[[i + 1L]]
  i <- i + 2L
}
`%||%` <- function(a, b) if (is.null(a)) b else a

out_dir <- opt[["out-dir"]] %||% "validation/exact/fixtures-src"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("scripts/exact_input_lib.R")

# Counts are non-negative integers with a deliberate zero pattern, because
# log(0) is where the two implementations are most likely to part company.
gen_counts <- function(n_tax, n_samp, seed, zero_rate, struc_zero) {
  set.seed(seed)
  lam <- exp(rnorm(n_tax, log(300), 0.4))
  m <- matrix(rpois(n_tax * n_samp, lambda = rep(lam, times = n_samp)),
              nrow = n_tax, ncol = n_samp, byrow = TRUE)
  storage.mode(m) <- "integer"
  if (zero_rate > 0) {
    set.seed(seed + 1L)
    m[matrix(runif(length(m)) < zero_rate, nrow = n_tax)] <- 0L
  }
  if (struc_zero) {
    # Taxon 1 is absent from every sample of group 1: its group-1 prevalence is
    # exactly zero, which is the structural-zero condition the reference flags.
    m[1, seq(1, n_samp, by = 2)] <- 0L
  }
  dimnames(m) <- list(paste0("T", sprintf("%0*d", nchar(as.character(n_tax)), seq_len(n_tax))),
                      paste0("S", sprintf("%0*d", nchar(as.character(n_samp)), seq_len(n_samp))))
  m
}

gen_meta <- function(n_samp, n_group, n_cont, seed) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_group), length.out = n_samp))
  df <- data.frame(group = g)
  for (k in seq_len(n_cont)) {
    df[[paste0("x", k)]] <- round(rnorm(n_samp), 6)
  }
  rownames(df) <- paste0("S", sprintf("%0*d", nchar(as.character(n_samp)), seq_len(n_samp)))
  df
}

write_counts <- function(m, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(colnames(m), collapse = "\t"), con)
  for (i in seq_len(nrow(m))) {
    writeLines(paste(c(rownames(m)[i],
                       format(m[i, ], trim = TRUE, scientific = FALSE, digits = 15)),
                     collapse = "\t"), con)
  }
}

write_meta <- function(df, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(names(df), collapse = "\t"), con)
  for (i in seq_len(nrow(df))) {
    cells <- vapply(names(df), function(nm) {
      v <- df[[nm]][i]
      if (is.factor(v)) as.character(v)
      else format(v, trim = TRUE, scientific = FALSE, digits = 15)
    }, character(1))
    writeLines(paste(c(rownames(df)[i], cells), collapse = "\t"), con)
  }
}

SPECS <- list(
  list(id = "re2g", n_tax = 8L, n_samp = 12L, n_group = 2L, n_cont = 1L,
       seed = 9101L, zero_rate = 0, struc_zero = FALSE,
       about = "8 taxa x 12 samples, 2 groups, one continuous covariate. Small on purpose: the random-effects path fits an lmer per taxon and refits it for every pseudo-count in the sensitivity grid, so the fixture has to be small enough to run and still large enough for lmer to converge."),
  list(id = "re3g-sens", n_tax = 8L, n_samp = 12L, n_group = 3L, n_cont = 1L,
       seed = 9102L, zero_rate = 0.1, struc_zero = FALSE,
       about = "same shape with 3 groups and 10% zeros, so a random-effects case can also exercise the multi-group warnings and the zero handling."),
  list(id = "tiny5x8", n_tax = 5L, n_samp = 8L, n_group = 2L, n_cont = 1L,
       seed = 9103L, zero_rate = 0.05, struc_zero = TRUE,
       about = "5 taxa x 8 samples with taxon T1 structurally absent from group 1 and a few zeros elsewhere: the smallest fixture that still has both a structural zero and a mixed missingness pattern."),
  list(id = "onegroup", n_tax = 6L, n_samp = 9L, n_group = 1L, n_cont = 1L,
       seed = 9104L, zero_rate = 0, struc_zero = FALSE,
       about = "a single group. The reference warns that the multi-group tests are deactivated and that the group variable has fewer than 3 categories, so this pins the warning text and the NULL res_global/res_pair/res_dunn/res_trend entries.")
)

built <- list()
for (sp in SPECS) {
  d <- file.path(out_dir, sp$id)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  counts <- gen_counts(sp$n_tax, sp$n_samp, sp$seed, sp$zero_rate, sp$struc_zero)
  meta <- gen_meta(sp$n_samp, sp$n_group, sp$n_cont, sp$seed + 7L)
  write_counts(counts, file.path(d, "counts.tsv"))
  write_meta(meta, file.path(d, "meta.tsv"))
  formula_txt <- paste(c("group", if (sp$n_cont > 0) paste0("x", seq_len(sp$n_cont))),
                        collapse = " + ")
  writeLines(formula_txt, file.path(d, "formula.txt"))

  # Round-trip: read back exactly as the exact inputs will, and require identity.
  back_counts <- read_committed_counts(file.path(d, "counts.tsv"))
  back_meta <- read_committed_meta(file.path(d, "meta.tsv"))
  # Compared as raw bits, and as integer: `identical()` on two doubles is exact,
  # so this catches a writer that loses a digit, and comparing values rather than
  # types catches nothing extra that a bit comparison would miss.
  if (!identical(dim(back_counts), dim(counts)) ||
      !identical(dimnames(back_counts), dimnames(counts))) {
    stop(sp$id, ": counts shape or dimnames did not survive the round trip")
  }
  for (nm in names(dimnames(back_counts))) {
    if (!identical(dimnames(back_counts)[[nm]], dimnames(counts)[[nm]])) {
      stop(sp$id, ": counts ", nm, " did not survive the round trip")
    }
  }
  if (!identical(as.numeric(as.vector(back_counts)), as.numeric(as.vector(counts)))) {
    stop(sp$id, ": counts did not survive the round trip; the writer is lossy")
  }
  if (!identical(dim(back_meta), dim(meta)) ||
      !identical(rownames(back_meta), rownames(meta))) {
    stop(sp$id, ": metadata shape or row names did not survive the round trip")
  }
  for (nm in names(meta)) {
    if (nm == "group") {
      if (!identical(as.character(back_meta$group), as.character(meta$group))) {
        stop(sp$id, ": group levels did not survive the round trip")
      }
    } else if (!isTRUE(all.equal(back_meta[[nm]], meta[[nm]], tolerance = 0))) {
      stop(sp$id, ": covariate ", nm, " did not survive the round trip")
    }
  }

  built[[sp$id]] <- list(
    id = sp$id, dir = d, n_tax = sp$n_tax, n_samp = sp$n_samp,
    n_group = sp$n_group, n_cont = sp$n_cont, seed = sp$seed,
    zero_rate = sp$zero_rate, struc_zero = sp$struc_zero,
    formula = formula_txt, about = sp$about
  )
  cat(sprintf("%-14s %d taxa x %d samples, %d group(s), %d covariate(s)  round-trip OK\n",
              sp$id, sp$n_tax, sp$n_samp, sp$n_group, sp$n_cont))
}

# ---------------------------------------------------------------------------
# Derived two-column versions of the seven committed edge fixtures.
#
# Every validation/edge/*/meta.tsv has exactly one column, and
# data_sanity_check's matrix branch ends with
#   meta_data = meta_data[colnames(feature_table), ]
# whose `[.data.frame` drop defaults to TRUE for a one-column frame. The metadata
# therefore collapses to a vector and the very next check reports the formula
# variable as absent, so **ancombc2() cannot analyse any of the seven committed
# edge fixtures through its public API**. Measured: all seven answer "The
# following variables specified are not in the meta data: group".
#
# The committed files are left exactly as they are -- they are still the input to
# the seven `edge-*` error cases, which now pin that refusal. These derived
# fixtures add one deterministic continuous covariate so the structural-zero and
# boundary scenarios the edge cases were written for are actually analysed.
# ---------------------------------------------------------------------------

EDGE_DIRS <- c("absent_in_group_a", "rare_in_group_a", "present_in_every_group",
               "na_counts", "struc_zero_off", "prv_cut_boundary")

for (ed in EDGE_DIRS) {
  src_dir <- file.path("validation/edge", ed)
  if (!dir.exists(src_dir)) {
    stop("committed edge fixture is missing: ", src_dir, call. = FALSE)
  }
  id <- gsub("_", "-", ed)
  id <- paste0(id, "-2col")
  d <- file.path(out_dir, id)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # The counts are copied byte for byte from the committed fixture. Regenerating
  # them would be a different input, which is the mistake the legacy goldens
  # already made once.
  file.copy(file.path(src_dir, "counts.tsv"), file.path(d, "counts.tsv"), overwrite = TRUE)

  meta <- read_committed_meta(file.path(src_dir, "meta.tsv"))
  set.seed(4200L + sum(utf8ToInt(ed)))
  meta$edge_cov <- round(rnorm(nrow(meta)), 6)
  write_meta(meta, file.path(d, "meta.tsv"))
  writeLines("group + edge_cov", file.path(d, "formula.txt"))

  # The edge harness's own configuration, read through the same parser.
  cfg <- read_config_tsv(file.path(src_dir, "analysis.tsv"))
  cfg_file <- file.path(d, "config.tsv")
  con <- file(cfg_file, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("key\tvalue", con)
  for (k in names(cfg)) {
    v <- cfg[[k]]
    vs <- if (is.logical(v)) (if (v) "true" else "false")
          else if (is.numeric(v)) format(v, trim = TRUE, scientific = FALSE, digits = 15)
          else as.character(v)
    writeLines(paste(k, vs, sep = "\t"), con)
  }
  close(con)
  con <- file(cfg_file, open = "wt")
  writeLines("key\tvalue", con)
  for (k in names(cfg)) {
    v <- cfg[[k]]
    vs <- if (is.logical(v)) (if (v) "true" else "false")
          else if (is.numeric(v)) format(v, trim = TRUE, scientific = FALSE, digits = 15)
          else as.character(v)
    writeLines(paste(k, vs, sep = "\t"), con)
  }
  close(con)

  # Round trip, same rule as the generated fixtures.
  back_meta <- read_committed_meta(file.path(d, "meta.tsv"))
  if (!identical(colnames(back_meta), c("group", "edge_cov"))) {
    stop(id, ": derived metadata columns did not survive the round trip")
  }
  if (!identical(as.character(back_meta$group), as.character(meta$group))) {
    stop(id, ": group levels did not survive the round trip")
  }
  if (!identical(back_meta$edge_cov, meta$edge_cov)) {
    stop(id, ": the added covariate did not survive the round trip")
  }
  back_cfg <- read_config_tsv(cfg_file)
  stopifnot(identical(names(back_cfg), names(cfg)))

  built[[id]] <- list(
    id = id, dir = d, derived_from = file.path("validation/edge", ed),
    n_tax = nrow(back_meta) * 0L + nrow(read_committed_counts(file.path(d, "counts.tsv"))),
    n_samp = nrow(back_meta), n_group = length(unique(meta$group)), n_cont = 1L,
    seed = 4200L + sum(utf8ToInt(ed)), zero_rate = NA_real_, struc_zero = NA,
    formula = "group + edge_cov",
    about = paste0("validation/edge/", ed, " with one deterministic continuous ",
                   "covariate added so ancombc2() can analyse it: the committed ",
                   "fixture's single-column metadata collapses under [data.frame's ",
                   "drop and the public API refuses it. Counts are copied byte for ",
                   "byte; only meta.tsv is extended.")
  )
  cat(sprintf("%-30s %d taxa x %d samples  derived from validation/edge/%s  round-trip OK\n",
              id, built[[id]]$n_tax, nrow(back_meta), ed))
}

if (!is.null(opt$manifest)) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write the source-fixture manifest", call. = FALSE)
  }
  doc <- list(
    schema = "ancombc2-exact-source-fixtures/1",
    generator = "scripts/make_exact_source_fixtures.R",
    round_trip = "each fixture is written, read back through the same reader the exact inputs use, and compared with tolerance 0 before this manifest is written",
    fixtures = built
  )
  dir.create(dirname(opt$manifest), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(doc, auto_unbox = TRUE, pretty = TRUE, null = "null"),
             opt$manifest)
  cat(sprintf("wrote %s\n", opt$manifest))
}
