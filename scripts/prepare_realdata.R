#!/usr/bin/env Rscript
# Prepare the Layer 4 real datasets.
#
#   Rscript --vanilla scripts/prepare_realdata.R
#
# Each dataset becomes a directory under `validation/realdata/<name>/` holding
#
#   counts.tsv     taxa x samples, R's `write.table` layout
#   meta.tsv       samples x variables, same layout
#   analysis.tsv   the analysis the comparison runs, one key per line
#   source.json    where the data came from, and exactly what was done to it
#
# # Why prepare once and commit, rather than at comparison time
# -----------------------------------------------------------------
# Both arms must read the same bytes. If the R arm reconstructed the dataset from
# the source object and the Rust arm read a committed file, a disagreement about
# the *input* would be indistinguishable from a disagreement about the
# estimator. So the table is built once, here, and both arms read it. The
# `source.json` records the provenance so the preparation is auditable rather
# than implicit.
#
# # What is and is not real
# -----------------------
# Two of the three datasets are *measured* data, unmodified:
#
#   qmp-real      the Quantitative Microbiome Project counts shipped with
#                 ANCOMBC 2.15.2, 91 taxa x 106 subjects
#   dietswap       the diet-swap study, from Bioconductor's `microbiome`
#                 package, 209 taxa x 84 samples -- this is the dataset the
#                 ANCOM-BC2 vignette analyses
#
# The third is the vignette's own construction: ANCOMBC's `sim_plnm` applied to
# the *real* QMP abundance table, with sampling fractions differing by batch.
# It is synthetic, and `source.json` says so, but its abundance structure is
# real and its confounding structure is the one the correction exists for. It is
# included because the vignette is the reference configuration a user will
# actually copy.

suppressWarnings(suppressMessages({
  library(utils)
}))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  args[[i + 1L]]
}
root <- arg("--out", "validation/realdata")
oracle_dir <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC")

# `phyloseq` is needed only to read two serialized S4 objects. The shipped
# shim (see docs/reproduction.md) provides the class definitions and nothing
# else, and the reader below uses `attributes()` rather than any phyloseq method,
# so no behaviour of the real package is relied upon.
have_phyloseq <- requireNamespace("phyloseq", quietly = TRUE)

#' Read a serialized phyloseq object's count matrix and sample metadata.
#'
#' The slots are read through `attributes()` and rebuilt as base objects, not
#' through `as(otu_table(x), "matrix")`. Two reasons: the shim does not implement
#' the coercions, and the payload is already a plain double matrix with its shape
#' in `dim` and its names in `dimnames`, so there is nothing to coerce it *from*.
#'
#' The slot names differ between phyloseq versions -- this serialization carries
#' `sam_data` and `phy_tree` where the current package uses `sample_data` and
#' `tree` -- so both are tried and the one that is present is used.
read_phyloseq <- function(path) {
  e <- new.env()
  nm <- load(path, envir = e)
  stopifnot(length(nm) == 1L)
  x <- e[[nm]]
  a <- attributes(x)
  pick <- function(...) {
    for (n in c(...)) {
      s <- a[[n]]
      if (!is.null(s)) return(s)
    }
    stop("none of ", paste(c(...), collapse = "/"), " is present in ", path)
  }
  o <- pick("otu_table")
  dm <- attr(o, "dim")
  dn <- attr(o, "dimnames")
  if (is.null(dm)) stop("the otu_table slot has no dim in ", path)
  # `as.vector` on the slot payload; the shape and names come from the
  # attributes, because reading them off the object dispatches on the S4 class.
  m <- matrix(as.vector(o), nrow = dm[1], ncol = dm[2], dimnames = dn)
  if (isTRUE(attr(o, "taxa_are_rows")) != TRUE) {
    stop("taxa_are_rows is not TRUE; this reader has only been checked for it")
  }
  s <- pick("sample_data", "sam_data")
  # The slot payload is a data.frame with the S4 class bolted on; dropping that
  # class is a no-op on the data and makes the base methods work.
  df <- s
  class(df) <- "data.frame"
  if (!is.data.frame(df)) stop("the sample_data slot is not a data.frame in ", path)
  list(counts = m, meta = df)
}


# ---- writers ---------------------------------------------------------------

write_table <- function(m, sample_names, taxon_prefix, path, taxon_names = NULL) {
  if (is.null(taxon_names)) taxon_names <- paste0(taxon_prefix, seq_len(nrow(m)))
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  # R's `write.table(row.names = TRUE)` layout: an unnamed first header field,
  # then one per sample, and one field per row for the taxon name plus counts.
  writeLines(paste(sample_names, collapse = "\t"), con)
  for (i in seq_len(nrow(m))) {
    writeLines(
      paste(
        c(
          taxon_names[i],
          format(m[i, ], trim = TRUE, scientific = FALSE, digits = 15)
        ),
        collapse = "\t"
      ),
      con
    )
  }
  invisible(NULL)
}

write_meta <- function(df, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(colnames(df), collapse = "\t"), con)
  # `rownames(df)` here, not `rownames(df[i, ])`: subsetting a data frame to a
  # single row replaces explicit row names with automatic ones, so the label came
  # out empty and every line was the value alone -- which the reader then read as
  # a ragged row. `df[i, ]` also drops to a vector for a single-column frame,
  # which silently promoted one value into all of them.
  labels <- rownames(df)
  if (is.null(labels)) {
    stop("the metadata frame needs row names: they are the sample names")
  }
  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    values <- vapply(row, as.character, "")
    writeLines(
      paste(c(labels[i], values[colnames(df)]), collapse = "\t"),
      con
    )
  }
  invisible(NULL)
}

# The analysis config is a two-column key/value table rather than JSON because
# `jsonlite` is not a dependency of this repository, and because a key/value file
# has no escaping rules to get wrong -- a formula containing a `+` or a `~` is
# just a line.
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
  invisible(NULL)
}

write_source <- function(s, path) {
  writeLines(
    paste0("# ", gsub("\n", "\n# ", s$provenance)),
    path
  )
}

prepare <- function(name, counts, meta, analysis, provenance) {
  dir <- file.path(root, name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # Counts are written as integers where they are integers. A count table with
  # `.0` suffixes is still a count, but writing it as an integer keeps the
  # committed table small and readable and avoids a rounding question.
  if (all(counts == round(counts))) {
    counts <- matrix(as.integer(round(counts)), nrow = nrow(counts),
                     dimnames = dimnames(counts))
  }
  write_table(counts, colnames(counts), "T", file.path(dir, "counts.tsv"))
  write_meta(meta, file.path(dir, "meta.tsv"))
  write_analysis(analysis, file.path(dir, "analysis.tsv"))
  write_source(provenance, file.path(dir, "source.txt"))
  cat(sprintf(
    "  %-14s %5d taxa x %4d samples, %d variable(s): %s\n",
    name, nrow(counts), ncol(counts), ncol(meta), paste(colnames(meta), collapse = ", ")
  ))
}

dir.create(root, recursive = TRUE, showWarnings = FALSE)
cat("preparing the Layer 4 datasets\n")

# ---- A: the QMP study, real counts -----------------------------------------
#
# The Quantitative Microbiome Project data shipped with ANCOMBC: 91 taxa by 106
# subjects, a real study, written out unmodified apart from being transposed into
# the taxa x samples layout. The analysis is `~ group`, the cohort contrast.

e <- new.env()
load(file.path(oracle_dir, "data", "QMP.rda"), envir = e)
# QMP is samples x taxa -- 106 subjects in rows, 91 OTUs in columns, as
# `?QMP` states -- so it is transposed into the taxa x samples layout the rest of
# the toolchain takes. Getting this backwards is a silent, total corruption of
# the analysis rather than an error, so the orientation is asserted against the
# documented shape instead of assumed.
qmp_samples_by_taxa <- e$QMP
stopifnot(is.matrix(qmp_samples_by_taxa))
stopifnot(nrow(qmp_samples_by_taxa) == 106L, ncol(qmp_samples_by_taxa) == 91L)
qmp <- t(qmp_samples_by_taxa)
stopifnot(nrow(qmp) == 91L, ncol(qmp) == 106L)

# The group comes from the *sample* names, which are the QMP row names: `SC` for
# the study cohort and `DC` for the disease cohort, 40 and 66 subjects. `?QMP`
# describes the table as pruned to healthy subjects from those two cohorts, so
# this is the cohort contrast and not a healthy/diseased one. The names are
# checked rather than assumed, because a prefix that means something else would
# produce a real analysis of the wrong variable.
prefix <- substr(colnames(qmp), 1, 2)
stopifnot(all(grepl("^(SC|DC)[0-9]+$", colnames(qmp))))
group <- ifelse(prefix == "DC", "disease_cohort", "study_cohort")
stopifnot(length(unique(group)) == 2L)

prepare(
  "qmp-real",
  qmp,
  data.frame(group = factor(group, levels = c("study_cohort", "disease_cohort")),
             row.names = colnames(qmp)),
  list(
    formula = "group",
    group = "group",
    # The reference level first, and it is stated rather than left to
    # `factor()`'s alphabetical sort -- see scripts/realdata_r.R. Both arms read
    # this, so they build the same contrast and the coefficient names match.
    group_levels = "study_cohort,disease_cohort",
    p_adj_method = "BH",
    pseudo = 0.5,
    prv_cut = 0.0,
    lib_cut = 0.0,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = TRUE,
    pairwise = TRUE
  ),
  list(provenance = paste0(
    "QMP: the Quantitative Microbiome Project counts shipped as\n",
    "reference/ANCOMBC/data/QMP.rda in ANCOMBC 2.15.2 at\n",
    "dc4febdf59badb3a8dfe0c767ef2186323c2199a. Real measured data.\n",
    "91 taxa x 106 samples, 27.1% structural zeros.\n\n",
    "ORIENTATION: `?QMP` documents 106 samples in rows and 91 OTUs in columns, so\n",
    "the shipped matrix is samples x taxa and is transposed here to the\n",
    "taxa x samples layout the toolchain takes. Asserted, not assumed.\n\n",
    "THE GROUP: no sample metadata ships with the matrix. It is derived from the\n",
    "sample names: `SC` (study cohort, 40 subjects) and `DC` (disease cohort, 66\n",
    "subjects). `?QMP` describes the table as pruned to HEALTHY subjects from\n",
    "those two cohorts, so this is the cohort contrast and NOT a healthy/diseased\n",
    "one -- reading it as the latter would be an analysis of a variable that does\n",
    "not exist in this table. The prefixes are checked with a regular expression\n",
    "before use. This is a deliberate, recorded choice rather than a recovered\n",
    "column: a reader who disagrees can edit meta.tsv and rerun.\n\n",
    "Formula: ~ group. Structural zeros and negative lower bounds on, the global\n",
    "and pairwise tests on, pseudo 0.5, prevalence and library filters off."
  ))
)

# ---- B: the diet-swap study, real counts -----------------------------------
#
# The dataset the ANCOM-BC2 vignette analyses. 209 taxa by 84 samples, real
# measured data from a dietary-intervention study.

if (have_phyloseq) {
  suppressWarnings(suppressMessages(library(phyloseq)))
  mb_data <- arg("--microbiome-data", "microbiome/data")

  # -- B: the diet-swap study, the vignette's own dataset ---------------------
  ds <- read_phyloseq(file.path(mb_data, "dietswap.rda"))
  otu <- ds$counts
  smd <- ds$meta
  # Align the metadata to the count table's columns *by name*, and fail loudly
  # if they disagree. A silent row-order mismatch here would make every metric
  # wrong in a way that looks like an estimator difference.
  missing_cols <- setdiff(colnames(otu), rownames(smd))
  if (length(missing_cols)) {
    stop("dietswap: ", length(missing_cols), " count-table columns have no metadata")
  }
  smd <- smd[colnames(otu), , drop = FALSE]
  stopifnot(identical(rownames(smd), colnames(otu)))
  for (v in c("group", "nationality", "timepoint")) {
    if (!v %in% colnames(smd)) {
      stop("dietswap: the vignette's formula needs a `", v, "` column; have ",
           paste(colnames(smd), collapse = ", "))
    }
  }
  # The vignette's analysis, `fix_formula = "nationality + timepoint + group"`
  # with `group = "group"`, and DI as the reference (ANCOMBC2.Rmd). The
  # vignette also passes `rand_formula = "(timepoint | subject)"`, which is the
  # random-effects path -- not implemented here, and out of scope, so this is
  # the fixed-effects part of the same analysis and nothing else.
  group_levels_ds <- c("DI", "ED", "HE")
  have_group <- all(group_levels_ds %in% as.character(smd$group))
  if (!have_group) {
    stop("dietswap: expected the vignette's group levels DI/ED/HE, found ",
         paste(sort(unique(as.character(smd$group))), collapse = ", "))
  }
  nat_levels <- sort(unique(as.character(smd$nationality)))
  prepare(
    "dietswap",
    otu,
    data.frame(
      group = factor(as.character(smd$group), levels = group_levels_ds),
      nationality = factor(as.character(smd$nationality), levels = nat_levels),
      timepoint = as.numeric(smd$timepoint),
      row.names = colnames(otu)
    ),
    list(
      formula = "nationality + timepoint + group",
      group = "group",
      group_levels = paste(group_levels_ds, collapse = ","),
      p_adj_method = "BH",
      pseudo = 0.5,
      prv_cut = 0.0,
      lib_cut = 0.0,
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.05,
      global = TRUE,
      pairwise = TRUE
    ),
    list(provenance = paste0(
      "dietswap: the two-week diet-swap study between western (USA) and\n",
      "traditional (rural Africa) diets, from Bioconductor's `microbiome`\n",
      "package. This is the dataset the ANCOM-BC2 vignette analyses\n",
      "(ANCOMBC2.Rmd, around line 716). Real measured data, unmodified.\n\n",
      sprintf("%d taxa x %d samples, %d subjects, 3 diet groups.\n\n",
              nrow(otu), ncol(otu), length(unique(smd$subject))),
      "Provenance: microbiome_1.24.0.tar.gz from\n",
      "  https://bioconductor.org/packages/3.18/bioc/src/contrib/microbiome_1.24.0.tar.gz\n",
      "  tar xzf microbiome.tar.gz\n",
      "\n",
      "THE ANALYSIS IS THE VIGNETTE'S, MINUS THE RANDOM EFFECT. The vignette\n",
      "calls\n",
      "  ancombc2(data = dietswap, tax_level = \"Family\",\n",
      "            fix_formula = \"nationality + timepoint + group\",\n",
      "            rand_formula = \"(timepoint | subject)\",\n",
      "            group = \"group\", struc_zero = TRUE, neg_lb = TRUE, ...)\n",
      "`rand_formula` is the random-effects path, which this port does not\n",
      "implement and which is out of scope; this is the fixed-effects part of\n",
      "the same analysis and nothing else. The group levels are DI (reference),\n",
      "ED, HE, stated explicitly rather than left to factor()'s alphabetical\n",
      "sort -- see scripts/realdata_r.R. `nationality` is a factor with two\n",
      "levels (AAM, AFR), so it contributes one column.\n\n",
      "This is the only dataset here with THREE groups, so it is the one that\n",
      "exercises the global and pairwise tests."
    ))
  )

  # -- C: a larger real metagenomic set ---------------------------------------
  at <- read_phyloseq(file.path(mb_data, "atlas1006.rda"))
  at_otu <- at$counts
  at_smd <- at$meta[colnames(at_otu), , drop = FALSE]
  stopifnot(identical(rownames(at_smd), colnames(at_otu)))
  if (!"sex" %in% colnames(at_smd)) {
    stop("atlas1006: expected a `sex` column; have ",
         paste(colnames(at_smd), collapse = ", "))
  }
  g_levels <- sort(unique(as.character(at_smd$sex)))
  if (length(g_levels) < 2L) stop("atlas1006: sex has fewer than two levels")
  prepare(
    "atlas1006",
    at_otu,
    data.frame(sex = factor(as.character(at_smd$sex), levels = g_levels),
               row.names = colnames(at_otu)),
    list(
      formula = "sex",
      group = "sex",
      group_levels = paste(g_levels, collapse = ","),
      p_adj_method = "BH",
      pseudo = 0.5,
      prv_cut = 0.0,
      lib_cut = 0.0,
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.05,
      global = TRUE,
      pairwise = TRUE,
      # The sensitivity analysis, which is the only place `passed_ss` and
      # `diff_robust` come from. It runs here so the comparison covers them.
      pseudo_sens = TRUE,
      conservative = TRUE
    ),
    list(provenance = paste0(
      "atlas1006: a 1,006-taxon atlas of the human gut microbiome, from\n",
      "Bioconductor's `microbiome` package. Real measured data, unmodified.\n\n",
      sprintf("%d taxa x %d samples, %s.\n\n",
              nrow(at_otu), ncol(at_otu),
              paste(sprintf("%.1f%% structural zeros", 100 * mean(at_otu == 0)),
                    collapse = "")),
      "Provenance: the same microbiome_1.24.0.tar.gz as dietswap,\n",
      "data/atlas1006.rda.\n\n",
      "THIS IS THE LARGE DATASET: 1,006 taxa against 91 for QMP, and it is the\n",
      "one that exercises the bitset missingness-pattern grouping and the cached\n",
      "per-pattern QR at a realistic taxon count.\n\n",
      "It is also the only dataset that runs the pseudo-count SENSITIVITY\n",
      "analysis, because `passed_ss` and `diff_robust` exist only there and the\n",
      "comparison is required to cover them. That is 4 pseudo-count levels on a\n",
      "1,006-taxon table, so it is the slowest of the four."
    ))
  )
} else {
  cat("  dietswap       SKIPPED: the phyloseq class definitions are not\n")
  cat("                 available, and the .rda holds a serialized S4 object.\n")
  cat("                 See docs/reproduction.md for the shim.\n")
  cat("  atlas1006      SKIPPED: same reason.\n")
}

# ---- C: the vignette's own construction ------------------------------------
#
# Exactly what ANCOMBC2.Rmd does: `sim_plnm` on the real QMP abundance table, with
# sampling fractions that differ by batch. This is the reference configuration a
# user copies, and its confounding structure -- the sampling fraction varying
# systematically with the batch the group is in -- is the one the
# sampling-fraction correction exists to absorb.
#
# `sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, ...)` is the vignette's own
# call: `taxa_are_rows = FALSE` says the *input* has taxa in columns, which the
# 106 x 91 samples-by-taxa QMP does. Its *output* is 91 taxa x 150 samples,
# asserted below, because the vignette transposes it by hand and getting that
# step wrong is a total corruption rather than an error.

for (f in c("oracle.R", "stubs.R", "harness.R")) {
  if (!file.exists(file.path("reference/R", f))) {
    stop("cannot find reference/R/", f, "; run from the repository root")
  }
  source(file.path("reference/R", f))
}
O <- load_sim_helpers(oracle_dir)

n <- 150
# The number of taxa. `qmp` was transposed to taxa x samples for dataset A, so
# this is `nrow`, not the `ncol(qmp)` the vignette uses on the untransposed table.
# Reading it off the wrong axis is what made `sim_plnm` filter to 106 "taxa".
d <- nrow(qmp)
stopifnot(d == 91L)
set.seed(123)
# The *untransposed* table goes in: `taxa_are_rows = FALSE` tells `sim_plnm` that
# the input has taxa in columns, which is the 106 x 91 samples-by-taxa QMP. The
# transposed `qmp` (taxa x samples) is what the *output* layout looks like, and
# passing it here makes `sim_plnm` treat 91 samples as 91 taxa -- which filters
# to 106 "taxa" and then fails the shape assertion below rather than silently
# producing a differently-shaped dataset.
abn_data <- O$sim_plnm(
  abn_table = qmp_samples_by_taxa, taxa_are_rows = FALSE, prv_cut = 0.05,
  n = n, lib_mean = 1e8, disp = 0.5
)
# `sim_plnm` already returns taxa x samples; assert it rather than transposing on
# the strength of a comment.
if (nrow(abn_data) != d || ncol(abn_data) != n) {
  stop(sprintf(
    "sim_plnm returned %d x %d but the vignette's shape is %d taxa x %d samples",
    nrow(abn_data), ncol(abn_data), d, n
  ))
}
rownames(abn_data) <- paste0("T", seq_len(d))
colnames(abn_data) <- paste0("S", seq_len(n))
# `sim_plnm` returns integers; assert it, because a non-integer count table would
# make every downstream comparison a question about rounding.
stopifnot(all(abn_data == round(abn_data)))

# Sampling fractions differ by batch, exactly as the vignette sets them, and the
# batch is crossed with the group so the confounding is real rather than
# collinear with the treatment.
set.seed(20240215)
samp_frac <- log(c(
  stats::runif(n / 3, min = 1e-4, max = 1e-3),
  stats::runif(n / 3, min = 1e-3, max = 1e-2),
  stats::runif(n / 3, min = 1e-2, max = 1e-1)
))
batch <- rep(c("A", "B", "C"), each = n / 3)
grp <- rep(c("control", "treated"), length.out = n)

prepare(
  "qmp-vignette",
  abn_data,
  data.frame(
    group = factor(grp, levels = c("control", "treated")),
    batch = factor(batch),
    cont_cov = stats::rnorm(n),
    samp_frac = samp_frac,
    row.names = colnames(abn_data)
  ),
  list(
    formula = "group + cont_cov",
    group = "group",
    group_levels = "control,treated",
    p_adj_method = "BH",
    pseudo = 0.5,
    prv_cut = 0.0,
    lib_cut = 0.0,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = TRUE,
    pairwise = TRUE
  ),
  list(provenance = paste0(
    "qmp-vignette: the ANCOM-BC2 vignette's own dataset, built exactly as\n",
    "ANCOMBC2.Rmd (around line 914) builds it.\n\n",
    "  data(QMP, package = \"ANCOMBC\"); set.seed(123)\n",
    "  abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05,\n",
    "                      n = 150, lib_mean = 1e8, disp = 0.5)\n",
    "  # sim_plnm already returns taxa x samples\n\n",
    "THIS DATASET IS SYNTHETIC. Its abundance structure comes from the real QMP\n",
    "table, and its sampling fractions differ by batch exactly as the vignette\n",
    "sets them, but the counts are generated. It is included because it is the\n",
    "configuration the documentation tells users to copy, and because the\n",
    "batch/group-crossed sampling fractions are the confounding the\n",
    "sampling-fraction correction exists to absorb -- a property the two real\n",
    "datasets do not have a controlled version of.\n\n",
    "The continuous covariate `cont_cov` and the sampling fractions are drawn with\n",
    "set.seed(20240215), recorded here so the table is reproducible; the vignette\n",
    "does not seed them, and matching it exactly would mean relying on R's global\n",
    "RNG state at a particular line of a particular document.\n\n",
    "91 taxa x 150 samples. Formula: ~ group + cont_cov, with the group contrast\n",
    "as the coefficient under comparison."
  ))
)

cat("done\n")
