#!/usr/bin/env Rscript
# Layer 4 -- real-data validation: the R arm.
#
#   Rscript --vanilla scripts/realdata_r.R --out results.r.jsonl
#   Rscript --vanilla scripts/realdata_r.R --dataset qmp-real --out one.jsonl
#
# Runs the pinned oracle on each prepared dataset under `validation/realdata/`
# and writes one JSONL row per dataset, in the schema
# `scripts/compare_realdata.py` reads. The Rust arm reads the *same* counts,
# metadata and analysis config from the *same* directory, so the comparison is
# between implementations on identical input, which is the only form in which a
# ">= 99.99% diff_abn agreement" statement means anything.
#
# Why the datasets are read as committed files
# ---------------------------------------------
# `scripts/prepare_realdata.R` builds them once from their sources and they are
# committed. Re-deriving them here would mean two independent reconstructions of
# the same table, and a disagreement about the *input* would be indistinguishable
# from a disagreement about the estimator. Each dataset's `source.txt` records
# exactly what was done to it.
#
# With no `--dataset`, every prepared dataset is analysed -- the same default as
# the Rust arm's `realdata` subcommand, so one invocation compares the same set
# on both sides.

suppressWarnings(suppressMessages({
  library(utils)
}))

# ---- arguments --------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  if (i == length(args)) stop("--", name, " needs a value")
  args[[i + 1L]]
}
dataset <- arg("--dataset")
data_root <- arg("--data", "validation/realdata")
out_path <- arg("--out")
if (is.null(out_path)) {
  stop("usage: realdata_r.R --out results.r.jsonl [--data DIR] [--dataset NAME]")
}

dirs <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)
dirs <- dirs[file.exists(file.path(dirs, "counts.tsv"))]
if (!length(dirs)) {
  stop("no prepared datasets under ", data_root,
       "\nrun scripts/prepare_realdata.R first")
}
datasets <- if (!is.null(dataset)) dataset else sort(basename(dirs))
for (d in datasets) {
  if (!dir.exists(file.path(data_root, d))) {
    stop("no such dataset directory: ", file.path(data_root, d),
         "\navailable: ", paste(sort(basename(dirs)), collapse = ", "))
  }
}

# ---- the pinned oracle ------------------------------------------------------

oracle_dir <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC")
harness_dir <- "reference/R"
for (f in c("oracle.R", "stubs.R", "harness.R")) {
  if (!file.exists(file.path(harness_dir, f))) {
    stop("cannot find ", file.path(harness_dir, f), "; run from the repository root")
  }
  source(file.path(harness_dir, f))
}
H <- tryCatch(
  load_harness(oracle_dir),
  error = function(e) {
    stop("cannot load the oracle from ", oracle_dir, ": ", conditionMessage(e),
         "\nSet ANCOMBC_ORACLE_DIR to the pinned ANCOMBC 2.15.2 checkout.")
  }
)

# ---- a minimal JSON writer --------------------------------------------------
#
# `jsonlite` is not a dependency of this repository and the row schema needs
# numbers, strings, booleans and one array of strings, so it is written by hand.
# `NA` becomes `null` so a quantity the oracle did not produce is absent rather
# than zero -- the difference between "not computed" and "computed as nothing",
# which is exactly what a parity comparison must not confuse.

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
# A JSON string. Taxon names are identifiers here, but the escaper handles
# backslashes and the three control characters anyway: an unescaped newline in a
# name would make the whole row unparseable, which loses the dataset rather than
# reporting it.
jstrs <- function(x) {
  x <- as.character(x)
  if (length(x) == 0L) return("[]")
  one <- function(s) {
    s <- gsub("\\", "\\\\", s, fixed = TRUE)
    s <- gsub("\"", "\\\"", s, fixed = TRUE)
    s <- gsub("\n", "\\n", s, fixed = TRUE)
    s <- gsub("\r", "\\r", s, fixed = TRUE)
    s <- gsub("\t", "\\t", s, fixed = TRUE)
    paste0("\"", s, "\"")
  }
  paste0("[", paste(vapply(x, one, ""), collapse = ","), "]")
}
# A vector of logicals as JSON booleans. `jvec` would write 0/1, and the
# comparison script requires booleans: a flag serialised as a number is a schema
# difference, and coercing it silently on the way in would let a genuine mismatch
# through as a match.
jbools <- function(x) {
  x <- as.logical(x)
  if (length(x) == 0L) return("[]")
  paste0("[", paste(ifelse(is.na(x), "null", ifelse(x, "true", "false")),
                    collapse = ","), "]")
}
# Flatten a message to one line and neutralise quotes: the oracle's condition
# messages carry newlines, and a newline inside a JSON string makes the row
# unparseable, which loses the dataset rather than reporting the failure.
flat <- function(s) gsub("[\"\r\n\t]", " ", trimws(as.character(s)))

# ---- the analysis config ----------------------------------------------------

# A two-column key/value table, written by `scripts/prepare_realdata.R` without
# a JSON dependency on either side. A key/value file also has no escaping rules
# to get wrong, so a formula containing `~` or `+` is just a line.
read_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")]
  lines <- lines[!startsWith(lines, "key")]
  cfg <- list()
  for (l in lines) {
    kv <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(kv) < 2L) next
    cfg[[trimws(kv[1])]] <- trimws(paste(kv[-1], collapse = "\t"))
  }
  if (is.null(cfg$formula)) stop(path, ": no `formula` key")
  cfg
}

# ---- one dataset ------------------------------------------------------------

run_one <- function(dataset, dir, cfg, emit) {
  counts <- as.matrix(utils::read.delim(file.path(dir, "counts.tsv"), row.names = 1,
                                        check.names = FALSE))
  # A matrix, not a data frame: the oracle's internal helpers dispatch on the
  # type, and a data frame makes `is.infinite()` fail with "default method not
  # implemented for type 'list'" deep inside the structural-zero scan. The
  # benchmark arm coerces for the same reason.
  storage.mode(counts) <- "double"
  counts[is.na(counts)] <- 0
  meta <- utils::read.delim(file.path(dir, "meta.tsv"), row.names = 1,
                            check.names = FALSE)

  # The group factor's level order, and so which level is the reference.
  #
  # This must be pinned rather than left to R's default. `factor(x)` sorts its
  # levels alphabetically, so the reference would be whichever sorts first, and
  # the two arms would build *different* contrasts: on a `disease_cohort` /
  # `study_cohort` variable R would drop `disease_cohort` and report
  # `groupstudy_cohort`, while a reader assuming the order of appearance would
  # report the opposite. The coefficients would then agree in name only by
  # accident and disagree in sign by construction.
  if (!is.null(cfg$group) && !is.null(cfg$group_levels)) {
    lv <- trimws(strsplit(as.character(cfg$group_levels), ",", fixed = TRUE)[[1]])
    meta[[cfg$group]] <- factor(as.character(meta[[cfg$group]]), levels = lv)
  } else if (!is.null(cfg$group)) {
    meta[[cfg$group]] <- factor(as.character(meta[[cfg$group]]))
  }

  t0 <- Sys.time()
  res <- try(suppressWarnings(suppressMessages(
    H$ref_run(
      counts,
      meta,
      as.character(cfg$formula),
      p_adj_method = if (is.null(cfg$p_adj_method)) "BH" else as.character(cfg$p_adj_method),
      pseudo = if (is.null(cfg$pseudo)) 0.5 else as.numeric(cfg$pseudo),
      pseudo_sens = identical(cfg$pseudo_sens, "true"),
      conservative = if (is.null(cfg$conservative)) TRUE else identical(cfg$conservative, "true"),
      prv_cut = if (is.null(cfg$prv_cut)) 0 else as.numeric(cfg$prv_cut),
      lib_cut = if (is.null(cfg$lib_cut)) 0 else as.numeric(cfg$lib_cut),
      s0_perc = if (is.null(cfg$s0_perc)) 0.05 else as.numeric(cfg$s0_perc),
      group = if (is.null(cfg$group)) NULL else as.character(cfg$group),
      struc_zero = identical(cfg$struc_zero, "true"),
      neg_lb = identical(cfg$neg_lb, "true"),
      alpha = if (is.null(cfg$alpha)) 0.05 else as.numeric(cfg$alpha),
      global = identical(cfg$global, "true"),
      pairwise = identical(cfg$pairwise, "true")
    )
  )), silent = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(res, "try-error")) {
    emit(paste0(
      '{"arm":"r","dataset":"', dataset,
      '","error":"', flat(res), '","elapsed_seconds":', jnum(elapsed), "}"
    ))
    cat("r arm on ", dataset, ": FAILED -- ", flat(res), "\n", sep = "")
    return(invisible(NULL))
  }

  taxa <- as.character(res$taxa_retained)
  coef_names <- colnames(res$beta)

  # The coefficient under comparison, by name: `<variable><level>` for the level
  # that is not the reference. Naming it explicitly rather than taking "the first
  # non-intercept column" or "the last column" is what makes the two arms
  # comparable on a design with a covariate -- a positional fallback compares
  # `cont_cov` on one side and the group on the other, and the signs then disagree
  # by construction while both arms report success.
  group_col <- NA_integer_
  if (!is.null(cfg$group)) {
    gv <- as.character(cfg$group)
    lv <- if (!is.null(cfg$group_levels)) {
      trimws(strsplit(as.character(cfg$group_levels), ",", fixed = TRUE)[[1]])
    } else {
      levels(meta[[gv]])
    }
    for (i in seq_along(lv)[-1]) {
      cand <- match(paste0(gv, lv[[i]]), coef_names)
      if (!is.na(cand)) {
        group_col <- cand
        break
      }
    }
  }
  if (is.na(group_col)) {
    non_int <- which(coef_names != "(Intercept)")
    if (!length(non_int)) {
      emit(paste0('{"arm":"r","dataset":"', dataset,
                  '","error":"the design has no coefficient to compare"}'))
      return(invisible(NULL))
    }
    group_col <- non_int[[1]]
  }

  # ---- the sensitivity flags ------------------------------------------------
  #
  # `ref_run` returns `ss_tab` and not the flags, and the oracle has TWO
  # definitions of `passed_ss`, selected by `conservative`:
  #
  #   conservative (ancombc2.R:556-564): `ss_prim_<coef>` is the FRACTION of the
  #   pseudo-count refits that flagged the taxon, and a taxon is robust when that
  #   fraction is 0 or 1 -- flagged by every refit or by none. The main run's
  #   p-value does not enter.
  #
  #   non-conservative (`flag_fun`, ancombc2.R:711): an agreement test against
  #   the main run's p-value,
  #   `(ss == 0 & p <= alpha) | (ss == 1 & p > alpha)`.
  #
  # Using one for the other is a silent 14% disagreement on `atlas1006`, and it is
  # invisible in `diff_robust` -- that column is `diff_abn & passed_ss`, and
  # `diff_abn` is FALSE for every affected taxon -- so `diff_robust` agreeing at
  # 100% while `passed_ss` agrees at 86% is the signature of exactly this bug.
  passed_ss <- rep(NA, length(taxa))
  diff_robust <- rep(NA, length(taxa))
  alpha_v <- if (is.null(cfg$alpha)) 0.05 else as.numeric(cfg$alpha)
  if (identical(cfg$pseudo_sens, "true") && !is.null(res$ss_tab)) {
    ss <- as.data.frame(res$ss_tab)
    want <- paste0("ss_prim_", coef_names[group_col])
    col <- if (want %in% colnames(ss)) {
      want
    } else {
      cand <- grep("^ss_prim_", colnames(ss), value = TRUE)
      if (!length(cand)) NULL else cand[1]
    }
    if (is.null(col)) {
      emit(paste0(
        '{"arm":"r","dataset":"', dataset, '","error":"',
        flat(paste0("ss_tab has no column for `", want, "`; it has ",
                    paste(colnames(ss), collapse = ", "))), '"}'
      ))
      return(invisible(NULL))
    }
    # `ss_tab` may be stored as a factor, and `as.numeric` on a factor is NA for
    # every value, so the labels are converted rather than the codes.
    raw_s <- as.character(ss[[col]])
    s <- suppressWarnings(as.numeric(raw_s))
    s[is.na(s)] <- suppressWarnings(as.numeric(as.integer(raw_s)))
    if (identical(cfg$conservative, "false")) {
      p_main <- as.numeric(res$p[, group_col])
      p_main[is.na(p_main)] <- 1
      s[is.na(s)] <- 1
      passed_ss <- (s == 0 & p_main <= alpha_v) | (s == 1 & p_main > alpha_v)
    } else {
      frac <- s
      frac[is.na(frac)] <- 0
      passed_ss <- frac == 0 | frac == 1
    }
    passed_ss <- as.logical(passed_ss)
    diff_robust <- as.logical(res$diff_abn[, group_col]) & passed_ss
  }

  emit(paste0(
    '{"arm":"r","dataset":"', dataset, '"',
    ',"n_taxa_in":', jnum(nrow(counts)),
    ',"n_samp_in":', jnum(ncol(counts)),
    ',"n_taxa_retained":', jnum(length(taxa)),
    ',"coefficient":"', coef_names[group_col], '"',
    ',"n_coefficients":', jnum(ncol(res$beta)),
    ',"elapsed_seconds":', jnum(elapsed),
    ',"taxa":', jstrs(taxa),
    ',"beta":', jvec(res$beta[, group_col]),
    ',"se":', jvec(res$se[, group_col]),
    ',"p":', jvec(res$p[, group_col]),
    ',"q":', jvec(res$q[, group_col]),
    ',"diff_abn":', jbools(res$diff_abn[, group_col]),
    ',"passed_ss":', jbools(passed_ss),
    ',"diff_robust":', jbools(diff_robust),
    ',"error":null}',
    sep = ""
  ))
  cat(sprintf("r arm on %s: %d taxa retained, %.1fs\n", dataset, length(taxa), elapsed))
  invisible(NULL)
}

# ---- the run ----------------------------------------------------------------

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
con <- file(out_path, open = "wt")
on.exit(close(con), add = TRUE)
emit <- function(line) writeLines(line, con)

for (dataset in datasets) {
  dir <- file.path(data_root, dataset)
  run_one(dataset, dir, read_config(file.path(dir, "analysis.tsv")), emit)
}
cat(sprintf("r arm: %d dataset(s) -> %s\n", length(datasets), out_path))
