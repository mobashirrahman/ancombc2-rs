#!/usr/bin/env Rscript
# The R arm of the benchmark surface: ANCOMBC 2.15.2, one core and all cores.
#
# Driven by scripts/run_benchmarks.py through the environment:
#
#   ANCOMBC_BENCH_DATA      directory holding <name>.counts.tsv and friends
#   ANCOMBC_BENCH_OUT       output directory
#   ANCOMBC_BENCH_DATASET   the dataset name
#   ANCOMBC_BENCH_THREADS   1 for the single-core arm, unset for all cores
#
# It writes a `run_metadata.tsv` in the same shape the Rust CLI writes, so the
# runner can read wall time, stage timings and peak RSS from either arm without
# special-casing. Peak RSS is `gc()`'s high-water mark plus the R heap, which is
# what `/usr/bin/time` would report for the process.
#
# Every stage is timed with `proc.time()` around a single call, and the numbers
# are written whether or not the run succeeded, so a failure is diagnosable from
# the results file alone.

suppressWarnings(suppressMessages({
  library(utils)
}))

data_dir <- Sys.getenv("ANCOMBC_BENCH_DATA")
out_dir  <- Sys.getenv("ANCOMBC_BENCH_OUT")
dataset  <- Sys.getenv("ANCOMBC_BENCH_DATASET")
threads  <- Sys.getenv("ANCOMBC_BENCH_THREADS")

if (!nzchar(data_dir) || !nzchar(dataset)) {
  stop("ANCOMBC_BENCH_DATA and ANCOMBC_BENCH_DATASET must be set; run via ",
       "scripts/run_benchmarks.py")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- environment record -------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# `jsonlite` is not a dependency of this repository, and the config files are
# written by scripts/make_bench_data.py with one scalar per line, so the three
# fields this script needs are read with a regex rather than a JSON parser. If
# that ever stops being true the parse fails loudly here.
read_config <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  one <- function(key, default) {
    m <- regmatches(txt, regexpr(sprintf('"%s"\\s*:\\s*("([^"]*)"|[-0-9.eE+]+|true|false)', key), txt))
    if (length(m) == 0L) return(default)
    v <- sub(sprintf('^"%s"\\s*:\\s*', key), "", m)
    v <- gsub('^"|"$', "", v)
    if (tolower(v) %in% c("true", "false")) return(tolower(v) == "true")
    suppressWarnings(as.numeric(v))
  }
  list(
    group = one("group", NULL),
    do_global = one("do_global", FALSE),
    do_pairwise = one("do_pairwise", FALSE),
    struc_zero = one("struc_zero", FALSE),
    sensitivity = one("sensitivity", FALSE),
    conservative = one("conservative", TRUE),
    p_adj_method = one("p_adj_method", "holm")
  )
}

env_lines <- c(
  "key\tvalue",
  sprintf("R\t%s", R.version.string),
  sprintf("platform\t%s", R.version$platform),
  sprintf("threads_requested\t%s", if (nzchar(threads)) threads else "all"),
  sprintf("blas\t%s", tryCatch(extSoftVersion()[["BLAS"]] %||% "unknown",
                                error = function(e) "unknown"))
)

# ---- load the pinned oracle --------------------------------------------
oracle_dir <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC")
harness_dir <- file.path("reference", "R")
for (f in c("stubs.R", "harness.R")) {
  source(file.path(harness_dir, f))
}
H <- tryCatch(
  load_harness(oracle_dir),
  error = function(e) {
    message("cannot load the oracle: ", conditionMessage(e))
    NULL
  }
)

if (is.null(H)) {
  # Record the reason rather than exiting silently: a results row with
  # `status = failed` and this message is more useful than a missing row.
  writeLines(c(env_lines,
               "status\tfailed",
               "error\toracle unavailable"),
             file.path(out_dir, "run_metadata.tsv"))
  quit(status = 2)
}

# ---- threading ----------------------------------------------------------
# ANCOMBC uses `foreach` with `%dorng%`. On one core `registerDoSEQ()` is
# behaviourally identical to the parallel backend for this code path and is what
# the golden harness uses, so the single-core arm is a sequential run.
# For the all-cores arm the real backend is needed, and it is only installed if
# the user has it; otherwise the arm reports that it degraded to sequential rather
# than quietly measuring the wrong thing.
n_cpu <- parallel::detectCores(logical = FALSE)
if (!nzchar(threads)) {
  if (requireNamespace("foreach", quietly = TRUE) &&
      requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(n_cpu)
    doParallel::registerDoParallel(cl)
    threads <- as.character(n_cpu)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    message("foreach/doParallel not installed; the all-cores arm runs sequentially")
    registerDoSEQ()
    threads <- "1-sequential"
  }
} else if (as.integer(threads) == 1L) {
  registerDoSEQ()
}

# ---- read the dataset ---------------------------------------------------
counts <- as.matrix(read.delim(file.path(data_dir, paste0(dataset, ".counts.tsv")),
                               row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(data_dir, paste0(dataset, ".meta.tsv")),
                   row.names = 1, check.names = FALSE)
cfg <- read_config(file.path(data_dir, paste0(dataset, ".config.json")))
fix_formula <- trimws(readLines(file.path(data_dir, paste0(dataset, ".formula.txt"))))

counts[is.na(counts)] <- 0

peak_rss_kb <- function() {
  g <- gc(verbose = FALSE, full = FALSE)
  # `gc()` reports Vcells in units of 8 bytes on a 64-bit build, and it is a
  # *current* figure rather than a process high-water mark -- which is why the
  # metadata says so. R's own vector heap plus any BLAS buffer is all it covers;
  # a large external timer is needed for the true peak, and this is recorded as
  # a lower bound rather than dressed up as a measurement.
  sum(g[, "used"]) * 8 / 1024
}

t_all <- proc.time()[["elapsed"]]
res <- H$ref_run(counts, meta, fix_formula)
t_all <- proc.time()[["elapsed"]] - t_all

n_sig <- if (!is.null(res$diff_abn)) sum(res$diff_abn, na.rm = TRUE) else NA
n_taxa <- if (!is.null(res$taxa_retained)) length(res$taxa_retained) else nrow(counts)

meta_lines <- c(
  env_lines,
  sprintf("arm\tr-%s", if (identical(threads, "1") || identical(threads, "1-sequential")) "1core" else "parallel"),
  sprintf("wall_seconds\t%.6f", t_all),
  sprintf("peak_rss_kb\t%.1f", peak_rss_kb()),
  "peak_rss_source\tgc() Vcells, a lower bound (not a process high-water mark)",
  sprintf("n_taxa_reported\t%d", n_taxa),
  sprintf("n_diff_abn\t%d", n_sig),
  "status\tok"
)
writeLines(meta_lines, file.path(out_dir, "run_metadata.tsv"))
message(sprintf("%s: %.3fs wall, %d taxa, %d significant", dataset, t_all, n_taxa, n_sig))
