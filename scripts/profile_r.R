#!/usr/bin/env Rscript
# The week-1 profiling gate: where does the reference's time actually go?
#
# The point is to decide *which* optimisation to make next, and to decide it on
# evidence rather than on intuition. The A/B/C decision rule is in
# docs/compatibility.md; this script produces the input to it.
#
# Usage:
#   Rscript scripts/profile_r.R <counts.tsv> [meta.tsv] [outdir]
#   Rscript scripts/profile_r.R validation/fixtures/fx03/counts.tsv
#
# Rprof writes a file of interval samples. The summary below is deliberately
# coarse -- one line per function, with the cumulative percentage -- because the
# decision rule asks which *stage* dominates, not which line of R interpreter
# overhead does. A finer breakdown would mostly show `[<-` and S3 dispatch, which
# is Outcome C.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("usage: Rscript scripts/profile_r.R <counts.tsv> [meta.tsv] [outdir]")
}
counts_path <- args[[1]]
meta_path <- if (length(args) >= 2L) args[[2]] else NULL
out_dir <- if (length(args) >= 3L) args[[3]] else "benchmarks/results/profile"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c("stubs.R", "harness.R")) {
  source(file.path("reference", "R", f))
}
H <- load_harness(Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC"))

counts <- as.matrix(read.delim(counts_path, row.names = 1, check.names = FALSE))
counts[is.na(counts)] <- 0

if (is.null(meta_path)) {
  # Derive a two-group metadata table from the column count, so a counts file
  # alone is enough to profile something.
  n <- ncol(counts)
  meta <- data.frame(group = factor(rep(seq_len(max(2L, n %/% 10L)), length.out = n)),
                     x1 = rnorm(n))
  message("no metadata given; using a synthetic two-group table with ", n, " samples")
} else {
  meta <- read.delim(meta_path, row.names = 1, check.names = FALSE)
  if ("group" %in% colnames(meta)) meta$group <- factor(meta$group)
}

# The formula is inferred from the metadata, since a counts file carries none.
terms <- setdiff(colnames(meta), "group")
fix_formula <- if (length(terms)) paste(c("group", terms), collapse = " + ") else "group"
message("profiling: ", nrow(counts), " x ", ncol(counts), ", formula ", fix_formula)

prof_file <- file.path(out_dir, "Rprof.out")
# `interval = 0.002` is fine enough to separate a stage that takes 5% of the run
# from one that takes 4%, and coarse enough that the output file stays readable.
Rprof(prof_file, interval = 0.002, line.profiling = FALSE, memory.profiling = FALSE)

t0 <- proc.time()[["elapsed"]]
res <- H$ref_run(counts, meta, fix_formula)
elapsed <- proc.time()[["elapsed"]] - t0

Rprof(NULL)
message(sprintf("profiled run: %.2fs", elapsed))

s <- summaryRprof(prof_file)
# Rprof only samples while R bytecode is executing, so a path dominated by
# compiled code (`lm.fit`, `crossprod`) yields few samples. An empty or tiny
# table is a *result*, not a failure: it means the attribution is too coarse to
# decide from, which is Outcome C. Guard it rather than crashing on it.
# `summaryRprof` puts the function names in the *row names*, not in a column, so
# they are lifted into an `fn` column here. Doing it once keeps the rest of the
# script from having to know that.
lift <- function(tab) {
  if (is.null(tab) || nrow(tab) == 0L) {
    return(data.frame(fn = character(0), total.pct = numeric(0),
                      self.pct = numeric(0), total.time = numeric(0),
                      self.time = numeric(0), stringsAsFactors = FALSE))
  }
  out <- data.frame(fn = rownames(tab), tab, row.names = NULL,
                    check.names = FALSE, stringsAsFactors = FALSE)
  out[order(-out$total.pct), , drop = FALSE]
}
n_samples <- if (is.null(s$by.self)) 0L else nrow(s$by.self)
by_total <- lift(s$by.total)
by_self <- lift(s$by.self)

top_n <- 25L
report <- c(
  sprintf("# Rprof of %s", counts_path),
  sprintf("# formula: %s", fix_formula),
  sprintf("# shape: %d taxa x %d samples", nrow(counts), ncol(counts)),
  sprintf("# wall: %.3f s", elapsed),
  sprintf("# samples: %d", n_samples),
  "",
  "## by total time (inclusive)",
  sprintf("%-42s %8s %8s %10s", "function", "total%", "self%", "seconds"),
  if (nrow(by_total)) {
    sprintf("%-42s %8.2f %8.2f %10.3f",
            head(by_total$fn, top_n), head(by_total$total.pct, top_n),
            head(by_total$self.pct, top_n), head(by_total$total.time, top_n))
  } else {
    "(no samples: the run was dominated by compiled code)"
  },
  "",
  "## by self time (exclusive)",
  sprintf("%-42s %8s %8s %10s", "function", "total%", "self%", "seconds"),
  if (nrow(by_self)) {
    sprintf("%-42s %8.2f %8.2f %10.3f",
            head(by_self$fn, top_n), head(by_self$total.pct, top_n),
            head(by_self$self.pct, top_n), head(by_self$self.time, top_n))
  } else {
    "(no samples: the run was dominated by compiled code)"
  }
)
writeLines(report, file.path(out_dir, "Rprof.summary.txt"))
cat(report, sep = "\n")

# --- the decision rule, applied -------------------------------------------
#
# Stage attribution needs a mapping from R functions to the algorithm's stages.
# The harness's `capture` list already separates the stages, so this uses it
# rather than guessing from function names.
capture_names <- function(capture) {
  if (is.null(capture)) character(0) else names(capture)
}
capture <- attr(res, "capture")
stage_of <- c(
  "log_transform_and_centre" = NA_real_,
  "iter_mle" = NA_real_,
  "bias_em" = NA_real_,
  "sampling_fractions" = NA_real_,
  "second_mle" = NA_real_,
  "sandwich" = NA_real_,
  "tests" = NA_real_,
  "sensitivity" = NA_real_
)
if (!is.null(capture)) {
  for (nm in names(stage_of)) {
    v <- capture[[nm]]
    if (is.null(v)) next
    t <- attr(v, "seconds")
    if (!is.null(t)) stage_of[[nm]] <- t
  }
}

# Fall back to the profile's own inclusive totals when the harness did not time
# the stages. Rprof qualifies the names (`"oracle$.bias_em"`, `"ref_core"`), so
# the match is a substring rather than a prefix.
if (all(is.na(stage_of)) && nrow(by_total) > 0L) {
  anchor <- c(
    iter_mle = ".iter_mle",
    bias_em = ".bias_em",
    sandwich = ".sandwich_vcov",
    tests = ".ancombc_global_F",
    sensitivity = ".ancombc2_sens_fit"
  )
  for (nm in names(anchor)) {
    # `vapply(X, f, ...)` passes `...` to `f` after its first argument, so the
    # function name is the *pattern* and the anchor is the string. Getting that
    # the wrong way round silently matches nothing.
    hit <- which(vapply(
      by_total$fn,
      function(one) grepl(anchor[[nm]], one, fixed = TRUE),
      logical(1)
    ))
    if (length(hit)) stage_of[[nm]] <- by_total$total.time[hit[1]]
  }
  # The two MLE stages are both inside `.iter_mle`, so its inclusive total covers
  # both; the sandwich is inside it too, so the stages overlap. The decision rule
  # only needs the *dominant* one, and an overlapping attribution is still sound
  # for that: the largest inclusive time is the largest cost.
  stage_of[["iter_mle"]] <- stage_of[["iter_mle"]]
}

# With overlapping attribution the stages do not sum to the run, so the share is
# taken against the run's own total rather than against their sum.
run_total <- max(by_total$total.time, 0)
if (is.null(run_total) || !is.finite(run_total) || run_total <= 0) {
  run_total <- elapsed
}
known <- sum(stage_of, na.rm = TRUE)
if (is.finite(known) && known > 0) {
  dominant <- names(which.max(ifelse(is.na(stage_of), -1, stage_of)))
  share <- max(stage_of, na.rm = TRUE) / run_total
  cat("\n## decision\n")
  cat(sprintf("run total %.2f s; %d stages attributed; %s is %.0f%% of the run\n",
              run_total, sum(!is.na(stage_of)), dominant, 100 * share))
  if (share >= 0.50) {
    cat(sprintf(paste("Outcome A: one stage (%.0f%%) dominates. Optimise %s, then",
                      "require its share below 35%% before declaring it done.\n"),
                100 * share, dominant))
  } else if (known > 0 && share < 0.35 && all(stage_of[!is.na(stage_of)] > 0)) {
    cat(paste("Outcome B: no stage dominates; the cost is structural, in the 20",
              "MLE iterations each rebuilding a full lm_fit_all. Fuse the",
              "iteration with the pattern solve before micro-optimising.\n"))
  } else {
    cat(paste("Outcome C: the attributed time is spread too thinly to call.",
              "Check the R-level overhead breakdown above; if it is [<- and S3",
              "dispatch, the algorithm is memory-bound and the next gain is in",
              "the layout.\n"))
  }
} else {
  cat("\n## decision\n")
  cat(paste("Outcome C: the stages could not be separated from this profile.",
            "The harness did not time them and the profile has no stage",
            "anchors. Re-run with a fixture large enough for the stages to",
            "differ measurably.\n"))
}
