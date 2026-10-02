#!/usr/bin/env Rscript
# Smoke test: prove the pinned oracle actually executes in this environment and
# that the instrumented mirror reproduces the oracle's own .ancombc2_core bit for
# bit. Nothing downstream is trustworthy until this passes.

suppressPackageStartupMessages(library(stats))
here <- Sys.getenv("ANCOMBC_REPO", unset = ".")
oracle_dir <- file.path(here, "reference", "ANCOMBC")
source(file.path(here, "reference", "R", "harness.R"))
source(file.path(here, "reference", "R", "fixtures.R"))

cat("R:", R.version.string, "\n")
h <- load_harness(oracle_dir)
cat("oracle loaded:", is.function(h$env$.iter_mle), "\n")

# The formula normaliser, checked on names that break it.
#
# `\s` inside a POSIX bracket expression is *literal* in R's default engine, so
# the class `[-+*\s]` is the five characters `-`, `+`, `*`, `\` and `s` -- and
# every `s` in a variable name used to be deleted. `sex` became `ex` and the
# analysis failed with "variables not in metadata", which reads like a metadata
# problem. It was invisible for every name the fixtures and the benchmark
# datasets use, and the real-data layer is what finally reached it. This check
# exists so it cannot come back.
cat("\n-- formula normalisation --")
formula_cases <- list(
  list("sex", "sex"),
  list("group", "group"),
  list("~sex", "sex"),
  list("group + cont_cov", "group + cont_cov"),
  list("samp_frac + sex", "samp_frac + sex"),
  list("group x1 x2", "group + x1 + x2"),
  list("nationality + timepoint + group", "nationality + timepoint + group")
)
formula_bad <- 0L
for (case in formula_cases) {
  got <- normalise_fix_formula(case[[1]])
  ok <- identical(got, case[[2]])
  if (!ok) formula_bad <- formula_bad + 1L
  cat(sprintf("\n  %-34s -> %-34s %s", case[[1]], got, if (ok) "ok" else "MISMATCH"))
}
if (!is.null(normalise_fix_formula(NULL))) {
  formula_bad <- formula_bad + 1L
  cat("\n  NULL did not stay NULL")
}
cat("\n")
if (formula_bad > 0L) {
  stop(formula_bad, " formula normalisation case(s) failed")
}
cat("all formula cases ok\n")

f <- gen_fixture(1)
cat("fixture 1:", nrow(f$counts), "taxa x", ncol(f$counts), "samples\n")
cat("formula:", fix_formula_for(f$spec), "\n")
print(table(f$meta$group))

ff <- fix_formula_for(f$spec)
cat("\n-- mirror self-check --\n")
t0 <- proc.time()[["elapsed"]]
ok <- h$ref_verify_mirror(f$counts, f$meta, fix_formula = ff,
                          p_adj_method = "holm", alpha = 0.05, prv_cut = 0.0)
cat("mirror identical to .ancombc2_core:", ok,
    sprintf("(%.2fs)", proc.time()[["elapsed"]] - t0), "\n")

cat("\n-- full ref_run --\n")
t0 <- proc.time()[["elapsed"]]
g <- h$ref_run(f$counts, f$meta, fix_formula = ff, p_adj_method = "holm",
               alpha = 0.05, prv_cut = 0.0)
cat(sprintf("ref_run done in %.2fs\n", proc.time()[["elapsed"]] - t0))
cat("taxa retained:", length(g$taxa_retained), "\n")
cat("fix_eff:", paste(g$fix_eff, collapse = ", "), "\n")
cat("delta_em:", sprintf("%.6f", g$delta_em), "\n")
cat("delta_wls:", sprintf("%.6f", g$delta_wls), "\n")
cat("s0:", sprintf("%.6f", g$s02), "\n")
print(head(g$res, 4))
cat("n diff_abn:", sum(g$diff_abn, na.rm = TRUE), "/", nrow(g$res), "\n")

cat("\n-- global + pairwise (3 groups) --\n")
f3 <- gen_fixture(3)
g3 <- h$ref_run(f3$counts, f3$meta, fix_formula = fix_formula_for(f3$spec),
                p_adj_method = "holm", alpha = 0.05, prv_cut = 0.0,
                global = TRUE, pairwise = TRUE, group = "group")
print(head(g3$res_global, 3))
cat("pairwise cols:", paste(grep("^lfc_", names(g3$res_pair), value = TRUE),
                            collapse = ", "), "\n")
print(head(g3$res_pair[, 1:4], 3))
cat("pairwise diff count:", sum(g3$res_pair$diff_global_group, na.rm = TRUE), "\n")

cat("\nOK\n")
