#!/usr/bin/env Rscript
# Guard `make r-test` against testing a stale build.
#
# `.Rprofile` puts the oracle's private library first on `.libPaths()`, because
# that is where `nloptr` lives. An `ancombc2rs` installed *there* therefore shadows
# the one this repository builds, and `r-test` will quietly exercise the stale copy
# -- which it did: the wrapper suite reported green for a month against a copy
# installed a build earlier, so new wrapper code was never actually under test.
#
# Two things are checked:
#
#   1. no stray copy of this package in the oracle library, removed if present;
#   2. the library `find.package()` resolves to is the one the build target writes.
#
# Both are cheap. Neither is optional if the number of green tests is to mean
# anything.

args <- commandArgs(trailingOnly = TRUE)
built <- args[[1]]
oracle_lib <- Sys.getenv("ANCOMBC_RLIB", Sys.getenv("R_LIBS_USER", ".rlib/deps"))

stray <- file.path(oracle_lib, "ancombc2rs")
if (dir.exists(stray)) {
  message("r-test: removing stale ancombc2rs from the oracle library: ", stray)
  unlink(stray, recursive = TRUE)
}

if (!nzchar(built)) stop("usage: check_r_lib.R <library-path-ancombc2rs-was-built-in>")
if (!requireNamespace("ancombc2rs", quietly = TRUE)) {
  stop("r-test: ancombc2rs is not installed; run `make r-install` first")
}

want <- normalizePath(file.path(built, "ancombc2rs"), mustWork = TRUE)
got <- normalizePath(find.package("ancombc2rs"), mustWork = TRUE)

if (!identical(want, got)) {
  stop(
    "r-test: loaded ancombc2rs from\n  ", got,
    "\nbut this target built\n  ", want,
    "\nThe tests would not be testing the build.",
    call. = FALSE
  )
}
message("r-test: testing ancombc2rs from ", got)
invisible(TRUE)