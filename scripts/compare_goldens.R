#!/usr/bin/env Rscript
# Compare two golden directories and report how far apart they are.
#
# Why this is not `diff -r`
# --------------------------
# Three things in a golden `.rds` are environment rather than contract:
#
#   * `session` records `sessionInfo()` and the oracle's path, so it differs
#     between a developer machine and CI;
#   * `spec` records the *fixture generator's* parameters, and it grew five fields
#     (`predictor`, `interact`, `balanced`, `pseudo`, `p_adj_method`) after the
#     four committed fixtures were generated, so a recomputation produces an `rds`
#     with a longer `spec` and the same data;
#   * a `capture` *environment*, attached by `reference/R/harness.R` as a
#     side-channel for internals the oracle does not return. `identical()` compares
#     environments by reference, so anything carrying one can never be `identical()`
#     to the same object built in another R process -- regardless of the data.
#
# And, for `fx01`..`fx04`, the numbers themselves cannot be bit-reproduced:
#
#   * the fixtures predate the `set.seed(spec$seed)` in `reference/R/fixtures.R`
#     and were produced from the ambient session state, so they cannot be
#     regenerated at all;
#   * even taking the committed fixture *files* as given, `write.table` renders a
#     double at 15 significant digits, so reading them back cannot reproduce the
#     in-memory doubles the goldens were computed from.
#
# That the oracle itself is deterministic was checked rather than assumed: two
# runs of `generate_goldens.R --from-committed 4` over the same text agree to
# **0** relative deviation. So every deviation below is the text round-trip, not
# oracle non-determinism.
#
# Measured floors, as a maximum relative deviation over the whole contract:
#
#   fx01  0          10 x 10, integers, representable exactly
#   fx02  4.13e-12   100 x 30
#   fx03  1.83e-11   1,000 x 100
#   fx04  9.55e-04   9,800 x 500, non-conservative sensitivity
#
# The `fx04` figure is the interesting one. `write.table` loses ~1e-15 relative
# per count, and that is amplified through the MLE over 9,800 taxa into `samp_frac`
# at 6.3e-08, then into `beta` at 2.0e-05 and `y_bias_crt` at 9.6e-04. It is the
# largest fixture by a factor of ten and the only one running the 50-refit
# non-conservative path, and it is a property of the committed fixture *text*, not
# of the algorithm. A single 1e-3 tolerance for all four would be far too loose for
# `fx01`..`fx03`, so the floors are per fixture.
#
# These are floors, not tolerances loosened to pass: each is within an order of
# magnitude of the measured value, and `--tol` overrides all of them. A deviation
# above a floor is drift and fails.
#
# Usage: compare_goldens.R DIR_A DIR_B [ids...] [--tol RELATIVE]

a_dir <- commandArgs(trailingOnly = TRUE)
# Per-fixture floors, as measured above. Anything not listed uses the smallest,
# which is the right default for a fixture that turns out to be exactly
# representable.
FLOORS <- c(fx01 = 1e-12, fx02 = 1e-10, fx03 = 1e-9, fx04 = 1e-3)
tol_override <- NA_real_
tol_i <- match("--tol", a_dir)
if (!is.na(tol_i)) {
  tol_override <- as.numeric(a_dir[tol_i + 1L])
  a_dir <- a_dir[-c(tol_i, tol_i + 1L)]
}
if (length(a_dir) < 2L) stop("usage: compare_goldens.R DIR_A DIR_B [ids...] [--tol RELATIVE]")
A <- a_dir[[1]]; B <- a_dir[[2]]
ids <- if (length(a_dir) > 2L) a_dir[-(1:2)] else
  sort(sub("^fx", "", list.files(A, pattern = "^fx[0-9]+$")))

strip_capture <- function(x) {
  if (is.environment(x)) return(NULL)
  attr(x, "capture") <- NULL
  if (is.list(x)) x[] <- lapply(x, strip_capture)
  x
}

# The compared quantities live in the golden's `.f64` files, not in
# `golden.rds`.
#
# `golden.rds` used to hold a second copy of every one of them. For `fx04` and
# the `shape-10000x500` cell those copies were 148 MB and 133 MB, and being the
# only two blobs in the repository over GitHub's 100 MB hard limit they made it
# unpushable. They are now written once, here, and `golden.rds` keeps only what
# has no file of its own -- see `slim_golden` in `reference/R/serialize.R`.
#
# This is a relocation, not a change of contract. `max_rel` below only ever
# considered numeric top-level fields, and every one of those is `.f64`-backed
# except `dof`, which the slim `golden.rds` still carries. The set of compared
# quantities, and the arithmetic, are identical to before -- the reported
# deviation is unchanged, which is the invariant worth preserving here.
#
# `manifest.rds` supplies the shapes. It is this repository's substitute for a
# JSON reader: `reference/R/serialize.R` writes JSON but has no reader on
# purpose, so the generator records the same manifest as an R object.
read_f64_payloads <- function(dir, manifest) {
  out <- list()
  for (nm in names(manifest)) {
    e <- manifest[[nm]]
    if (!is.list(e) || is.null(e$file) || !identical(e$kind, "f64")) next
    nc <- if (is.null(e$nc) || is.na(e$nc)) NA_integer_ else as.integer(e$nc)
    nr <- as.integer(e$nr)
    v <- readBin(file.path(dir, e$file), what = "double",
                 n = nr * (if (is.na(nc)) 1L else nc), size = 8, endian = "little")
    out[[nm]] <- if (is.na(nc)) v else matrix(v, nrow = nr, ncol = nc)
  }
  out
}

# The `.f64` payloads plus the slim residual, which is how the whole numeric
# contract is reassembled.
golden_quantities <- function(dir) {
  rds <- readRDS(file.path(dir, "golden.rds"))
  mf_path <- file.path(dir, "manifest.rds")
  if (!file.exists(mf_path)) {
    stop(dir, ": no manifest.rds, so the .f64 shapes are unknown. ",
         "Regenerate this golden; a manifest.rds holding the summary data.frame ",
         "means it predates the slim golden.rds and must be rebuilt.")
  }
  c(read_f64_payloads(dir, readRDS(mf_path)), strip_capture(rds$golden))
}

# Largest relative deviation between two numerics, per top-level field.
max_rel <- function(a, b) {
  worst <- 0
  for (k in union(names(a), names(b))) {
    x <- a[[k]]; y <- b[[k]]
    if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) next
    d <- max(abs(x - y) / pmax(abs(x), 1e-300), na.rm = TRUE)
    if (is.finite(d) && d > worst) worst <- d
  }
  worst
}

fail <- FALSE
for (id in ids) {
  id <- if (startsWith(id, "fx")) id else paste0("fx", id)
  pa <- file.path(A, id, "golden.rds")
  pb <- file.path(B, id, "golden.rds")
  if (!file.exists(pa) || !file.exists(pb)) {
    cat(sprintf("%s: MISSING\n", id))
    fail <- TRUE
    next
  }
  a <- readRDS(pa); b <- readRDS(pb)
  ga <- golden_quantities(dirname(pa)); gb <- golden_quantities(dirname(pb))
  d <- max_rel(ga, gb)
  c_ok <- identical(a$config, b$config)
  # A field present in one spec and not the other is a generator change, which the
  # header explains; report it so it is never silently tolerated.
  spec_extra <- setdiff(names(b$spec), names(a$spec))
  tol <- if (is.na(tol_override)) FLOORS[[id]] else tol_override
  if (is.null(tol)) tol <- min(FLOORS)
  ok <- is.finite(d) && d <= tol && c_ok
  cat(sprintf("%s: max relative deviation %.2e (floor %.0e) config %s%s -> %s\n",
              id, d, tol, if (c_ok) "same" else "DIFFERS",
              if (length(spec_extra))
                sprintf(" (spec gained: %s)", paste(spec_extra, collapse = ", "))
              else "",
              if (ok) "within the text round-trip floor" else "DRIFT"))
  if (!ok) fail <- TRUE
}
if (fail) quit(status = 1L)