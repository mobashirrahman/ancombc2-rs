# Run the bridge's acceptance checks as part of `R CMD check`.
#
# `r/ANCOMBC/tests/` also holds upstream's own `testthat.R`, which `R CMD check`
# runs. Both are plain scripts with no framework dependency, because the bridge
# checks are about a `.Call` boundary and a test framework's fixtures would obscure
# exactly what is being asserted.
library(ANCOMBC)

if (!ANCOMBC:::rb_available()) {
  stop("the typed native bridge is not present in this installation; ",
       "the transport checks cannot be skipped, only failed",
       call. = FALSE)
}

src <- system.file("tests", "bridge_selftest.R", package = "ANCOMBC")
if (!nzchar(src)) {
  # `tests/` is not installed, so run the file from the source tree when there is
  # one, and otherwise say so rather than passing silently.
  cand <- file.path(dirname(dirname(getwd())), "bridge_selftest.R")
  stop("bridge_selftest.R was not installed with the package and no source tree ",
       "is available at ", cand, "; the transport checks did not run",
       call. = FALSE)
}

status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(src)))
if (!identical(status, 0L)) {
  stop("the typed native transport selftest failed; see its output above",
       call. = FALSE)
}

# The preprocessing stages need no fixture -- the script builds its own and compares
# against expressions transcribed from the pinned R -- so it runs here whenever the
# bridge is present.
stage_src <- system.file("tests", "stage_selftest.R", package = "ANCOMBC")
if (!nzchar(stage_src)) {
  stage_src <- file.path("..", "scripts", "check_preprocess_stages.R")
}
if (file.exists(stage_src)) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(stage_src)))
  if (!identical(status, 0L)) {
    stop("the preprocessing stage comparison failed; see its output above",
         call. = FALSE)
  }
} else {
  cat("SKIP  the preprocessing stage checks: ", stage_src, " was not found\n", sep = "")
}

# The output transport needs a payload fixture captured from the pinned original,
# which a source check has no way to produce. When one is present it is run; when
# it is not, that is said out loud rather than left as a silent gap.
out_src <- system.file("tests", "output_selftest.R", package = "ANCOMBC")
payload <- Sys.getenv("ANCOMBC_PAYLOAD_RDS", "validation/exact/payloads/tiny-defaults.rds")
if (!nzchar(out_src)) {
  cat("SKIP  the output transport checks: tests/ is not installed\n")
} else if (!file.exists(payload)) {
  cat(sprintf("SKIP  the output transport checks: no payload fixture at %s\n", payload))
  cat("      They are run by `make r-output-selftest`, which needs the pinned\n",
      "      original. Set ANCOMBC_PAYLOAD_RDS to run them here.\n")
} else {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(out_src)))
  if (!identical(status, 0L)) {
    stop("the output transport selftest failed; see its output above", call. = FALSE)
  }
}
