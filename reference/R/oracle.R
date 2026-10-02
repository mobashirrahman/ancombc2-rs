# Load the pinned ANCOMBC oracle and expose an ancombc2-rs-compatible entry point.
#
# The oracle's R files are sourced rather than the package being installed:
#   * the package declares R >= 4.5.0, and the interpreter available here is 4.3.3
#   * the golden contract needs the internal (dot-prefixed) functions
#     .iter_mle, .lm_fit_all, .sandwich_vcov and .bias_em, none of which are
#     exported
#
# See reference/env/ORACLE.md.

ORACLE_SHA <- "dc4febdf59badb3a8dfe0c767ef2186323c2199a"
ORACLE_VERSION <- "2.15.2"
ORACLE_DIR <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC")

verify_oracle <- function(dir = ORACLE_DIR) {
  have_git <- dir.exists(file.path(dir, ".git"))
  if (have_git) {
    sha <- system2("git", c("-C", shQuote(dir), "rev-parse", "HEAD"),
                   stdout = TRUE, stderr = FALSE)
    if (length(sha) != 1 || !identical(trimws(sha), ORACLE_SHA)) {
      stop(sprintf("oracle commit mismatch: expected %s, found %s",
                   ORACLE_SHA, paste(sha, collapse = " ")))
    }
  }
  desc <- read.dcf(file.path(dir, "DESCRIPTION"))
  ver <- trimws(unname(desc[1L, "Version"]))
  if (!identical(ver, ORACLE_VERSION)) {
    stop(sprintf("oracle version mismatch: expected %s, found %s",
                 ORACLE_VERSION, ver))
  }
  invisible(TRUE)
}

# Sequential stubs for the foreach/doRNG machinery. With n_cl = 1 the real
# package calls foreach::registerDoSEQ(), for which %dopar% is a plain
# lapply/iteration, so a sequential stub is behaviourally identical for the
# fixed-effects path.
stub_foreach <- function() {
  if (!exists(".stubbed_foreach", envir = globalenv())) {
    assign("foreach", structure(list(), class = "stub_foreach"), envir = globalenv())
    assign("%dopar%", function(obj, ex) obj, envir = globalenv())
    assign("%dorng%", function(obj, ex) obj, envir = globalenv())
    assign("registerDoSEQ", function(...) invisible(NULL), envir = globalenv())
    assign(".stubbed_foreach", TRUE, envir = globalenv())
  }
  invisible(TRUE)
}

# The oracle uses `%dorng%` as an infix operator applied to the *result* of a
# call, e.g. `foreach(i = 1:3) %dorng% { ... }`. Sourcing the file with a
# stubbed operator in globalenv() is enough, but the operator is looked up in
# the calling frame of the loop, so we install it globally.
stub_foreach()

# The files the fixed-effects path needs. `sim_data.R` is deliberately *not*
# here: it holds the data generators, which the golden contract does not use, and
# sourcing them would put a second RNG consumer into the harness's environment.
# `load_sim_helpers()` adds them for the callers that want them -- the vignette
# reconstruction in scripts/prepare_realdata.R.
ORACLE_R_FILES <- c("utils.R", "ancombc_prep.R", "ancombc_bias_correct.R",
                    "ancombc_mult.R")
ORACLE_SIM_FILES <- c("sim_data.R")

load_oracle <- function(dir = ORACLE_DIR, files = ORACLE_R_FILES) {
  verify_oracle(dir)
  env <- new.env(parent = globalenv())
  for (f in files) {
    sys.source(file.path(dir, "R", f), envir = env)
  }
  env
}

#' The oracle's data generators, for reconstructing the vignette's dataset.
#'
#' `sim_plnm` is what `ANCOMBC2.Rmd` uses to build its worked example, and it is
#' the only way to reproduce that dataset faithfully -- regenerating it any other
#' way would be a different dataset with the same name.
load_sim_helpers <- function(dir = ORACLE_DIR) {
  load_oracle(dir, files = c(ORACLE_R_FILES, ORACLE_SIM_FILES))
}
