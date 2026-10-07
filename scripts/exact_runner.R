#!/usr/bin/env Rscript
# One arm of an exact original-vs-candidate comparison (IMPROVED_PLAN.md S02).
#
#   Rscript --vanilla scripts/exact_runner.R \
#     --library <dir>          the ONLY library this arm may resolve ANCOMBC from
#     --arm original|candidate
#     --input <input.rds>      the argument object + initial state, saved once
#     --out <dir>              where the capture is written
#     [--forbid-path <dir>]    a path the resolved package must NOT live under
#
# What this runner deliberately does not do:
#
#   * it does not source reference/R/harness.R, reference/R/oracle.R or
#     reference/R/stubs.R;
#   * it does not install a foreach / %dorng% / registerDoSEQ stub, so the
#     installed package's real parallel backend is what runs;
#   * it does not read a golden, an oracle result, or anything the candidate
#     produced;
#   * it does not modify the result. The result is serialized once, by the
#     contract function, and that is the only copy anything downstream sees.
#
# Both arms load the *same* input file. That is the whole point: the input is a
# saved R object, not a text file both sides might parse differently.
#
# The arm's package identity is written out as canonical bytes
# (`pkg_identity.txt`); scripts/run_exact.py hashes it with the same function
# scripts/verify_profile.py uses, so there is one definition of "same package".
# This script stays on base R so a missing package cannot stop the capture.

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
  key <- sub("^--", "", a)
  if (i == length(args)) stop("option --", key, " needs a value", call. = FALSE)
  opt[[key]] <- args[[i + 1L]]
  i <- i + 2L
}

`%||%` <- function(a, b) if (is.null(a)) b else a

need <- c("library", "arm", "input", "out")
absent <- need[vapply(need, function(k) is.null(opt[[k]]), logical(1))]
if (length(absent)) {
  stop("missing required option(s): ", paste0("--", absent, collapse = " "),
       call. = FALSE)
}
arm <- opt[["arm"]]
if (!arm %in% c("original", "candidate")) {
  stop("--arm must be 'original' or 'candidate', not ", sQuote(arm), call. = FALSE)
}
lib <- normalizePath(opt[["library"]], mustWork = FALSE)
if (!dir.exists(lib)) stop("library does not exist: ", lib, call. = FALSE)
if (!file.exists(opt[["input"]])) {
  stop("input does not exist: ", opt[["input"]], call. = FALSE)
}
out_dir <- opt[["out"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Resolve the package, and refuse to run if this arm could be the other arm.
# ---------------------------------------------------------------------------

pkg_dir <- find.package("ANCOMBC", quiet = TRUE)
if (!length(pkg_dir)) {
  stop("ANCOMBC is not installed in ", lib,
       "; refusing to run. An arm with no package is not a result.", call. = FALSE)
}
pkg_dir <- normalizePath(pkg_dir[[1]])
if (!startsWith(pkg_dir, lib)) {
  # find.package() searched other libraries too. That is exactly how an arm ends
  # up measuring the other arm, so it is a hard failure, not a warning.
  stop("ANCOMBC resolved to ", pkg_dir, ", which is outside the arm's own ",
       "library ", lib, "; refusing to run", call. = FALSE)
}
if (!is.null(opt[["forbid-path"]])) {
  fp <- normalizePath(opt[["forbid-path"]], mustWork = FALSE)
  if (nzchar(fp) && startsWith(pkg_dir, fp)) {
    stop("ANCOMBC resolved under the forbidden path ", fp, call. = FALSE)
  }
}

# The package's own identity is NOT digested here: R has no base sha256, and a
# second definition of "same package" in two languages is exactly how two
# identities drift apart. The resolved path is written to meta.tsv and
# scripts/run_exact.py digests that directory with the same function
# scripts/verify_profile.py uses.

descf <- read.dcf(file.path(pkg_dir, "DESCRIPTION"))
pkg_version <- trimws(unname(descf[1, "Version"]))

# A candidate must announce itself. Until S05 there is no replacement package,
# and silently measuring the original under the candidate's name is the failure
# mode this guard exists to prevent.
# `R CMD INSTALL` copies an inst/ directory's contents to the installed package
# root, so the marker is looked for there and not at the source-tree path. A
# guard that looked in the wrong place would let an unmarked package through.
marker <- file.path(pkg_dir, "REPLACEMENT_PROVENANCE.md")
is_replacement <- file.exists(marker)
if (arm == "candidate" && !is_replacement) {
  stop("arm=candidate but the package at ", pkg_dir,
       " carries no REPLACEMENT_PROVENANCE.md marker; refusing to run the ",
       "original as the candidate", call. = FALSE)
}
if (arm == "original" && is_replacement) {
  stop("arm=original but the package at ", pkg_dir, " is marked as a replacement; ",
       "refusing to run", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Initial state. Both arms get it from the same saved object, so a difference in
# RNG state cannot be mistaken for a difference in the implementation.
# ---------------------------------------------------------------------------

input <- readRDS(opt[["input"]])
if (!is.list(input) || !identical(input$schema, "ancombc2-exact-input/1")) {
  stop("input is not an ancombc2-exact-input/1 object", call. = FALSE)
}
for (k in c("case_id", "fn", "args")) {
  if (is.null(input[[k]])) stop("input is missing ", sQuote(k), call. = FALSE)
}
state <- input$state %||% list()
if (!is.null(state$rng_kind)) {
  suppressWarnings(do.call(RNGkind, as.list(state$rng_kind)))
}
if (!is.null(state$sample_kind)) {
  suppressWarnings(RNGkind(sample.kind = state$sample_kind))
}
if (!is.null(state$seed)) set.seed(state$seed)
if (!is.null(state$options)) do.call(options, state$options)

# Observable state that a call is allowed to perturb. Recorded before and after,
# so a mutation is a diff rather than a mystery.
OBSERVED_OPTIONS <- c("digits", "OutDec", "scipen", "warn", "stringsAsFactors",
                      "useFancyQuotes", "contrasts", "showErrorCalls",
                      "nwarnings", "keep.source", "max.print", "expressions",
                      "width")
observe <- function() {
  list(
    options = options()[OBSERVED_OPTIONS],
    globalenv_names = sort(ls(envir = globalenv(), all.names = TRUE)),
    search_path = .libPaths(),
    rng_kind = RNGkind(),
    dot_random_seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    },
    loaded_namespaces = sort(loadedNamespaces())
  )
}
before <- observe()

# ---------------------------------------------------------------------------
# Conditions. Recorded by a *calling* handler that does not muffle, so default
# warning/message/error behaviour is exactly what the installed package would get
# without the runner. Their order is part of the comparison.
# ---------------------------------------------------------------------------

cond_log <- list()
record <- function(c) {
  cond_log[[length(cond_log) + 1L]] <<- list(
    class = class(c),
    message = conditionMessage(c),
    call = if (is.null(conditionCall(c))) NA_character_
           else paste(deparse(conditionCall(c)), collapse = " "),
    type = if (inherits(c, "error")) "error"
           else if (inherits(c, "warning")) "warning"
           else if (inherits(c, "message")) "message"
           else "condition"
  )
  invisible(NULL)
}

fn <- get(input$fn, envir = asNamespace("ANCOMBC"))
if (!exists(input$fn, envir = asNamespace("ANCOMBC"), inherits = FALSE)) {
  stop("input$fn is not an internal of the installed package: ", input$fn,
       call. = FALSE)
}
call_args <- input$args
# `call_style = "positional"` means the arguments were recorded positionally, so
# the names must be dropped before do.call() or the two are combined into a
# different, shifted call. "named" keeps them.
if (identical(input$call_style %||% "named", "positional")) names(call_args) <- NULL

outcome <- "success"
err_obj <- NULL
t0 <- proc.time()[["elapsed"]]
result <- withCallingHandlers(
  tryCatch(do.call(fn, call_args),
           error = function(e) {
             outcome <<- "error"
             err_obj <<- e
             NULL
           }),
  warning = record,
  message = record
)
t1 <- proc.time()[["elapsed"]]

after <- observe()

# ---------------------------------------------------------------------------
# Serialisation. The release contract, verbatim: one function, both arms, no
# normalisation on the way in or out.
# ---------------------------------------------------------------------------

result_bytes <- function(x) serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3)

if (identical(outcome, "success")) {
  writeBin(result_bytes(result), file.path(out_dir, "result.bytes"))
} else {
  # An error is a comparable outcome, not a missing result: its condition goes
  # through the same contract so its bytes can be compared too.
  writeBin(result_bytes(list(
    message = conditionMessage(err_obj),
    class = class(err_obj),
    call = if (is.null(conditionCall(err_obj))) NULL else deparse(conditionCall(err_obj))
  )), file.path(out_dir, "error.bytes"))
}

saveRDS(list(conditions = cond_log, before = before, after = after,
             outcome = outcome),
        file.path(out_dir, "state.rds"), version = 3)

writeLines(
  c(sprintf("schema=ancombc2-exact-run/1"),
    sprintf("arm=%s", arm),
    sprintf("case_id=%s", input$case_id),
    sprintf("fn=%s", input$fn),
    sprintf("call_style=%s", input$call_style %||% "named"),
    sprintf("outcome=%s", outcome),
    sprintf("pkg_path=%s", pkg_dir),
    sprintf("pkg_version=%s", pkg_version),
    sprintf("pkg_is_replacement=%s", is_replacement),
    sprintf("n_exports=%d", length(getNamespaceExports("ANCOMBC"))),
    sprintf("export_names=%s", paste(sort(getNamespaceExports("ANCOMBC")),
                                     collapse = ",")),
    sprintf("n_conditions=%d", length(cond_log)),
    sprintf("r_version=%s", R.version.string),
    sprintf("r_platform=%s", R.version$platform),
    sprintf("blas=%s", tryCatch(extSoftVersion()[["BLAS"]],
                                error = function(e) "<absent>"))),
  con = file.path(out_dir, "meta.tsv")
)

# Exit 3 means "the analysis raised an error", which S03 compares as an outcome.
# Any other non-zero status is a runner failure.
quit(save = "no", status = if (identical(outcome, "error")) 3L else 0L)
