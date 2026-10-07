#!/usr/bin/env Rscript
# Build every exact input RDS named by a fixture manifest (IMPROVED_PLAN.md S04).
#
#   Rscript --vanilla scripts/make_exact_inputs.R \
#     --spec validation/exact/fixtures.json \
#     --out-dir validation/exact/inputs \
#     [--index validation/exact/inputs/built.json] \
#     [--only id1,id2] [--profile linux-r453-openblas]
#
# Both arms read these RDS files and nothing else. The legacy goldens under
# validation/golden are NOT touched: this writes only under --out-dir.
#
# This script builds and records *what* it built. It does not compute digests:
# R has no base sha256, and two implementations of "the digest" is two things to
# keep in step. scripts/check_exact_inputs.py owns hashing, writes the manifest,
# and is also the regeneration check -- it rebuilds into a temporary directory
# and diffs the digests.

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
  k <- sub("^--", "", a)
  if (i == length(args)) stop("option --", k, " needs a value", call. = FALSE)
  opt[[k]] <- args[[i + 1L]]
  i <- i + 2L
}
`%||%` <- function(a, b) if (is.null(a)) b else a

repo_root <- normalizePath(opt[["repo"]] %||% ".", mustWork = TRUE)
spec_path <- normalizePath(opt$spec, mustWork = TRUE)
out_dir <- normalizePath(opt[["out-dir"]], mustWork = FALSE)
if (!nzchar(out_dir) || out_dir == "") stop("--out-dir is required", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required to read the fixture manifest", call. = FALSE)
}
spec <- jsonlite::fromJSON(spec_path, simplifyVector = FALSE)
if (!identical(spec$schema, "ancombc2-exact-fixtures/1")) {
  stop("fixture manifest schema must be ancombc2-exact-fixtures/1, got ",
       spec$schema %||% "<absent>", call. = FALSE)
}
cases <- spec$cases
if (!length(cases)) stop("fixture manifest names no cases", call. = FALSE)

only <- if (!is.null(opt$only)) trimws(strsplit(opt$only, ",", fixed = TRUE)[[1]]) else NULL
if (!is.null(only)) {
  unknown <- setdiff(only, vapply(cases, function(cs) cs$id, character(1)))
  if (length(unknown)) {
    stop("--only names cases not in the manifest: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  cases <- Filter(function(cs) cs$id %in% only, cases)
}

source(file.path(repo_root, "scripts", "exact_input_lib.R"))

.relpath <- function(path, root) {
  p <- normalizePath(path, mustWork = FALSE)
  r <- normalizePath(root, mustWork = TRUE)
  if (startsWith(p, paste0(r, "/"))) substring(p, nchar(r) + 2L) else p
}

ids <- vapply(cases, function(cs) cs$id, character(1))
if (anyDuplicated(ids)) {
  stop("duplicate case id(s) in the manifest: ",
       paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
}

index <- list()
for (cs in cases) {
  for (k in c("id", "expect")) {
    if (is.null(cs[[k]])) stop("case ", cs$id %||% "<no id>", ": missing ", k, call. = FALSE)
  }
  if (is.null(cs$expect$outcome)) {
    stop("case ", cs$id, ": expect.outcome is required (success or error); ",
         "a case with no declared outcome is a case nobody checks", call. = FALSE)
  }
  res <- build_exact_input(cs, repo_root)
  inp <- res$input
  out <- file.path(out_dir, paste0(cs$id, ".rds"))
  saveRDS(inp, out, version = 3)

  index[[length(index) + 1L]] <- list(
    id = cs$id,
    input = out,
    input_bytes = file.info(out)$size,
    fn = inp$fn,
    call_style = inp$call_style,
    expect = cs$expect,
    tags = as.character(cs$tags %||% character(0)),
    source = list(
      dir = cs$source$dir,
      counts = res$counts_path,
      meta = res$meta_path,
      config = if (is.null(res$cfg_path)) NA_character_ else res$cfg_path
    ),
    n_tax = res$n_tax,
    n_samp = res$n_samp,
    overrides = as.character(cs$overrides %||% character(0)),
    null_args = as.character(cs$null_args %||% character(0)),
    threads = as.integer(cs$threads %||% 1L),
    seed = inp$state$seed
  )
  cat(sprintf("%-34s %-8s n_tax=%-6d n_samp=%-6d %d bytes\n", cs$id,
              cs$expect$outcome, res$n_tax, res$n_samp, file.info(out)$size))
}

if (!is.null(opt$index)) {
  doc <- list(
    schema = "ancombc2-exact-input-index/1",
    spec = .relpath(spec_path, repo_root),
    out_dir = .relpath(out_dir, repo_root),
    profile = opt$profile %||% NA_character_,
    builder = "scripts/make_exact_inputs.R",
    hashed_by = "scripts/check_exact_inputs.py (this file carries no digests, ",
    hashing_note = "because R has no base sha256 and two implementations of ",
    cases = index
  )
  dir.create(dirname(opt$index), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite::toJSON(doc, auto_unbox = TRUE, pretty = TRUE, null = "null"),
             opt$index)
  cat(sprintf("wrote %s (%d cases)\n", opt$index, length(index)))
}
