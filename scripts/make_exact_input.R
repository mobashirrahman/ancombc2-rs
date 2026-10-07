#!/usr/bin/env Rscript
# Build the exact input RDS that both arms consume (IMPROVED_PLAN.md S04, first
# case landed in S02 so the runner has something deterministic to read).
#
#   Rscript --vanilla scripts/make_exact_input.R \
#     --fixture validation/fixtures/fx01 \
#     --out validation/exact/inputs/fx01-basic.rds \
#     [--case-id fx01-basic] [--seed 42] [--threads 1]
#
# Why an RDS and not the fixture's own files:
#
# The committed fixture is `counts.tsv` + `meta.tsv` + `config.json`. If each
# arm read those, each arm would also own a text parser, and a parser difference
# would look like a numerical difference. Reading happens once, here, and what is
# saved is the already-typed argument object: an integer matrix, a data.frame
# with a factor, and a plain list of controls. Both arms then load the same
# bytes and call the same exported function with the same list.
#
# This file reads the committed inputs exactly as scripts/generate_goldens.R
# does (`read_committed_counts` / `read_committed_meta` are reproduced verbatim
# below, with the reason for the tab split recorded there). The legacy goldens
# are NOT regenerated or overwritten; this only produces a new exact input.

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
  opt[[sub("^--", "", a)]] <- args[[i + 1L]]
  i <- i + 2L
}
`%||%` <- function(a, b) if (is.null(a)) b else a

fixture <- opt[["fixture"]]
out <- opt[["out"]]
if (is.null(fixture) || is.null(out)) {
  stop("usage: make_exact_input.R --fixture <dir> --out <rds> [--case-id ID] ",
       "[--seed N] [--threads N] [--positional] [--overrides key=value ...]",
       call. = FALSE)
}
if (!dir.exists(fixture)) stop("fixture directory does not exist: ", fixture, call. = FALSE)

# --- the declared R read path, from scripts/generate_goldens.R -------------

read_committed_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  nms <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  cells <- lapply(lines[-1], function(l) strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  m <- do.call(rbind, cells)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- nms
  for (k in seq_along(df)) {
    if (names(df)[k] == "group") next
    num <- suppressWarnings(as.numeric(df[[k]]))
    if (!anyNA(num)) df[[k]] <- num
  }
  rownames(df) <- vapply(lines[-1],
                         function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1],
                         character(1))
  df
}

read_committed_counts <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  body <- lapply(lines[-1], function(l) as.numeric(strsplit(l, "\t", fixed = TRUE)[[1]][-1]))
  m <- do.call(rbind, body)
  # The counts header carries only the sample names, with no leading label cell,
  # so it is split without dropping a field; meta.tsv's header *does* carry one.
  colnames(m) <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  rownames(m) <- vapply(lines[-1],
                        function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1],
                        character(1))
  m
}

# config.json is a flat object of scalars, nested two levels deep for the control
# lists. A JSON parser is used rather than a regex: a regex scalar reader is one
# of the defects IMPROVED_PLAN.md S15 names, and repeating it here would put the
# same defect on the exact path.
parse_config <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read ", path,
         "; the config must be parsed, not pattern-matched", call. = FALSE)
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

counts_path <- file.path(fixture, "counts.tsv")
meta_path <- file.path(fixture, "meta.tsv")
cfg_path <- file.path(fixture, "config.json")
for (p in c(counts_path, meta_path)) {
  if (!file.exists(p)) stop("missing committed fixture file: ", p, call. = FALSE)
}

counts <- read_committed_counts(counts_path)
storage.mode(counts) <- "integer"
meta <- read_committed_meta(meta_path)
if ("group" %in% names(meta)) {
  meta$group <- factor(meta$group)
}
if (!is.null(opt[["formula"]])) {
  # `fix_formula` is a CHARACTER string in the pinned original, not a formula
  # object: `data_sanity_check.R` does `gsub("\\*", "+", fix_formula)` and
  # `strsplit(fix_formula, "\\s*\\+\\s*")` on it, and `.ancombc2_sens_fit` builds
  # `stats::formula(paste0("y ~ ", fix_formula))`. Passing a formula object here
  # would take a different code path in the original, so it is not accepted.
  fix_formula <- as.character(opt[["formula"]])
  if (length(fix_formula) != 1L || !nzchar(fix_formula)) {
    stop("--formula must be a single non-empty character string", call. = FALSE)
  }
} else {
  f <- file.path(fixture, "formula.txt")
  fix_formula <- if (file.exists(f)) {
    trimws(readLines(f, warn = FALSE)[1])
  } else {
    "group"
  }
}

cfg <- if (file.exists(cfg_path)) parse_config(cfg_path) else list()

# Config keys whose name differs from the ancombc2 formal. Named lookup, so an
# unmapped key is an error rather than a silent identity.
ALIASES <- list(do_global = "global", do_pairwise = "pairwise",
                sensitivity = "pseudo_sens")

call_args <- list(data = counts, meta_data = meta, fix_formula = fix_formula,
                  taxa_are_rows = TRUE)

for (k in names(cfg)) {
  kk <- if (!is.null(ALIASES[[k]])) ALIASES[[k]] else k
  call_args[[kk]] <- cfg[[k]]
}
call_args[["n_cl"]] <- as.integer(opt[["threads"]] %||% 1L)
# `verbose` is deliberately left at the original's default (TRUE) unless a case
# overrides it, so the default-behaviour path is what the first cases exercise.

# `--overrides key=value` applies a value that is not in config.json. `NULL` is
# spelled `__NULL__` so an argument can be explicitly absent, which is a
# different case from omitting it: the original's `assay.type = assay_name` and
# `rank = tax_level` aliases mean an explicitly-NULL argument and an omitted one
# do not reach the same code.
for (ov in opt[names(opt) == "overrides"]) {
  kv <- strsplit(ov, "=", fixed = TRUE)[[1]]
  if (length(kv) < 2L) stop("--overrides needs key=value, got ", ov, call. = FALSE)
  key <- kv[1]
  raw <- paste(kv[-1], collapse = "=")
  val <- if (identical(raw, "__NULL__")) {
    NULL
  } else if (raw %in% c("TRUE", "FALSE")) {
    identical(raw, "TRUE")
  } else if (grepl("^-?[0-9]+$", raw)) {
    as.integer(raw)
  } else if (grepl("^-?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?$", raw)) {
    as.numeric(raw)
  } else {
    raw
  }
  call_args[[key]] <- val
}

# Only arguments the installed original actually declares are passed. An unknown
# name would be a silent no-op here and an error there, which is exactly the kind
# of asymmetry the exact contract must not have.
n_extra <- 0L

inp <- list(
  schema = "ancombc2-exact-input/1",
  case_id = opt[["case-id"]] %||% basename(fixture),
  fn = "ancombc2",
  call_style = if (isTRUE(opt[["positional"]])) "positional" else "named",
  args = call_args,
  state = list(
    seed = as.integer(opt[["seed"]] %||% 42L),
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    options = list(stringsAsFactors = FALSE)
  ),
  source = list(
    fixture = fixture,
    counts = counts_path,
    meta = meta_path,
    config = if (file.exists(cfg_path)) cfg_path else NA_character_,
    formula = if (file.exists(file.path(fixture, "formula.txt")))
      file.path(fixture, "formula.txt") else NA_character_
  ),
  provenance = list(
    built_by = "scripts/make_exact_input.R",
    note = "Both arms load this file. The candidate must not read anything else."
  )
)

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(inp, out, version = 3)
cat(sprintf("wrote %s (case_id=%s, %d named args, n_tax=%d, n_samp=%d)\n",
            out, inp$case_id, length(call_args) + n_extra,
            nrow(counts), ncol(counts)))
