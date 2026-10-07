# Shared input construction for the exact path (IMPROVED_PLAN.md S04).
#
# Sourced by scripts/make_exact_input.R (one case) and scripts/make_exact_inputs.R
# (the whole manifest). One implementation, because two implementations of "read
# the committed fixture" would eventually disagree and the disagreement would look
# like a numerical difference.

# The declared R read path, reproduced from scripts/generate_goldens.R. The tab
# split is not incidental: `read.table(..., row.names = 1)` consumes `group` as
# the row-name label and shifts every value by one, which fails downstream with
# a mangled covariate name rather than at the point of the mistake.
read_committed_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  nms <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  cells <- lapply(lines[-1], function(l) strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  m <- do.call(rbind, cells)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- nms
  # `group` holds *labels*. In the committed fixtures those labels are the
  # strings "1" and "2", which parse as numbers, and converting them would make
  # `model.matrix` treat the group factor as a numeric slope instead of contrasts.
  for (k in seq_along(df)) {
    if (names(df)[k] == "group") next
    num <- suppressWarnings(as.numeric(df[[k]]))
    if (!anyNA(num)) df[[k]] <- num
  }
  # `unname()` is load-bearing. `vapply(lines[-1], ...)` inherits each line as
  # the result's `names`, so without it the row-name vector arrives carrying a
  # names attribute whose values are the whole tab-separated data rows. That is a
  # real difference in the object, not a cosmetic one: it serializes into the
  # input RDS, and `identical()` on a data.frame with such row names is FALSE
  # against one without them. The legacy reader in scripts/generate_goldens.R has
  # this defect; this copy does not.
  rownames(df) <- unname(vapply(lines[-1],
                                function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1],
                                character(1)))
  df
}

# One count cell, read explicitly.
#
# `as.numeric("NA")` returns NA *with* the warning "NAs introduced by coercion".
# That is harmless in a one-off read and unacceptable here: warnings are part of
# the compared condition set, so a warning raised by the input reader would sit
# in the capture as if `ancombc2` had produced it. The non-finite tokens are
# therefore recognised as tokens, and a token that is neither numeric nor one of
# them is a named read error rather than a silent NA -- a silently-missing count
# is exactly the kind of defect this whole tree exists to prevent.
parse_count_cell <- function(tok) {
  if (is.na(tok)) return(NA_real_)
  if (tok %in% c("", "NA")) return(NA_real_)
  if (tok == "NaN") return(NaN)
  if (tok == "Inf") return(Inf)
  if (tok == "-Inf") return(-Inf)
  v <- suppressWarnings(as.numeric(tok))
  if (is.na(v)) {
    stop("count file holds a non-numeric cell: ", sQuote(tok), call. = FALSE)
  }
  v
}

read_committed_counts <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  body <- lapply(lines[-1], function(l) {
    vapply(strsplit(l, "\t", fixed = TRUE)[[1]][-1], parse_count_cell, numeric(1))
  })
  m <- do.call(rbind, body)
  # The counts header carries only the sample names, with no leading label cell,
  # so it is split without dropping a field; meta.tsv's header *does* carry one.
  colnames(m) <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  # See read_committed_meta: without unname() the row names inherit the line text
  # as their names attribute.
  rownames(m) <- unname(vapply(lines[-1],
                               function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1],
                               character(1)))
  m
}

# Two config formats exist in this repository and both are read, neither
# pattern-matched: validation/fixtures/*/config.json (JSON) and
# validation/edge/*/analysis.tsv (key<TAB>value). The regex scalar reader that
# scripts/bench_r.R uses is one of the defects IMPROVED_PLAN.md S15 names, and it
# does not belong anywhere near the exact path.
require_jsonlite <- function(what) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to ", what,
         "; configs must be parsed, not pattern-matched", call. = FALSE)
  }
}

read_config_json <- function(path) {
  require_jsonlite("read a JSON config")
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

read_config_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(list())
  head <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  # `validation/edge/*/analysis.tsv` starts with a `key<TAB>value` header. Without
  # skipping it the config gains an entry called `key` whose value is the literal
  # string "value", every one of those seven cases is handed an unknown argument,
  # and the reference answers `unused argument (key = "value")`. Measured: all
  # seven edge cases produced exactly that before this line existed.
  if (identical(tolower(head[1]), "key") && length(head) >= 2L &&
      identical(tolower(head[2]), "value")) {
    lines <- lines[-1]
  }
  kv <- strsplit(lines, "\t", fixed = TRUE)
  kv <- kv[lengths(kv) >= 2L]
  out <- list()
  for (p in kv) {
    key <- p[1]
    val <- paste(p[-1], collapse = "\t")
    out[[key]] <- if (val %in% c("true", "false")) {
      identical(val, "true")
    } else if (grepl("^-?[0-9]+$", val)) {
      as.integer(val)
    } else if (grepl("^-?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?$", val)) {
      as.numeric(val)
    } else {
      val
    }
  }
  out
}

read_config <- function(path) {
  if (!file.exists(path)) return(list())
  if (grepl("[.]json$", path)) read_config_json(path) else read_config_tsv(path)
}

# Config keys whose name differs from the ancombc2 formal. Named lookup, so an
# unmapped key stays the key itself rather than vanishing.
CONFIG_ALIASES <- list(do_global = "global", do_pairwise = "pairwise",
                       sensitivity = "pseudo_sens", formula = "fix_formula")

# An R literal for a structured argument value, written in the manifest as
#   {"__rlit__": "list(tol = 0.01, max_iter = 20L, verbose = FALSE)"}
#
# JSON has no NULL and no nested list literal, but `iter_control`, `em_control`,
# `lme_control`, `mdfdr_control` and `trend_control` are all R lists, and
# `data`/`meta_data` for the secom_* exports are lists of matrices. Only
# `parse()` of the given text is used: no `eval` of arbitrary code with a side
# effect, and nothing here reads a file or reaches the network.
#
# It is a tagged object rather than a string prefix because it has to nest:
# `secom_dist(data = list(matrix(...)), meta_data = list(data.frame(...)))` is a
# list containing lists. Two earlier spellings were measured failing --
# `R:<expr>` collided with `secom_dist`'s own argument named `R` ("object 'R' not
# found"), and stripping a single `R:` prefix left `list(R:matrix(...))`, which R
# parses as the `:` operator and rejects with "object 'Rexpr' not found". A
# reserved key has neither problem and recurses.
R_LITERAL_KEY <- "__rlit__"

# The override spelling that removes an argument rather than passing it a value.
# See parse_override below.
OMIT_SENTINEL <- "__OMITTED__"

as_r_literal <- function(x) {
  if (is.list(x) && length(x) == 1L && !is.null(names(x)) && names(x) == R_LITERAL_KEY) {
    src <- x[[1L]]
    if (!is.character(src) || length(src) != 1L) {
      stop(R_LITERAL_KEY, " must hold a single character string of R source",
           call. = FALSE)
    }
    exprs <- parse(text = src)
    if (length(exprs) != 1L) {
      stop(R_LITERAL_KEY, " must be a single R expression, got ", length(exprs),
           call. = FALSE)
    }
    return(eval(exprs, envir = baseenv()))
  }
  x
}

# Interpret one `--overrides key=value` token.
#
#   key=__OMITTED__   remove the argument entirely. Measured: passing the string
#                     through instead reaches the package as an unknown argument
#                     and the reference answers
#                     `unused argument (sensitivity = "__OMITTED__")`.
#   key=__NULL__      keep the key with an explicit NULL.
#
# The two are different cases, and the original's `assay.type = assay_name` and
# `rank = tax_level` aliases make the difference observable, so both spellings
# exist rather than one overloading the other.
parse_override <- function(raw) {
  kv <- strsplit(raw, "=", fixed = TRUE)[[1]]
  if (length(kv) < 2L) stop("--overrides needs key=value, got ", raw, call. = FALSE)
  key <- kv[1]
  val <- paste(kv[-1], collapse = "=")
  if (identical(val, OMIT_SENTINEL)) return(list(key = key, value = NULL, omit = TRUE))
  if (identical(val, "__NULL__")) return(list(key = key, value = NULL, omit = FALSE))
  if (val %in% c("TRUE", "FALSE")) {
    return(list(key = key, value = identical(val, "TRUE"), omit = FALSE))
  }
  if (grepl("^-?[0-9]+$", val)) {
    return(list(key = key, value = as.integer(val), omit = FALSE))
  }
  if (grepl("^-?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?$", val)) {
    return(list(key = key, value = as.numeric(val), omit = FALSE))
  }
  list(key = key, value = val, omit = FALSE)
}

# `fix_formula` is a CHARACTER string in the pinned original, not a formula
# object: data_sanity_check.R does `gsub("\\*", "+", fix_formula)` and
# `strsplit(fix_formula, "\\s*\\+\\s*")` on it, and .ancombc2_sens_fit builds
# `stats::formula(paste0("y ~ ", fix_formula))`. A formula object would take a
# different path in the original, so it is not accepted here either.
resolve_fix_formula <- function(case, fixture_dir) {
  if (!is.null(case$formula)) {
    f <- case$formula
    if (length(f) != 1L || !is.character(f) || !nzchar(f)) {
      stop("case formula must be a single non-empty character string", call. = FALSE)
    }
    return(f)
  }
  p <- file.path(fixture_dir, "formula.txt")
  if (file.exists(p)) return(trimws(readLines(p, warn = FALSE)[1]))
  NULL
}

build_exact_input <- function(case, repo_root) {
  src <- case$source
  if (is.null(src) || is.null(src$dir)) stop("case ", case$id, ": source.dir is required")
  fixture_dir <- if (grepl("^/", src$dir)) src$dir else file.path(repo_root, src$dir)

  counts_path <- file.path(fixture_dir, "counts.tsv")
  meta_path <- file.path(fixture_dir, "meta.tsv")
  for (p in c(counts_path, meta_path)) {
    if (!file.exists(p)) stop("case ", case$id, ": missing ", p, call. = FALSE)
  }

  counts <- read_committed_counts(counts_path)
  # Counts are integers in the microbiome sense, so they stay integer and a
  # double-vs-integer difference is a candidate defect rather than an artefact of
  # how the input was loaded. NA is *kept*: validation/edge/na_counts exists to
  # pin how the original treats an NA count, and scripts/bench_r.R's
  # `counts[is.na(counts)] <- 0` would erase the very case that case is for. The
  # matrix falls back to double only if a non-integer value is genuinely present.
  if (all(is.na(counts) | counts == trunc(counts))) {
    storage.mode(counts) <- "integer"
  }
  meta <- read_committed_meta(meta_path)
  if (!is.null(case$group_column) && !is.null(meta[[case$group_column]])) {
    meta[[case$group_column]] <- factor(meta[[case$group_column]])
  } else if ("group" %in% names(meta)) {
    meta$group <- factor(meta$group)
  }

  cfg_path <- if (!is.null(src$config)) {
    if (grepl("^/", src$config)) src$config else file.path(fixture_dir, src$config)
  } else if (file.exists(file.path(fixture_dir, "config.json"))) {
    file.path(fixture_dir, "config.json")
  } else {
    NULL
  }
  cfg <- if (is.null(cfg_path)) list() else read_config(cfg_path)

  fix_formula <- resolve_fix_formula(case, fixture_dir)

  # `args_from = "none"` builds the argument list from `args`/`overrides` alone.
  # It exists for the other exports: ancombc(), ancom(), secom_*() and sim_plnm()
  # have their own signatures, and handing them ancombc2's arguments produces
  # `unused argument (fix_formula = "group + x1")` -- measured on all six before
  # this switch existed.
  mode <- case$args_from %||% "fixture"
  from_fixture <- !identical(mode, "none")
  # secom_dist()/secom_linear() index data[[i]] and meta_data[[i]]: they take a
  # *list* of abundance matrices, one per rank. `fixture-list` wraps the same
  # fixture counts and metadata in a one-element list rather than inventing a
  # second input path for them.
  as_list <- identical(mode, "fixture-list")

  call_args <- list()
  if (from_fixture) {
    call_args[["data"]] <- if (as_list) list(counts) else counts
    call_args[["meta_data"]] <- if (as_list) list(meta) else meta
    if (!as_list) {
      call_args[["taxa_are_rows"]] <- TRUE
      if (!is.null(fix_formula)) call_args[["fix_formula"]] <- fix_formula
    }
  }

  if (from_fixture && !as_list) {
    for (k in names(cfg)) {
      kk <- if (!is.null(CONFIG_ALIASES[[k]])) CONFIG_ALIASES[[k]] else k
      # `group_levels` is an edge-case-harness directive, not an ancombc2 argument.
      if (kk == "group_levels") next
      call_args[[kk]] <- cfg[[k]]
    }

    call_args[["n_cl"]] <- as.integer(case$threads %||% 1L)
    # `verbose` stays at the original's default (TRUE) unless a case overrides it,
    # so the default-behaviour path is exercised from the first case.
    if (!is.null(case$verbose)) call_args[["verbose"]] <- isTRUE(case$verbose)
  }

  for (ov in case$overrides %||% character(0)) {
    p <- parse_override(ov)
    if (isTRUE(p$omit)) {
      call_args[[p$key]] <- NULL
    } else {
      call_args[[p$key]] <- p$value
    }
  }
  for (k in names(case$args %||% list())) {
    v <- case$args[[k]]
    # "__OMITTED__" means the argument is not passed at all. It is spelled out
    # rather than implied by leaving the key out of `args`, because a case that
    # wants to *state* that an argument is absent should be visible in the
    # manifest instead of being invisible in it.
    if (identical(v, OMIT_SENTINEL)) {
      call_args[[k]] <- NULL
      next
    }
    v <- as_r_literal(v)
    if (identical(v, quote(expr = NULL)) || (is.character(v) && identical(v, "__NULL__"))) {
      stop("case ", case$id, ": argument ", k,
           " cannot be set to NULL here. `args.<name> = NULL` removes the ",
           "argument; to pass an explicitly NULL use a dedicated input built for ",
           "it, because an R list cannot hold both.", call. = FALSE)
    }
    call_args[[k]] <- v
  }

  # Dropping a NULL argument is not the same as passing it. `call_args[[k]] <-
  # NULL` removes the element; a case that wants an explicitly-NULL argument must
  # say so, which is what `null_args` is for.
  for (k in case$null_args %||% character(0)) call_args[[k]] <- NULL

  style <- case$call_style %||% "named"
  if (identical(style, "positional")) {
    # Positional cases name the arguments in formal order and the names are
    # dropped at call time. Anything not supplied takes the original's default,
    # which is the whole point of the case.
    call_args <- call_args[order(match(names(call_args), ancombc2_formal_order()))]
    call_args <- call_args[!vapply(call_args, is.null, logical(1))]
    names(call_args) <- NULL
  }

  state <- list(
    seed = as.integer(case$state$seed %||% 42L),
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    options = list(stringsAsFactors = FALSE)
  )
  if (!is.null(case$state$options)) state$options <- case$state$options

  list(
    input = list(
      schema = "ancombc2-exact-input/1",
      case_id = case$id,
      fn = case$fn %||% "ancombc2",
      call_style = style,
      args = call_args,
      state = state,
      source = list(
        kind = src$kind %||% "committed",
        dir = src$dir,
        counts = src$dir,
        meta = file.path(src$dir, "meta.tsv"),
        config = if (is.null(cfg_path)) NA_character_ else cfg_path,
        formula = if (file.exists(file.path(fixture_dir, "formula.txt")))
          file.path(src$dir, "formula.txt") else NA_character_
      ),
      provenance = list(
        built_by = "scripts/make_exact_inputs.R",
        note = "Both arms load this file. The candidate must not read anything else."
      )
    ),
    n_tax = nrow(counts), n_samp = ncol(counts),
    counts_path = counts_path, meta_path = meta_path, cfg_path = cfg_path
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# The formal order of the pinned ancombc2, recorded here rather than read from the
# package so a positional case is pinned to the *documented* signature. It is
# asserted against the installed package by check_exact_inputs.py.
ancombc2_formal_order <- function() {
  c("data", "taxa_are_rows", "assay.type", "assay_name", "rank", "tax_level",
    "aggregate_data", "meta_data", "fix_formula", "rand_formula", "p_adj_method",
    "pseudo", "pseudo_sens", "conservative", "prv_cut", "lib_cut", "s0_perc",
    "group", "struc_zero", "neg_lb", "alpha", "n_cl", "verbose", "global",
    "pairwise", "dunnet", "trend", "iter_control", "em_control", "lme_control",
    "mdfdr_control", "trend_control")
}
