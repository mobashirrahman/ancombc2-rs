# Minimal JSON writer/reader. The R installation here has no jsonlite, and the
# only consumer that matters is the Rust harness, so a self-contained writer is
# preferable to another dependency. The reader exists so the test generator can
# read back what it wrote (a round-trip check, not a general parser).

json_escape <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  s
}

json_num <- function(x) {
  if (length(x) != 1L) {
    # vapply() hands over one element at a time; anything else is a bug
    stop("json_num expects a scalar, got length ", length(x))
  }
  # An absent dimension is NA_integer_, not a double: it must serialise as the
  # JSON null, and as.numeric(NA_integer_) would print "NA", which is not JSON.
  if (is.na(x) && !is.nan(x)) return("null")
  x <- as.numeric(x)
  if (is.nan(x)) return("null")
  if (is.infinite(x)) return(if (x > 0) "1e999" else "-1e999")
  # 17 significant digits round-trips an f64 exactly
  sprintf("%.17g", x)
}

# A matrix is serialised column-major (R's native order) with an explicit
# nrow/ncol, which is also what the Rust reader expects.
to_json <- function(x) {
  if (is.null(x)) return("null")
  # Named scalars (the manifest is full of them) must stay scalars, not become
  # one-element arrays.
  if (is.atomic(x) && length(x) == 1L) {
    if (is.character(x)) return(sprintf("\"%s\"", json_escape(x)))
    if (is.logical(x)) return(if (x) "true" else "false")
    return(json_num(x))
  }
  if (is.list(x)) {
    nms <- names(x)
    # An *unnamed* list is a JSON array, not an object. It used to serialise as
    # `{}` -- `names()` is NULL, the `vapply` produced nothing, and every element
    # was silently dropped. That is not a formatting nicety: the fixture matrix
    # manifest is a list of cell objects, so it was written as `"cells":{}` and
    # the Rust side could not see a single cell. The golden contract has no
    # unnamed lists (verified: no manifest entry is an empty object), so this
    # cannot change a committed golden.
    if (is.null(nms)) {
      parts <- vapply(x, to_json, character(1))
      if (length(parts) == 0L) return("[]")
      return(paste0("[", paste(parts, collapse = ","), "]"))
    }
    parts <- vapply(nms, function(n) {
      sprintf("\"%s\":%s", json_escape(n), to_json(x[[n]]))
    }, character(1))
    if (length(parts) == 0L) return("{}")
    return(paste0("{", paste(parts, collapse = ","), "}"))
  }
  if (is.data.frame(x)) {
    # Column by column, preserving each column's own type: mixing a factor and a
    # numeric column into a matrix would coerce both to character.
    cols <- vapply(names(x), function(n) {
      sprintf("\"%s\":%s", json_escape(n), to_json(x[[n]]))
    }, character(1))
    return(paste0("{", paste(cols, collapse = ","), "}"))
  }
  # A factor is a character vector with a level attribute; as.character() drops
  # the attribute, so this must be tested before the numeric branch below.
  if (is.factor(x)) {
    return(to_json(as.character(x)))
  }
  if (is.matrix(x)) {
    if (is.character(x)) {
      return(paste0("{\"nr\":", nrow(x), ",\"nc\":", ncol(x),
                    ",\"v\":[", paste(sprintf("\"%s\"", json_escape(as.vector(x))),
                                       collapse = ","), "]}"))
    }
    return(paste0("{\"nr\":", nrow(x), ",\"nc\":", ncol(x),
                  ",\"v\":[", paste(vapply(as.vector(x), json_num, character(1)),
                                    collapse = ","), "]}"))
  }
  if (is.logical(x)) {
    return(paste0("[", paste(ifelse(as.logical(x), "true", "false"),
                                  collapse = ","), "]"))
  }
  if (is.numeric(x) || is.integer(x)) {
    return(paste0("[", paste(vapply(as.numeric(x), json_num, character(1)),
                             collapse = ","), "]"))
  }
  if (is.character(x)) {
    return(paste0("[", paste(sprintf("\"%s\"", json_escape(x)), collapse = ","), "]"))
  }
  "null"
}

write_json <- function(path, x) {
  writeLines(to_json(x), path)
  invisible(path)
}

write_raw_f64 <- function(path, x) {
  # canonical little-endian IEEE-754 doubles, C order
  con <- file(path, "wb")
  on.exit(close(con))
  v <- as.vector(as.matrix(x))
  v[is.nan(v)] <- NaN
  writeBin(as.numeric(v), con, size = 8, endian = "little")
  invisible(path)
}

# Quantities written to canonical f64 + JSON. Character / logical quantities are
# written to JSON only.
CANON_F64 <- c(
  "O1", "O2", "x", "y1", "y2",
  "beta_star", "theta", "var1", "beta_corr_stage1", "samp_frac",
  "y_bias_crt", "beta", "var_hat", "s02", "var_final", "se",
  "W", "p", "q", "delta_em", "delta_wls", "var_delta"
)
# vcov is a list of p x p matrices per taxon
CANON_VCOV <- c("vcov1", "vcov")
CANON_TEXT <- c("taxa_retained", "samples_retained", "fix_eff")
CANON_BOOL <- c("zero_ind")

# Nested list quantities (tables) go to JSON as-is.
write_canonical <- function(g, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()

  # Per-stage wall time, the contract's "per-stage timings" quantity.
  #
  # Written as JSON rather than as `.f64`, and deliberately *not* compared for
  # equality by the parity suite: wall-clock is not reproducible, so asserting on
  # it would assert that the oracle ran at the same speed. It is in the golden so
  # that both sides name the same stages and a stage which stopped being timed is
  # visible; `crates/ancombc2-core/tests/golden/mod.rs` checks exactly that, and
  # checks that the values are finite and non-negative.
  # The convergence trace: the MLE's epsilon at each iteration, plus the
  # tolerance and cap it ran under. Written as JSON because it is a small
  # vector of iteration bookkeeping rather than a matrix, and because the Rust
  # side reads it as JSON rather than as a canonical `.f64` blob -- it is
  # metadata about the run, not a result quantity.
  if (!is.null(g$convergence_trace) &&
      length(g$convergence_trace$epsilons) > 0) {
    write_json(file.path(outdir, "convergence_trace.json"), list(
      note = paste("Per-iteration epsilon of the first .iter_mle, captured from",
                   "the oracle's own verbose trace. Compared by level A: the",
                   "Rust run must take the same number of iterations to the same",
                   "epsilon, or one of them stopped early."),
      epsilons = as.list(g$convergence_trace$epsilons),
      iterations = as.list(g$convergence_trace$iterations),
      tol = g$convergence_trace$tol,
      max_iter = g$convergence_trace$max_iter
    ))
  }

  # The EM mixture parameters: the contract's "EM mixture parameters".
  #
  # Written as JSON with the terms named, rather than as canonical `.f64`, for the
  # same reason as the convergence trace: these are a small amount of iteration
  # bookkeeping about the bias fit, not a result matrix. The `.f64` side of the
  # contract already carries `delta_em`, and the Rust check reads this file to
  # confirm that its own mixture is the one the oracle fitted.
  # The missingness pattern assignment: the contract's Level A quantity of that
  # name. Written as JSON rather than `.f64` because it is an assignment of taxa to
  # patterns, not a measurement -- `group` is an integer pattern id, `key` is the
  # 0/1 row the oracle grouped on, and both are compared exactly.
  if (!is.null(g$pattern_assignment) && !is.null(g$pattern_assignment$group)) {
    pa <- g$pattern_assignment
    stopifnot(length(pa$group) == pa$reported_set_n_taxa)
    write_json(file.path(outdir, "pattern_assignment.json"), list(
      note = paste("Missingness pattern assignment, read off .lm_fit_all as it",
                   "exits. `group` is the 1-based pattern id in the reference's",
                   "first-appearance order (factor(keys, levels = unique(keys))),",
                   "`key` the 0/1 usable-sample row the grouping is derived from,",
                   "and `n_used` its row sum. Compared exactly: the assignment decides",
                   "which taxa share one QR, so a different assignment is a",
                   "different factorisation and a different beta. Verified on the",
                   "oracle side against two things it returns: each fitted taxon's",
                   "dof, which is n_used - p for the group it was fitted in, and the",
                   "literal zeros .lm_fit_all writes at the samples that group did",
                   "not use."),
      reported_set_n_taxa = pa$reported_set_n_taxa,
      n_p = pa$n_p,
      n_groups = pa$n_groups,
      group = as.list(pa$group),
      key = as.list(pa$key),
      n_used = as.list(pa$n_used),
      unfittable = as.list(pa$unfittable)
    ))
  }

  if (!is.null(g$em_mixture) && length(g$em_mixture$terms) > 0) {
    m <- g$em_mixture
    write_json(file.path(outdir, "em_mixture.json"), list(
      note = paste("Final three-component Gaussian mixture of .bias_em, read",
                   "off the oracle's own function as it exits. `pi` columns are",
                   "the component weights, `l` the component offsets (l1, l2),",
                   "`kappa` the component variance offsets (kappa1, kappa2),",
                   "and `delta` is the same bias the .f64 contract carries as",
                   "delta_em -- compared here to tie the mixture to it."),
      terms = as.list(m$terms),
      pi = lapply(seq_along(m$terms), function(i) as.list(m$pi[i, ])),
      delta = as.list(m$delta),
      l = lapply(seq_along(m$terms), function(i) as.list(m$l[i, ])),
      kappa = lapply(seq_along(m$terms), function(i) as.list(m$kappa[i, ])),
      iterations = as.list(m$iterations),
      epsilon = as.list(m$epsilon),
      tol = m$tol,
      max_iter = m$max_iter
    ))
  }

  if (!is.null(g$stage_seconds) && length(g$stage_seconds) > 0) {
    st <- g$stage_seconds
    write_json(file.path(outdir, "stage_seconds.json"), list(
      note = paste("Wall-clock seconds per stage, recorded from the oracle.",
                   "Not a compared quantity: wall time is not reproducible and",
                   "parity is asserted on the stage names and their being finite",
                   "and non-negative, never on the values."),
      seconds = as.list(st)
    ))
  }

  for (nm in CANON_F64) {
    v <- g[[nm]]
    if (is.null(v)) next
    nmv <- if (is.matrix(v)) sprintf("%s.f64", nm) else sprintf("%s.v.f64", nm)
    write_raw_f64(file.path(outdir, nmv), v)
    manifest[[nm]] <- list(
      file = nmv, kind = "f64",
      nr = if (is.matrix(v)) nrow(v) else length(v),
      nc = if (is.matrix(v)) ncol(v) else NA_integer_)
  }
  for (nm in CANON_VCOV) {
    v <- g[[nm]]
    if (is.null(v)) next
    k <- length(v)
    p <- if (k > 0) nrow(v[[1]]) else 0L
    flat <- unlist(lapply(v, as.vector), use.names = FALSE)
    write_raw_f64(file.path(outdir, sprintf("%s.f64", nm)), flat)
    manifest[[nm]] <- list(file = sprintf("%s.f64", nm), kind = "vcov",
                           n_tax = k, p = p)
  }
  for (nm in CANON_TEXT) {
    v <- g[[nm]]
    if (is.null(v)) next
    manifest[[nm]] <- list(kind = "text", value = as.character(v))
  }
  for (nm in CANON_BOOL) {
    v <- g[[nm]]
    if (is.null(v)) next
    write_json(file.path(outdir, sprintf("%s.json", nm)), as.matrix(v))
    manifest[[nm]] <- list(file = sprintf("%s.json", nm), kind = "logical_matrix",
                           nr = nrow(v), nc = ncol(v),
                           colnames = colnames(v), rownames = rownames(v))
  }
  for (nm in c("res", "res_global", "res_pair", "ss_tab", "diff_abn")) {
    v <- g[[nm]]
    if (is.null(v)) next
    write_json(file.path(outdir, sprintf("%s.json", nm)), v)
    manifest[[nm]] <- list(file = sprintf("%s.json", nm), kind = "table",
                           nr = nrow(v), ncol = ncol(v),
                           colnames = colnames(v), rownames = rownames(v))
  }
  if (!is.null(g$sens_pseudo)) {
    write_raw_f64(file.path(outdir, "sens_pseudo.f64"), g$sens_pseudo)
    manifest[["sens_pseudo"]] <- list(file = "sens_pseudo.f64", kind = "f64",
                                      nr = length(g$sens_pseudo), nc = NA_integer_)
  }
  write_json(file.path(outdir, "manifest.json"), manifest)
  # The manifest itself, not a summary. It is the only place the `.f64` shapes
  # (`nr`, `nc`) are recorded outside `manifest.json`, and this repository has a
  # JSON writer but deliberately no JSON reader -- so `scripts/compare_goldens.R`
  # reads `manifest.rds` to learn how to reshape the `.f64` payloads. It is a few
  # hundred bytes per golden.
  manifest
}

# Every quantity `write_canonical` puts a file in the golden's directory, and
# which is therefore recoverable without consulting `golden.rds`.
CANON_JSON_TABLES <- c("res", "res_global", "res_pair", "ss_tab", "diff_abn")
FILE_BACKED <- c(CANON_F64, CANON_VCOV, CANON_BOOL, CANON_JSON_TABLES,
                 "sens_pseudo", "stage_seconds",
                 "taxa_retained", "samples_retained", "fix_eff")

# `golden.rds` is a comparison artifact, not the payload store.
#
# It used to hold every quantity, which duplicated the `.f64`/`.json` files
# sitting beside it in the same directory. For `fx04` and the `shape-10000x500`
# cell that duplication was 148 MB and 133 MB -- and those two blobs are the
# only ones in the repository over GitHub's 100 MB hard limit, so it could not
# be pushed at all. Deleting the files would have broken the parity suite, which
# reads them; duplicating them in LFS would have hidden the duplication behind a
# 310 MB dependency. So they are written once.
#
# What stays is what the directory cannot supply: the numeric quantities with no
# file of their own. In practice that is `dof` and nothing else, since every
# other numeric quantity is in `CANON_F64`.
#
# No compared quantity is lost. `scripts/compare_goldens.R` only ever compared
# numeric top-level fields -- 24 of 41 for `fx04` -- and 23 of those 24 are
# `.f64`-backed; the 24th, `dof`, is what this function keeps. The non-numeric
# remainder was already skipped by that comparison, and the tables remain on
# disk as `.json`.
#
# `sens_ss_list` is the one quantity that leaves the repository entirely: it is
# assigned in `reference/R/harness.R`, read by nothing here, was never part of
# the compared contract, and retaining it cost 90 MB of `fx04`.
slim_golden <- function(g) {
  keep <- g[setdiff(names(g), FILE_BACKED)]
  Filter(is.numeric, keep)
}

write_session <- function(path) {
  writeLines(capture_session_lines(), path)
}

capture_session_lines <- function() {
  utils::capture.output(utils::sessionInfo())
}
