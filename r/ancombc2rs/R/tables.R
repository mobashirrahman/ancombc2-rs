#' The compatibility target of the linked core
#'
#' @return A string naming the ANCOMBC version, commit, R version and seed that
#'   the core reproduces, or `NA_character_` when the core is not linked.
#' @export
ancombc2_compat_target <- function() {
  if (!ancombc2_available()) return(NA_character_)
  .Call(C_ancombc2_rs_version)
}

#' Whether the Rust core is linked in
#'
#' @return `TRUE` when the shared object loaded and its entry points resolved,
#'   `FALSE` otherwise. A wrapper that silently reported "not available" without
#'   a reason would be worse than one that fails at install time, so the
#'   registration in `src/init.c` is strict: a renamed Rust symbol is a load
#'   error, not a runtime surprise.
#' @export
ancombc2_available <- function() {
  !is.null(tryCatch(getNativeSymbolInfo("C_ancombc2_rs_run", PACKAGE = "ancombc2rs"),
                    error = function(e) NULL))
}

call_core <- function(request) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("the `jsonlite` package is required", call. = FALSE)
  }
  # `digits = NA` is the load-bearing argument: the default is 4 significant
  # digits, which would round the design matrix and the counts on the way across
  # and put the R path on a different input from the CLI's.
  text <- jsonlite::toJSON(request, auto_unbox = TRUE, digits = NA, null = "null")
  raw <- .Call(C_ancombc2_rs_run, as.character(text))
  if (is.null(raw) || length(raw) != 1L || is.na(raw) || !nzchar(raw)) {
    stop("the ancombc2-rs core returned nothing", call. = FALSE)
  }
  parsed <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
  if (!is.null(parsed$error)) {
    stop(sprintf("ancombc2-rs: %s", parsed$error), call. = FALSE)
  }
  # `simplifyVector = FALSE` keeps NA out of the numbers, which is what the golden
  # contract needs, but it also leaves the string fields as lists. Those are
  # scalars-as-strings and a character vector would be wrong for a caller, so they
  # are unwrapped here rather than left for `primary_table` to trip over.
  for (k in c("fix_eff", "taxa", "taxa_bias", "warnings")) {
    if (!is.null(parsed[[k]])) parsed[[k]] <- unlist(parsed[[k]], use.names = FALSE)
  }
  if (!is.null(parsed$samples)) parsed$samples <- unlist(parsed$samples, use.names = FALSE)
  parsed$res <- primary_table(parsed)
  parsed
}

# `stats::model.matrix` drops rows with a missing covariate by default. The
# reference switches `na.action` off for the design and back on afterwards, keeps
# the row, and lets the core report an incomplete design as an unidentifiable
# covariate -- rather than a silently shorter table that still fits.
with_local_na_pass <- function(expr) {
  old <- options(na.action = "na.pass")
  on.exit(options(old), add = TRUE)
  force(expr)
}

# The input types the reference's `data_sanity_check` accepts: a matrix, a
# data.frame, a phyloseq object, or a (Tree)SummarizedExperiment. For the two S4
# types the reference reaches for `microbiome::abundances()` and
# `microbiome::meta()` -- or, for a TSE, `mia::convertToPhyloseq()` first -- and
# requires those packages to be installed.
#
# Those packages are *Suggests*, not *Imports*, so the wrapper must not fail to
# load when they are absent. Each accessor therefore has a fallback that reaches
# the same slots through the defining package (`phyloseq::otu_table`,
# `SummarizedExperiment::assay`), and failing that, reads the S4 slots directly.
# The fallback is not a different computation -- it is the same feature table by
# another route -- and `test-coerce.R` asserts the three agree on an object built
# to have the phyloseq class.
#
# What is deliberately *not* reproduced is `tax_level` aggregation. Aggregating
# taxa is a `microbiome`/`mia` operation, not a numerical one, and silently
# passing the argument through without aggregating would produce results for a
# different model than the user asked for. It is rejected with a message naming
# the packages instead.
is_phyloseq <- function(x) inherits(x, "phyloseq")
is_tse <- function(x) inherits(x, c("TreeSummarizedExperiment", "SummarizedExperiment"))

require_or_stop <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "the '%s' package is needed to %s but is not installed. Install it, or pass \
a count matrix and the sample metadata separately.",
      pkg, why
    ), call. = FALSE)
  }
  invisible(TRUE)
}

# `phyloseq` keeps the feature table in the `@otu_table` slot and the sample
# metadata in `@sam_data`. Those slots are readable with `attr()` on an S4 object
# without the class's own package being installed, and they hold exactly what
# `microbiome::abundances()` and `microbiome::meta()` return for a phyloseq input
# -- so reading them is the same computation by a shorter route, not an
# approximation of one.
#
# The slot is therefore tried *first*, and `microbiome` is the fallback for an
# object whose slots are not named this way. That ordering matters for more than
# tidiness: reaching for `phyloseq::otu_table()` requires the package to be
# installed, and on a machine where something else has registered a class by that
# name the `::` lookup resolves against the wrong namespace and fails with a
# message about a package that is not installed at all.
slot_or_null <- function(obj, name) {
  v <- attr(obj, name, exact = TRUE)
  if (is.null(v)) NULL else v
}

abundances_of <- function(obj) {
  otu <- slot_or_null(obj, "otu_table")
  if (!is.null(otu)) return(as.matrix(otu))
  if (requireNamespace("microbiome", quietly = TRUE)) {
    return(as.matrix(microbiome::abundances(obj)))
  }
  NULL
}

meta_of <- function(obj) {
  sd <- slot_or_null(obj, "sam_data")
  if (!is.null(sd)) return(as.data.frame(sd))
  if (requireNamespace("microbiome", quietly = TRUE)) {
    return(microbiome::meta(obj))
  }
  NULL
}

# `x` has already been converted to a phyloseq object by the time this runs, so
# both phyloseq and TSE inputs share one implementation, as they do in the
# reference.
coerce_input <- function(data, meta_data, tax_level, assay_type) {
  if (!is.null(tax_level)) {
    stop(paste0(
      "`tax_level` aggregation is not implemented: it is a 'microbiome'/'mia' ",
      "operation on the taxonomy, not a numerical one, and silently ignoring it ",
      "would return results for a different model than the one requested. ",
      "Aggregate the feature table yourself and pass the result as `data`."
    ), call. = FALSE)
  }
  if (is_phyloseq(data)) {
    counts <- abundances_of(data)
    if (is.null(counts)) stop("could not read the feature table", call. = FALSE)
    meta <- if (is.null(meta_data)) meta_of(data) else meta_data
    return(list(counts = counts, meta = meta))
  }
  if (is_tse(data)) {
    obj <- data
    if (requireNamespace("mia", quietly = TRUE)) {
      obj <- mia::convertToPhyloseq(data, assay.type = assay_type)
    }
    counts <- abundances_of(obj)
    if (is.null(counts)) {
      # No `microbiome`, no `phyloseq`: read the assay and the colData directly.
      # `assay.type` selects among the assays, defaulting to the first, which is
      # what `SummarizedExperiment::assay(se, assay.type)` does with a name.
      require_or_stop("SummarizedExperiment", "read a (Tree)SummarizedExperiment")
      counts <- as.matrix(SummarizedExperiment::assay(data, assay_type))
      if (is.null(meta_data)) {
        meta <- as.data.frame(SummarizedExperiment::colData(data))
      } else {
        meta <- meta_data
      }
      return(list(counts = counts, meta = meta))
    }
    meta <- if (is.null(meta_data)) meta_of(obj) else meta_data
    return(list(counts = counts, meta = meta))
  }
  NULL
}

as_count_matrix <- function(data) {
  if (is.data.frame(data)) {
    rn <- rownames(data)
    data <- as.matrix(data)
    rownames(data) <- rn
  }
  if (!is.matrix(data)) {
    stop("`data` must be a matrix, a data.frame, a phyloseq object, or a \
(Tree)SummarizedExperiment", call. = FALSE)
  }
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("sample_", seq_len(ncol(data)))
  }
  if (is.null(rownames(data))) {
    rownames(data) <- paste0("taxon_", seq_len(nrow(data)))
  }
  storage.mode(data) <- "double"
  data
}

as_meta_data <- function(meta_data, samples) {
  if (!is.data.frame(meta_data)) {
    stop("`meta_data` must be a data.frame with samples as rows", call. = FALSE)
  }
  rownames(meta_data) <- samples
  meta_data
}

check_group_requirement <- function(group, global, pairwise, struc_zero) {
  if (is.null(group) && any(c(global, pairwise, struc_zero))) {
    # The reference's own wording, so a user migrating sees the same message.
    stop(paste0("Group variable is required for the multi-group comparison\n",
                "`group` is `NULL` while some of the arguments ",
                "(`global`, `pairwise`) are `TRUE`"), call. = FALSE)
  }
  invisible(TRUE)
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (length(a) == 1L && is.na(a)) return(b)
  a
}

#' @export
print.ancombc2_result <- function(x, ...) {
  cat("<ancombc2-rs result>\n")
  cat("  target: ", ancombc2_compat_target(), "\n", sep = "")
  cat("  taxa retained: ", length(x$taxa), "\n", sep = "")
  cat("  coefficients: ", paste(x$fix_eff, collapse = ", "), "\n", sep = "")
  if (length(x$warnings)) {
    cat("  warnings:\n")
    for (w in x$warnings) cat("    - ", gsub("\n", " ", w), "\n", sep = "")
  }
  cat("  $res holds the primary table; see `?ancombc2`.\n")
  invisible(x)
}

#' @export
summary.ancombc2_result <- function(object, ...) {
  p <- object$p
  # `unlist()` drops `NULL`, and the core serialises an unfittable coefficient as
  # JSON `null` -- so a taxon with no estimate would silently vanish from the
  # table and every column after it would be attributed to the wrong taxon. The
  # positions are preserved as `NA` instead. See the note on `mat()` in
  # `ancombc2.R` for the full account.
  is_mat <- function(name) {
    v <- object[[name]]
    if (is.null(v)) return(NULL)
    n <- length(object$taxa)
    k <- length(object$fix_eff)
    flat <- vapply(v, function(x) {
      if (is.null(x)) NA_real_ else as.numeric(x)[1L]
    }, numeric(1))
    if (length(flat) != n * k) {
      stop("the core returned ", length(flat), " `", name, "` value(s) for ",
           n, " taxa x ", k, " fix_eff = ", n * k,
           "; the table would be mis-shaped", call. = FALSE)
    }
    matrix(flat, nrow = n, byrow = TRUE,
           dimnames = list(object$taxa, object$fix_eff))
  }
  structure(list(
    n_taxa = length(object$taxa),
    fix_eff = object$fix_eff,
    beta = is_mat("beta"),
    se = is_mat("se"),
    p = is_mat("p"),
    q = is_mat("q"),
    diff_abn = is_mat("diff_abn"),
    ml_iterations = object$ml_iterations,
    warnings = object$warnings,
    res = object$res
  ), class = "summary.ancombc2_result")
}

#' @export
print.summary.ancombc2_result <- function(x, ...) {
  cat("ancombc2-rs summary\n")
  cat("  taxa retained: ", x$n_taxa, "\n", sep = "")
  cat("  MLE iterations: ", x$ml_iterations, "\n", sep = "")
  if (!is.null(x$diff_abn)) {
    calls <- colSums(x$diff_abn, na.rm = TRUE)
    for (i in seq_along(calls)) {
      cat(sprintf("  %-20s %d significant\n", x$fix_eff[i], calls[i]))
    }
  }
  if (length(x$warnings)) {
    cat("  warnings:\n")
    for (w in x$warnings) cat("    - ", gsub("\n", " ", w), "\n", sep = "")
  }
  invisible(x)
}
