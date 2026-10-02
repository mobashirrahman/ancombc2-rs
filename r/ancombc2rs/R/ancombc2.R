#' ANCOM-BC2 differential abundance, backed by Rust
#'
#' `ancombc2()` is a wrapper for the `ancombc2-rs` core. The analysis runs in
#' Rust; this function's job is to build the design matrix exactly as
#' `ANCOMBC::ancombc2()` does, marshal the inputs across a C ABI, and return the
#' result tables with the reference's names and column order.
#'
#' Building the design here rather than in Rust is deliberate. The design matrix
#' is `model.matrix(formula, meta_data)`, and the column *order and names* are part
#' of the compatibility contract -- the global and pairwise tests identify group
#' columns by name. Re-implementing `model.matrix` in Rust would be a second
#' source of divergence for no gain.
#'
#' @param data A count matrix, taxa as rows and samples as columns, with an ID
#'   column in the corner. A `data.frame`, a `phyloseq` object, or a
#'   (Tree)`SummarizedExperiment` are also accepted, exactly as in the reference.
#'   For the two S4 types the feature table is read from the `@otu_table` slot and
#'   the sample metadata from `@sam_data` -- the values `microbiome::abundances()`
#'   and `microbiome::meta()` would return -- so `phyloseq` and `microbiome` are
#'   not required to be installed.
#' @param taxa_are_rows Whether taxa are in the rows of a matrix or data.frame.
#'   Ignored for the S4 types, which carry their own orientation. The reference's
#'   default is `TRUE`.
#' @param assay.type Which assay to read from a (Tree)`SummarizedExperiment`.
#' @param tax_level Aggregate counts to this taxonomic level before the analysis.
#'   **Not implemented**: aggregation is a `microbiome`/`mia` operation rather than
#'   a numerical one, and silently ignoring it would return results for a different
#'   model than the one requested. Aggregate the table yourself and pass it as
#'   `data`.
#' @param meta_data A `data.frame` of sample metadata, samples as rows. Optional
#'   for a `phyloseq` or (Tree)`SummarizedExperiment` input, in which case it is
#'   taken from the object; required for a matrix or data.frame, which has no such
#'   source.
#' @param fix_formula A one-sided formula naming the fixed effects, e.g.
#'   `~ group + age`.
#' @param group The name of the grouping column in `meta_data`, as a string.
#'   Required for `global`, `pairwise` and `struc_zero`, as in the reference.
#' @param p_adj_method One of `"holm"`, `"hochberg"`, `"hommel"`,
#'   `"bonferroni"`, `"BH"`, `"BY"`, `"none"`.
#' @param pseudo Added to every count before the log transform.
#' @param pseudo_sens Run the pseudo-count sensitivity analysis.
#' @param conservative Rerun the whole algorithm per pseudo-count. The default,
#'   `TRUE`, is the reference's. `FALSE` estimates the sampling fractions once and
#'   refits only the inference 50 times.
#' @param prv_cut,lib_cut Prevalence and library-size cutoffs.
#' @param s0_perc Quantile of the SE distribution used for the `s0` regulariser.
#' @param struc_zero Detect structural zeros.
#' @param neg_lb Classify structural zeros by the asymptotic lower bound.
#' @param alpha Significance level.
#' @param global,pairwise Run the multi-group tests.
#' @param mdfdr_control A list with `fwer_ctrl_method`, the family-wise method
#'   for the pairwise test.
#' @param iter_control,em_control Lists with `tol` and `max_iter`.
#' @param compat `"ancombc2-2.15"` to reproduce the reference including its
#'   quirks, or `"strict"` for the mathematically intended behaviour where the two
#'   differ. See `vignette("reference-behaviour")` and
#'   `docs/reference_behavior.md`.
#' @param ... Ignored, with a warning. Present so that a caller passing a
#'   reference argument the core does not implement gets told, rather than having
#'   it silently dropped by `...`.
#'
#' @return An object of class `ancombc2_result`: a list with `res` (the primary
#'   table), `res_global`, `res_pair`, `zero_ind`, `samp_frac`, `delta_em`,
#'   `delta_wls`, `var_delta`, and the `warnings` the core recorded.
#'
#' @examples
#' \dontrun{
#' res <- ancombc2(data = counts, meta_data = meta, fix_formula = ~ group + age,
#'                 group = "group")
#' head(res$res)
#' }
#' @export
ancombc2 <- function(data, meta_data = NULL, fix_formula, group = NULL,
                     taxa_are_rows = TRUE,
                     assay.type = c("counts", "abundance", "relabundance"),
                     tax_level = NULL,
                     p_adj_method = c("holm", "hochberg", "hommel",
                                      "bonferroni", "BH", "BY", "none"),
                     pseudo = 0, pseudo_sens = FALSE, conservative = TRUE,
                     prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                     struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05,
                     global = FALSE, pairwise = FALSE,
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                     iter_control = list(tol = 0.01, max_iter = 20),
                     em_control = list(tol = 1e-05, max_iter = 100),
                     compat = c("ancombc2-2.15", "strict"),
                     ...) {
  p_adj_method <- match.arg(p_adj_method)
  compat <- match.arg(compat)
  extra <- list(...)
  if (length(extra)) {
    # Naming the ignored arguments is the whole point: a caller who passed
    # `rand_formula` and got a result has been told nothing useful.
    warning(sprintf(
      "ignoring unsupported argument(s): %s. The fixed-effects path is implemented; \
random effects, the trend test and Dunnett's test are not.",
      paste(sQuote(names(extra)), collapse = ", ")
    ), call. = FALSE)
  }
  if (missing(data) || missing(fix_formula)) {
    stop("`data` and `fix_formula` are both required", call. = FALSE)
  }
  assay.type <- match.arg(assay.type)
  if (!ancombc2_available()) {
    stop(
      "the ancombc2-rs core is not linked into this package. Build it with \
`cargo build --release -p ancombc2-ffi` and re-install, or see r/ancombc2rs/README.md.",
      call. = FALSE
    )
  }

  # phyloseq and (Tree)SummarizedExperiment inputs carry their own sample
  # metadata, so `meta_data` is optional for them and, as in the reference, is
  # taken from the object when it is not supplied. A matrix or data.frame has no
  # such source and still needs it.
  coerced <- coerce_input(data, meta_data, tax_level, assay.type)
  if (is.null(coerced)) {
    if (is.null(meta_data)) {
      stop("Missing sample metadata. Please provide the sample metadata in \
'data.frame' format.", call. = FALSE)
    }
    counts <- as_count_matrix(data)
    if (!taxa_are_rows) counts <- t(as.matrix(counts))
  } else {
    counts <- as_count_matrix(coerced$counts)
    meta_data <- coerced$meta
  }
  if (!is.data.frame(meta_data)) {
    stop("`meta_data` must be a data.frame with samples as rows", call. = FALSE)
  }
  if (!all(colnames(counts) %in% rownames(meta_data))) {
    stop(paste0("Sample names do not match between the feature table and sample metadata.",
                "\nPlease ensure the column names of the feature abundance matrix and the row names of the metadata data frame are consistent.",
                sep = "\n"), call. = FALSE)
  }
  meta <- as_meta_data(meta_data, colnames(counts))
  check_group_requirement(group, global, pairwise, struc_zero)

  # `data_sanity_check` validates the group *before* the design is built, and its
  # three thresholds are not interchangeable: fewer than two categories is a hard
  # error, but fewer than three only *deactivates* the multi-group comparisons.
  # The check has to live here, before `model.matrix`, because a one-level factor
  # would otherwise fail inside R with a contrasts error that says nothing about
  # ANCOM-BC2.
  group_levels <- if (is.null(group)) 0L else {
    g <- meta[[group]]
    length(unique(g[!is.na(g)]))
  }
  if (!is.null(group) && group_levels < 2L) {
    stop("The group variable should have >= 2 categories", call. = FALSE)
  }

  # `options(na.action = "na.pass")` around the model matrix, as
  # `.ancombc2_prep` does: a missing covariate must produce an NA row that the
  # core reports, not a silently shorter design.
  design <- with_local_na_pass({
    stats::model.matrix(stats::reformulate(attr(stats::terms(fix_formula), "term.labels"),
                                           intercept = TRUE), meta)
  })
  fix_eff <- colnames(design)

  group_index <- NULL
  if (!is.null(group)) {
    g <- meta[[group]]
    if (group_levels < 3L && (global || pairwise)) {
      # The core deactivates the comparisons and records the warning; this
      # surfaces it at the point the user asked for it. Emitting it twice would
      # be noise, so the wrapper warns and the core's copy is the one that ends
      # up in `res$warnings`.
      warning(paste0("The group variable has < 3 categories\n",
                     "The multi-group comparisons (global/pairwise) will be deactivated"),
              call. = FALSE)
    }
    if (is.numeric(g)) {
      # `data_sanity_check` coerces a non-numeric group to a factor and leaves a
      # numeric one alone; either way the *level order* is the sorted numeric
      # order, so the dummy columns come out the same.
      g <- factor(g, levels = sort(unique(g[!is.na(g)])))
    } else {
      g <- factor(g)
    }
    group_index <- as.integer(g) - 1L
  }

  request <- list(
    counts = as.numeric(t(as.matrix(counts))),
    n_taxa = nrow(counts),
    n_samp = ncol(counts),
    design = as.numeric(design),
    fix_eff = fix_eff,
    taxon_names = rownames(counts),
    sample_names = colnames(counts),
    group = group_index,
    group_name = group,
    p_adj_method = p_adj_method,
    pseudo = pseudo,
    pseudo_sens = pseudo_sens,
    conservative = conservative,
    prv_cut = prv_cut,
    lib_cut = lib_cut,
    s0_perc = s0_perc,
    struc_zero = struc_zero,
    neg_lb = neg_lb,
    alpha = alpha,
    global = global,
    pairwise = pairwise,
    mdfdr_fwer = mdfdr_control$fwer_ctrl_method %||% "holm",
    iter_tol = iter_control$tol %||% 0.01,
    iter_max = as.integer(iter_control$max_iter %||% 20L),
    em_tol = em_control$tol %||% 1e-05,
    em_max = as.integer(em_control$max_iter %||% 100L),
    compat = compat
  )

  out <- call_core(request)
  # The deactivation is the reference's behaviour and the core's, and the core
  # records it in `warnings`; the check above only reports it earlier so the design
  # matrix can still be built from a two-level factor.
  structure(out, class = "ancombc2_result")
}

# The FFI returns the columns; the reference's *table* is assembled here, because
# its column order is part of the contract and is easier to read in R than to
# thread through JSON. `.ancombc2_prep` builds it with
# `cbind(taxon, beta_prim, se_prim, W_prim, p_prim, q_prim, diff_prim)`, so all the
# `lfc_` columns come first, then all the `se_`, and so on.
primary_table <- function(parsed) {
  n <- length(parsed$taxa)
  k <- length(parsed$fix_eff)
  # `unlist()` is wrong here. The core serialises an `f64::NAN` as JSON `null`,
  # which `jsonlite` hands back as `NULL`, and `unlist()` **drops** `NULL`
  # elements -- so a taxon whose coefficients could not be estimated silently
  # disappears and the reshape comes up short. That is what produced
  # "data length [39] is not a sub-multiple ... of the number of rows [14]" on a
  # 14-taxon fixture with `p = 3`, where 42 values were expected: the one taxon
  # observed in a single group level has no estimable coefficient, its three
  # entries are `null`, and `unlist()` took them out.
  #
  # A mis-shaped result table with a reshape warning is worse than a failure -- it
  # reads as a table, and every column after the missing taxon is attributed to the
  # wrong one. So the positions are preserved here, as `NA`, and the length is
  # checked.
  mat <- function(name) {
    v <- parsed[[name]]
    if (is.null(v)) {
      stop("the core returned no `", name, "`", call. = FALSE)
    }
    flat <- vapply(v, function(x) {
      if (is.null(x)) NA_real_ else as.numeric(x)[1L]
    }, numeric(1))
    if (length(flat) != n * k) {
      stop("the core returned ", length(flat), " `", name,
           "` value(s) for ", n, " taxa x ", k, " fix_eff = ", n * k,
           "; the table would be mis-shaped", call. = FALSE)
    }
    matrix(flat, nrow = n, byrow = TRUE)
  }
  parts <- list(taxon = parsed$taxa)
  for (col in c("beta", "se", "w", "p", "q")) {
    v <- mat(col)
    for (j in seq_len(k)) {
      parts[[paste0(prefix_of(col), "_", parsed$fix_eff[j])]] <- v[, j]
    }
  }
  d <- mat("diff_abn")
  for (j in seq_len(k)) {
    parts[[paste0("diff_", parsed$fix_eff[j])]] <- as.logical(d[, j])
  }
  if (!is.null(parsed$sensitivity$passed_ss)) {
    ps <- matrix(unlist(parsed$sensitivity$passed_ss), nrow = n, byrow = TRUE)
    dr <- matrix(unlist(parsed$sensitivity$diff_robust), nrow = n, byrow = TRUE)
    for (j in seq_len(k)) {
      parts[[paste0("passed_ss_", parsed$fix_eff[j])]] <- as.logical(ps[, j])
      parts[[paste0("diff_robust_", parsed$fix_eff[j])]] <- as.logical(dr[, j])
    }
  }
  as.data.frame(parts, stringsAsFactors = FALSE, check.names = FALSE)
}

prefix_of <- function(col) {
  switch(col, beta = "lfc", se = "se", w = "W", p = "p", q = "q", col)
}
