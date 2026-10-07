# The R side of the typed native bridge.
#
# This file is deliberately the *only* place that knows the bridge exists. The
# retained upstream R code is byte-identical to the pinned commit and must stay
# that way; the bridge is reached from here, not by editing `R/ancombc2.R`.
#
# Nothing here serialises. There is no JSON, no string encoding of a number, and
# no `null` for a missing value: every number crosses as a number in its own type,
# with its bits intact. See `docs/exact_seam_inventory.md` section 6.

# Probe selectors, mirrored in crates/ancombc2-rbridge/src/abi.rs and
# r/ANCOMBC/src/init.c. Named constants rather than bare integers so a mismatch is
# a comparison failure with a readable name rather than a wrong buffer.
.rb_echo <- c(data_int = 0L, data_real = 1L, aggregate = 2L, design = 3L,
             design_complete = 4L, group_index = 5L)

#' Is the typed native bridge present in this installation?
#'
#' Not an error: a source installation built with
#' `ANCOMBC_RBRIDGE_SKIP_BUILD=1` and no prebuilt library has a working package
#' with no bridge, and that is a state worth being able to ask about rather than
#' one worth crashing in.
rb_available <- function() {
  # Deliberately no trailing `n = 0L`. `PACKAGE` is a formal of `.Call`, so
  # anything after it is passed to the C function as a real argument, and
  # `C_ancombc2_rb_oracle_sha` takes none.
  #
  # The test is that the call *returns a string*, not that it returns a function.
  # A registered routine is reached through a NativeSymbolInfo and `.Call` returns
  # its value -- a character scalar here. Asking `is.function()` of that was the
  # bug this replaced, and it reported a working bridge as absent.
  #
  # The call is the real probe: an installation built without the Rust library
  # fails to load `ANCOMBC.so` at all, and one that loads but cannot reach the
  # bridge raises here. Either way the answer is "not available", not a crash.
  is.character(tryCatch(.Call(C_ancombc2_rb_oracle_sha, PACKAGE = "ANCOMBC"),
                        error = function(e) NULL))
}

#' The pinned original this bridge reproduces.
rb_oracle_sha <- function() {
  .Call(C_ancombc2_rb_oracle_sha, PACKAGE = "ANCOMBC")
}

#' The bridge's own version string.
rb_version <- function() {
  .Call(C_ancombc2_rb_version, PACKAGE = "ANCOMBC")
}

#' Check that `x` is a numeric matrix, giving it a dim if it is a bare vector.
#'
#' **The type is left alone.** The count matrices are not widened: `NA_integer_` is
#' `INT_MIN` and the returned `feature_table` is `O2` with whatever type it
#' arrived as, so widening here would be a schema change dressed as a convenience.
#' An earlier version of this file called `as_double_matrix` on `data` and every
#' integer case then arrived at the bridge as a double.
as_numeric_matrix <- function(x, field) {
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix, not %s", field,
                 paste(class(x), collapse = "/")), call. = FALSE)
  }
  if (length(dim(x)) == 1L) dim(x) <- c(length(x), 1L)
  if (length(dim(x)) != 2L) {
    stop(sprintf("`%s` must have exactly two dimensions", field), call. = FALSE)
  }
  x
}

#' Widen to double. Used only for the design.
#'
#' `.ancombc2_core` calls `model.matrix`, which is always double, so an integer
#' design would be widened by the original anyway; doing it here rather than in
#' Rust keeps `REAL()` on the R side of the boundary.
as_double_matrix <- function(x, field) {
  x <- as_numeric_matrix(x, field)
  if (is.integer(x)) storage.mode(x) <- "double"
  x
}

#' Run one transport probe.
#'
#' Returns a list with the echoed buffer, the four name vectors as they came back
#' out of the bridge, and the `(nrow, ncol)` of the count matrix. Nothing is
#' interpreted: this is the transport, and S06's acceptance is that it loses
#' nothing.
#'
#' @param which one of `data_int`, `data_real`, `aggregate`, `design`,
#'   `design_complete`, `group_index`.
#' @param data,aggregate_data count matrices. Integer is preserved.
#' @param design the design matrix, `n_samp x p`.
#' @param design_complete `complete.cases(design)`, a logical of length `n_samp`.
#'   Sent as raw 0/1 bytes rather than as `LGLSXP`, because a mask is a byte.
#' @param group_index integer of length `n_samp`: a **1-based** level index, or
#'   `0` for a sample with no group. `0` rather than `NA` because `NA` would have
#'   to survive as `INT_MIN` and `0` reads unambiguously in a diagnostic.
#' @param group_labels `levels(group)`.
#' @param controls a named list of the numeric controls.
#' @param flags a named list of the logical switches.
rb_probe <- function(which,
                     data,
                     aggregate_data = data,
                     design,
                     design_complete = rep(TRUE, nrow(design)),
                     group_index = NULL,
                     group_labels = character(0),
                     fix_eff = colnames(design),
                     taxon_names = rownames(data),
                     sample_names = colnames(data),
                     controls = list(pseudo = 0, s0_perc = 0.05, prv_cut = 0.1,
                                     lib_cut = 0, alpha = 0.05,
                                     iter_tol = 0.01, iter_max = 20L,
                                     em_tol = 1e-5, em_max = 100L, mdfdr_b = 100L),
                     flags = list(global = FALSE, pairwise = FALSE, dunnet = FALSE,
                                  trend = FALSE, pseudo_sens = TRUE,
                                  conservative = TRUE, struc_zero = FALSE,
                                  neg_lb = FALSE)) {
  which <- match.arg(which, names(.rb_echo))
  selector <- .rb_echo[[which]]

  data <- as_numeric_matrix(data, "data")
  aggregate_data <- as_numeric_matrix(aggregate_data, "aggregate_data")
  design <- as_double_matrix(design, "design")

  # The integer probe reads `data` as INTSXP. Refusing here keeps the error on the
  # R side, where the argument name still means something; the bridge refuses it
  # again, and a probe selected against the wrong table type must not silently
  # answer a question nobody asked.
  if (which == "data_int" && !is.integer(data)) {
    stop("`data` must be an integer matrix for the integer probe, not ",
         paste(class(data), collapse = "/"), call. = FALSE)
  }
  if (which == "data_real" && is.integer(data)) {
    stop("`data` must be a double matrix for the real probe, not integer",
         call. = FALSE)
  }

  n_samp <- nrow(design)
  if (length(design_complete) != n_samp) {
    stop(sprintf("`design_complete` has length %d but the design has %d rows",
                 length(design_complete), n_samp), call. = FALSE)
  }
  if (anyNA(design_complete)) {
    stop("`design_complete` must not contain NA", call. = FALSE)
  }

  if (is.null(group_index)) group_index <- integer(n_samp)
  group_index <- as.integer(group_index)
  if (length(group_index) != n_samp) {
    stop(sprintf("`group_index` has length %d but there are %d samples",
                 length(group_index), n_samp), call. = FALSE)
  }
  if (is.null(fix_eff)) fix_eff <- character(ncol(design))
  if (is.null(taxon_names)) taxon_names <- rownames(data)
  if (is.null(taxon_names)) taxon_names <- as.character(seq_len(nrow(data)))
  if (is.null(sample_names)) sample_names <- colnames(data)
  if (is.null(sample_names)) sample_names <- as.character(seq_len(ncol(data)))

  num <- function(name, default) {
    v <- controls[[name]]
    if (is.null(v)) default else as.numeric(v)[1L]
  }
  int <- function(name, default) {
    v <- controls[[name]]
    if (is.null(v)) default else as.integer(v)[1L]
  }
  flg <- function(name) {
    v <- flags[[name]]
    if (is.null(v)) FALSE else as.logical(v)[1L]
  }

  # Capacity: enough for whichever probe is selected.
  cap <- switch(which,
                data_int = nrow(data) * ncol(data),
                data_real = nrow(data) * ncol(data),
                aggregate = nrow(aggregate_data) * ncol(aggregate_data),
                design = nrow(design) * ncol(design),
                design_complete = n_samp,
                group_index = n_samp,
                stop("unreachable"))

  args <- list(
    data, aggregate_data, design,
    as.raw(as.integer(design_complete)),
    group_index,
    as.character(group_labels),
    as.character(fix_eff),
    as.character(taxon_names),
    as.character(sample_names),
    num("pseudo", 0), num("s0_perc", 0.05), num("prv_cut", 0.1),
    num("lib_cut", 0), num("alpha", 0.05),
    num("iter_tol", 0.01), int("iter_max", 20L),
    num("em_tol", 1e-5), int("em_max", 100L), int("mdfdr_b", 100L),
    flg("global"), flg("pairwise"), flg("dunnet"), flg("trend"),
    flg("pseudo_sens"), flg("conservative"), flg("struc_zero"), flg("neg_lb"),
    selector, cap)

  out <- .Call(C_ancombc2_rb_probe, args, PACKAGE = "ANCOMBC")

  list(
    which = which,
    payload = out[[1L]],
    group_labels = out[[2L]],
    fix_eff = out[[3L]],
    taxon_names = out[[4L]],
    sample_names = out[[5L]],
    dim = out[[6L]]
  )
}

#' The raw 64 bits of each element of a double vector, as hex.
#'
#' `identical(a, b)` is TRUE for `-0` and `0`, and `NA_real_` and `NaN` print the
#' same, so a transport claim cannot be checked with either. This is what it is
#' checked with: one 16-character word per element, in order.
#'
#' A first version returned the whole concatenated byte buffer, which made
#' `rb_bits(x)[1]` the first *byte* rather than the first element -- and three
#' checks passed or failed for the wrong reason. One word per element is the only
#' shape where an index means what it looks like it means.
rb_bits <- function(x) {
  b <- writeBin(as.double(x), raw(), size = 8L, endian = "big")
  n <- length(b) %/% 8L
  if (n == 0L) return(character(0))
  m <- matrix(as.integer(b), nrow = 8L)
  apply(m, 2L, function(col) paste(sprintf("%02x", col), collapse = ""))
}

#' The raw 32-bit integers behind a double vector.
rb_int_bits <- function(x) {
  as.integer(x)
}

# ---------------------------------------------------------------------------
# The output transport
# ---------------------------------------------------------------------------
#
# `rb_emit_payload()` is the mirror image of `rb_probe()`: where the probe copies
# R's inputs into the bridge and brings typed values back, this copies the core's
# payloads out of R and brings them back again, bit for bit, in the types and
# shapes `crates/ancombc2-rbridge/src/output.rs` declares.
#
# It computes nothing. Before any stage's output can be trusted, the path it
# travels has to be shown to lose nothing -- and the only way to show that is to
# put a known payload through it and compare the whole result's bytes.

#' Flatten a list of matrices the way the output wire carries them.
#'
#' The wire is **row-major within each taxon**, while R stores matrices
#' column-major, so `t()` appears here and in exactly one other place:
#' `src/init.c`'s list rebuild. Keeping both visible is deliberate -- a symmetric
#' covariance transposed is a covariance transposed, and would pass every symmetry
#' check while producing wrong standard errors for `p >= 3`.
rb_flatten_vcov <- function(vcov_hat) {
  if (!is.list(vcov_hat)) {
    stop("`vcov_hat` must be a list of matrices, not ", class(vcov_hat)[1L])
  }
  if (length(vcov_hat) == 0L) {
    return(numeric(0))
  }
  p <- nrow(vcov_hat[[1L]])
  vapply(vcov_hat, function(m) {
    if (!is.matrix(m) || nrow(m) != p || ncol(m) != p) {
      stop("every element of `vcov_hat` must be a ", p, " x ", p, " matrix")
    }
    as.vector(t(m))
  }, numeric(p * p))
}

#' The `dimnames` the original gives each covariance matrix.
#'
#' `.sandwich_vcov` sets both margins to the fixed-effect names, and `assemble.R`
#' puts them back. Captured separately rather than reconstructed from `fix_eff`,
#' because reconstructing it is exactly the step that could go wrong and look
#' right.
rb_vcov_dimnames <- function(vcov_hat) {
  lapply(vcov_hat, dimnames)
}

#' Round-trip the core's payloads through the native bridge.
#'
#' Returns a list with `payload` (named by the plan), `plan_names`, `types`,
#' `rtype`, `dims`, `taxa`, `taxa_bias`, `taxon_names`, `sample_names`, `fix_eff`,
#' `dof_is_int`, and `vcov_dimnames`.
#'
#' `taxa` and `taxa_bias` are indices into the *input* table, not into `O2`.
#' `n_tax` is `nrow(O2)`, so `taxa_bias` can be longer than `taxa`: `struc_zero`
#' drops rows from `O2` that `O1` keeps, because the flagged taxa still estimate
#' the sampling fractions. Collapsing the two would silently report the wrong
#' number of taxa whenever a structural zero was flagged.
rb_emit_payload <- function(int) {
  n_tax <- nrow(int$O2)
  n_samp <- ncol(int$O2)
  p <- ncol(int$x)
  ybc <- int$y_bias_crt
  if (is.data.frame(ybc)) {
    # A data.frame on the wire would be a second, differently-shaped
    # representation of a matrix; `as.matrix` here and `data.frame()` in
    # `assemble.R` are the pair that cancels out.
    ybc <- as.matrix(ybc)
  }
  taxa <- if (is.null(int$taxa)) seq_len(n_tax) else as.integer(int$taxa)
  taxa_bias <- if (is.null(int$taxa_bias)) taxa else as.integer(int$taxa_bias)
  fix_eff <- if (is.null(int$fix_eff)) colnames(int$x) else as.character(int$fix_eff)
  # `O2` usually has rownames, but "usually" is not a contract, and `pack_names`
  # refuses NULL rather than inventing `"1"`, `"2"`. Named fallbacks instead, so a
  # missing name is a visible value and not a silent one.
  taxon_names <- rownames(int$O2)
  if (is.null(taxon_names)) {
    taxon_names <- if (is.null(int$tax_name)) as.character(seq_len(n_tax)) else as.character(int$tax_name)
  }
  sample_names <- colnames(int$O2)
  if (is.null(sample_names)) {
    sample_names <- as.character(seq_len(n_samp))
  }
  if (length(taxon_names) != n_tax || length(sample_names) != n_samp) {
    stop("`taxon_names` has ", length(taxon_names), " entries and `sample_names` has ",
         length(sample_names), ", but O2 is ", n_tax, " x ", n_samp)
  }

  args <- list(
    n_tax = as.double(n_tax),
    n_samp = as.double(n_samp),
    p = as.double(p),
    n_taxa_bias = as.double(length(taxa_bias)),
    beta_hat = int$beta_hat,
    var_hat = int$var_hat,
    dof = int$dof,
    vcov_hat = rb_flatten_vcov(int$vcov_hat),
    y_bias_crt = ybc,
    theta_hat = int$theta_hat,
    delta_em = int$delta_em,
    delta_wls = int$delta_wls,
    var_delta = int$var_delta,
    s02 = int$s02,
    O2 = int$O2,
    taxa = taxa,
    taxa_bias = taxa_bias,
    taxon_names = taxon_names,
    sample_names = sample_names,
    fix_eff = fix_eff
  )
  out <- .Call(C_ancombc2_rb_emit, args, PACKAGE = "ANCOMBC")
  out <- rb_restore_names(out)
  out
}

#' Put the names back on the transported payloads.
#'
#' The wire carries values and shapes; names are reconstructed here, in one place,
#' from the three name vectors that *did* travel. `serialize()` compares names and
#' encoding flags as firmly as it compares doubles, so an unnamed `samp_frac` is as
#' wrong as a perturbed one.
#'
#' The rules are not uniform -- `theta_hat` is named by the samples and `s02` by the
#' fixed effects, while `delta_em`, `delta_wls` and `var_delta` carry no names at
#' all -- so they are written out per entry rather than inferred. Inferring them
#' from a shape would be exactly the sort of "almost right" that costs an afternoon.
rb_restore_names <- function(out) {
  tax <- out$taxon_names
  samp <- out$sample_names
  fix <- out$fix_eff
  dn <- function(a, b) list(a, b)
  pl <- out$payload
  dimnames(pl$beta_hat) <- dn(tax, fix)
  dimnames(pl$var_hat) <- dn(tax, fix)
  dimnames(pl$dof) <- dn(tax, fix)
  dimnames(pl$O2) <- dn(tax, samp)
  dimnames(pl$y_bias_crt) <- dn(tax, samp)
  names(pl$theta_hat) <- samp
  names(pl$s02) <- fix
  # `s02` comes from `apply(var_hat, 2, quantile)`, so `apply` gave it the column
  # names. `delta_em`, `delta_wls` and `var_delta` are plain vectors upstream and
  # stay plain: giving them names would be a difference the byte comparison sees.
  dn_vcov <- out$vcov_dimnames
  if (is.null(dn_vcov)) {
    dn_vcov <- rep(list(dn(fix, fix)), length(pl$vcov_hat))
  }
  for (i in seq_along(pl$vcov_hat)) {
    dimnames(pl$vcov_hat[[i]]) <- dn_vcov[[i]]
  }
  out$payload <- pl
  out
}

# ---------------------------------------------------------------------------
# The preprocessing stage probe (IMPROVED_PLAN.md S08)
# ---------------------------------------------------------------------------
#
# `rb_preprocess_stages()` runs every preprocessing stage in the bridge and hands
# each stage's array back. `scripts/check_preprocess_stages.R` then runs the
# *reference's own expressions* on the same input and compares the bytes.
#
# The split matters: R supplies the semantics, Rust supplies the arithmetic, and the
# comparison says which stage differs. A whole-result comparison would say only that
# something, somewhere, differs.

#' Run the preprocessing stages in the bridge and return each stage's array.
#'
#' @param data `feature_table`, `n_tax x n_samp`.
#' @param aggregate `feature_table_aggregate`, the same shape.
#' @param meta_data The full `meta_data`, for `group`.
#' @param group The grouping column name, or `NULL` when `struc_zero` is `FALSE`.
#' @param struc_zero,neg_lb,prv_cut,lib_cut,pseudo The reference's controls.
#' @return A named list of the nineteen stages, plus `n_tax`, `n_samp`, `n_groups`
#'   and `has_group`.
rb_preprocess_stages <- function(data, aggregate, meta_data = NULL, group = NULL,
                                 struc_zero = FALSE, neg_lb = FALSE,
                                 prv_cut = 0, lib_cut = 0, pseudo = 0) {
  data <- as.matrix(data)
  aggregate <- as.matrix(aggregate)
  if (!identical(dim(data), dim(aggregate))) {
    stop("`feature_table` is ", paste(dim(data), collapse = " x "),
         " and `feature_table_aggregate` is ",
         paste(dim(aggregate), collapse = " x "),
         "; the reference requires the same shape")
  }
  n_tax <- nrow(data)
  n_samp <- ncol(data)

  # `G` is the `n_samp x n_group` 0/1 indicator, column-major as R stores it, built
  # exactly as `.get_struc_zero` builds it: a sample whose group label is missing has
  # an all-zero row and belongs to no group.
  n_groups <- 0L
  gind <- NULL
  if (!is.null(group)) {
    if (is.null(meta_data)) {
      stop("`group` is given but `meta_data` is not")
    }
    grp <- factor(meta_data[, group])
    grp_int <- as.integer(grp)
    grp_ok <- !is.na(grp_int)
    n_groups <- nlevels(grp)
    gind <- matrix(0, nrow = n_samp, ncol = n_groups)
    gind[cbind(which(grp_ok), grp_int[grp_ok])] <- 1
  }

  # The four wrappers go in as arguments rather than being looked up by name from C.
  # See `eval_na_rm_reduction` in `src/init.c` for why: `Rf_findNamespace` is not in
  # the public API, and building `rowMeans(x, na.rm = TRUE)` by hand in C has two
  # failure modes that both report as a missing function.
  args <- list(
    data, aggregate, gind,
    as.double(n_tax), as.double(n_samp), as.double(n_groups),
    as.double(struc_zero), as.double(prv_cut), as.double(lib_cut), as.double(pseudo),
    as.double(neg_lb),
    rb_row_means_na_rm, rb_col_means_na_rm, rb_col_sums_na_rm, rb_row_sums_na_rm
  )
  .Call(C_ancombc2_rb_preprocess_probe, args, PACKAGE = "ANCOMBC")
}

# ---------------------------------------------------------------------------
# R's reductions, as one-argument wrappers
# ---------------------------------------------------------------------------
#
# `src/init.c` calls these by name. Each exists so that the call it constructs from C
# is `f(x)` -- positional, one formal -- rather than `rowMeans(x, na.rm = TRUE)`,
# which needs a named argument built by hand in C. Getting that hand-built form wrong
# fails as "could not find function na.rm" rather than as a name mismatch, which is a
# poor way to spend an afternoon.
#
# They are the *only* place the reductions are written, so there is one spelling of
# `na.rm = TRUE` to be right about.

rb_row_means_na_rm <- function(x) base::rowMeans(x, na.rm = TRUE)
rb_col_means_na_rm <- function(x) base::colMeans(x, na.rm = TRUE)
rb_col_sums_na_rm <- function(x) base::colSums(x, na.rm = TRUE)
rb_row_sums_na_rm <- function(x) base::rowSums(x, na.rm = TRUE)

# ---------------------------------------------------------------------------
# The least-squares probe (IMPROVED_PLAN.md S09)
# ---------------------------------------------------------------------------
#
# `rb_fit_probe()` runs the pipeline's *own* `lm_fit_all` on the same `x` and `Y` the
# reference's `.lm_fit_all` receives, and hands back `beta`, `fitted` and `dof`. It
# computes the fit the estimator will compute -- not a copy made for the test -- so a
# difference here is a difference in the code that will run.
#
# `scripts/check_fit_stages.R` compares those three against `.lm_fit_all` transcribed
# from the pinned `ancombc_bias_correct.R:11-84`, and against `stats::lm.fit` directly
# for a single group.

#' Run the pipeline's `lm_fit_all` and return its three outputs.
#'
#' @param x The design, `n_samp x p`.
#' @param y The theta-adjusted response, `n_taxa x n_samp`.
#' @param observed Logical, `n_taxa x n_samp`: `TRUE` where the cell is usable.
#' @return A list of `beta` (`n_taxa x p`), `fitted` (`n_taxa x n_samp`) and `dof`
#'   (`n_taxa`).
rb_fit_probe <- function(x, y, observed) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  observed <- as.matrix(observed)
  n_samp <- nrow(x)
  p <- ncol(x)
  n_taxa <- nrow(y)
  if (ncol(y) != n_samp) {
    stop("`y` has ", ncol(y), " columns and `x` has ", n_samp, " rows")
  }
  if (!identical(dim(observed), dim(y))) {
    stop("`observed` is ", paste(dim(observed), collapse = " x "),
         " and `y` is ", paste(dim(y), collapse = " x "))
  }
  args <- list(x, y, observed, as.double(n_samp), as.double(p), as.double(n_taxa))
  .Call(C_ancombc2_rb_fit_probe, args, PACKAGE = "ANCOMBC")
}
