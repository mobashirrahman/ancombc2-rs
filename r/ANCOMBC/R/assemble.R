# Rebuilding the original's result from transported payloads.
#
# This is the other half of the output transport. `rb_emit_payload()` brings the
# core's payloads back out of the bridge; this turns them into the list the
# original's `.ancombc2_core()` returns.
#
# Every expression below is transcribed from the pinned
# `reference/ANCOMBC/R/ancombc_prep.R` step 4, in its order, with its operators
# and its `data.frame` construction intact. It is not a reimplementation: the
# rounding, the `is.na` masking, the `p_hat[is.na(p_hat)] = 1` substitution and the
# column-naming order are all load-bearing for the byte comparison, and every one
# of them is reproduced rather than improved on. Section references are to the
# pinned file.
#
# The reference is read at *development* time. This file never opens
# `reference/` or `validation/exact/payloads/`: a candidate run may not read either
# one, and a function that only works when a directory happens to exist is not
# part of the replacement.

#' Assemble the original's core result from transported payloads.
#'
#' @param pay The list `rb_emit_payload()` returns.
#' @param p_adj_method Passed through to [p.adjust()]; `ancombc_prep.R` lines
#'   253 and 264.
#' @param alpha The significance level, `ancombc_prep.R` line 254.
#' @param tax_name The reported taxon names, used for `res$taxon`.
#' @param global,pairwise,dunnet,trend The four test switches. Only `FALSE` is
#'   transported yet; see the note below.
#' @return The ten-element list `.ancombc2_core()` returns.
ancombc2_assemble_core <- function(pay, p_adj_method = "BH", alpha = 0.05,
                                   tax_name = NULL, global = FALSE, pairwise = FALSE,
                                   dunnet = FALSE, trend = FALSE) {
  pl <- pay$payload
  need <- c("beta_hat", "var_hat", "dof", "theta_hat", "delta_em", "delta_wls",
            "y_bias_crt", "O2")
  missing <- setdiff(need, names(pl))
  if (length(missing)) {
    stop("the transported payload is missing ", paste(missing, collapse = ", "))
  }

  # Ancombc-prep.R lines 226-241, *already applied*.
  #
  # The transported `var_hat` is the value the original finished with, not the one
  # it started with: `t(t(var_hat) + s02)` and `var_hat[is.na(beta_hat)] = NA` have
  # both run on the original's side, and `vcov_hat`'s diagonal has been overwritten
  # with the same numbers. So `se_hat <- sqrt(var_hat)` and nothing else.
  #
  # Adding `s02` here was the first version of this file, and it was wrong in a way
  # that looked right: it inflated every standard error, changed every p-value, and
  # would still have been "nearly right" on any case where `s0_perc = 0`. The
  # transported payload is a capture of the core's *outputs*, not of its inputs.
  var_hat <- pl$var_hat
  beta_hat <- pl$beta_hat
  dof <- pl$dof
  se_hat <- sqrt(var_hat)

  # Ancombc-prep.R lines 239-256. The primary test.
  W <- beta_hat / se_hat
  p_hat <- 2 * stats::pt(abs(W), df = dof, lower.tail = FALSE)
  p_hat[is.na(p_hat)] <- 1
  q_hat <- apply(p_hat, 2, function(x) stats::p.adjust(x, method = p_adj_method))
  diff_abn <- q_hat <= alpha & !is.na(q_hat)

  beta_prim <- data.frame(beta_hat, check.names = FALSE)
  se_prim <- data.frame(se_hat, check.names = FALSE)
  W_prim <- data.frame(W, check.names = FALSE)
  p_prim <- data.frame(p_hat, check.names = FALSE)
  q_prim <- data.frame(q_hat, check.names = FALSE)
  diff_prim <- data.frame(diff_abn, check.names = FALSE)
  colnames(beta_prim) <- paste0("lfc_", colnames(beta_hat))
  colnames(se_prim) <- paste0("se_", colnames(se_hat))
  colnames(W_prim) <- paste0("W_", colnames(W))
  colnames(p_prim) <- paste0("p_", colnames(p_hat))
  colnames(q_prim) <- paste0("q_", colnames(p_hat))
  colnames(diff_prim) <- paste0("diff_", colnames(p_hat))
  if (is.null(tax_name)) {
    tax_name <- pay$taxon_names
  }
  res <- do.call("cbind", list(data.frame(taxon = tax_name), beta_prim, se_prim,
                               W_prim, p_prim, q_prim, diff_prim))
  rownames(res) <- NULL

  # `bias_correct_log_table` is a data.frame in the original's return value and a
  # matrix on the wire: `as.matrix` on the way out (bridge.R), `data.frame` here.
  ybc <- pl$y_bias_crt
  # `data.frame` on a matrix takes both margins from its `dimnames`, and the
  # transport restored them, so this needs no explicit naming. The original's
  # `bias_correct_log_table` is indexed by taxon, not by row number, so unlike
  # `res` it keeps its rownames -- `res` is the one the original sets to NULL at
  # line 266.
  bias_correct_log_table <- data.frame(ybc, check.names = FALSE)

  # A re-assertion: `rb_restore_names` has already set both margins. Kept because
  # this function should not depend on that having happened, and because it costs
  # nothing -- assigning a name that is already there is not a silent change.
  O2 <- pl$O2
  rownames(O2) <- pay$taxon_names
  colnames(O2) <- pay$sample_names

  # Ancombc2.R lines 783-805. The four optional tables are `NULL` when their
  # switch is off, and `NULL` is a value the comparison sees, so it has to be the
  # same `NULL` and not an absent element.
  unsupported <- c(global = global, pairwise = pairwise, dunnet = dunnet, trend = trend)
  if (any(unsupported)) {
    stop("the global, pairwise, Dunnett and trend assemblies are not transported ",
         "yet; requested: ",
         paste(names(unsupported)[unsupported], collapse = ", "),
         ". This is a missing transport, not a missing result.")
  }

  list(
    feature_table = O2,
    bias_correct_log_table = bias_correct_log_table,
    samp_frac = pl$theta_hat,
    delta_em = pl$delta_em,
    delta_wls = pl$delta_wls,
    res = res,
    res_global = NULL,
    res_pair = NULL,
    res_dunn = NULL,
    res_trend = NULL
  )
}
