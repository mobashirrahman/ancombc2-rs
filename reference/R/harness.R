# Instrumented reference runner for ancombc2-rs.
#
#   harness <- load_harness()
#   g <- harness$ref_run(counts, meta, fix_formula = "group", ...)
#
# `ref_run` reproduces ANCOMBC 2.15.2's fixed-effects path. It mirrors
# `.ancombc2_core` statement by statement so that intermediate quantities the
# package does not return can be captured, and it asserts at every call site that
# the instrumented mirror agrees with the oracle's own functions. The final
# tables are additionally cross-checked against an unmodified call to
# `.ancombc2_core`; see `ref_verify_mirror()`.
#
# Sources the pinned oracle's R files rather than installing the package; see
# reference/env/ORACLE.md for why.

load_harness <- function(oracle_dir = NULL) {
  if (is.null(oracle_dir)) {
    # `Sys.getenv`'s `unset` must be a character; `NULL` is a type error, so the
    # documented `load_harness()` call -- with the directory coming from the
    # environment -- failed outright with "wrong type for argument". Every script
    # here happens to pass `oracle_dir` explicitly, which is why it went unnoticed.
    oracle_dir <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "")
    if (!nzchar(oracle_dir)) {
      stop("set ANCOMBC_ORACLE_DIR or pass oracle_dir=")
    }
  }
  here <- environment()
  # Where the harness's own R files live. By default this is deduced from the
  # oracle's location, assuming the repository layout (`<root>/reference/ANCOMBC`
  # beside `<root>/reference/R`). That assumption breaks when the oracle is
  # checked out elsewhere -- notably `/opt/ANCOMBC` in the benchmark image, for
  # which `..` resolves to `/opt/R` and the `source()` fails with "No such file
  # or directory", which is how the container's R arms failed twice with nothing
  # in the log. `ANCOMBC_REFERENCE_R` overrides it explicitly.
  ref_dir <- Sys.getenv("ANCOMBC_REFERENCE_R", unset = "")
  if (!nzchar(ref_dir)) {
    ref_dir <- file.path(oracle_dir, "..", "R")
  }
  source(file.path(ref_dir, "oracle.R"), local = TRUE)
  source(file.path(ref_dir, "stubs.R"), local = TRUE)

  verify_oracle(oracle_dir)
  # foreach / %dorng% / registerDoSEQ are resolved in globalenv() by the sourced
  # oracle code, so they must live there.
  assign("foreach", foreach_stub, envir = globalenv())
  assign("%dorng%", `%dorng%`, envir = globalenv())
  assign("registerDoSEQ", registerDoSEQ, envir = globalenv())

  O <- new.env(parent = globalenv())
  for (f in c("utils.R", "ancombc_prep.R", "ancombc_bias_correct.R",
              "ancombc_mult.R")) {
    sys.source(file.path(oracle_dir, "R", f), envir = O)
  }

  for (fn in c(".data_core", ".get_struc_zero", ".lm_fit_all",
               ".sandwich_vcov", ".iter_mle", ".bias_em", ".ancombc2_core",
               ".ancombc_global_F", ".ancombc_pair", ".mdfdr", ".combn_fun",
               ".combn_fun2", ".var_diff")) {
    if (!exists(fn, envir = O, inherits = FALSE)) {
      stop("oracle function missing: ", fn)
    }
  }

  list(
    env = O,
    ref_run = function(...) ref_run(..., oracle = O),
    ref_verify_mirror = function(...) ref_verify_mirror(..., oracle = O)
  )
}

# ---------------------------------------------------------------------------
# data_sanity_check, restricted to the matrix / data.frame input path
# ---------------------------------------------------------------------------

# R's `terms()` accepts whitespace as a term separator as well as `+`, and a
# formula stored in a file is often written either way. `formula()` itself does
# *not* accept the whitespace form, so the separator is normalised once, here,
# rather than in every caller.
#
# The class is `[-+*]` plus a *POSIX* `[[:space:]]`, deliberately not `\s`.
# R's default engine is TRE, where a backslash inside a bracket expression is an
# ordinary literal -- so the class `[-+*\s]` is the five characters `-`, `+`,
# `*`, `\` and **`s`**, and every `s` in every variable name is silently deleted.
# `sex` becomes `ex` and the analysis then fails with "variables not in
# metadata", which reads like a metadata problem rather than a regex one. It is
# invisible for `group`, `x1`, `cont_cov` and every other name the fixtures and
# the benchmark datasets happen to use, which is exactly why it survived until a
# real dataset with an `s` in a variable name reached it.
# Expand `a * b` into `a + b + a:b` the way `terms()` does, before the operators
# are stripped.
#
# The strip below removes `*` and keeps only the variable names, which silently
# deleted every interaction: `group + x4 * x1` normalised to
# `group + x1 + x2 + x3 + x4 + x1`, and the design matrix had no interaction
# column at all. A user asking for `*` got main effects only, with no warning --
# and the two are different models, so results were reported for a model that was
# not the one requested.
#
# Two details are load-bearing. `*` binds tighter than `+`, so the expansion is
# per `+`-separated term. And `terms()` writes an interaction label with its
# variables in sorted order, so `x10 * x1` is `x1:x10` and not `x10:x1`; the
# column *name* is part of the Level A contract, so getting this backwards would
# be a parity failure rather than a cosmetic difference.
expand_interactions <- function(fix_formula) {
  parts <- strsplit(fix_formula, split = "+", fixed = TRUE)[[1]]
  out <- character(0)
  for (part in parts) {
    if (!grepl("*", part, fixed = TRUE)) {
      out <- c(out, part)
      next
    }
    f <- trimws(strsplit(part, split = "*", fixed = TRUE)[[1]])
    f <- f[nzchar(f)]
    if (length(f) < 2L) {
      out <- c(out, part)
      next
    }
    out <- c(out, f)
    if (length(f) == 2L) {
      out <- c(out, paste(sort(f), collapse = ":"))
    } else {
      # `terms()` orders the interactions of three or more factors by the
      # increasing length of the term, and each label is sorted.
      k <- 2L
      while (k <= length(f)) {
        cmb <- utils::combn(f, k)
        for (j in seq_len(ncol(cmb))) {
          out <- c(out, paste(sort(cmb[, j]), collapse = ":"))
        }
        k <- k + 1L
      }
    }
  }
  paste(out, collapse = " + ")
}

normalise_fix_formula <- function(fix_formula) {
  if (is.null(fix_formula)) return(NULL)
  fix_formula <- expand_interactions(fix_formula)
  # Two passes rather than one character class. A single class mixing the
  # operators and whitespace is easy to get subtly wrong -- see the note above
  # about `\s` -- and `[*+-]` is unambiguous only because `*`, `+`, `,` and `-`
  # happen to be contiguous in ASCII, which is not a property anyone should have
  # to know to read this line.
  s <- sub("~", " ", fix_formula)
  s <- gsub("[*+-]", " ", s)
  s <- trimws(gsub("[[:space:]]+", " ", s))
  parts <- strsplit(s, split = "[[:space:]]+")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) return(NULL)
  # `terms()` keeps each term once, at its first occurrence, so the expansion of
  # `a * b` -- which repeats `a` and `b` -- must not leave the repeats behind.
  parts <- parts[!duplicated(parts)]
  paste(parts, collapse = " + ")
}

ref_sanity_check <- function(counts, meta_data, fix_formula, group = NULL,
                             struc_zero = FALSE, global = FALSE,
                             pairwise = FALSE) {
  if (!is.null(group) && any(c(global, pairwise))) {
    if (!group %in% colnames(meta_data)) {
      stop("group variable not in metadata: ", group)
    }
  }
  if (!is.null(fix_formula)) {
    # R's `terms()` accepts whitespace as a term separator as well as `+`, and
    # `*` as an interaction, so the variable list has to be built the same way or
    # a formula written as `group x1 x2` -- which is what the benchmark data
    # generator writes, matching the style of an R formula in a file -- would be
    # reported as one unknown variable.
    #
    # This used to re-implement the split with the same `[-+*\s]` character
    # class as `normalise_fix_formula`, and therefore had the same bug: the class
    # is `-`, `+`, `*`, `\` and `s` in TRE, so `sex` arrived here as `ex` and the
    # check reported a variable that was not in the metadata. Two copies of a
    # regex is two places for it to be wrong, so this now calls the one function
    # that has the tested behaviour.
    #
    # `normalise_fix_formula` returns the terms *joined* with " + ", so the
    # joined string is split again here. `strsplit` on a fixed separator and not
    # on whitespace, because a variable name may itself contain no space but a
    # term must not be re-split into characters.
    vars <- trimws(strsplit(normalise_fix_formula(fix_formula), split = " + ",
                            fixed = TRUE)[[1]])
    vars <- vars[nzchar(vars)]
    # An interaction term is `a:b`, and neither `a:b` nor the string is a column
    # of the metadata -- `model.matrix` builds it from `a` and `b`. So the check
    # is over the *variables* each term is made of, not over the terms.
    # Otherwise the first interaction a fixture ever declared was rejected as
    # "variables not in metadata: x1:x10", which is a false alarm about a term
    # the design matrix does construct.
    vars <- unique(unlist(strsplit(vars, split = ":", fixed = TRUE)))
    vars <- vars[nzchar(vars)]
    miss <- setdiff(vars, colnames(meta_data))
    if (length(miss) > 0) stop("variables not in metadata: ",
                              paste(miss, collapse = ", "))
  }
  if (!all(colnames(counts) %in% rownames(meta_data))) {
    stop("sample names do not match between feature table and metadata")
  }
  meta_data <- meta_data[colnames(counts), , drop = FALSE]
  meta_data[] <- lapply(meta_data, function(x) if (is.factor(x)) factor(x) else x)

  if (is.null(group)) {
    if (any(c(global, pairwise))) stop("group required for multi-group comparison")
    if (struc_zero) stop("group required for structural zeros")
  } else {
    if (!is.numeric(meta_data[, group])) {
      meta_data[, group] <- as.factor(meta_data[, group])
      n_level <- nlevels(meta_data[, group])
      if (n_level < 2) stop("group must have >= 2 categories")
      if (n_level < 3) {
        global <- FALSE
        pairwise <- FALSE
      }
      size_per_group <- table(meta_data[, group])
      if (any(size_per_group < 2)) {
        stop("sample size per group must be >= 2; small groups: ",
             paste(names(size_per_group)[which(size_per_group < 2)],
                   collapse = ", "))
      }
    } else {
      warning("group variable contains numerical values")
    }
  }
  list(feature_table = counts, feature_table_aggregate = counts,
       meta_data = meta_data, global = global, pairwise = pairwise)
}

# ---------------------------------------------------------------------------
# Instrumented mirror of .ancombc2_core (fixed-effects branch only)
# ---------------------------------------------------------------------------

ref_core <- function(data, aggregate_data, meta_data, fix_formula,
                     p_adj_method = "holm", pseudo = 0, s0_perc = 0.05,
                     group = NULL, alpha = 0.05, global = FALSE,
                     pairwise = FALSE,
                     iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                     oracle = NULL) {
  capture <- new.env(parent = emptyenv())

  O1 <- data + pseudo
  O2 <- aggregate_data + pseudo
  n_tax <- nrow(O2)
  tax_name <- rownames(O2)

  old_opt <- options(na.action = "na.pass")
  on.exit(options(old_opt), add = TRUE)
  x <- stats::model.matrix(stats::as.formula(paste0("~", fix_formula)),
                           data = meta_data)
  options(old_opt)
  fix_eff <- colnames(x)
  n_fix_eff <- length(fix_eff)

  o1 <- log(O1)
  o1[is.infinite(o1)] <- NA
  y1 <- o1 - rowMeans(o1, na.rm = TRUE)

  capture$y1 <- y1
  capture$x <- x
  capture$fix_eff <- fix_eff
  capture$O1 <- O1
  capture$O2 <- O2

  # ---- stage 1: iterative MLE with theta estimated from zero -------------
  #
  # `verbose = TRUE` regardless of what the caller asked for, because the
  # per-iteration epsilon is the contract's "convergence trace" and `.iter_mle`
  # does not return it -- it only *prints* it, one `message()` per iteration:
  #
  #   ML iteration = 1, epsilon = 1.2e+01
  #
  # The messages are intercepted and muffled, so the trace is captured as data
  # without the run writing to stderr. This is the same reason `ref_core` exists
  # rather than calling the package's `.ancombc2`: the quantities the contract
  # names are not all in what the package returns, and ORACLE.md records that the
  # harness sources the pinned R files precisely so the `.`-prefixed internals
  # can be instrumented.
  #
  # `ref_run` pins `verbose = FALSE` on the refits so a 50-run sensitivity sweep
  # does not print 50 traces; here it is forced on for the main run only.
  iter_trace <- list()
  collect_trace <- function(m) {
    txt <- conditionMessage(m)
    if (!grepl("^ML iteration = ", txt)) return(invisible(NULL))
    it <- as.integer(sub("^ML iteration = ([0-9]+).*$", "\\1", txt))
    ep <- as.numeric(sub(".*epsilon = ", "", txt))
    iter_trace[[length(iter_trace) + 1L]] <<-
      list(iteration = it, epsilon = ep)
    invokeRestart("muffleMessage")
  }
  para1 <- withCallingHandlers(
    oracle$.iter_mle(x = x, y = y1, meta_data = meta_data,
                     formula = fix_formula, theta = NULL,
                     tol = iter_control$tol,
                     max_iter = iter_control$max_iter,
                     verbose = TRUE),
    message = function(m) collect_trace(m))
  beta1 <- para1$beta
  var_hat1 <- para1$var_hat
  capture$beta_star <- beta1
  capture$theta <- para1$theta
  capture$vcov1 <- para1$vcov_hat
  capture$var1 <- var_hat1
  capture$dof1 <- para1$dof
  # The convergence trace: epsilon at each iteration, in order, plus the tolerance
  # and cap it was run under so the trace is interpretable on its own.
  capture$convergence_trace <- list(
    epsilons = vapply(iter_trace, function(z) z$epsilon, numeric(1)),
    iterations = vapply(iter_trace, function(z) z$iteration, numeric(1)),
    tol = iter_control$tol,
    max_iter = iter_control$max_iter
  )

  # ---- stage 2: E-M bias estimation --------------------------------------
  # The mixture parameters are the contract's "EM mixture parameters".
  #
  # `.bias_em` computes the three-component Gaussian mixture's weights and
  # component offsets/variances on every iteration but returns only
  # `c(delta_em, delta_wls, var_delta)`, so the mixture is not in what the
  # package hands back.
  #
  # These are captured by tracing the *exit* of the oracle's own `.bias_em` and
  # reading the loop's final `*_new` values -- the values it actually used, not a
  # re-derivation. Re-implementing the routine here instead would be a ~150-line
  # transcription of a Nelder-Mead EM fit into the one number the whole method
  # rests on, where a transcription slip would be a silent parity risk that no
  # assertion could see. Reading them off the function that produced them cannot
  # drift from the oracle.
  #
  # The trace is installed on the harness's own `oracle` env and removed in the
  # `on.exit` below, so a repeated run in the same session cannot accumulate
  # handlers or leak the environment it writes into.
  em_mix <- new.env(parent = emptyenv())
  em_mix$calls <- list()
  # Assigned into the oracle env, which is `.bias_em`'s enclosing environment, so
  # the traced code finds it by ordinary scoping.
  assign("ancombc2rs_em_mix", em_mix, envir = oracle)
  on.exit({
    suppressMessages(try(
      untrace(".bias_em", where = oracle), silent = TRUE))
    if (exists("ancombc2rs_em_mix", envir = oracle, inherits = FALSE)) {
      rm("ancombc2rs_em_mix", envir = oracle)
    }
  }, add = TRUE)
  suppressMessages(trace(
    what = ".bias_em", where = oracle, print = FALSE,
    exit = quote({
      ancombc2rs_em_mix$calls <- c(
        ancombc2rs_em_mix$calls,
        list(list(pi = c(pi0_new, pi1_new, pi2_new),
                  delta = delta_new, l1 = l1_new, l2 = l2_new,
                  kappa1 = kappa1_new, kappa2 = kappa2_new,
                  iterations = iterNum,
                  # The E-M's own final step size: `epsilon` is the Euclidean norm
                  # of the parameter change on the last iteration. It is recorded
                  # because it is the only measure of *how far from converged* a
                  # fit that stopped at `max_iter` actually is, and therefore the
                  # only defensible basis for a tolerance on that fit's parameters.
                  epsilon = epsilon,
                  tol = tol, max_iter = max_iter)))
    })))

  n_em_before <- length(em_mix$calls)
  bias1 <- foreach(i = seq_len(ncol(beta1)), .combine = rbind) %dorng% {
    oracle$.bias_em(beta = beta1[, i], var_hat = var_hat1[, i],
                    tol = em_control$tol, max_iter = em_control$max_iter)
  }
  bias1 <- as.data.frame(bias1, row.names = fix_eff, check.names = FALSE)
  colnames(bias1) <- c("delta_em", "delta_wls", "var_delta")
  delta_em <- bias1$delta_em
  delta_wls <- bias1$delta_wls
  var_delta <- bias1$var_delta
  capture$delta_em <- delta_em
  capture$delta_wls <- delta_wls
  capture$var_delta <- var_delta

  # The EM mixture parameters, one entry per fixed effect, in the order
  # `.bias_em` was called.
  #
  # Only the calls made by *this* stage are kept. A later `.bias_em` call would
  # otherwise be appended to the same list -- the non-conservative path refits
  # bias for each of its 50 sensitivity pseudo-counts, and a naive read of the
  # list would attribute one of those refits to the main run.
  em_calls <- em_mix$calls[seq.int(n_em_before + 1L, length(em_mix$calls))]
  if (length(em_calls) != ncol(beta1)) {
    stop("expected ", ncol(beta1), " .bias_em trace(s) for the main run, got ",
         length(em_calls))
  }
  # The traced `delta` is `.bias_em`'s own `delta_new`, i.e. the `delta_em` it
  # returned. Asserting that ties the captured mixture to the already-golden
  # `delta_em`, so a trace that silently captured a different call cannot pass as
  # this one.
  if (!isTRUE(all.equal(vapply(em_calls, function(z) z$delta, numeric(1)),
                        as.numeric(delta_em), tolerance = 1e-12))) {
    stop("the captured EM mixture does not correspond to delta_em")
  }
  capture$em_mixture <- list(
    terms = fix_eff,
    pi = do.call(rbind, lapply(em_calls, function(z) z$pi)),
    delta = vapply(em_calls, function(z) z$delta, numeric(1)),
    l = do.call(rbind, lapply(em_calls, function(z) c(z$l1, z$l2))),
    kappa = do.call(rbind, lapply(em_calls, function(z) c(z$kappa1, z$kappa2))),
    iterations = vapply(em_calls, function(z) z$iterations, numeric(1)),
    epsilon = vapply(em_calls, function(z) z$epsilon, numeric(1)),
    tol = em_control$tol,
    max_iter = em_control$max_iter
  )
  # Remove the trace now that the main run's parameters are captured, so the
  # sensitivity refits below cannot add to them even accidentally.
  # `untrace()` takes no `print` argument; passing one errors with "unused
  # argument", and an earlier version of this line wrapped the call in
  # `try(..., silent = TRUE)`, so the error was swallowed and the trace stayed
  # installed -- surfacing much later as "object 'ancombc2rs_em_mix' not found"
  # from an unrelated `.bias_em` call. Hence no `try`: if the trace cannot be
  # removed, this should fail here rather than somewhere downstream.
  suppressMessages(untrace(".bias_em", where = oracle))

  # ---- stage 3: bias correction and sampling fractions --------------------
  beta1 <- t(t(beta1) - delta_em)
  n_samp <- nrow(x)
  theta_hat <- matrix(NA, nrow = nrow(y1), ncol = ncol(y1))
  for (i in seq_len(nrow(y1))) {
    theta_hat[i, ] <- y1[i, ] -
      base::rowSums(x * rep(beta1[i, ], each = n_samp), na.rm = TRUE)
  }
  theta_hat <- colMeans(theta_hat, na.rm = TRUE)
  names(theta_hat) <- colnames(y1)
  capture$beta_corr_stage1 <- beta1
  capture$samp_frac <- theta_hat

  # ---- stage 4: final MLE with theta fixed -------------------------------
  o2 <- log(O2)
  o2[is.infinite(o2)] <- NA
  y2 <- o2 - rowMeans(o2, na.rm = TRUE)
  # The result carries a data.frame (R's `t(t(y) - theta)`), but the golden
  # contract needs a matrix so its shape is explicit in the manifest.
  y_bias_crt <- as.matrix(t(t(y2) - theta_hat))
  capture$y2 <- y2
  capture$y_bias_crt <- y_bias_crt


  # The raw capture map is handed out unresolved; `ref_run` knows `nrow(O1)` and
  # `nrow(O2)`, which are the authoritative taxon counts of the two sets (the
  # tables the MLEs are actually fitted on), and resolves it there.

  para2 <- oracle$.iter_mle(x = x, y = y2, meta_data = meta_data,
                            formula = fix_formula, theta = theta_hat,
                            tol = iter_control$tol,
                            max_iter = iter_control$max_iter,
                            verbose = iter_control$verbose)
  beta_hat <- para2$beta
  var_hat <- para2$var_hat
  dof <- para2$dof
  capture$beta <- beta_hat
  capture$var_hat <- var_hat
  capture$dof <- dof
  capture$vcov2 <- para2$vcov_hat

  # The missingness pattern assignment for the reported set: which taxa share a
  # usable-sample pattern, and in what order the patterns are numbered.
  #
  # This is the contract's "missingness pattern assignment" and a Level A quantity
  # -- the assignment decides which taxa are solved together by one QR, so a
  # different assignment is a different factorisation and a different `beta` for
  # every taxon in the affected patterns.
  #
  # An earlier version traced the oracle's `.lm_fit_all` and read `groups` off its
  # exit. That was abandoned: `trace()` re-binds the traced function, and the
  # traced call did not reliably correspond to the stage being measured -- a
  # direct probe of `.iter_mle` reported its own `theta` argument as non-`NULL`
  # even on the call that was passed `NULL`, so the exit expression was reading a
  # different binding than the one the stage ran with. The assignment is derived
  # here instead, from the rule `.lm_fit_all` states in its own source:
  #
  #   use   = is.finite(Ymat) & matrix(complete.cases(x), byrow = TRUE)
  #   keys  = do.call(paste0, asplit(use * 1L, 2L))
  #   groups = split(seq_len(n_tax), factor(keys, levels = unique(keys)))
  #
  # -- and then **verified against two observables the oracle returns**, so this is
  # not a transcription taken on trust:
  #
  #   1. `para2$dof` is `n_used - p` for every taxon that was fitted, so it says
  #      each taxon's pattern size directly. An assignment that disagreed would
  #      disagree with the oracle's own degrees of freedom.
  #   2. `.lm_fit_all` writes a literal `0` into `fitted` at every sample its
  #      group's fit did *not* use, so the zero pattern of `para2$fitted` is the
  #      pattern itself, seen from the other side.
  n_tax <- nrow(y2)
  n_p <- ncol(x)
  x_ok <- stats::complete.cases(x)
  y2_ok <- is.finite(y2)
  if (all(x_ok) && all(y2_ok)) {
    use_mat <- matrix(TRUE, nrow = n_tax, ncol = ncol(y2))
  } else {
    use_mat <- y2_ok & matrix(x_ok, nrow = n_tax, ncol = ncol(y2),
                              byrow = TRUE)
  }
  pat_keys <- do.call(paste0, asplit(use_mat * 1L, 2L))
  pat_groups <- split(seq_len(n_tax), factor(pat_keys, levels = unique(pat_keys)))
  pat_gid <- integer(n_tax)
  for (k in seq_along(pat_groups)) pat_gid[pat_groups[[k]]] <- k
  pat_n_used <- as.integer(rowSums(use_mat))

  # Verification 1: the oracle's own residual degrees of freedom.
  # `.iter_mle` returns `dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff,
  # byrow = FALSE)` -- the per-taxon vector replicated once per coefficient. Taking
  # `as.numeric()` of that yields `n_tax * p` values, so a flat read silently
  # recycles the taxon index and produces NAs rather than a mismatch. One column is
  # the per-taxon vector.
  dof_mat <- para2$dof
  if (is.matrix(dof_mat)) {
    if (nrow(dof_mat) != n_tax || ncol(dof_mat) != n_p) {
      stop("dof is ", nrow(dof_mat), " x ", ncol(dof_mat), ", expected ",
           n_tax, " x ", n_p)
    }
    oracle_dof <- as.numeric(dof_mat[, 1])
  } else {
    oracle_dof <- as.numeric(dof_mat)
  }
  if (length(oracle_dof) != n_tax) {
    stop("per-taxon dof has ", length(oracle_dof), " entries for ", n_tax, " taxa")
  }
  # `NA` counts as unfittable too: `.lm_fit_all` leaves a taxon's `dof` at 999
  # when its per-taxon `lm` fails, and `NA` when the fit produced no degrees of
  # freedom at all. Indexing with an `NA` mask would silently drop taxa from the
  # check below, which is the opposite of what a verification should do.
  unfittable <- is.na(oracle_dof) | oracle_dof >= 999
  fitted_ok <- !unfittable
  if (any(fitted_ok)) {
    # `.lm_fit_all` reports `dof` as `nrow(xr) - rank`, where `rank` is the rank of
    # the group's *sub*-design, not `p`. For a full-rank group those agree; for a
    # rank-deficient one they do not, and assuming `rank == p` makes the check fail
    # on exactly the taxa that are hardest to get right -- which is how `fx04`'s
    # 200 rank-deficient taxa were the first thing to disagree. So the rank is
    # computed per pattern, as the oracle computes it.
    pat_rank <- integer(length(pat_groups))
    for (k in seq_along(pat_groups)) {
      t0 <- pat_groups[[k]][1L]
      xr <- x[which(use_mat[t0, ]), , drop = FALSE]
      pat_rank[k] <- if (nrow(xr) == 0L) 0L else qr(xr)$rank
    }
    implied <- numeric(n_tax)
    implied[fitted_ok] <- pat_n_used[fitted_ok] - pat_rank[pat_gid[fitted_ok]]
    if (!identical(implied[fitted_ok], oracle_dof[fitted_ok])) {
      bad <- which(fitted_ok & implied != oracle_dof)[1L]
      stop(paste0("the derived pattern assignment disagrees with the oracle's ",
                  "own dof at taxon ", bad, ": n_used = ", pat_n_used[bad],
                  ", pattern rank = ", pat_rank[pat_gid[bad]], ", so dof should be ",
                  implied[bad], ", and the oracle reports ", oracle_dof[bad],
                  ". The grouping is not what was fitted."))
    }
  }

  # Verification 2: the literal zeros `.lm_fit_all` writes outside each group's
  # samples. This is a subset relation rather than equality -- a fitted value can
  # coincide with zero -- so it is checked as "every sample the group did not use
  # is a literal zero", which is the property that matters and cannot hold for a
  # wrong assignment.
  fitted_mat <- para2$fitted
  for (t in seq_len(n_tax)) {
    if (unfittable[t]) next
    unused <- which(!use_mat[t, ])
    if (length(unused) == 0) next
    col <- fitted_mat[t, unused]
    if (!all(col == 0)) {
      stop("taxon ", t, ": samples the pattern says were unused (",
           paste(unused, collapse = ","), ") are not literal zeros in the ",
           "oracle's fitted values, so the derived grouping is not the one the ",
           "oracle fitted")
    }
  }

  capture$pattern_assignment <- list(
    reported_set_n_taxa = n_tax,
    n_p = n_p,
    n_groups = length(pat_groups),
    group = as.integer(pat_gid),
    key = pat_keys,
    n_used = pat_n_used,
    unfittable = as.logical(unfittable)
  )


  # ---- stage 5: variance of delta, regularization, SEs -------------------
  var_hat <- sweep(var_hat, 2, var_delta, "+") +
    2 * sqrt(sweep(var_hat, 2, var_delta, "*"))
  if (is.null(s0_perc)) {
    s02 <- 0
  } else {
    s02 <- apply(var_hat, 2, function(x) stats::quantile(x, s0_perc, na.rm = TRUE))
  }
  var_hat <- t(t(var_hat) + s02)
  var_hat[is.na(beta_hat)] <- NA
  se_hat <- sqrt(var_hat)
  vcov_hat <- lapply(seq_len(n_tax), function(i) {
    diag(para2$vcov_hat[[i]]) <- var_hat[i, ]
    para2$vcov_hat[[i]]
  })
  capture$s02 <- s02
  capture$var_final <- var_hat
  capture$se <- se_hat
  capture$vcov <- vcov_hat

  # ---- stage 6: primary inference ----------------------------------------
  W <- beta_hat / se_hat
  p_hat <- 2 * stats::pt(abs(W), df = dof, lower.tail = FALSE)
  p_hat[is.na(p_hat)] <- 1
  q_hat <- apply(p_hat, 2, function(x) stats::p.adjust(x, method = p_adj_method))
  diff_abn <- q_hat <= alpha & !is.na(q_hat)
  capture$W <- W
  capture$p <- p_hat
  capture$q <- q_hat
  capture$diff_abn <- diff_abn

  beta_prim <- as.data.frame(beta_hat, check.names = FALSE)
  se_prim <- as.data.frame(se_hat, check.names = FALSE)
  W_prim <- as.data.frame(W, check.names = FALSE)
  p_prim <- as.data.frame(p_hat, check.names = FALSE)
  q_prim <- as.data.frame(q_hat, check.names = FALSE)
  diff_prim <- as.data.frame(diff_abn, check.names = FALSE)
  colnames(beta_prim) <- paste0("lfc_", colnames(beta_hat))
  colnames(se_prim) <- paste0("se_", colnames(se_hat))
  colnames(W_prim) <- paste0("W_", colnames(W))
  colnames(p_prim) <- paste0("p_", colnames(p_hat))
  colnames(q_prim) <- paste0("q_", colnames(q_hat))
  colnames(diff_prim) <- paste0("diff_", colnames(diff_abn))
  res <- cbind(data.frame(taxon = tax_name), beta_prim, se_prim, W_prim,
               p_prim, q_prim, diff_prim)
  rownames(res) <- NULL

  res_global <- NULL
  if (global) {
    res_global <- oracle$.ancombc_global_F(x = x, group = group,
                                          beta_hat = beta_hat,
                                          vcov_hat = vcov_hat, dof = dof,
                                          p_adj_method = p_adj_method,
                                          alpha = alpha)
    rownames(res_global) <- NULL
  }

  res_pair <- NULL
  if (pairwise) {
    rp <- oracle$.ancombc_pair(x = x, group = group, beta_hat = beta_hat,
                              var_hat = var_hat, vcov_hat = vcov_hat, dof = dof,
                              fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                              alpha = alpha, full_model = NULL,
                              fix_formula = fix_formula, rand_formula = NULL,
                              control = NULL, y = y_bias_crt,
                              meta_data = meta_data)
    beta_pair <- as.data.frame(rp$beta, check.names = FALSE)
    se_pair <- as.data.frame(rp$se, check.names = FALSE)
    W_pair <- as.data.frame(rp$W, check.names = FALSE)
    p_pair <- as.data.frame(rp$p_val, check.names = FALSE)
    q_pair <- as.data.frame(rp$q_val, check.names = FALSE)
    diff_pair <- as.data.frame(rp$diff_abn, check.names = FALSE)
    colnames(beta_pair) <- paste0("lfc_", colnames(beta_pair))
    colnames(se_pair) <- paste0("se_", colnames(se_pair))
    colnames(W_pair) <- paste0("W_", colnames(W_pair))
    colnames(p_pair) <- paste0("p_", colnames(p_pair))
    colnames(q_pair) <- paste0("q_", colnames(q_pair))
    colnames(diff_pair) <- paste0("diff_", colnames(diff_pair))
    res_pair <- cbind(data.frame(taxon = tax_name), beta_pair, se_pair, W_pair,
                      p_pair, q_pair, diff_pair)
    rownames(res_pair) <- NULL
    capture$pair_raw <- rp
  }

  out <- list(feature_table = O2, bias_correct_log_table = y_bias_crt,
              samp_frac = theta_hat, delta_em = delta_em,
              delta_wls = delta_wls, res = res, res_global = res_global,
              res_pair = res_pair, fix_eff = fix_eff, x = x, meta_data = meta_data)
  attr(out, "capture") <- capture
  out
}

# ---------------------------------------------------------------------------
# Top-level run: preprocessing + core (+ sensitivity analysis)
# ---------------------------------------------------------------------------

ref_run <- function(counts, meta_data, fix_formula, oracle = NULL,
                    p_adj_method = "holm", pseudo = 0, pseudo_sens = FALSE,
                    conservative = TRUE, prv_cut = 0.10, lib_cut = 0,
                    s0_perc = 0.05, group = NULL, struc_zero = FALSE,
                    neg_lb = FALSE, alpha = 0.05, global = FALSE,
                    pairwise = FALSE,
                    iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
                    em_control = list(tol = 1e-5, max_iter = 100),
                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)) {
  fix_formula <- normalise_fix_formula(fix_formula)

  # Per-stage wall time, for the golden contract's "per-stage timings" quantity.
  #
  # These are *recorded, not compared*: wall-clock is not reproducible and a
  # parity test that asserted on it would be asserting that the oracle ran at the
  # same speed, which is not a property of the contract. What the Rust side
  # asserts is that both sides name the same stages and that the timings are
  # finite and non-negative -- so a stage that silently stopped being timed is
  # caught, while a faster or slower machine is not.
  #
  # The stage names match `StageTimings` in `crates/ancombc2-core/src/pipeline.rs`
  # so the two sides can be listed against each other.
  stage_seconds <- list()
  stage_mark <- function(name) {
    list(name = name, t0 = proc.time()[["elapsed"]])
  }
  stage_done <- function(mark) {
    stage_seconds[[mark$name]] <<- proc.time()[["elapsed"]] - mark$t0
  }

  t0 <- stage_mark("sanity_check")
  check <- ref_sanity_check(counts, meta_data, fix_formula = fix_formula,
                            group = group, struc_zero = struc_zero,
                            global = global, pairwise = pairwise)
  stage_done(t0)
  ft <- check$feature_table
  fta <- check$feature_table_aggregate
  meta <- check$meta_data
  global <- check$global
  pairwise <- check$pairwise

  zero_ind <- NULL
  t0 <- stage_mark("structural_zeros")
  if (struc_zero) {
    zero_ind <- oracle$.get_struc_zero(data = fta, meta_data = meta,
                                       group = group, neg_lb = neg_lb)
    tax_keep <- which(apply(zero_ind[, -1], 1, function(x) all(x == FALSE)))
  } else {
    tax_keep <- seq(nrow(fta))
  }

  stage_done(t0)

  t0 <- stage_mark("preprocess")
  core1 <- oracle$.data_core(data = ft, meta_data = meta, prv_cut = prv_cut,
                              lib_cut = lib_cut, tax_keep = NULL,
                              samp_keep = NULL)
  O1 <- core1$feature_table
  samp_keep <- colnames(O1)

  core2 <- oracle$.data_core(data = fta, meta_data = meta, prv_cut = prv_cut,
                             lib_cut = lib_cut, tax_keep = tax_keep,
                             samp_keep = samp_keep)
  O2 <- core2$feature_table
  meta <- core2$meta_data
  stage_done(t0)

  t0 <- stage_mark("core")
  main <- ref_core(data = O1, aggregate_data = O2, meta_data = meta,
                   fix_formula = fix_formula, p_adj_method = p_adj_method,
                   pseudo = pseudo, s0_perc = s0_perc, group = group,
                   alpha = alpha, global = global, pairwise = pairwise,
                   iter_control = iter_control, em_control = em_control,
                   mdfdr_control = mdfdr_control, oracle = oracle)
  capture <- attr(main, "capture")

  stage_done(t0)

  golden <- list(
    taxa_retained = rownames(O2),
    samples_retained = colnames(O2),
    O1 = O1, O2 = O2,
    fix_eff = main$fix_eff, x = main$x,
    zero_ind = zero_ind,
    y1 = capture$y1, y2 = capture$y2,
    beta_star = capture$beta_star, theta = capture$theta,
    convergence_trace = capture$convergence_trace,
    pattern_assignment = capture$pattern_assignment,
    em_mixture = capture$em_mixture,
    vcov1 = capture$vcov1, var1 = capture$var1, dof1 = capture$dof1,
    delta_em = capture$delta_em, delta_wls = capture$delta_wls,
    var_delta = capture$var_delta,
    beta_corr_stage1 = capture$beta_corr_stage1,
    samp_frac = capture$samp_frac,
    y_bias_crt = capture$y_bias_crt,
    beta = capture$beta, var_hat = capture$var_hat, dof = capture$dof,
    s02 = capture$s02, var_final = capture$var_final, se = capture$se,
    vcov = capture$vcov,
    W = capture$W, p = capture$p, q = capture$q, diff_abn = capture$diff_abn,
    res = main$res, res_global = main$res_global, res_pair = main$res_pair
  )

  if (pseudo_sens && conservative) {
    pseudo_list <- c(0.1, 0.5, 1)
    ic <- iter_control
    ic$verbose <- FALSE
    t0 <- stage_mark("sensitivity")
    ss_list <- lapply(pseudo_list, function(pc) {
      ref_core(data = O1, aggregate_data = O2, meta_data = meta,
               fix_formula = fix_formula, p_adj_method = p_adj_method,
               pseudo = pc, s0_perc = s0_perc, group = group, alpha = alpha,
               global = global, pairwise = pairwise, iter_control = ic,
               em_control = em_control, mdfdr_control = mdfdr_control,
               oracle = oracle)
    })
    all_runs <- c(list(main), ss_list)
    all_pseudo <- c(0, pseudo_list)

    # The oracle stacks the runs into a taxa x effect x run array and reduces over
    # the first two margins, so each run contributes one observation per cell.
    qcols <- grep("^q_", colnames(all_runs[[1]]$res), value = TRUE)
    qmat <- array(unlist(lapply(all_runs, function(r) as.matrix(r$res[, qcols]))),
                  c(nrow(all_runs[[1]]$res), length(qcols), length(all_runs)))
    ss_prim <- apply(qmat, c(1, 2), function(x) sum(x > alpha) / length(all_pseudo))
    colnames(ss_prim) <- gsub("^q_", "ss_prim_", qcols)
    ss_prim_log <- (ss_prim == 0 | ss_prim == 1)
    colnames(ss_prim_log) <- gsub("ss_prim_", "passed_ss_", colnames(ss_prim))
    res <- cbind(all_runs[[1]]$res, ss_prim_log)
    # The suffix list is taken once, before any diff_robust_ column exists:
    # re-matching "^diff_" inside the loop would also match the columns this loop
    # is creating.
    diff_suffixes <- sub("^diff_", "", grep("^diff_", names(res), value = TRUE))
    for (suffix in diff_suffixes) {
      dc <- paste0("diff_", suffix); pc <- paste0("passed_ss_", suffix)
      if (!dc %in% names(res) || !pc %in% names(res)) {
        stop(sprintf("diff_robust loop: missing column (diff=%s present=%s, passed=%s present=%s); have %s",
                     dc, dc %in% names(res), pc, pc %in% names(res),
                     paste(names(res), collapse = ",")))
      }
      res[[paste0("diff_robust_", suffix)]] <- res[[dc]] & res[[pc]]
    }

    res_global <- NULL
    if (global) {
      qg <- array(unlist(lapply(all_runs, function(r) r$res_global[, "q_val", drop = FALSE])),
                  c(nrow(all_runs[[1]]$res_global), 1L, length(all_runs)))
      ss_g <- apply(qg, c(1, 2), function(x) sum(x > alpha) / length(all_pseudo))
      colnames(ss_g) <- "ss_global"
      log_g <- (ss_g == 0 | ss_g == 1)
      colnames(log_g) <- "passed_ss"
      res_global <- cbind(all_runs[[1]]$res_global, log_g)
      res_global[["diff_robust_abn"]] <-
        res_global[["diff_abn"]] & res_global[["passed_ss"]]
    }

    res_pair <- NULL
    if (pairwise) {
      pcols <- grep("^q_", colnames(all_runs[[1]]$res_pair), value = TRUE)
      qp <- array(unlist(lapply(all_runs,
                                function(r) as.matrix(r$res_pair[, pcols]))),
                  c(nrow(all_runs[[1]]$res_pair), length(pcols), length(all_runs)))
      ss_p <- apply(qp, c(1, 2), function(x) sum(x > alpha) / length(all_pseudo))
      colnames(ss_p) <- gsub("^q_", "ss_pair_", pcols)
      log_p <- (ss_p == 0 | ss_p == 1)
      colnames(log_p) <- gsub("ss_pair_", "passed_ss_", colnames(ss_p))
      res_pair <- cbind(all_runs[[1]]$res_pair, log_p)
      for (suffix in sub("^diff_", "", grep("^diff_", names(res_pair), value = TRUE))) {
        res_pair[[paste0("diff_robust_", suffix)]] <-
          res_pair[[paste0("diff_", suffix)]] & res_pair[[paste0("passed_ss_", suffix)]]
      }
    }

    ss_tab <- data.frame(taxon = rownames(O2), ss_prim, check.names = FALSE)
    if (global) ss_tab$ss_global <- ss_g[, 1]
    if (pairwise) ss_tab <- cbind(ss_tab, ss_p)

    golden$res <- res
    golden$res_global <- res_global
    golden$res_pair <- res_pair
    golden$ss_tab <- ss_tab
    golden$sens_pseudo <- all_pseudo
    golden$sens_runs <- all_runs
    stage_done(t0)
  }

  if (pseudo_sens && !conservative) {
    t0 <- stage_mark("sensitivity")
    golden <- ref_run_nonconservative(golden, O1, O2, meta, fix_formula, main,
                                      p_adj_method, alpha, group, global,
                                      pairwise, oracle)
    stage_done(t0)
  }

  golden$stage_seconds <- stage_seconds
  golden
}

# ---------------------------------------------------------------------------
# Non-conservative sensitivity analysis (pseudo added to zeros of the
# bias-corrected table; 50 OLS refits of the inference step)
# ---------------------------------------------------------------------------

ref_run_nonconservative <- function(golden, O1, O2, meta, fix_formula, main,
                                    p_adj_method, alpha, group, global,
                                    pairwise, oracle) {
  pseudo_list <- seq(0.01, 0.5, 0.01)
  n_tax <- nrow(O2)
  fix_eff <- gsub("^q_", "", grep("^q_", colnames(main$res), value = TRUE))
  if (is.null(group)) {
    group_eff <- NULL
  } else {
    group_eff <- fix_eff[grepl(group, fix_eff) & !grepl(":", fix_eff)]
  }
  samp_frac <- main$samp_frac

  ss_list <- foreach(pseudo_count = pseudo_list) %dorng% {
    O <- as.matrix(O2)
    O[O == 0] <- pseudo_count
    o <- log(O)
    y <- o - rowMeans(o)
    y_bias_crt <- t(t(y) - samp_frac)
    Y <- as.data.frame(t(y_bias_crt), check.names = FALSE)

    p_list <- lapply(Y, function(ycol) {
      df <- data.frame(y = ycol, meta)
      lm_fit <- stats::lm(stats::as.formula(paste0("y ~ ", fix_formula)), data = df)
      coef_tab <- summary(lm_fit)$coefficients
      p_val <- rep(NA_real_, length(fix_eff))
      names(p_val) <- fix_eff
      p_val[rownames(coef_tab)] <- coef_tab[, "Pr(>|t|)"]

      if (pairwise) {
        beta_hat <- stats::coef(lm_fit)[group_eff]
        vcov_fit <- stats::vcov(lm_fit)
        vcov_hat <- matrix(NA_real_, nrow = length(fix_eff), ncol = length(fix_eff),
                           dimnames = list(fix_eff, fix_eff))
        vcov_hat[rownames(vcov_fit), colnames(vcov_fit)] <- vcov_fit
        dof <- lm_fit$df.residual
        combn_mat <- utils::combn(group_eff, 2)
        pair_p <- vapply(seq_len(ncol(combn_mat)), function(i) {
          id1 <- combn_mat[2, i]
          id2 <- combn_mat[1, i]
          beta_diff <- beta_hat[id1] - beta_hat[id2]
          var_diff <- vcov_hat[id1, id1] + vcov_hat[id2, id2] - 2 * vcov_hat[id1, id2]
          2 * stats::pt(abs(beta_diff / sqrt(var_diff)), df = dof, lower.tail = FALSE)
        }, FUN.VALUE = double(1))
        names(pair_p) <- paste(combn_mat[2, ], combn_mat[1, ], sep = "_")
        p_val <- c(p_val, pair_p)
      }
      if (global) {
        anova_tab <- stats::anova(lm_fit)
        p_val <- c(p_val, global = anova_tab$`Pr(>F)`[rownames(anova_tab) == group])
      }
      p_val
    })
    p_hat <- do.call(rbind, p_list)
    p_hat[is.na(p_hat)] <- 1
    p_hat
  }

  ss_tab_fun <- function(sens_col) {
    ss <- vapply(seq_along(sens_col), function(j) {
      p_pseudo <- vapply(ss_list, function(p) p[, sens_col[j]], FUN.VALUE = double(n_tax))
      rowMeans(p_pseudo > alpha)
    }, FUN.VALUE = double(n_tax))
    matrix(ss, nrow = n_tax, ncol = length(sens_col))
  }
  passed_fun <- function(ss_tab, p_main) {
    p_main <- as.matrix(p_main)
    p_main[is.na(p_main)] <- 1
    (ss_tab == 0 & p_main <= alpha) | (ss_tab == 1 & p_main > alpha)
  }
  flag_fun <- function(res_tab, ss_tab, p_main) {
    ss_log <- passed_fun(ss_tab, p_main)
    colnames(ss_log) <- sub("^ss_(prim|pair|dunn)_", "passed_ss_", colnames(ss_tab))
    res_tab <- cbind(res_tab, ss_log)
    for (suffix in sub("^diff_", "", grep("^diff_", names(res_tab), value = TRUE))) {
      res_tab[[paste0("diff_robust_", suffix)]] <-
        res_tab[[paste0("diff_", suffix)]] & res_tab[[paste0("passed_ss_", suffix)]]
    }
    res_tab
  }

  ss_prim <- ss_tab_fun(fix_eff)
  colnames(ss_prim) <- paste0("ss_prim_", fix_eff)
  res <- flag_fun(main$res, ss_prim, main$res[, paste0("p_", fix_eff)])

  res_global <- NULL
  if (global) {
    ss_g <- ss_tab_fun("global")
    colnames(ss_g) <- "ss_global"
    log_g <- passed_fun(ss_g, main$res_global[, "p_val"])
    colnames(log_g) <- "passed_ss"
    res_global <- cbind(main$res_global, log_g)
    res_global[["diff_robust_abn"]] <-
      res_global[["diff_abn"]] & res_global[["passed_ss"]]
  }

  res_pair <- NULL
  if (pairwise) {
    pair_col <- gsub("^q_", "", grep("^q_", colnames(main$res_pair), value = TRUE))
    ss_p <- ss_tab_fun(pair_col)
    colnames(ss_p) <- paste0("ss_pair_", pair_col)
    res_pair <- flag_fun(main$res_pair, ss_p, main$res_pair[, paste0("p_", pair_col)])
  }

  ss_tab <- data.frame(taxon = rownames(O2), ss_prim, check.names = FALSE)
  if (global) ss_tab$ss_global <- ss_g[, 1]
  if (pairwise) ss_tab <- cbind(ss_tab, ss_p)

  golden$res <- res
  golden$res_global <- res_global
  golden$res_pair <- res_pair
  golden$ss_tab <- ss_tab
  golden$sens_pseudo <- pseudo_list
  golden$sens_ss_list <- ss_list
  golden
}

# ---------------------------------------------------------------------------
# Guard: the instrumented mirror must agree with the oracle's own .ancombc2_core
# ---------------------------------------------------------------------------

ref_verify_mirror <- function(counts, meta_data, fix_formula, oracle = NULL,
                              p_adj_method = "holm", alpha = 0.05,
                              prv_cut = 0.10, lib_cut = 0, ...) {
  check <- ref_sanity_check(counts, meta_data, fix_formula = fix_formula)
  O1 <- oracle$.data_core(check$feature_table, check$meta_data,
                          prv_cut = prv_cut, lib_cut = lib_cut,
                          tax_keep = NULL, samp_keep = NULL)$feature_table
  o <- oracle$.ancombc2_core(data = O1, aggregate_data = O1,
                             meta_data = check$meta_data,
                             fix_formula = fix_formula, rand_formula = NULL,
                             p_adj_method = p_adj_method, alpha = alpha,
                             verbose = FALSE, global = FALSE, pairwise = FALSE, ...)
  m <- ref_core(data = O1, aggregate_data = O1, meta_data = check$meta_data,
                fix_formula = fix_formula, p_adj_method = p_adj_method,
                alpha = alpha, oracle = oracle, ...)
  # The mirror stores the bias-corrected table as a matrix where the package
  # returns a data.frame, so the numeric content is compared after coercing.
  same <- isTRUE(all.equal(o$res, m$res, tolerance = 0)) &&
    isTRUE(all.equal(o$samp_frac, m$samp_frac, tolerance = 0)) &&
    isTRUE(all.equal(o$delta_em, m$delta_em, tolerance = 0)) &&
    isTRUE(all.equal(o$delta_wls, m$delta_wls, tolerance = 0)) &&
    isTRUE(all.equal(as.matrix(o$bias_correct_log_table),
                     m$bias_correct_log_table, tolerance = 0))
  if (!same) {
    print(all.equal(o$res, m$res))
    stop("instrumented mirror disagrees with the oracle's .ancombc2_core")
  }
  TRUE
}
