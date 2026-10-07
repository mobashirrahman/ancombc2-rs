#!/usr/bin/env Rscript
# S09: does one least-squares fit agree with the reference, bit for bit?
#
# The claim under test is the narrowest one in the plan: *for the same design and the
# same theta-adjusted response, the pipeline's `lm_fit_all` produces R's `beta`,
# `fitted` and `dof` exactly.* One fit, one group, before any iteration and before the
# trace.
#
# The R side is transcribed from the pinned `ancombc_bias_correct.R:11-84`, and the
# single-group comparison is against `stats::lm.fit` directly -- which is the primitive
# the reference calls, so there is nothing between the two numbers being compared.
#
# What this does NOT establish: anything about the MLE's iterations, the EM step, the
# sandwich, or the assembled result. Those come later and each has its own instrument.

suppressWarnings(suppressMessages(library(ANCOMBC)))

fails <- 0L
checks <- 0L
diffs <- new.env(parent = emptyenv())

ok <- function(label, cond, why = "") {
  checks <<- checks + 1L
  if (isTRUE(cond)) {
    cat(sprintf("ok    %s\n", label))
  } else {
    fails <<- fails + 1L
    cat(sprintf("FAIL  %s%s\n", label, if (nzchar(why)) paste0(" -- ", why) else ""))
    assign(label, why, envir = diffs)
  }
}

# `==` on doubles is the wrong comparison twice over: `NaN != NaN`, and a difference
# of one ulp passes nothing useful. Bits, and `serialize()`.
same_bits <- function(a, b) {
  if (length(a) != length(b)) return(FALSE)
  ba <- writeBin(as.double(a), raw())
  bb <- writeBin(as.double(b), raw())
  identical(ba, bb)
}
same_bytes <- function(a, b) {
  identical(serialize(a, NULL, ascii = FALSE, xdr = TRUE, version = 3),
            serialize(b, NULL, ascii = FALSE, xdr = TRUE, version = 3))
}

# Where a stage differs: the first few differing positions with both bit patterns.
# "differs" alone is not actionable when the array is 5000 x 3.
# The bit pattern of one value, 16 hex digits.
#
# Four unsigned 16-bit reads rather than one 8-byte read: `readBin` rejects
# `size = 8L` for "integer", and reading 4 bytes and zero-extending is what made
# every pattern print as `0x00000000........`. Little-endian, so the halves go in
# reverse order. `rev(w)` rather than reading backwards, because the halves are
# already in the wrong order by being little-endian.
hex1 <- function(v) {
  if (is.integer(v)) {
    return(sprintf("%08x", as.numeric(readBin(writeBin(v, raw()), "integer", 4L,
                                             signed = FALSE))))
  }
  r8 <- writeBin(as.double(v), raw())
  w <- vapply(0:3, function(k) as.numeric(
    readBin(r8[(2L * k + 1L):(2L * k + 2L)], "integer", n = 1L, size = 2L,
            signed = FALSE)), numeric(1L))
  sprintf("%04x%04x%04x%04x", w[4], w[3], w[2], w[1])
}

# Where a stage differs: the first few differing positions with both bit patterns.
# "differs" alone is not actionable when the array is 5000 x 3.
first_diff <- function(got, want, k = 4L) {
  g <- as.vector(got); w <- as.vector(want)
  if (length(g) != length(w)) {
    return(paste0("lengths ", length(g), " vs ", length(w)))
  }
  if (!is.null(dim(want)) && !identical(dim(got), dim(want))) {
    return(paste0("shape ", paste(dim(got), collapse = "x"), " vs ",
                  paste(dim(want), collapse = "x")))
  }
  bg <- vapply(seq_along(g), function(i) hex1(g[i]), character(1L))
  bw <- vapply(seq_along(w), function(i) hex1(w[i]), character(1L))
  at <- which(bg != bw)
  if (length(at) == 0L) return("bits agree but serialize() does not")
  paste0(head(paste0(at, ": got 0x", bg[at], " want 0x", bw[at]), k), collapse = "; ")
}

# ---------------------------------------------------------------------------
# `.lm_fit_all`, transcribed from ancombc_bias_correct.R:11-84
# ---------------------------------------------------------------------------
ref_lm_fit_all <- function(x, Ymat, meta_data = NULL, tformula = NULL) {
  # `meta_data` and `tformula` exist only for the per-taxon `lm()` fallback, which is
  # where the reference's surprises live. A fixture that triggers the fallback
  # without them is a fixture this transcription cannot judge, and says so.
  n_tax <- nrow(Ymat); n_samp <- ncol(Ymat); p <- ncol(x)
  fix_eff <- colnames(x)
  tax_id <- rownames(Ymat); samp_id <- colnames(Ymat)
  x_ok <- stats::complete.cases(x)
  beta <- matrix(NA_real_, nrow = n_tax, ncol = p, dimnames = list(tax_id, fix_eff))
  fitted <- matrix(NA_real_, nrow = n_tax, ncol = n_samp, dimnames = list(tax_id, samp_id))
  dof <- rep(999L, n_tax); names(dof) <- tax_id

  fit_one <- function(i) {
    if (is.null(tformula)) {
      stop("this fixture reaches the per-taxon lm() fallback and has no formula, so ",
           "this transcription cannot judge it; supply meta_data and tformula",
           call. = FALSE)
    }
    df <- data.frame(y_crt = Ymat[i, ], meta_data)
    fit <- suppressWarnings(try(stats::lm(tformula, data = df), silent = TRUE))
    if (inherits(fit, "lm")) {
      bi <- rep(0, p); ci <- stats::coef(fit)
      bi[match(names(ci), fix_eff)] <- ci
      beta[i, ] <<- bi
      fi <- rep(0, n_samp); fv <- stats::fitted(fit)
      fi[match(names(fv), samp_id)] <- fv
      fitted[i, ] <<- fi
      dof[i] <<- fit$df.residual
    }
  }

  if (all(x_ok) && all(is.finite(Ymat))) {
    use <- matrix(TRUE, nrow = n_tax, ncol = n_samp)
    groups <- list(seq_len(n_tax))
  } else {
    use <- is.finite(Ymat) & matrix(x_ok, nrow = n_tax, ncol = n_samp, byrow = TRUE)
    keys <- do.call(paste0, asplit(use * 1L, 2L))
    groups <- split(seq_len(n_tax), factor(keys, levels = unique(keys)))
  }

  for (idx in groups) {
    rows <- use[idx[1L], ]
    if (!any(rows)) { for (i in idx) fit_one(i); next }
    xr <- x[rows, , drop = FALSE]
    Yr <- t(Ymat[idx, rows, drop = FALSE])
    fit <- stats::lm.fit(xr, Yr)
    if (fit$rank < ncol(xr)) { for (i in idx) fit_one(i); next }
    co <- fit$coefficients
    if (is.null(dim(co))) co <- matrix(co, ncol = 1L)
    beta[idx, ] <- t(co)
    fv <- fit$fitted.values
    if (is.null(dim(fv))) fv <- matrix(fv, ncol = 1L)
    fitted[idx, rows] <- t(fv)
    fitted[idx, !rows] <- 0
    dof[idx] <- fit$df.residual
  }
  list(beta = beta, fitted = fitted, dof = dof)
}

# ---------------------------------------------------------------------------
# The fixtures
# ---------------------------------------------------------------------------
#
# Full-rank and well-conditioned first, then the cases where the plan says the
# reference does something surprising: a scaled column (where pivoting changes the
# answer), a missing count (which changes a taxon's group), an incomplete design row
# (which changes every group's rows), and a rank-deficient design (which takes the
# per-taxon fallback).

fixture <- function(name, x, y, obs = NULL) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  y <- as.matrix(y); storage.mode(y) <- "double"
  n_tax <- nrow(y); n_samp <- nrow(x)
  if (is.null(obs)) {
    obs <- matrix(is.finite(y), nrow = n_tax, ncol = n_samp)
  }
  obs <- matrix(as.logical(obs), nrow = n_tax, ncol = n_samp)
  list(name = name, x = x, y = y, observed = obs,
       n_samp = n_samp, p = ncol(x), n_taxa = n_tax)
}

set.seed(11)
n <- 12
x_int <- cbind(1, rep(c(-1, 1), each = n / 2))
x_fac <- cbind(1, rep(c(0, 1), each = n / 2), rep(c(0, 1), times = n / 2))
y_a <- matrix(rnorm(6 * n), 6, n) + 1        # finite everywhere
y_na <- y_a; y_na[2, 3] <- NA                 # one missing count
y_big <- y_a; y_big[1, ] <- y_big[1, ] * 1e12 # one taxon on a much larger scale

fixtures <- list(
  fixture("full-rank", x_int, y_a),
  fixture("two-factors", x_fac, y_a),
  fixture("scaled-column", cbind(1, rep(c(0.5, -1.25, 2, 0.75), length.out = n)), y_a),
  fixture("missing-one", x_int, y_na),
  fixture("missing-all-in-a-taxon", x_int, {
    z <- y_a; z[3, ] <- NA; z
  }),
  fixture("huge-response", x_int, y_big)
)

# An incomplete design row: `stats::complete.cases(x)` is FALSE for it, so every
# taxon's group loses that sample. `lm.fit` is never handed an NA.
x_incomplete <- x_int
x_incomplete[5, 2] <- NA
fixtures <- c(fixtures, list(
  fixture("incomplete-design-row", x_incomplete, y_a),
  fixture("incomplete-row-and-missing", x_incomplete, y_na)
))

ok("the bridge is present", ANCOMBC:::rb_available())

# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

for (f in fixtures) {
  cat(sprintf("\n--- %s  (x %d x %d, y %d x %d)\n", f$name, f$n_samp, f$p,
              f$n_taxa, f$n_samp))
  got <- try(ANCOMBC:::rb_fit_probe(f$x, f$y, f$observed), silent = TRUE)
  if (inherits(got, "try-error")) {
    ok(paste0(f$name, ": the probe runs"), FALSE, as.character(got))
    next
  }
  ok(paste0(f$name, ": the probe runs"), TRUE)

  # Which rows and taxa each pattern group got. `dof` is `n_used - rank`, so a `dof`
  # difference is either a row-count difference (visible here) or a rank difference.
  g <- got$groups
  ng <- g[4L]
  # One (rows, taxa) pair per group, and only `ng` of them: the buffer is sized for
  # `n_taxa` groups, so reading past `ng` walks into unwritten memory.
  pairs <- if (ng > 0L) {
    g[8L + seq_len(2L * ng) - 1L]
  } else numeric(0L)
  cat(sprintf("  dims: n_samp %d, p %d, n_taxa %d, n_groups %d, raw_mask %d of %d, x_ok %d, x NaN %s\n",
              g[1L], g[2L], g[3L], ng, g[5L], g[1L] * g[3L], g[6L], g[7L]))
  cat(sprintf("  groups: %d %s\n", ng,
              if (ng > 0L) paste(sprintf("[%d rows x %d taxa]",
                                        pairs[seq.int(1L, 2L * ng, by = 2L)],
                                        pairs[seq.int(2L, 2L * ng, by = 2L)]),
                                collapse = " ") else ""))

  # A two-column design is `(Intercept)` + one covariate, so the per-taxon `lm()`
  # fallback can be given the formula it would have had. The probe returns plain
  # matrices (names are restored by `rb_restore_names` in the assembly), so the
  # reference's dimnames are dropped before the byte comparison: a name is not a
  # number, and comparing it here would fail the fit for a transport detail.
  xr_ <- f$x
  meta_ <- NULL; tf_ <- NULL
  if (ncol(xr_) == 2L) {
    colnames(xr_) <- c("(Intercept)", "x1")
    meta_ <- data.frame(x1 = xr_[, 2L]); tf_ <- y_crt ~ x1
  }
  want <- ref_lm_fit_all(xr_, f$y, meta_, tf_)
  want$beta <- unname(want$beta); want$fitted <- unname(want$fitted)
  want$dof <- unname(want$dof)

  ok(paste0(f$name, ": beta is byte-identical"),
     same_bytes(got$beta, want$beta), first_diff(got$beta, want$beta))
  ok(paste0(f$name, ": fitted is byte-identical"),
     same_bytes(got$fitted, want$fitted), first_diff(got$fitted, want$fitted))
  # `dof` is `INTSXP` on the reference side and `f64` here, so compare the values:
  # 10 as an integer and 10 as a double are the same degree of freedom, and a
  # container comparison would report a failure that says nothing about the fit.
  ok(paste0(f$name, ": dof agrees"),
     isTRUE(all.equal(got$dof, want$dof)), first_diff(got$dof, want$dof))

  # Per-taxon, because a whole-array comparison cannot say *which* taxon, and
  # because `.lm_fit_all` is per taxon for everything except the shared QR.
  for (i in seq_len(f$n_taxa)) {
    lbl <- sprintf("%s: taxon %d", f$name, i)
    ok(paste(lbl, "beta"),
       same_bits(got$beta[i, ], want$beta[i, ]), first_diff(got$beta[i, ], want$beta[i, ]))
    ok(paste(lbl, "fitted"),
       same_bits(got$fitted[i, ], want$fitted[i, ]),
       first_diff(got$fitted[i, ], want$fitted[i, ]))
    ok(paste(lbl, "dof"),
       got$dof[i] == want$dof[i],
       paste("got", got$dof[i], "want", want$dof[i]))
  }
}

# ---------------------------------------------------------------------------
# One group against `stats::lm.fit` itself
# ---------------------------------------------------------------------------
#
# `.lm_fit_all` is a thin wrapper over one `lm.fit` per group, so comparing a single
# group against `lm.fit` directly removes the wrapper from the question. It is also
# where the primitive's own behaviour -- `fitted.values = y - residuals`, `NA` for an
# aliased coefficient, `dof = n - rank` -- is visible without anything on top of it.

cat("\n--- one group against stats::lm.fit directly\n")
x1 <- cbind(1, rep(c(-1, 1), each = 3), rep(c(0.5, 1.5), times = 3))
y1 <- matrix(c(2, -1, 0.5, 3, -2, 1, 4, -3, 2.5, -4, 3.5, 1.5), 4, 6)
obs1 <- matrix(TRUE, 4, 6)
g <- ANCOMBC:::rb_fit_probe(x1, y1, obs1)
direct <- stats::lm.fit(x1, t(y1))

ok("one group: beta matches lm.fit's, bit for bit",
   same_bits(as.vector(g$beta), as.vector(t(direct$coefficients))),
   first_diff(as.vector(g$beta), as.vector(t(direct$coefficients))))
ok("one group: dof matches lm.fit's",
   all(g$dof == direct$df.residual),
   paste("got", paste(g$dof, collapse = ","), "want", direct$df.residual))
ok("one group: fitted matches lm.fit's, bit for bit",
   same_bits(as.vector(g$fitted), as.vector(t(direct$fitted.values))),
   first_diff(as.vector(g$fitted), as.vector(t(direct$fitted.values))))
ok("lm.fit's own fitted.values is not x %*% coef on this input",
   !same_bits(as.vector(direct$fitted.values), as.vector(x1 %*% direct$coefficients)),
   "it is, so the fitted path needs no care here")
ok("lm.fit's own fitted.values is not y - (y - x %*% coef) either",
   !same_bits(as.vector(direct$fitted.values),
              as.vector(t(y1) - (t(y1) - x1 %*% direct$coefficients))),
   "it is, so the residuals come from the solver")

# A rank-deficient group, which takes `.lm_fit_all`'s per-taxon `lm()` fallback.
#
# The fallback exists for one situation: a taxon whose usable samples miss a whole level
# of the group factor, so a treatment contrast is all zero over its rows. `lm` then
# *re-levels* (`drop.unused.levels`), the dropped contrast has no name in `coef()`, and
# `fit_one` leaves a literal 0 there. The probe treats every non-intercept column as a
# group contrast, so this is the case it models -- a numeric duplicated column is not,
# and comparing against it would be comparing two different models.
g3 <- factor(rep(c("a", "b", "c"), each = 4))
x3 <- stats::model.matrix(~ g3); storage.mode(x3) <- "double"
set.seed(5)
y3 <- matrix(rnorm(4 * 12), 4, 12)
y3[2, g3 == "c"] <- NA      # taxon 2 never sees level c
y3[3, g3 == "b"] <- NA      # taxon 3 never sees level b
obs3 <- is.finite(y3)
gd <- try(ANCOMBC:::rb_fit_probe(x3, y3, obs3), silent = TRUE)
ok("a missing factor level does not crash the probe", !inherits(gd, "try-error"),
   if (inherits(gd, "try-error")) as.character(gd) else "")
if (!inherits(gd, "try-error")) {
  y3n <- y3; colnames(y3n) <- as.character(seq_len(ncol(y3)))
  x3n <- x3; colnames(x3n) <- c("(Intercept)", "g3b", "g3c")
  wd <- ref_lm_fit_all(x3n, y3n, data.frame(g3 = g3), y_crt ~ g3)
  wd$beta <- unname(wd$beta); wd$fitted <- unname(wd$fitted); wd$dof <- unname(wd$dof)
  ok("missing level: beta is byte-identical to the per-taxon lm()'s",
     same_bytes(gd$beta, wd$beta), first_diff(gd$beta, wd$beta))
  ok("missing level: fitted is byte-identical to the per-taxon lm()'s",
     same_bytes(gd$fitted, wd$fitted), first_diff(gd$fitted, wd$fitted))
  ok("missing level: dof is n_used - rank of the re-levelled model",
     isTRUE(all.equal(gd$dof, wd$dof)),
     paste("got", paste(gd$dof, collapse = ","), "want", paste(wd$dof, collapse = ",")))
}

# ---------------------------------------------------------------------------
# What R's primitive actually does, recorded so the port can be checked against it
# ---------------------------------------------------------------------------

cat("\n--- the primitive's own behaviour\n")
xdup <- cbind(1, c(1, 2, 3, 4), c(1, 2, 3, 4))
fd <- stats::lm.fit(xdup, c(1, 2, 3, 4))
ok("lm.fit's rank is 2 on a duplicated column", fd$rank == 2L, fd$rank)
ok("lm.fit's aliased coefficient is NA_real_, not zero",
   is.na(fd$coefficients[3]) && identical(unname(fd$coefficients[3]), NA_real_),
   paste("got", fd$coefficients[3]))
ok("lm.fit's df.residual is n - rank", fd$df.residual == nrow(xdup) - fd$rank,
   fd$df.residual)
ok("lm.fit's rank does not depend on y",
   {
     rs <- vapply(1:8, function(k) { set.seed(k); stats::lm.fit(x1, matrix(rnorm(6 * 4), 6, 4))$rank }, 0L)
     length(unique(rs)) == 1L
   },
   "the .iter_mle smoke test draws rnorm, so this decides whether it is deterministic")

xu <- cbind(1, c(1, 2), c(2, 5))
fu <- stats::lm.fit(xu, c(3, 7))
ok("lm.fit's df.residual is 0 when underdetermined", fu$df.residual == 0L,
   fu$df.residual)

cat(sprintf("\n%d/%d checks passed\n", checks - fails, checks))
if (fails > 0L) {
  cat("\nFailing checks:\n")
  for (nm in sort(ls(diffs))) cat(sprintf("  %s: %s\n", nm, get(nm, envir = diffs)))
  quit(status = 1L)
}
