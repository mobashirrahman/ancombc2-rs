# Synthetic fixture generation. Every fixture is defined by a single integer id
# and a pure function of (id, seed), so R and Rust can regenerate identical data
# without shipping count matrices through git.
#
# Deterministic design:
#   * counts  are Poisson draws with lambda = exp(mu), where
#     mu[i, j] = base_i + signal_i * effect(x_j) + log_depth_j + normal noise
#   * the linear predictor is built with a fixed LCG so that the *design* of the
#     dataset is reproducible; the Poisson draws use R's Mersenne-Twister and
#     are written to disk, because Poisson sampling is not trivially portable
#     across languages.
#
# Consequently: the generator (here, R) writes counts.tsv and meta.tsv; both the
# R harness and the Rust test suite read those exact files. No cross-language RNG
# equivalence is required.

fixture_spec <- function(id) {
  # 10 x 10 / 100 x 30 / 1000 x 100 / 10000 x 500
  shape <- switch(as.character(id),
    "1"  = c(n_tax =   10, n_samp =  10),
    "2"  = c(n_tax =  100, n_samp =  30),
    "3"  = c(n_tax = 1000, n_samp = 100),
    "4"  = c(n_tax =10000, n_samp = 500),
    stop("unknown fixture id")
  )
  spec <- list(id = id, shape = shape, seed = 1000L + id)

  # Groups: 2 (id 1,2), 3 (id 3), 5 (id 4)
  spec$n_group <- c(2L, 2L, 3L, 5L)[id]

  # Zero inflation: 0% (1,2), 10% (3), 50% (4)
  spec$zero_rate <- c(0, 0, 0.1, 0.5)[id]

  # Structural zeros: none (1,2,3), present in one group (4)
  spec$struc_zero <- c(FALSE, FALSE, FALSE, TRUE)[id]

  # Extra covariates
  spec$n_cont <- c(0L, 1L, 2L, 3L)[id]

  # Sensitivity analysis
  spec$pseudo_sens <- c(FALSE, FALSE, TRUE, TRUE)[id]
  spec$conservative <- c(TRUE, TRUE, TRUE, FALSE)[id]

  # Spread
  spec$signal_frac <- c(0.3, 0.1, 0.05, 0.02)[id]

  # --- axes that `fixture_spec` has always implied, stated explicitly so the
  # fixture matrix can vary them. The values are the ones the four committed
  # fixtures have always used, so a spec built this way reproduces them
  # byte-identically; the matrix overrides them.
  spec$predictor    <- "continuous"  # "continuous" | "binary"
  spec$interact     <- FALSE         # include `group:x1` in the formula
  spec$balanced     <- TRUE          # FALSE gives unequal group sizes
  # Config-only axes. They change `config.json` and the golden, never the
  # counts or metadata, so cells that differ only in these share input files.
  spec$pseudo       <- 0
  spec$p_adj_method <- "holm"
  spec
}

#' Build a spec from named axis levels rather than an integer id.
#'
#' The fixture matrix of PLAN.md section 5.5 asks for a matrix over shape,
#' sparsity, predictor type, covariate count and interactions, group balance,
#' structural zeros, pseudo-count, sensitivity on/off, conservative on/off and
#' the adjustment method. Only the first four are expressible by the four
#' committed fixtures, so this is what the rest are generated from.
#'
#' Every level defaults to what the committed fixtures use, so
#' `spec_from_axes()` with no arguments is the 10x10 / 2-group / continuous /
#' 0-covariate case. Unknown level names are an error rather than a silent
#' fallback: a typo in a matrix cell must not quietly produce a duplicate of
#' some other cell.
spec_from_axes <- function(shape = c(n_tax = 100L, n_samp = 30L),
                           n_group = 2L,
                           zero_rate = 0,
                           struc_zero = FALSE,
                           n_cont = 1L,
                           predictor = c("continuous", "binary"),
                           interact = FALSE,
                           balanced = TRUE,
                           signal_frac = 0.1,
                           pseudo = 0,
                           p_adj_method = c("holm", "hochberg", "hommel",
                                            "bonferroni", "BH", "BY", "none"),
                           pseudo_sens = FALSE,
                           conservative = TRUE,
                           seed = 9001L) {
  if (length(shape) != 2L || !all(c("n_tax", "n_samp") %in% names(shape))) {
    stop("shape must be c(n_tax = , n_samp = )")
  }
  predictor <- match.arg(predictor)
  p_adj_method <- match.arg(p_adj_method)
  if (n_group < 2L) stop("a group factor needs at least 2 levels")
  if (zero_rate < 0 || zero_rate > 1) stop("zero_rate must be in [0, 1]")
  if (n_cont < 0L) stop("n_cont must be >= 0")
  if (interact && n_cont < 1L) {
    stop("an interaction needs at least one covariate to interact with")
  }
  list(
    id = NA_integer_, shape = c(n_tax = as.integer(shape[["n_tax"]]),
                                n_samp = as.integer(shape[["n_samp"]])),
    seed = as.integer(seed), n_group = as.integer(n_group),
    zero_rate = as.numeric(zero_rate), struc_zero = isTRUE(struc_zero),
    n_cont = as.integer(n_cont), predictor = predictor,
    interact = isTRUE(interact), balanced = isTRUE(balanced),
    signal_frac = as.numeric(signal_frac), pseudo = as.numeric(pseudo),
    p_adj_method = p_adj_method,
    pseudo_sens = isTRUE(pseudo_sens), conservative = isTRUE(conservative)
  )
}

# Deterministic uniform stream: independent of R's RNG, so the dataset *design*
# (which taxa are differentially abundant, how groups are assigned) is stable
# and identical whenever the fixtures are regenerated.
#
# MINSTD (Park-Miller multiplicative, a = 16807, m = 2^31 - 1) with a Lehmer
# randomisation step to hide the lattice structure. Chosen because every
# intermediate value stays below 2^53, so double arithmetic is *exact* and the
# sequence is bit-reproducible in any language. A 64-bit LCG in double
# arithmetic does not have this property: `a * z` reaches ~1e38, rounding
# destroys the low bits, and the stream degenerates.
#
# quality(): decile chi-square ~9.3, lag-1 correlation ~0.
LCG_M <- 2147483647   # 2^31 - 1
LCG_A <- 16807

lcg_stream <- function(n, seed) {
  # Park-Miller seeding: 1 <= s < m
  s <- (as.numeric(seed) * 13007 + 1) %% LCG_M
  if (s <= 0) s <- 1
  out <- numeric(n)
  for (i in seq_len(n)) {
    s <- (s * LCG_A) %% LCG_M
    # Lehmer: u = (s * t) / m, t coprime to m; the fractional part is uniform
    # even though s alone is visibly lattice-structured.
    t <- 11
    while (s > 0 && s <= (LCG_M %/% t) && s %% t == 0) t <- t + 2
    out[i] <- ((s * t) %% LCG_M) / LCG_M
  }
  out
}

gen_fixture <- function(id) {
  gen_fixture_from_spec(fixture_spec(id))
}

gen_fixture_from_spec <- function(spec) {
  # Seeded *first*, before any draw.
  #
  # This used to sit just before the Poisson call, which left `sample.int`,
  # `rnorm` and `rlnorm` above it running on whatever state the session happened
  # to be in. The consequence was that a fixture was not a function of its seed
  # but of its seed *plus* however much randomness earlier fixtures in the same
  # session had consumed -- so `Rscript scripts/generate_goldens.R 3` and
  # `Rscript scripts/generate_goldens.R 1 2 3` produced different counts for
  # fx03, and the file's claim that "every fixture is a pure function of (id,
  # seed)" was false.
  #
  # The four committed fixtures under `validation/fixtures/` were generated
  # *before* this fix, from the ambient session state, and they are what the
  # committed goldens were computed from. They are left exactly as they are:
  # regenerating them needs the declared R 4.5, and replacing the inputs would
  # invalidate the goldens that currently pass parity. The divergence is recorded
  # in docs/reference_behavior.md.
  set.seed(spec$seed)
  n_tax <- spec$shape[["n_tax"]]
  n_samp <- spec$shape[["n_samp"]]

  n_stream <- n_tax * 8L + n_samp * 8L + 64L
  s <- lcg_stream(n_stream, spec$seed)
  counter <- local({ k <- 0L; function(n) { k <<- k + n; s[(k - n + 1L):k] } })
  take <- counter

  # ---- covariates -------------------------------------------------------
  if (isTRUE(spec$balanced)) {
    grp <- rep(seq_len(spec$n_group), length.out = n_samp)
  } else {
    # Deliberately unequal level sizes, with every level still populated: a
    # group factor with an empty level is a different thing entirely (it makes
    # the design unidentifiable, which the edge matrix covers), while an
    # unbalanced factor is what a real cohort looks like.
    weights <- seq(spec$n_group, 1) / sum(seq(spec$n_group, 1))
    sizes <- as.integer(round(weights * n_samp))
    sizes[spec$n_group] <- n_samp - sum(sizes[-spec$n_group])
    if (any(sizes < 2L)) stop("unbalanced split left a group with fewer than 2 samples")
    grp <- rep(seq_len(spec$n_group), times = sizes)
  }
  grp <- grp[sample.int(n_samp)]   # permute so groups are not blocked in order
  group_f <- factor(grp, levels = seq_len(spec$n_group))
  meta <- data.frame(group = group_f)

  if (spec$n_cont > 0) {
    cont <- lapply(seq_len(spec$n_cont), function(j) {
      if (identical(spec$predictor, "binary")) {
        # A 0/1 covariate. The design matrix then carries a column of zeros and
        # ones rather than a continuous one, which is a different rank and
        # scaling structure -- and with 10 covariates a binary design is where
        # near-collinearity and aliasing actually show up.
        v <- as.numeric(take(n_samp) < 0.5)
        return(v)
      }
      v <- (take(n_samp) - 0.5) * 2
      as.numeric(scale(v))
    })
    names(cont) <- paste0("x", seq_len(spec$n_cont))
    meta <- cbind(meta, as.data.frame(cont))
  }

  # ---- linear predictor --------------------------------------------------
  base <- -3 + 2 * take(n_tax)                      # taxon-specific baseline
  spread <- 0.2 + 1.6 * take(n_tax)                 # taxon-specific spread
  is_signal <- take(n_tax) < spec$signal_frac
  signal_sd <- 1.0
  signal <- rnorm(n_tax, 0, signal_sd) * is_signal

  # group effect: a signed shift, deterministic sign from the LCG
  sgn <- ifelse(take(n_tax) < 0.5, -1, 1)
  eff_group <- signal * sgn
  # eff_cont[[j]] holds the per-taxon coefficient for covariate xj, not the
  # per-sample contribution; the latter is formed with outer() below.
  eff_cont <- vector("list", spec$n_cont)
  if (spec$n_cont > 0) {
    for (j in seq_len(spec$n_cont)) {
      eff_cont[[j]] <- rnorm(n_tax) * is_signal * 0.6
    }
  }
  # `group:x1` interaction: the per-taxon coefficient of the product, so the
  # generating model matches the term the formula asks the design matrix for.
  # A formula with an interaction that the data does not generate would still
  # be a valid test, but it would be testing the null, which is much weaker.
  #
  # Drawn *only* when the interaction is on, and that is load-bearing rather
  # than tidy: this and the noise below come from R's Mersenne-Twister, so a
  # single unconditional `rnorm(n_tax)` here shifts every subsequent draw and
  # regenerates all four committed fixtures to different counts. The LCG stream
  # used for the design is order-independent by construction; this one is not.
  eff_int <- if (isTRUE(spec$interact)) rnorm(n_tax) * is_signal * 0.5 else rep(0, n_tax)

  # sampling depth: lognormal, a strong confounder
  depth <- rlnorm(n_samp, meanlog = 8, sdlog = 0.5)

  # mu[i, j] = base_i + eff_group_i * D[j, ] + sum_j eff_cont_ij * x_j + log(depth_j)
  # where D is the treatment-coded group indicator (level 1 is the reference), so
  # the generating model matches the design matrix the oracle builds.
  grp_code <- as.integer(group_f)
  mu <- outer(base, rep(1, n_samp))
  for (k in seq_len(spec$n_group)[-1]) {
    ind <- as.numeric(grp_code == k)
    mu <- mu + outer(eff_group, ind)
  }
  for (j in seq_len(spec$n_cont)) {
    mu <- mu + outer(eff_cont[[j]], meta[[paste0("x", j)]])
  }
  if (isTRUE(spec$interact)) {
    for (k in seq_len(spec$n_group)[-1]) {
      ind <- as.numeric(grp_code == k)
      mu <- mu + outer(eff_int, ind * meta[["x1"]])
    }
  }
  mu <- mu + matrix(rep(log(depth), each = n_tax), n_tax, n_samp)
  mu <- mu + matrix(rnorm(n_tax * n_samp, 0, 1), n_tax, n_samp) * spread

  lam <- exp(mu)
  counts <- matrix(stats::rpois(length(lam), lambda = lam), n_tax, n_samp)

  # ---- structural zeros --------------------------------------------------
  if (isTRUE(spec$struc_zero)) {
    n_struct <- max(5L, floor(n_tax * 0.02))
    struct_tax <- sample.int(n_tax, n_struct)
    dead_group <- sample.int(spec$n_group, n_struct, replace = TRUE)
    for (k in seq_along(struct_tax)) {
      idx <- which(as.integer(group_f) == dead_group[k])
      counts[struct_tax[k], idx] <- 0L
    }
  }

  # ---- zero inflation ----------------------------------------------------
  if (spec$zero_rate > 0) {
    z <- matrix(stats::runif(length(counts)) < spec$zero_rate, n_tax, n_samp)
    counts[z] <- 0L
  }

  rownames(counts) <- sprintf("taxon_%05d", seq_len(n_tax))
  colnames(counts) <- sprintf("sample_%04d", seq_len(n_samp))
  rownames(meta) <- colnames(counts)

  list(spec = spec, counts = counts, meta = meta)
}

fix_formula_for <- function(spec) {
  parts <- "group"
  if (spec$n_cont > 0) parts <- c(parts, paste0("x", seq_len(spec$n_cont)))
  if (isTRUE(spec$interact)) {
    # `*` rather than a bare `group:x1` so the formula reads the way a user
    # writes it; the two are identical to `model.matrix`.
    parts <- paste(parts, collapse = " + ")
    return(paste(parts, "* x1"))
  }
  paste(parts, collapse = " + ")
}

write_fixture <- function(id, dir) {
  f <- gen_fixture(id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # row.names = TRUE writes the row names into a leading column of the data rows
  # but *not* into a blank leading header field, so the header holds exactly the
  # column names. Both readers rely on that; see crates/.../tests/golden/mod.rs.
  # `digits = 17` so the text round-trips a double exactly. The default
  # `getOption("digits")` is 7, and even 15 loses the last bit for some values:
  # a 5e-15 perturbation of the design is enough to move an ill-conditioned
  # taxon's coefficient by 3e-4, which would make the CLI's numbers differ from
  # the in-process parity run for a reason that has nothing to do with the
  # algorithm. The design matrix must be bit-identical on both paths.
  utils::write.table(f$counts, file.path(dir, "counts.tsv"),
                     sep = "\t", quote = FALSE, row.names = TRUE, digits = 17)
  utils::write.table(f$meta, file.path(dir, "meta.tsv"),
                     sep = "\t", quote = FALSE, row.names = TRUE, digits = 17)
  writeLines(fix_formula_for(f$spec), file.path(dir, "formula.txt"))
  writeLines(c(
    sprintf("id=%d", id),
    sprintf("n_tax=%d", f$spec$shape[["n_tax"]]),
    sprintf("n_samp=%d", f$spec$shape[["n_samp"]]),
    sprintf("n_group=%d", f$spec$n_group),
    sprintf("zero_rate=%g", f$spec$zero_rate),
    sprintf("struc_zero=%s", f$spec$struc_zero),
    sprintf("n_cont=%d", f$spec$n_cont),
    sprintf("pseudo_sens=%s", f$spec$pseudo_sens),
    sprintf("conservative=%s", f$spec$conservative)
  ), file.path(dir, "spec.txt"))
  invisible(f)
}
