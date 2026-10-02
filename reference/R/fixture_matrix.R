# The golden fixture matrix of PLAN.md section 5.5.
#
# PLAN.md asks for "a golden fixture matrix over shape 10x10/100x30/1000x100/
# 10000x500, sparsity zero-free/10/50/90%, continuous/binary/3-group/5-group
# predictors, 1/5/10 covariates and interactions, balanced/unbalanced groups,
# structural zeros present/absent, pseudo 0/0.1/0.5/1.0, sensitivity on/off,
# conservative on/off, and all 7 adjust methods".
#
# The full cross product of those axes is 4*4*4*4*2*2*4*2*2*7 = 245,760 cells,
# which is neither runnable nor useful: almost every one of them is a
# combination of axes the pipeline treats independently. What is actually
# required is that every axis level is exercised and that the interesting
# interactions are covered, so the matrix is built as
#
#   * one **core** cell per axis, at a middle level of every other axis, and
#   * a **sweep** over each axis's levels, one factor at a time, and
#   * a small set of **interaction** cells for the combinations that are not
#     independent in the implementation.
#
# The core cell is what a regression would be caught by; the sweeps are what
# proves each level works; the interactions are the ones with a real reason to
# interact. Every cell is a pure function of its name and seed (see
# `fixtures.R`), so the whole matrix regenerates identically.
#
# Cells that differ only in `pseudo` or `p_adj_method` share their counts and
# metadata: those are configuration, not data, so `fixture_matrix` groups them
# and `scripts/generate_matrix_goldens.R` writes one input directory per input
# group and one golden per cell.

# `core` is the reference cell: 100x30, 2 groups, continuous, 1 covariate,
# balanced, no structural zeros, no zero inflation, pseudo 0, holm, no
# sensitivity. Every sweep below varies exactly one thing from it.
core_axes <- list(
  shape = c(n_tax = 100L, n_samp = 30L),
  n_group = 2L, zero_rate = 0, struc_zero = FALSE, n_cont = 1L,
  predictor = "continuous", interact = FALSE, balanced = TRUE,
  signal_frac = 0.1, pseudo = 0, p_adj_method = "holm",
  pseudo_sens = FALSE, conservative = TRUE
)

# How much of a cell may be in the rank-deficient indeterminacy class.
#
# The committed golden fixtures are held to 1% -- see `MAX_RANK_DEFICIENT_SHARE`
# in the Rust harness, which is the Level D contract for them. A matrix cell is a
# different thing: a 90%-zero table with a five-level factor and six columns
# *should* produce a large class, because a taxon observed in a handful of samples
# has a singular sub-design by construction rather than by accident.
#
# So the ceiling is declared per cell, from the cell's axes and not from its
# result -- a budget fitted to the answer would be no budget at all. Cells that are
# sparse, or that carry many covariates against few samples, get a ceiling; the
# rest keep the 1% default and a `modifyList` typo would show up as a cell
# inheriting the wrong one.
DEFAULT_RANK_DEFICIENT_CAP <- 0.01

#' The declared ceiling for a cell, derived from its axes.
rank_deficient_cap_for <- function(axes) {
  # Sparsity and width are the two axes that drive a taxon into a singular
  # sub-design. 90% zeros and `p = 6` against 30 samples is the extreme corner;
  # everything else stays at the default.
  cap <- DEFAULT_RANK_DEFICIENT_CAP
  if (axes$zero_rate >= 0.5) cap <- max(cap, 0.30)
  if (axes$n_cont >= 10L) cap <- max(cap, 0.10)
  if (axes$n_cont >= 10L && axes$zero_rate >= 0.5) cap <- max(cap, 0.40)
  # A table with fewer than 10 samples per design column cannot determine a
  # coefficient for most taxa, whatever its sparsity.
  if (axes$shape[["n_samp"]] < 4L * (axes$n_cont + axes$n_group)) {
    cap <- max(cap, 0.30)
  }
  cap
}

# Build one cell: a name, the axis levels, and the seed.
#
# The seed comes from `seed_group` when given, otherwise from the cell name --
# so the matrix does not have to be reseeded by hand when a cell is added.
#
# `seed_group` exists for the config-only axes (pseudo-count and adjustment
# method). Those change `config.json` and the golden but not the data, so a
# sweep over them is only a clean test of the configuration if every cell in the
# sweep analyses *the same table*. Giving each its own seed would quietly turn a
# configuration sweep into a data sweep as well, and a failure would not say
# which of the two broke.
cell <- function(name, ..., seed_group = name) {
  axes <- utils::modifyList(core_axes, list(...))
  list(name = name, axes = axes,
       rank_deficient_cap = rank_deficient_cap_for(axes),
       seed = fixture_matrix_seed(seed_group))
}

# A stable seed per cell name: `sum(utf8ToInt(name))` is too weak on its own, so
# it is combined with a fixed multiplier and reduced mod the MINSTD prime the
# generator uses. Deterministic, order-independent, and a new cell cannot
# accidentally reuse a neighbour's data.
fixture_matrix_seed <- function(name) {
  s <- sum(utf8ToInt(name))
  as.integer((s * 7919L + 104729L) %% 2147483647L)
}

# The spec for a cell: every axis it declares, with the cell's own seed.
spec_of <- function(cell) {
  a <- cell$axes
  spec_from_axes(shape = a$shape, n_group = a$n_group, zero_rate = a$zero_rate,
                 struc_zero = a$struc_zero, n_cont = a$n_cont,
                 predictor = a$predictor, interact = a$interact,
                 balanced = a$balanced, signal_frac = a$signal_frac,
                 pseudo = a$pseudo, p_adj_method = a$p_adj_method,
                 pseudo_sens = a$pseudo_sens, conservative = a$conservative,
                 seed = cell$seed)
}

fixture_matrix_cells <- function() {
  cells <- list()

  # --- shape -------------------------------------------------------------
  # Each shape also has to be big enough for the feature to be visible: a 10x10
  # table cannot exercise a 5-group factor or 10 covariates, and pretending
  # otherwise would put a cell in the matrix that silently tests nothing. It can
  # carry one covariate, so the shape sweep varies *only* the shape and every
  # other axis is held at the core level -- otherwise the sweep would be two
  # sweeps at once and a failure could not be attributed.
  cells[["shape-10x10"]] <- cell("shape-10x10", shape = c(n_tax = 10L, n_samp = 10L))
  cells[["shape-100x30"]] <- cell("shape-100x30")
  cells[["shape-1000x100"]] <- cell("shape-1000x100",
                                     shape = c(n_tax = 1000L, n_samp = 100L),
                                     signal_frac = 0.05)
  cells[["shape-10000x500"]] <- cell("shape-10000x500",
                                      shape = c(n_tax = 10000L, n_samp = 500L),
                                      signal_frac = 0.02)

  # --- sparsity ----------------------------------------------------------
  for (z in c(0, 0.1, 0.5, 0.9)) {
    cells[[sprintf("sparsity-%03d", round(z * 100))]] <-
      cell(sprintf("sparsity-%03d", round(z * 100)), zero_rate = z)
  }

  # --- predictor type ----------------------------------------------------
  # A "3-group" or "5-group" predictor here means the group factor itself is
  # that many levels, which is what switches the global and pairwise tests on
  # (`global = n_group >= 3` in the config). Those cells are the only ones in
  # the matrix that exercise either test.
  cells[["predictor-continuous"]] <- cell("predictor-continuous")
  cells[["predictor-binary"]] <- cell("predictor-binary", predictor = "binary")
  cells[["predictor-3group"]] <- cell("predictor-3group", n_group = 3L)
  cells[["predictor-5group"]] <- cell("predictor-5group", n_group = 5L)

  # --- covariate count and interactions ----------------------------------
  for (k in c(1L, 5L, 10L)) {
    cells[[sprintf("covariates-%d", k)]] <- cell(sprintf("covariates-%d", k), n_cont = k)
  }
  cells[["covariates-10-interaction"]] <-
    cell("covariates-10-interaction", n_cont = 10L, interact = TRUE)

  # --- group balance -----------------------------------------------------
  cells[["groups-balanced"]] <- cell("groups-balanced", balanced = TRUE)
  cells[["groups-unbalanced"]] <- cell("groups-unbalanced", balanced = FALSE)

  # --- structural zeros --------------------------------------------------
  cells[["struczero-absent"]] <- cell("struczero-absent", struc_zero = FALSE)
  cells[["struczero-present"]] <- cell("struczero-present", struc_zero = TRUE)

  # --- pseudo-count ------------------------------------------------------
  # Configuration, not data: these four share one counts/meta pair.
  for (p in c(0, 0.1, 0.5, 1.0)) {
    nm <- sprintf("pseudo-%s", format(p, trim = TRUE))
    cells[[nm]] <- cell(nm, pseudo = p, seed_group = "config-sweep")
  }

  # --- sensitivity -------------------------------------------------------
  # `conservative` is only meaningful with `pseudo_sens` on; a cell with
  # `pseudo_sens = FALSE` and `conservative = FALSE` is the same run twice, so
  # the on/off pair is taken with the flag that selects the branch.
  cells[["sensitivity-off"]] <- cell("sensitivity-off", pseudo_sens = FALSE)
  cells[["sensitivity-on-conservative"]] <-
    cell("sensitivity-on-conservative", pseudo_sens = TRUE, conservative = TRUE)
  cells[["sensitivity-on-nonconservative"]] <-
    cell("sensitivity-on-nonconservative", pseudo_sens = TRUE, conservative = FALSE)

  # --- adjustment method -------------------------------------------------
  # All seven, including `none`, which is the one that makes the pipeline's
  # choice of method unobservable and so tests that it is actually threaded
  # through rather than defaulted.
  for (m in c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none")) {
    nm <- sprintf("adjust-%s", m)
    cells[[nm]] <- cell(nm, p_adj_method = m, seed_group = "config-sweep")
  }

  # --- interactions ------------------------------------------------------
  # The combinations with a reason to interact rather than being independent:
  # a 5-group factor with 10 binary covariates is the near-collinear design
  # where aliasing and rank deficiency show up; zero inflation *and* structural
  # zeros is where the two zero mechanisms are told apart; and a 90%-zero table
  # with a 5-group factor is where the prevalence filter bites hardest.
  cells[["int-5group-10cov-binary"]] <-
    cell("int-5group-10cov-binary", n_group = 5L, n_cont = 10L, predictor = "binary")
  cells[["int-struczero-sparsity"]] <-
    cell("int-struczero-sparsity", struc_zero = TRUE, zero_rate = 0.5)
  cells[["int-sparsity90-5group"]] <-
    cell("int-sparsity90-5group", zero_rate = 0.9, n_group = 5L)
  cells[["int-nonconservative-binary"]] <-
    cell("int-nonconservative-binary", predictor = "binary",
         pseudo_sens = TRUE, conservative = FALSE)

  cells
}

# Cells that share their counts.tsv and meta.tsv because they differ only in
# configuration. Grouped so the generator writes each input set once; the golden
# is still per cell, because the configuration changes the answer.
fixture_matrix_input_groups <- function(cells) {
  key_of <- function(c) {
    # Everything that shapes the *data*, ignoring the config-only axes.
    a <- c$axes
    paste0(
      a$shape[["n_tax"]], "x", a$shape[["n_samp"]], "|",
      a$n_group, "|", format(a$zero_rate, digits = 15), "|", a$struc_zero, "|",
      a$n_cont, "|", a$predictor, "|", a$interact, "|", a$balanced, "|",
      a$signal_frac, "|", c$seed
    )
  }
  groups <- list()
  for (nm in names(cells)) {
    k <- key_of(cells[[nm]])
    groups[[k]] <- c(groups[[k]], nm)
  }
  groups
}
