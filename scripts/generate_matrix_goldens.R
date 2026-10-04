#!/usr/bin/env Rscript
# Generate the golden fixture matrix of PLAN.md section 5.5: every cell declared
# in reference/R/fixture_matrix.R, with its inputs and the full golden contract
# in both .rds and the canonical little-endian f64 + JSON form, so the Rust
# harness never needs R.
#
#   Rscript scripts/generate_matrix_goldens.R [cell_name...]
#
# Cells that differ only in configuration share one input directory; the golden
# is written per cell, because the configuration changes the answer.
#
# Re-running must reproduce byte-identical canonical output. `make matrix-drift`
# checks that, and a difference means the oracle moved.

suppressPackageStartupMessages(library(stats))
here <- Sys.getenv("ANCOMBC_REPO", unset = ".")
source(file.path(here, "reference", "R", "harness.R"))
source(file.path(here, "reference", "R", "fixtures.R"))
source(file.path(here, "reference", "R", "fixture_matrix.R"))
source(file.path(here, "reference", "R", "serialize.R"))

# The *full* matrix, always, before any subsetting: `cells.json` is the
# declaration the Rust coverage test reads, so a partial regeneration must not
# shrink it. Regenerating one cell to debug it would otherwise leave a one-cell
# manifest behind, and the coverage test would then report the matrix as missing
# nine axis levels rather than reporting the real problem.
all_cells <- fixture_matrix_cells()

wanted <- commandArgs(trailingOnly = TRUE)
oracle_dir <- file.path(here, "reference", "ANCOMBC")
h <- load_harness(oracle_dir)

cells <- all_cells
if (length(wanted) > 0L) {
  unknown <- setdiff(wanted, names(cells))
  if (length(unknown)) stop("unknown cell(s): ", paste(unknown, collapse = ", "))
  cells <- cells[wanted]
}
groups <- fixture_matrix_input_groups(cells)

fix_root <- file.path(here, "validation", "matrix", "fixtures")
gold_root <- file.path(here, "validation", "matrix", "golden")
dir.create(fix_root, recursive = TRUE, showWarnings = FALSE)
dir.create(gold_root, recursive = TRUE, showWarnings = FALSE)

capture_session <- function() {
  c(utils::capture.output(utils::sessionInfo()), "", "ANCOMBC oracle:",
    "  version 2.15.2",
    "  commit   dc4febdf59badb3a8dfe0c767ef2186323c2199a",
    sprintf("  path     %s", normalizePath(oracle_dir)))
}

# --- inputs, one directory per input group ----------------------------------
# The directory is named for the *first* cell in the group, so a cell's inputs
# are where a reader would look for them, and `cells.json` records the mapping
# so nothing depends on that convention holding.
cell_inputs <- list()
for (k in seq_along(groups)) {
  names_in_group <- groups[[k]]
  first <- names_in_group[[1]]
  spec <- spec_of(cells[[first]])
  f <- gen_fixture_from_spec(spec)
  ff <- fix_formula_for(spec)
  dir <- file.path(fix_root, first)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(f$counts, file.path(dir, "counts.tsv"),
                     sep = "\t", quote = FALSE, row.names = TRUE)
  utils::write.table(f$meta, file.path(dir, "meta.tsv"),
                     sep = "\t", quote = FALSE, row.names = TRUE)
  writeLines(ff, file.path(dir, "formula.txt"))
  writeLines(c(
    sprintf("cell=%s", first),
    sprintf("n_tax=%d", nrow(f$counts)),
    sprintf("n_samp=%d", ncol(f$counts)),
    sprintf("n_group=%d", spec$n_group),
    sprintf("zero_rate=%g", spec$zero_rate),
    sprintf("struc_zero=%s", spec$struc_zero),
    sprintf("n_cont=%d", spec$n_cont),
    sprintf("predictor=%s", spec$predictor),
    sprintf("interact=%s", spec$interact),
    sprintf("balanced=%s", spec$balanced),
    sprintf("seed=%d", spec$seed),
    sprintf("zero_fraction=%.6f", mean(f$counts == 0)),
    sprintf("group_sizes=%s", paste(sprintf("%d", table(f$meta$group)), collapse = ","))
  ), file.path(dir, "spec.txt"))
  for (nm in names_in_group) cell_inputs[[nm]] <- first
  cat(sprintf("inputs %-30s %5d x %-5d %-14s zeros=%.3f\n", first, nrow(f$counts),
              ncol(f$counts), ff, mean(f$counts == 0)))
}

# --- the manifest, so the Rust side needs no convention knowledge -----------
# The `inputs` directory of a cell is the input *group* it belongs to, which is a
# property of the full matrix: a cell generated on its own still shares a table
# with its sweep-mates, and the Rust test checks that they do.
all_groups <- fixture_matrix_input_groups(all_cells)
all_inputs <- list()
for (k in seq_along(all_groups)) {
  first <- all_groups[[k]][[1]]
  for (nm in all_groups[[k]]) all_inputs[[nm]] <- first
}
manifest <- list(
  target = "ancombc2-rs v0.1 equivalent to ANCOMBC 2.15.2 at dc4febdf59badb3a8dfe0c767ef2186323c2199a, R >= 4.5.0, seed 42",
  generator = "scripts/generate_matrix_goldens.R",
  n_cells = length(all_cells),
  # Which cells this invocation actually wrote. The Rust test reports a missing
  # golden against this list, so a partial run is visible rather than looking
  # like a smaller matrix.
  generated = names(cells),
  cells = lapply(names(all_cells), function(nm) {
    a <- all_cells[[nm]]$axes
    list(name = nm, inputs = all_inputs[[nm]], seed = all_cells[[nm]]$seed,
         rank_deficient_cap = all_cells[[nm]]$rank_deficient_cap,
         n_tax = a$shape[["n_tax"]], n_samp = a$shape[["n_samp"]],
         n_group = a$n_group, zero_rate = a$zero_rate,
         struc_zero = a$struc_zero, n_cont = a$n_cont,
         predictor = a$predictor, interact = a$interact,
         balanced = a$balanced, pseudo = a$pseudo,
         p_adj_method = a$p_adj_method,
         pseudo_sens = a$pseudo_sens, conservative = a$conservative)
  })
)
write_json(file.path(here, "validation", "matrix", "cells.json"), manifest)

# --- goldens, one per cell --------------------------------------------------
for (nm in names(cells)) {
  t0 <- proc.time()[["elapsed"]]
  a <- cells[[nm]]$axes
  fdir <- file.path(fix_root, cell_inputs[[nm]])
  # `as.matrix` on the counts is load-bearing: `read.delim` returns a
  # data.frame, and the oracle indexes it with `is.infinite`, which has no method
  # for a list column. The committed generator never hits this because it passes
  # the in-memory matrix straight through; this one re-reads the file, so it has
  # to convert. The metadata stays a data.frame, because `data_sanity_check`
  # coerces the group column to a factor itself.
  f <- list(counts = as.matrix(utils::read.delim(file.path(fdir, "counts.tsv"),
                                                 row.names = 1, check.names = FALSE)),
           meta = utils::read.delim(file.path(fdir, "meta.tsv"), row.names = 1,
                                    check.names = FALSE))
  # The group column must be a *factor* with its levels in the order the file
  # implies, or `model.matrix` builds one numeric column named `group` instead of
  # `group2`, `group3`, ... and every multi-group result is wrong: the global
  # test loses its numerator, and the pairwise test then calls
  # `combn(x, 2)` on a length-1 vector and dies with "n < m".
  #
  # The committed generator never hits this because it passes `f$meta` straight
  # from `gen_fixture`, where `group` is already a factor. This script re-reads
  # the file to prove the on-disk inputs are what the goldens were computed from,
  # and `read.delim` types the column as integer.
  f$meta$group <- factor(f$meta$group, levels = sort(unique(f$meta$group)))
  ff <- readLines(file.path(fdir, "formula.txt"), warn = FALSE)
  ff <- paste(ff, collapse = " ")

  cfg <- list(
    fix_formula = ff,
    p_adj_method = a$p_adj_method,
    pseudo = a$pseudo,
    pseudo_sens = a$pseudo_sens,
    conservative = a$conservative,
    prv_cut = 0.10,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "group",
    struc_zero = a$struc_zero,
    neg_lb = FALSE,
    alpha = 0.05,
    # The reference's own rule: the multi-group tests need somewhere to compare.
    global = a$n_group >= 3,
    pairwise = a$n_group >= 3,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 100),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
  )
  write_json(file.path(fdir, sprintf("config-%s.json", nm)), cfg)

  # The mirror self-check runs before anything is recorded, exactly as
  # generate_goldens.R does for the four committed fixtures.
  ok <- h$ref_verify_mirror(f$counts, f$meta, fix_formula = ff,
                            p_adj_method = cfg$p_adj_method, alpha = cfg$alpha,
                            prv_cut = cfg$prv_cut, lib_cut = cfg$lib_cut)
  if (!isTRUE(ok)) stop(sprintf("%s: the oracle's own mirror check failed", nm))

  g <- h$ref_run(f$counts, f$meta, fix_formula = ff,
                 p_adj_method = cfg$p_adj_method, pseudo = cfg$pseudo,
                 pseudo_sens = cfg$pseudo_sens, conservative = cfg$conservative,
                 prv_cut = cfg$prv_cut, lib_cut = cfg$lib_cut,
                 s0_perc = cfg$s0_perc, group = cfg$group,
                 struc_zero = cfg$struc_zero, neg_lb = cfg$neg_lb,
                 alpha = cfg$alpha, global = cfg$global, pairwise = cfg$pairwise,
                 iter_control = cfg$iter_control, em_control = cfg$em_control,
                 mdfdr_control = cfg$mdfdr_control)

  out <- file.path(gold_root, nm)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(golden = g, config = cfg, session = capture_session()),
          file.path(out, "golden.rds"))
  m <- write_canonical(g, out)
  saveRDS(m, file.path(out, "manifest.rds"))
  cat(sprintf("%-30s %5d x %-5d %3d fix_eff sens=%-5s %d quantities %.1fs\n",
              nm, length(g$taxa_retained), length(g$samples_retained),
              length(g$fix_eff), cfg$pseudo_sens, nrow(m),
              proc.time()[["elapsed"]] - t0))
}

cat(sprintf("\n%d cell(s) across %d input group(s) -> %s\n",
            length(cells), length(groups), gold_root))
