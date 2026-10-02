#!/usr/bin/env Rscript
# Generate golden fixtures: counts/meta inputs, the full 26-quantity golden
# contract, in both .rds and a canonical little-endian f64 + JSON form so the
# Rust test harness never needs R.
#
#   Rscript scripts/generate_goldens.R [fixture_ids...]
#
# Re-running must reproduce byte-identical canonical output. If it does not, the
# ANCOMBC oracle moved and the compatibility target must be re-frozen.

suppressPackageStartupMessages(library(stats))
here <- Sys.getenv("ANCOMBC_REPO", unset = ".")
source(file.path(here, "reference", "R", "harness.R"))
source(file.path(here, "reference", "R", "fixtures.R"))
source(file.path(here, "reference", "R", "serialize.R"))

# The committed fixture files, read back exactly as `write.table` wrote them.
#
# `read.table(..., row.names = 1)` is wrong for these files, and wrong in a way
# that does not fail: the header holds *every* value-column name, while the data
# rows carry a leading, unlabelled row-name field. So `fx02/meta.tsv` is
#
#     group<TAB>x1
#     sample_0001<TAB>1<TAB>-1.5923679284651
#
# -- two header fields, three data fields. `row.names = 1` consumes `group` as
# the row-name label, leaves `x1` as the only column, and then the values shift by
# one: the group level ends up in a column called `x1` and the covariate ends up
# in a column with no name. The oracle then fails with
# "Estimation failed for the following covariates: x11.67025515850448", which is
# the mangled name and the first value run together.
#
# So these are read by splitting on tabs directly, with the row-name field
# dropped from each row and the header used verbatim.
read_committed_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  nms <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  cells <- lapply(lines[-1], function(l) strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  m <- do.call(rbind, cells)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- nms
  # Numeric where every value parses, character otherwise -- the group factor is
  # labels, and `colClasses = "character"` would make the covariates labels too.
  # Numeric where every value parses, except `group`. `group` holds *labels*; in
  # the committed fixtures those labels are the strings "1" and "2", which parse
  # as numbers, and converting them would make `model.matrix` treat the group
  # factor as a numeric slope instead of a set of contrasts. The generator builds
  # it as a factor, and `write.table` wrote its levels out, so it comes back as
  # text and must stay text.
  for (k in seq_along(df)) {
    if (names(df)[k] == "group") next
    num <- suppressWarnings(as.numeric(df[[k]]))
    if (!anyNA(num)) df[[k]] <- num
  }
  # The row name is the *first* field. (Taking the remainder instead yields the
  # covariate values, and R then rejects them as duplicate row names.)
  rownames(df) <- vapply(lines[-1],
                         function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1],
                         character(1))
  df
}

read_committed_counts <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  # `write.table` quotes nothing here, so the first field is the taxon id and
  # every field after it is a count, in column order.
  body <- lapply(lines[-1], function(l) as.numeric(strsplit(l, "\t", fixed = TRUE)[[1]][-1]))
  m <- do.call(rbind, body)
  # The header holds every column name and the rows carry a leading taxon id, the
  # same convention as `meta.tsv`. Without these, `colnames(counts)` is `NULL` and
  # the harness's `meta_data[colnames(counts), ]` silently yields zero rows.
  colnames(m) <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  rownames(m) <- vapply(lines[-1],
                        function(l) strsplit(l, "\t", fixed = TRUE)[[1]][1],
                        character(1))
  m
}

argv <- commandArgs(trailingOnly = TRUE)
# `--from-committed` recomputes the golden contract from the fixture files already
# in `validation/fixtures/`, without generating the fixture.
#
# This exists because the four committed fixtures cannot be regenerated: they were
# produced before `fixtures.R` gained its `set.seed(spec$seed)`, from the ambient
# session state, and `reference/R/fixtures.R` records that they are left alone
# because replacing them would invalidate goldens that currently pass parity. So
# asking the generator to regenerate `fx01`..`fx04` produces *different counts* and
# the resulting diff is not drift -- it is a different input.
#
# That made a golden-drift check for the four fixtures impossible to express: any
# run that regenerated them would both destroy the fixtures and report a false
# failure. This mode is the check that is actually meaningful, which is the one
# the property we care about: *these goldens are what the pinned oracle produces
# for these committed inputs.*
from_committed <- "--from-committed" %in% argv
ids <- setdiff(argv, "--from-committed")
if (length(ids) == 0L) ids <- as.character(1:4)
ids <- as.integer(ids)

oracle_dir <- file.path(here, "reference", "ANCOMBC")
h <- load_harness(oracle_dir)
gold_dir <- file.path(here, "validation", "golden")
fix_dir <- file.path(here, "validation", "fixtures")
dir.create(gold_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fix_dir, recursive = TRUE, showWarnings = FALSE)

capture_session <- function() {
  si <- utils::capture.output(utils::sessionInfo())
  c(si, "", "ANCOMBC oracle:",
    sprintf("  version 2.15.2"),
    sprintf("  commit   dc4febdf59badb3a8dfe0c767ef2186323c2199a"),
    sprintf("  path     %s", normalizePath(oracle_dir)))
}

for (id in ids) {
  t0 <- proc.time()[["elapsed"]]
  spec <- fixture_spec(id)
  dir <- file.path(fix_dir, sprintf("fx%02d", id))

  if (from_committed) {
    # Read the committed fixture back rather than regenerating it, and require it
    # to be there: a missing input must fail here, not silently fall through to
    # generating a fresh one, which is the mistake this mode exists to prevent.
    if (!file.exists(file.path(dir, "counts.tsv"))) {
      stop("no committed fixture for fx", id, " at ", dir,
           "; rerun without --from-committed to generate it")
    }
    f <- list(
      counts = read_committed_counts(file.path(dir, "counts.tsv")),
      meta = read_committed_meta(file.path(dir, "meta.tsv"))
    )
    ff <- trimws(readLines(file.path(dir, "formula.txt"), warn = FALSE))
  } else {
    f <- gen_fixture(id)
    ff <- fix_formula_for(spec)
    # --- inputs -----------------------------------------------------------
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.table(f$counts, file.path(dir, "counts.tsv"),
                       sep = "\t", quote = FALSE, row.names = TRUE)
    utils::write.table(f$meta, file.path(dir, "meta.tsv"),
                       sep = "\t", quote = FALSE, row.names = TRUE)
    writeLines(ff, file.path(dir, "formula.txt"))
  }

  # --- config -------------------------------------------------------------
  cfg <- list(
    fix_formula = ff,
    p_adj_method = "holm",
    pseudo = 0,
    pseudo_sens = isTRUE(spec$pseudo_sens),
    conservative = isTRUE(spec$conservative),
    prv_cut = 0.10,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "group",
    struc_zero = isTRUE(spec$struc_zero),
    neg_lb = FALSE,
    alpha = 0.05,
    global = spec$n_group >= 3,
    pairwise = spec$n_group >= 3,
    iter_control = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-5, max_iter = 100),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
  )
  write_json(file.path(dir, "config.json"), cfg)

  # --- the mirror self-check must pass before anything is recorded --------
  stopifnot(h$ref_verify_mirror(f$counts, f$meta, fix_formula = ff,
                                p_adj_method = cfg$p_adj_method,
                                alpha = cfg$alpha, prv_cut = cfg$prv_cut,
                                lib_cut = cfg$lib_cut))

  # --- golden contract ----------------------------------------------------
  g <- h$ref_run(f$counts, f$meta, fix_formula = ff,
                 p_adj_method = cfg$p_adj_method, pseudo = cfg$pseudo,
                 pseudo_sens = cfg$pseudo_sens, conservative = cfg$conservative,
                 prv_cut = cfg$prv_cut, lib_cut = cfg$lib_cut,
                 s0_perc = cfg$s0_perc, group = cfg$group,
                 struc_zero = cfg$struc_zero, neg_lb = cfg$neg_lb,
                 alpha = cfg$alpha, global = cfg$global, pairwise = cfg$pairwise,
                 iter_control = cfg$iter_control, em_control = cfg$em_control,
                 mdfdr_control = cfg$mdfdr_control)

  out <- file.path(gold_dir, sprintf("fx%02d", id))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(golden = g, config = cfg, spec = spec, session = capture_session()),
          file.path(out, "golden.rds"))
  manifest <- write_canonical(g, out)
  saveRDS(manifest, file.path(out, "manifest.rds"))

  cat(sprintf("fx%02d: %d taxa x %d samples, %d fix_eff, sens=%s(%s), %d quantities, %.1fs\n",
              id, length(g$taxa_retained), length(g$samples_retained),
              length(g$fix_eff), cfg$pseudo_sens,
              if (cfg$conservative) "conservative" else "nonconservative",
              nrow(manifest), proc.time()[["elapsed"]] - t0))
}

cat("\nwrote goldens to", gold_dir, "\n")
