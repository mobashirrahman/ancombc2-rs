#!/usr/bin/env Rscript
# S06 acceptance: the typed native transport loses nothing.
#
#   Rscript --vanilla r/ANCOMBC/tests/bridge_selftest.R
#
# Run against the *installed* package, so what is exercised is the compiled
# `init.c` and the compiled Rust bridge, not the source files.
#
# Every check here is about the transport and nothing else. No stage computes
# anything: if one of these fails, the boundary is wrong, and no amount of work on
# the numerical side can be trusted until it is fixed. That is the whole reason
# the probe exists separately from the seam.

suppressWarnings(suppressMessages({
  library(ANCOMBC)
}))

failures <- 0L
ok <- function(cond, name, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(cond)) "PASS" else "FAIL", name,
              if (!isTRUE(cond) && nzchar(detail)) paste0("\n         ", detail) else ""))
  if (!isTRUE(cond)) failures <<- failures + 1L
}

get_or <- function(name, default) {
  ns <- asNamespace("ANCOMBC")
  if (exists(name, envir = ns, inherits = FALSE)) get(name, envir = ns) else default
}

rb_probe <- get_or("rb_probe", NULL)
rb_bits <- get_or("rb_bits", NULL)

cat("typed native transport selftest\n")

if (is.null(rb_probe)) {
  cat("  [FAIL] the package exposes no `rb_probe`\n")
  cat("  [FAIL] nothing about the transport can be checked\n")
  quit(save = "no", status = 1L)
}

ok(rb_bits_ok <- !is.null(rb_bits),
   "the package exposes `rb_bits`, which is how transports are compared")
if (!rb_bits_ok) quit(save = "no", status = 1L)

cat("\n## the bridge identifies its target\n")
ok(grepl("dc4febdf59badb3a8dfe0c767ef2186323c2199a", ANCOMBC:::rb_oracle_sha(), fixed = TRUE),
   "the bridge names the pinned commit it is written against",
   ANCOMBC:::rb_oracle_sha())
ok(nzchar(ANCOMBC:::rb_version()), "the bridge reports a version",
   ANCOMBC:::rb_version())

# ---------------------------------------------------------------------------
# A fixture with every awkward value in it.
# ---------------------------------------------------------------------------

n_tax <- 4L; n_samp <- 6L; p <- 3L

# Counts: plain integers, and NA_integer_.
taxon_names <- sprintf("T%02d", seq_len(n_tax))
sample_names <- sprintf("S%02d", seq_len(n_samp))
counts <- matrix(c(10L, 0L, NA_integer_, 5L, 6L, 7L,
                   3L, 3L, 1L, 0L, 2L, 9L,
                   4L, 4L, 1L, 1L, 1L, 8L,
                   12L, 11L, 13L, 14L, 15L, 16L),
                 nrow = n_tax, ncol = n_samp, byrow = TRUE,
                 dimnames = list(taxon_names, sample_names))
counts_int <- counts

# Real counts carrying every non-finite shape R can hold.
counts_real <- matrix(
  c(0, -0, 1, 2, 3, 4,
   NA_real_, NaN, Inf, -Inf, 5, 6,
   1.5, 2.5, 3.5, 4.5, 5.5, 6.5,
   -1, -2, -3, -4, -5, -6),
  nrow = n_tax, ncol = n_samp, byrow = TRUE,
  dimnames = list(taxon_names, sample_names))

# Design, column-major: an intercept, a two-level group contrast, a covariate.
design <- cbind(`(Intercept)` = 1, g2 = c(0, 0, 0, 1, 1, 1), x1 = c(1, 2, 3, 4, 5, 6))
rownames(design) <- sample_names

fix_eff <- colnames(design)
group_index <- c(1L, 1L, 1L, 2L, 2L, 2L)
group_labels <- c("a", "b")
design_complete <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

base <- list(data = counts_int, aggregate_data = counts_int, design = design,
             design_complete = design_complete, group_index = group_index,
             group_labels = group_labels, fix_eff = fix_eff,
             taxon_names = taxon_names, sample_names = sample_names)

cat("\n## integer counts, including NA_integer_\n")
r <- rb_probe("data_int", data = counts_int, aggregate_data = counts_int,
              design = design, design_complete = design_complete,
              group_index = group_index, group_labels = group_labels,
              fix_eff = fix_eff, taxon_names = taxon_names,
              sample_names = sample_names)
ok(identical(r$payload, as.vector(counts_int)),
   "the integer count matrix comes back element for element",
   "a mismatch here means the payload was widened or reordered")
ok(identical(r$payload[is.na(counts_int)], NA_integer_),
   "NA_integer_ comes back as NA_integer_ and not as a number")
ok(r$payload[is.na(counts_int)] == NA_integer_ ||
     identical(intToBits(r$payload[is.na(counts_int)]), intToBits(NA_integer_)),
   "NA_integer_ is still INT_MIN after the round trip")
ok(length(r$payload) == n_tax * n_samp, "the payload length is rows * columns",
   sprintf("%d vs %d", length(r$payload), n_tax * n_samp))
ok(is.integer(r$payload), "the payload came back as an integer vector, not a double")

cat("\n## real counts: NA_real_, NaN, Inf, -Inf and -0.0\n")
r <- rb_probe("data_real", data = counts_real, aggregate_data = counts_real,
              design = design, design_complete = design_complete,
              group_index = group_index, group_labels = group_labels,
              fix_eff = fix_eff, taxon_names = taxon_names,
              sample_names = sample_names)
ok(identical(rb_bits(r$payload), rb_bits(as.vector(counts_real))),
   "every one of the 64 bits of every cell survives",
   "NA_real_ and NaN are different NaNs and -0.0 is not 0.0; identical() cannot see any of that")
pos0 <- which(!is.na(counts_real) & !is.nan(counts_real) &
                rb_bits(as.vector(counts_real)) == rb_bits(0))
neg0 <- which(rb_bits(as.vector(counts_real)) == rb_bits(-0))
ok(length(pos0) == 1L && identical(rb_bits(r$payload)[pos0], rb_bits(0)),
   "+0.0 keeps its own bits")
ok(length(neg0) == 1L &&
     identical(rb_bits(r$payload)[neg0], rb_bits(-0)),
   "-0.0 keeps its own bits and is not turned into +0.0",
   sprintf("sprintf and identical() both report these as '0'; the fixture has %d",
           length(neg0)))
ok(length(pos0) == 1L && length(neg0) == 1L &&
     !identical(rb_bits(r$payload)[pos0], rb_bits(r$payload)[neg0]),
   "-0.0 and +0.0 are distinguishable after the round trip")
na_idx <- which(is.na(counts_real) & !is.nan(counts_real))
nan_idx <- which(is.nan(counts_real))
ok(length(na_idx) == 1L &&
     identical(rb_bits(r$payload)[na_idx], rb_bits(NA_real_)),
   "NA_real_ arrives with R's own NaN payload",
   sprintf("expected %s, got %s", rb_bits(NA_real_), rb_bits(r$payload)[na_idx]))
ok(length(nan_idx) == 1L &&
     !identical(rb_bits(r$payload)[nan_idx], rb_bits(NA_real_)),
   "a real NaN is not turned into NA_real_",
   "this is the distinction a JSON bridge loses")
inf_idx <- which(is.infinite(counts_real))
ok(length(inf_idx) == 2L &&
     identical(rb_bits(r$payload)[inf_idx], rb_bits(as.vector(counts_real)[inf_idx])),
   "+Inf and -Inf both keep their own bits")
ok(is.double(r$payload), "the payload came back as a double vector")

cat("\n## dimensions, names and masks\n")
ok(identical(r$dim, c(n_tax, n_samp)), "the count matrix's shape is reported back")
ok(identical(r$taxon_names, taxon_names), "taxon names survive intact")
ok(identical(r$sample_names, sample_names), "sample names survive intact")
ok(identical(r$fix_eff, fix_eff), "coefficient names survive intact")
ok(identical(r$group_labels, group_labels), "group labels survive intact")

rd <- rb_probe("design", data = counts_int, aggregate_data = counts_int,
               design = design, design_complete = design_complete,
               group_index = group_index, group_labels = group_labels,
               fix_eff = fix_eff, taxon_names = taxon_names,
               sample_names = sample_names)
ok(identical(rd$payload, as.vector(design)),
   "the design comes back element for element, column-major",
   "a row-major round trip would give a different but plausible design")

# An incomplete design row: complete.cases(x) is FALSE for it, and the mask must
# carry that.
dc2 <- design; dc2[3, 3] <- NA_real_
rcomp <- rb_probe("design_complete", data = counts_int, aggregate_data = counts_int,
                  design = dc2, design_complete = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
                  group_index = group_index, group_labels = group_labels,
                  fix_eff = fix_eff, taxon_names = taxon_names,
                  sample_names = sample_names)
ok(is.raw(rcomp$payload), "the mask came back as a raw vector")
ok(identical(as.integer(rcomp$payload), c(1L, 1L, 0L, 1L, 1L, 1L)),
   "the completeness mask survives element for element")

cat("\n## the group index\n")
rg <- rb_probe("group_index", data = counts_int, aggregate_data = counts_int,
               design = design, design_complete = design_complete,
               group_index = group_index, group_labels = group_labels,
               fix_eff = fix_eff, taxon_names = taxon_names,
               sample_names = sample_names)
ok(identical(rg$payload, group_index), "the 1-based group index survives")

rabsent <- rb_probe("group_index", data = counts_int, aggregate_data = counts_int,
                    design = design, design_complete = design_complete,
                    group_index = c(0L, 1L, 1L, 2L, 2L, 2L),
                    group_labels = group_labels,
                    fix_eff = fix_eff, taxon_names = taxon_names,
                    sample_names = sample_names)
ok(identical(rabsent$payload[1], 0L),
   "0 means 'no group' and is not confused with a level index")

cat("\n## the aggregate table\n")
agg <- counts_int + 100L
ra <- rb_probe("aggregate", data = counts_int, aggregate_data = agg,
               design = design, design_complete = design_complete,
               group_index = group_index, group_labels = group_labels,
               fix_eff = fix_eff, taxon_names = taxon_names,
               sample_names = sample_names)
ok(identical(ra$payload, as.vector(agg)),
   "the aggregate table survives, in its own integer type",
   sprintf("payload type %s, input type %s", class(ra$payload)[1], class(agg)[1]))

cat("\n## lifetime: the returned buffer is R's, not a view of the input\n")
# Mutating the input after the call must not change what came back. If the bridge
# had returned a view, this would change it.
snap <- rb_probe("data_int", data = counts_int, aggregate_data = counts_int,
                 design = design, design_complete = design_complete,
                 group_index = group_index, group_labels = group_labels,
                 fix_eff = fix_eff, taxon_names = taxon_names,
                 sample_names = sample_names)$payload
scratch <- counts_int
scratch[1, 1] <- 999L
ok(identical(snap, as.vector(counts_int)),
   "the returned buffer does not alias the input, so a later change to the input cannot alter it")
ok(counts_int[1, 1] == 10L, "and the change was made to a copy, not the fixture")

cat("\n## rejections: each must be an R error naming the argument\n")
expect_error <- function(expr, pattern, name) {
  e <- tryCatch({ force(expr); NULL }, error = function(e) e)
  ok(!is.null(e) && grepl(pattern, conditionMessage(e), fixed = FALSE),
     name,
     if (is.null(e)) "no error was raised" else conditionMessage(e))
}

expect_error(rb_probe("data_int", data = counts_real, aggregate_data = counts_real,
                      design = design, design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "integer matrix", "a real count table is refused by the integer probe")

expect_error(rb_probe("design_complete", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = c(TRUE, TRUE, NA, TRUE, TRUE, TRUE),
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "design_complete", "an NA in the mask is refused")

expect_error(rb_probe("design_complete", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = c(TRUE, TRUE),
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "length", "a mask of the wrong length is refused")

expect_error(rb_probe("group_index", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = design_complete,
                      group_index = c(1L, 2L), group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "group_index", "a group index of the wrong length is refused")

expect_error(rb_probe("design", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = design_complete,
                      group_index = c(1L, 1L, 1L, 2L, 2L, 9L),
                      group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "group index", "a group index past the last level is refused")

expect_error(rb_probe("design", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = c("only-one"),
                      taxon_names = taxon_names, sample_names = sample_names),
             "fix_eff", "a coefficient-name count that disagrees with the design is refused")

expect_error(rb_probe("design", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff,
                      taxon_names = taxon_names[1:2],
                      sample_names = sample_names),
             "taxon_names", "a taxon-name count that disagrees with the table is refused")

expect_error(rb_probe("design", data = counts_int,
                      aggregate_data = counts_int[1:2, , drop = FALSE],
                      design = design, design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "shapes disagree|aggregate_data", "an aggregate table of the wrong shape is refused")

expect_error(rb_probe("design", data = counts_int, aggregate_data = counts_int,
                      design = rbind(design, design), design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "shapes disagree|design", "a design with the wrong number of rows is refused")

expect_error(rb_probe("design", data = counts_int, aggregate_data = counts_int,
                      design = design, design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names,
                      controls = list(s0_perc = NaN)),
             "s0_perc|must be finite", "a non-finite control is refused, named")

expect_error(rb_probe("data_int", data = "not a matrix", aggregate_data = counts_int,
                      design = design, design_complete = design_complete,
                      group_index = group_index, group_labels = group_labels,
                      fix_eff = fix_eff, taxon_names = taxon_names,
                      sample_names = sample_names),
             "numeric matrix", "a character count table is refused by type")

cat(sprintf("\n  %d check(s) failed\n", failures))
if (failures > 0L) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
