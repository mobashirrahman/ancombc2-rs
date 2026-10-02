# Input coercion: matrices, data.frames, phyloseq and (Tree)SummarizedExperiment.
#
# `phyloseq`, `microbiome` and `mia` are all Suggests. The coercion path has to be
# *tested* without them, or it is only tested on the machines that happen to have
# the microbiome stack installed -- which is not this one. So these tests build
# objects that carry the right S4 classes and slots, with no dependency on the
# packages that would normally provide them, and assert that the wrapper reads the
# same feature table and metadata out of them that it reads from a plain matrix.
#
# The assertion that matters is not "it does not error" but "it produces the same
# analysis as passing the matrix", because that is the property the reference's own
# coercion has: `microbiome::abundances()` is just a different route to the same
# counts.

skip_without_core <- function() {
  testthat::skip_if_not(ancombc2_available(), "ancombc2-rs core not linked")
}

# A stand-in for a phyloseq object, built with base R only.
#
# `phyloseq` stores its feature table in `@otu_table` and its sample metadata in
# `@sam_data`. Declaring an S4 class of the same name and slot layout gives an
# object that `inherits(x, "phyloseq")` accepts and whose slots are reachable with
# `attr()`, which is exactly the last tier of the wrapper's fallback chain -- so
# this test exercises the real code path on a machine where `phyloseq`,
# `microbiome` and `mia` are all absent.
#
# The class is only declared when the real package is not installed, so this file
# does not shadow the genuine thing on a machine that has it.
fake_phyloseq <- function(counts, meta) {
  # `getClassDef()` *throws* for an undefined class rather than returning NULL,
  # so the guard is a tryCatch. `isVirtualClass()` is not a substitute: it answers
  # a different question and is what made this fail.
  if (is.null(tryCatch(methods::getClassDef("phyloseq"),
                       error = function(e) NULL))) {
    methods::setClass(
      "phyloseq",
      representation(otu_table = "ANY", sam_data = "ANY")
    )
  }
  otu <- counts
  class(otu) <- c("otu_table", "matrix", "array")
  sd <- meta
  class(sd) <- c("sample_data", "data.frame")
  methods::new("phyloseq", otu_table = otu, sam_data = sd)
}

# 12 samples, 14 taxa, two groups of six.
#
# Sized away from the edge on purpose. With `~ group + sex` the design has three
# columns, and `.iter_mle` groups taxa by missingness pattern and fits each group
# on the samples it actually observed, so a table with only three or four samples
# can hand a group `n < p` and trip the QR's `n >= p` assertion -- a property of
# the fixture, not of the coercion under test.
tiny_counts <- function() {
  set.seed(42)
  m <- matrix(rnbinom(14 * 12, mu = 60, size = 4), nrow = 14)
  # A taxon that is a structural zero in group `a`, so the coercion is also
  # exercised on a table that has zeros rather than only dense counts.
  m[1, 1:6] <- 0
  m[2, 7:12] <- 0
  dimnames(m) <- list(paste0("taxon", 1:14), paste0("sample", 1:12))
  m
}

tiny_meta <- function() {
  data.frame(
    group = factor(rep(c("a", "b"), each = 6)),
    sex = factor(rep(c("f", "m"), 6)),
    row.names = paste0("sample", 1:12)
  )
}

test_that("a phyloseq object is read the same as the matrix it holds", {
  skip_without_core()
  counts <- tiny_counts()
  meta <- tiny_meta()

  direct <- ancombc2(data = counts, meta_data = meta,
                     fix_formula = ~ group + sex, group = "group")
  via_pseq <- ancombc2(data = fake_phyloseq(counts, meta),
                       fix_formula = ~ group + sex, group = "group")

  # Same retained taxa, same order, and every estimate to the last bit: the
  # coercion changed how the table was *read*, nothing else.
  expect_equal(rownames(direct$res), rownames(via_pseq$res))
  expect_equal(direct$res$taxa, via_pseq$res$taxa)

  # Guard against vacuous equality. The result's estimates live at the *top
  # level* of the object (`$beta`, `$se`, ...), not in `$res`, which is the
  # reference's wide `lfc_*`/`se_*`/`p_*` table. An earlier version of this test
  # compared `res$beta` against `res$beta` -- two `NULL`s -- and passed without
  # comparing anything at all. So the lengths are asserted before the equality.
  for (nm in c("beta", "se", "p", "q", "w", "theta", "delta_em", "delta_wls")) {
    expect_gt(length(direct[[nm]]), 0L)
    expect_gt(length(via_pseq[[nm]]), 0L)
  }
  for (nm in c("beta", "se", "p", "q", "w", "theta", "delta_em", "delta_wls")) {
    expect_equal(direct[[nm]], via_pseq[[nm]], info = paste("quantity", nm))
  }
  expect_equal(length(direct$diff_abn), 42L)
  expect_identical(direct$diff_abn, via_pseq$diff_abn)
  # And the wide table itself, which is what a user actually looks at.
  expect_equal(direct$res, via_pseq$res)
})

test_that("metadata supplied alongside a phyloseq object wins", {
  skip_without_core()
  counts <- tiny_counts()
  meta <- tiny_meta()
  # A different metadata frame, so the test distinguishes "read from the object"
  # from "used what was passed".
  meta2 <- meta
  meta2$sex <- factor(rep(c("f", "f", "m", "m"), 3))

  fml <- ~ group + sex
  # Passing `meta` explicitly must give exactly the analysis that `meta` gives,
  # and omitting it must give exactly the analysis `meta2` gives -- compared
  # against the plain-matrix route, which is the reference's own path. That
  # identifies *which* metadata was used, rather than merely asserting the two
  # happen to differ.
  explicit <- ancombc2(data = fake_phyloseq(counts, meta2), meta_data = meta,
                       fix_formula = fml, group = "group")
  from_object <- ancombc2(data = fake_phyloseq(counts, meta2),
                          fix_formula = fml, group = "group")
  via_meta <- ancombc2(data = counts, meta_data = meta,
                       fix_formula = fml, group = "group")
  via_meta2 <- ancombc2(data = counts, meta_data = meta2,
                        fix_formula = fml, group = "group")

  expect_equal(explicit$beta, via_meta$beta)
  expect_equal(from_object$beta, via_meta2$beta)
  expect_equal(explicit$res, via_meta$res)
  # ...and the two metadata frames really are different analyses, so the test is
  # not vacuous: `sex` is nested within `group` here and crossed with it above.
  expect_gt(length(via_meta$beta), 0L)
  expect_true(any(abs(unlist(via_meta$beta) - unlist(via_meta2$beta)) > 1e-6))
})

test_that("taxa_are_rows = FALSE transposes a plain matrix, as the reference does", {
  skip_without_core()
  counts <- tiny_counts()
  meta <- tiny_meta()
  rows <- ancombc2(data = counts, meta_data = meta,
                   fix_formula = ~ group + sex, group = "group")
  cols <- ancombc2(data = t(counts), meta_data = meta,
                   fix_formula = ~ group + sex, group = "group",
                   taxa_are_rows = FALSE)
  expect_equal(rows$res$taxa, cols$res$taxa)
  expect_equal(rows$res, cols$res)
  expect_equal(rows$beta, cols$beta)
})

test_that("sample names that do not match the metadata are rejected", {
  skip_without_core()
  counts <- tiny_counts()
  meta <- tiny_meta()
  meta <- meta[1:3, , drop = FALSE]
  expect_error(
    ancombc2(data = counts, meta_data = meta,
             fix_formula = ~ group, group = "group"),
    "Sample names do not match"
  )
})

test_that("a matrix with no metadata says the reference's message", {
  skip_without_core()
  expect_error(
    ancombc2(data = tiny_counts(), fix_formula = ~ group, group = "group"),
    "Missing sample metadata"
  )
})

test_that("tax_level is rejected rather than silently ignored", {
  # Aggregation is a `microbiome`/`mia` operation. Accepting the argument and
  # ignoring it would return results for a different model than the one asked
  # for, which is the same class of defect as the `*`-in-a-formula bug this
  # project already had to fix.
  counts <- tiny_counts()
  meta <- tiny_meta()
  expect_error(
    ancombc2(data = fake_phyloseq(counts, meta), fix_formula = ~ group,
             group = "group", tax_level = "Genus"),
    "not implemented"
  )
})

test_that("an unsupported input type names the types that are supported", {
  skip_without_core()
  expect_error(
    ancombc2(data = 1:10, meta_data = tiny_meta(), fix_formula = ~ group,
             group = "group"),
    "phyloseq"
  )
})