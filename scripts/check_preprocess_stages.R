#!/usr/bin/env Rscript
# S08: are the preprocessing stage arrays exact?
#
# The claim under test is narrow and stated exactly: *for the same input, the
# preprocessing stages the bridge computes are bit-identical to the stages R's own
# expressions produce.* Nothing here claims the bridge has been shown to compute
# anything the original does not -- that is S09 onwards. What it establishes is that
# when the arithmetic does run, it lands on the same bits at every step.
#
# R's side is transcribed from the pinned original, line by line:
#
#   .data_core        reference/ANCOMBC/R/ancombc_prep.R:2-44
#   .get_struc_zero   reference/ANCOMBC/R/ancombc_prep.R:47-84
#   O/o/y             reference/ANCOMBC/R/ancombc_prep.R:110-130, 212-217
#
# Those expressions are *run here*, not imported from the installed original. That is
# deliberate: importing `.ancombc2_core` would let a divergence hide inside the
# package under test, and the whole point is to compare the bridge against
# transcribed semantics.
#
# Usage: Rscript scripts/check_preprocess_stages.R [--fixtures DIR] [--json PATH]

suppressWarnings(suppressMessages(library(ANCOMBC)))

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- which(args == name)
  if (length(i) == 0L) default else args[i[1L] + 1L]
}

fails <- 0L
checks <- 0L
# An environment rather than a list: `<<-` from inside the loop depends on where
# the loop body is evaluated, and an environment is unambiguous about it.
stage_report <- new.env(parent = emptyenv())

ok <- function(label, cond, why = "") {
  checks <<- checks + 1L
  if (isTRUE(cond)) {
    cat(sprintf("ok    %s\n", label))
  } else {
    fails <<- fails + 1L
    cat(sprintf("FAIL  %s%s\n", label, if (nzchar(why)) paste0(" -- ", why) else ""))
  }
}

# `serialize()` of each, byte for byte. Not `all.equal`, not a tolerance: a
# `NA_real_` and a `NaN` are `identical()`-equal to nothing and differ in the last
# two bytes, and that difference is the one this whole file exists to find.
same_bytes <- function(x, y) {
  identical(serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3),
            serialize(y, NULL, ascii = FALSE, xdr = TRUE, version = 3))
}

# Where a stage differs: the first few differing positions with both bit patterns,
# because "differs" alone is not actionable when the array is 5000 x 10.
where_differs <- function(got, want, k = 4L) {
  if (identical(dim(got), dim(want)) && length(got) == length(want)) {
    # The bit pattern of one element, as 16 hex digits. `readBin` on the raw bytes
    # rather than a format string, because `%a` prints `0x1.8p+1023` for both
    # `NA_real_` and a computed `NaN` and this file needs to tell them apart.
    bits1 <- function(z) {
      if (is.integer(z)) {
        sprintf("%08x", as.numeric(readBin(writeBin(z, raw()), "integer", 4L,
                                          signed = FALSE)))
      } else if (is.double(z)) {
        # Four unsigned 16-bit reads rather than two signed 32-bit ones: `signed =
        # FALSE` is rejected for a 4-byte read, and `as.integer(x %% 2^32)` is NA for
        # anything above 2^31, which is most of the exponent half. Little-endian, so
        # the first two bytes are the least significant -- printed last.
        r8 <- writeBin(z, raw())
        # `readBin` has no offset argument -- it always reads from the start of what
        # it is given -- so each half has to be sliced out first. Reading all four from
        # `r8` gives the same 16 bits four times, which is how this produced
        # `0xab73ab73ab73ab73`.
        w <- vapply(0:3, function(k) as.numeric(
          readBin(r8[(2L * k + 1L):(2L * k + 2L)], "integer", n = 1L, size = 2L,
                  signed = FALSE)), numeric(1L))
        sprintf("%04x%04x%04x%04x", w[4], w[3], w[2], w[1])
      } else {
        format(z)
      }
    }
    bits <- function(v) vapply(v, bits1, character(1L), USE.NAMES = FALSE)
    gb <- bits(got); wb <- bits(want)
    at <- which(gb != wb)
    if (length(at) == 0L) {
      # Every value agrees to the bit and the shapes match, so the difference is an
      # attribute -- a name, a dim, a class, a type. `str()` of both is the shortest
      # way to see which.
      return(paste0("attributes differ:\n        got  ",
                    paste(capture.output(str(got)), collapse = "\n        got  "),
                    "\n        want ",
                    paste(capture.output(str(want)), collapse = "\n        want ")))
    }
    head(paste0(at, ": got 0x", gb[at], " want 0x", wb[at]), k)
  } else {
    paste0("shape: got ", paste(dim(got), collapse = "x"),
           " want ", paste(dim(want), collapse = "x"),
           " (lengths ", length(got), " vs ", length(want), ")")
  }
}

# ---------------------------------------------------------------------------
# The reference's expressions, transcribed
# ---------------------------------------------------------------------------

# .data_core (ancombc_prep.R:2-44), with `tax_keep`/`samp_keep` optional.
ref_data_core <- function(data, prv_cut, lib_cut, tax_keep = NULL, samp_keep = NULL) {
  feature_table <- data
  prevalence_fun <- function(x) {
    x <- as.matrix(x)
    rowSums(x != 0, na.rm = TRUE) / rowSums(!is.na(x))
  }
  if (is.null(tax_keep)) {
    prevalence <- prevalence_fun(feature_table)
    tax_keep <- which(prevalence >= prv_cut)
  } else if (length(tax_keep) == 0) {
    stop("All taxa contain structural zeros", call. = FALSE)
  } else {
    feature_table <- feature_table[tax_keep, , drop = FALSE]
    prevalence <- prevalence_fun(feature_table)
    tax_keep <- which(prevalence >= prv_cut)
  }
  if (length(tax_keep) > 0) {
    feature_table <- feature_table[tax_keep, , drop = FALSE]
  } else {
    stop("No taxa remain under the current cutoff", call. = FALSE)
  }
  if (is.null(samp_keep)) {
    lib_size <- colSums(feature_table, na.rm = TRUE)
    samp_keep <- which(lib_size >= lib_cut)
  }
  if (length(samp_keep) > 0) {
    feature_table <- feature_table[, samp_keep, drop = FALSE]
  } else {
    stop("No samples remain under the current cutoff", call. = FALSE)
  }
  list(feature_table = feature_table, tax_keep = tax_keep, samp_keep = samp_keep,
       prevalence = prevalence_fun(feature_table))
}

# The library size the reference actually computes: `colSums` of the *taxon*-filtered
# table, before the sample axis is subsetted (ancombc_prep.R:29-31).
ref_lib_size <- function(data, tax_keep) colSums(data[tax_keep, , drop = FALSE], na.rm = TRUE)

# .get_struc_zero (ancombc_prep.R:47-84), split into its three arrays.
ref_struc_zero <- function(data, meta_data, group, neg_lb) {
  feature_table <- data
  group_data <- factor(meta_data[, group])
  feature_mat <- as.matrix(feature_table)
  present_table <- feature_mat
  present_table[is.na(present_table)] <- 0
  present_table[present_table != 0] <- 1
  n_tax <- nrow(feature_table)
  n_group <- nlevels(group_data)
  grp_int <- as.integer(group_data)
  grp_ok <- !is.na(grp_int)
  G <- matrix(0, nrow = ncol(feature_table), ncol = n_group)
  G[cbind(which(grp_ok), grp_int[grp_ok])] <- 1
  n_g <- colSums(G)
  p_hat <- sweep(present_table %*% G, 2, n_g, "/")
  samp_size <- ((!is.na(feature_mat)) * 1) %*% G
  p_hat_lo <- p_hat - 1.96 * sqrt(p_hat * (1 - p_hat) / samp_size)
  zero_ind <- (p_hat == 0)
  if (neg_lb) zero_ind[p_hat_lo <= 0] <- TRUE
  list(p_hat = p_hat, samp_size = samp_size, p_hat_lo = p_hat_lo,
       zero_ind = zero_ind,
       tax_keep = which(apply(zero_ind, 1, function(x) all(x == FALSE))))
}

# `.data_core`'s local `prevalence`, which it computes and does not return
# (ancombc_prep.R:9). Kept here so the bridge's `prevalence1` has something to be
# compared against.
prevalence_of <- function(x) {
  x <- as.matrix(x)
  rowSums(x != 0, na.rm = TRUE) / rowSums(!is.na(x))
}

# Step 1 of .ancombc2_core (ancombc_prep.R:110-130): `O1`, `o1`, `y1`, and the row
# means as their own stage -- the reduction is where R's accumulator width shows.
ref_log_center <- function(counts, pseudo) {
  O <- counts + pseudo
  o <- log(O)
  o[is.infinite(o)] <- NA
  means <- rowMeans(o, na.rm = TRUE)
  y <- o - means
  list(O = O, log = o, means = means, y = y)
}

# ---------------------------------------------------------------------------
# The fixtures
# ---------------------------------------------------------------------------
#
# S08's acceptance names four: tiny, sparse, missing-value and boundary. Each is
# built here rather than read from a file, so the shape of the test is visible in
# one place and a fixture cannot drift from the claim it is meant to support.

make_fixture <- function(name, data, group = NULL, struc_zero = FALSE,
                         neg_lb = FALSE, prv_cut = 0, lib_cut = 0, pseudo = 0,
                         aggregate = NULL) {
  data <- as.matrix(data)
  storage.mode(data) <- "double"
  if (is.null(aggregate)) aggregate <- data else {
    aggregate <- as.matrix(aggregate)
    storage.mode(aggregate) <- "double"
  }
  n_tax <- nrow(data)
  n_samp <- ncol(data)
  # The grouping column, built here rather than in a second pass. Two alternating
  # levels, so a group boundary lands in the middle of every table, and both levels
  # are declared explicitly so `factor()` cannot drop an unused one.
  meta <- NULL
  if (!is.null(group)) {
    labels <- rep(c("a", "b"), length.out = n_samp)
    meta <- data.frame(.group = factor(labels, levels = c("a", "b")))
    group <- ".group"
  }
  list(name = name, data = data, aggregate = aggregate, meta_data = meta,
       group = group, struc_zero = struc_zero, neg_lb = neg_lb,
       prv_cut = prv_cut, lib_cut = lib_cut, pseudo = pseudo,
       n_tax = n_tax, n_samp = n_samp)
}

fixtures <- list(
  make_fixture("tiny", rbind(c(4, 0, 9, 1), c(0, 7, 0, 3), c(2, 2, 2, 2)),
               group = "g"),
  make_fixture("sparse", rbind(c(0, 0, 1000, 0), c(0, 0, 0, 5), c(7, 0, 0, 0),
                               c(0, 3, 0, 0), c(9, 9, 9, 9)),
               group = "g", prv_cut = 0.3, lib_cut = 5),
  make_fixture("missing-value",
               rbind(c(NA_real_, 4, 9, 1), c(2, NA_real_, 2, 2), c(0, 0, 0, NA_real_),
                     c(5, 5, NA_real_, 5)),
               group = "g", pseudo = 0.5),
  make_fixture("missing-value-int",
               matrix(as.integer(c(3, NA_integer_, 1, 4, 0, 2, 0, 5, 7, 7, 7, 1)),
                      nrow = 3, byrow = TRUE),
               group = "g", struc_zero = TRUE, pseudo = 0),
  # Boundaries: prv_cut exactly at a prevalence, lib_cut exactly at a library size,
  # pseudo 0 so `log(0)` is `-Inf`, and a single sample and a single group.
  make_fixture("boundary-prevalence", rbind(c(1, 0, 0, 0), c(5, 5, 5, 5), c(0, 0, 0, 3)),
               group = "g", prv_cut = 0.25),
  make_fixture("boundary-library", rbind(c(1, 1, 1, 10), c(0, 0, 0, 0), c(4, 4, 4, 4)),
               group = "g", lib_cut = 10),
  make_fixture("boundary-zero-pseudo", rbind(c(0, 0, 0, 0), c(1, 2, 3, 4)),
               group = "g", pseudo = 0, struc_zero = TRUE),
  make_fixture("one-sample", rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)), group = "g"),
  make_fixture("one-group", rbind(c(1, 0, 3, 0), c(2, 2, 2, 2), c(0, 0, 0, 1)),
               group = "g", struc_zero = TRUE),
  make_fixture("no-group", rbind(c(1, 0, 3, 0), c(2, 2, 2, 2), c(0, 0, 0, 1)),
               struc_zero = FALSE, prv_cut = 0.3),
  # The screen reads the *aggregate* table (ancombc2.R:453-462), so a fixture where
  # the two tables differ is the only one that can catch it reading the wrong one.
  # Here taxon 1 is present in every sample of `feature_table` but absent from group b
  # in the aggregate, so screening `feature_table` would keep it.
  make_fixture("aggregate-differs",
               rbind(c(5, 5, 5, 5), c(5, 5, 5, 5), c(1, 2, 3, 4)),
               aggregate = rbind(c(5, 5, 5, 5), c(5, 0, 5, 5), c(1, 2, 3, 4)),
               group = "g", struc_zero = TRUE),
  make_fixture("neg-lb", rbind(c(5, 0, 0, 0), c(1, 2, 3, 4)),
               group = "g", struc_zero = TRUE, neg_lb = TRUE),
  make_fixture("neg-lb-off", rbind(c(5, 0, 0, 0), c(1, 2, 3, 4)),
               group = "g", struc_zero = TRUE, neg_lb = FALSE)
)

# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------

ok("the bridge is present", ANCOMBC:::rb_available())

for (f in fixtures) {
  cat(sprintf("\n--- %s  (%d x %d, struc_zero=%s, neg_lb=%s, prv_cut=%g, lib_cut=%g, pseudo=%g)\n",
              f$name, f$n_tax, f$n_samp, f$struc_zero, f$neg_lb, f$prv_cut,
              f$lib_cut, f$pseudo))
  
  # R's side. Computed here, from the transcribed expressions.
  zero <- if (f$struc_zero) ref_struc_zero(f$aggregate, f$meta_data, f$group, f$neg_lb) else NULL
  tax_keep <- if (f$struc_zero) zero$tax_keep else seq_len(f$n_tax)
  first <- ref_data_core(f$data, f$prv_cut, f$lib_cut, tax_keep = NULL, samp_keep = NULL)
  lib_size <- ref_lib_size(f$data, first$tax_keep)
  # `ancombc2.R:475-476` passes `samp_keep = colnames(O1)`, i.e. *names*. Transcribed
  # to indices here, which is the same set of columns because the aggregate table
  # keeps `feature_table`'s column names and order.
  second <- ref_data_core(f$aggregate, f$prv_cut, f$lib_cut,
                          tax_keep = tax_keep, samp_keep = first$samp_keep)
  # The second pass's retained taxa in `aggregate`'s coordinates: take the subset-relative
  # indices `.data_core` reports and map them back through the structural-zero subset.
  agg_absolute_keep <- if (f$struc_zero) sort(zero$tax_keep)[second$tax_keep] else second$tax_keep
  lc1 <- ref_log_center(first$feature_table, f$pseudo)
  lc2 <- ref_log_center(second$feature_table, f$pseudo)

  # The bridge's side.
  got <- try(ANCOMBC:::rb_preprocess_stages(
    data = f$data, aggregate = f$aggregate, meta_data = f$meta_data, group = f$group,
    struc_zero = f$struc_zero, neg_lb = f$neg_lb, prv_cut = f$prv_cut,
    lib_cut = f$lib_cut, pseudo = f$pseudo), silent = TRUE)

  if (inherits(got, "try-error")) {
    ok(paste0(f$name, ": the probe runs"), FALSE, as.character(got))
    next
  }
  ok(paste0(f$name, ": the probe runs"), TRUE)

  stages <- list(
    # `.data_core`, pass 1.
    tax_keep1 = first$tax_keep,
    lib_size1 = lib_size,
    samp_keep1 = first$samp_keep,
    O1 = lc1$O,
    log1 = lc1$log,
    means1 = lc1$means,
    y1 = lc1$y,
    # `.data_core`, pass 2.
    # `.data_core` subsets and *then* `which()`es, so its own `tax_keep` counts within
    # the subset. `ref_data_core` reproduces that; the bridge reports the same number
    # and also the absolute set, which is what a Rust caller needs.
    tax_keep2 = second$tax_keep,
    tax_keep2_absolute = agg_absolute_keep,
    samp_keep2 = second$samp_keep,
    O2 = lc2$O,
    log2 = lc2$log,
    means2 = lc2$means,
    y2 = lc2$y,
    # `.get_struc_zero`.
    zero_keep = if (f$struc_zero) as.integer(zero$tax_keep) else seq_len(f$n_tax),
    group_prevalence = if (f$struc_zero) zero$p_hat else numeric(0),
    group_size = if (f$struc_zero) zero$samp_size else numeric(0),
    group_lower = if (f$struc_zero) zero$p_hat_lo else numeric(0),
    # `storage.mode<-` rather than `as.integer()`, which drops the `dim`. The bridge
    # returns the flags as an `n_tax x n_group` matrix and so does the reference's
    # `(p_hat == 0)`, so the shapes are part of what is being compared.
    zero_ind = if (f$struc_zero) {
      zi <- zero$zero_ind
      storage.mode(zi) <- "integer"
      zi
    } else {
      integer(0)
    }
  )
  for (nm in names(stages)) {
    want <- stages[[nm]]
    have <- got[[nm]]
    label <- paste0(f$name, ": ", nm)
    if (is.null(have)) {
      ok(label, FALSE, "the bridge did not return this stage")
      next
    }
    if (nm %in% c("tax_keep1", "samp_keep1", "tax_keep2", "tax_keep2_absolute",
                  "samp_keep2", "zero_keep")) {
      # Integer *indices*: 1-based on both sides, so `identical` is the right check
      # and it is the "retained identifiers" half of S08's acceptance.
      ok(label, identical(as.integer(have), as.integer(want)),
         paste0("got ", paste(have, collapse = ","), " want ", paste(want, collapse = ",")))
    } else if (length(want) == 0L && length(have) == 0L) {
      ok(label, TRUE)
    } else {
      same <- same_bytes(have, want)
      if (!same) assign(label, where_differs(have, want), envir = stage_report)
      ok(label, same, if (same) "" else paste(stage_report[[label]], collapse = "; "))
    }
  }
}

# `prevalence1` needs its own comparison because `.data_core` keeps it local.
for (f in fixtures) {
  have <- try(ANCOMBC:::rb_preprocess_stages(
    data = f$data, aggregate = f$aggregate, meta_data = f$meta_data, group = f$group,
    struc_zero = FALSE, neg_lb = FALSE, prv_cut = f$prv_cut, lib_cut = f$lib_cut,
    pseudo = f$pseudo), silent = TRUE)
  if (inherits(have, "try-error")) next
  want <- prevalence_of(f$data)
  ok(paste0(f$name, ": prevalence1"), same_bytes(have$prevalence1, want),
     if (same_bytes(have$prevalence1, want)) "" else
       paste(where_differs(have$prevalence1, want), collapse = "; "))
}

cat(sprintf("\n%d/%d stage comparisons passed\n", checks - fails, checks))
if (fails > 0L) {
  cat("\nFailing stages:\n")
  for (nm in sort(ls(stage_report))) {
    cat(sprintf("  %s: %s\n", nm, paste(get(nm, envir = stage_report), collapse = "; ")))
  }
  quit(status = 1L)
}
