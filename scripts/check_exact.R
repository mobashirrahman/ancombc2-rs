#!/usr/bin/env Rscript
# Strict exact-result comparator (IMPROVED_PLAN.md S03).
#
#   Rscript --vanilla scripts/check_exact.R --left <capture-dir> --right <capture-dir>
#   Rscript --vanilla scripts/check_exact.R --left-obj a.rds --right-obj b.rds
#   Rscript --vanilla scripts/check_exact.R --selftest
#   Rscript --vanilla scripts/check_exact.R --cases validation/exact/cases.json \
#          --case-set required [--observed observed.json]
#
# The acceptance rule is one expression, applied to whole objects and to nothing
# smaller:
#
#   identical(serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3),
#             serialize(y, NULL, ascii = FALSE, xdr = TRUE, version = 3))
#
# Everything else in this file exists to *name* a difference once one has been
# found: which field, which attribute, which element, and for a number its raw
# 64 bits. No diagnostic here can accept a result. There is no tolerance
# parameter, no `all.equal` fallback and no normalization, because a comparator
# with an escape hatch is a comparator that will eventually use it.
#
# Raw bytes, not a text rendering: `NA_real_` and `NaN` print the same, `-0` and
# `0` print the same, and `1L` and `1` compare equal with `==`.

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

opt <- list()
FLAGS <- c("selftest")
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
  k <- sub("^--", "", a)
  if (k %in% FLAGS) {
    opt[[k]] <- TRUE
    i <- i + 1L
    next
  }
  if (i == length(args)) stop("option --", k, " needs a value", call. = FALSE)
  opt[[k]] <- args[[i + 1L]]
  i <- i + 2L
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# The contract
# ---------------------------------------------------------------------------

result_bytes <- function(x) serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3)

is_byte_identical <- function(a, b) {
  ba <- result_bytes(a)
  bb <- result_bytes(b)
  identical(ba, bb)
}

# Raw representation of one double, big-endian, as hex. This is what actually
# distinguishes -0 from 0 and NA_real_ from NaN.
hex64 <- function(x) {
  if (length(x) != 1L) return(paste(vapply(x, hex64, character(1)), collapse = " "))
  r <- writeBin(as.double(x), raw(), size = 8L, endian = "big")
  paste(sprintf("%02x", as.integer(r)), collapse = "")
}

class_tag <- function(x) {
  paste(class(x), collapse = "/")
}

# ---------------------------------------------------------------------------
# Recursive first-difference diagnostic. Never decides acceptance.
# ---------------------------------------------------------------------------

MAX_REPORTED <- 40L

diff_objects <- function(a, b, path = "$", out = list()) {
  if (length(out) >= MAX_REPORTED) return(out)

  note <- function(kind, detail) {
    out[[length(out) + 1L]] <<- list(path = path, kind = kind, detail = detail)
  }

  # NULL is checked before class, because `class(NULL)` is the string "NULL" and
  # comparing it to `class(list())` would report the symptom ("class NULL vs
  # list") instead of the cause (a field that is present and empty rather than
  # absent). IMPROVED_PLAN.md S07 requires absent-test NULL entries to survive.
  if (!identical(is.null(a), is.null(b))) {
    note("null", sprintf("%s is NULL, %s is %s", class_tag(a), class_tag(b),
                         if (is.null(b)) "NULL" else describe(b)))
    return(out)
  }

  if (!identical(class_tag(a), class_tag(b))) {
    note("class", sprintf("class %s vs %s", class_tag(a), class_tag(b)))
    return(out)
  }

  # Attributes first: a changed dim or a changed dimname changes the reading of
  # every element, and reporting 400 element diffs that are all consequences of
  # one changed attribute wastes the report.
  an <- attributes(a)
  bn <- attributes(b)
  akeys <- sort(names(an))
  bkeys <- sort(names(bn))
  if (!identical(akeys, bkeys)) {
    note("attribute-set", sprintf("attributes %s vs %s",
                                  paste(akeys, collapse = ","),
                                  paste(bkeys, collapse = ",")))
    return(out)
  }
  for (k in akeys) {
    if (identical(an[[k]], bn[[k]])) next
    if (k %in% c("names", "dim", "dimnames", "class", "levels", "row.names")) {
      note(paste0("attribute:", k), paste(describe(an[[k]]), collapse = " | "))
    } else {
      note(paste0("attribute:", k),
           sprintf("%s vs %s", describe(an[[k]]), describe(bn[[k]])))
    }
  }

  if (is.null(a) || is.null(b)) {
    if (!identical(is.null(a), is.null(b))) note("null", "one side is NULL")
    return(out)
  }

  if (is.raw(a)) {
    if (!identical(as.integer(a), as.integer(b))) {
      i <- which(as.integer(a) != as.integer(b))
      note("raw", sprintf("%d byte(s) differ, first at %d: %s vs %s", length(i),
                          i[1], hex64(as.numeric(a[i[1]])), hex64(as.numeric(b[i[1]]))))
    }
    return(out)
  }

  if (is.atomic(a)) {
    if (is.double(a) && is.double(b)) {
      na <- is.na(a) & !is.nan(a)
      nb <- is.na(b) & !is.nan(b)
      if (!identical(na, nb)) {
        i <- which(xor(na, nb))
        note("na-kind", sprintf("%d cell(s) differ in NA vs NaN, first at %d", length(i), i[1]))
        return(out)
      }
      if (!identical(is.nan(a), is.nan(b))) {
        i <- which(xor(is.nan(a), is.nan(b)))
        note("nan", sprintf("%d cell(s) differ in NaN, first at %d", length(i), i[1]))
        return(out)
      }
      finite <- !is.nan(a)
      if (any(finite)) {
        bits_a <- hex64(a[finite]); bits_b <- hex64(b[finite])
        diff <- which(bits_a != bits_b)
        if (length(diff)) {
          k <- diff[1]
          note("double-bits",
               sprintf("index %d: %.17g [%s] vs %.17g [%s] (%d of %d doubles differ in raw bits)",
                       k, a[finite][k], bits_a[k], b[finite][k], bits_b[k],
                       length(diff), length(bits_a)))
        }
      }
      return(out)
    }
    if (is.integer(a) && is.double(b)) {
      note("type", sprintf("integer vs double at %s", path))
      return(out)
    }
    if (!identical(a, b)) {
      # Name the position. "the vector differs" in a 10,000-element field is a
      # restatement of the failure, not a diagnosis; the element index and its
      # two values are the diagnosis.
      idx <- which(!mapply(identical, a, b))
      idx <- idx[!vapply(idx, function(j) identical(a[[j]], b[[j]]), logical(1))]
      if (length(idx)) {
        j <- idx[1]
        note("value", sprintf("index %d (%d of %d differ): %s vs %s; first %s",
                              j, length(idx), length(a), describe(a[j]), describe(b[j]),
                              if (length(idx) > 1L) "(and more)" else ""))
      } else {
        note("value", sprintf("%s vs %s", describe(a), describe(b)))
      }
    }
    return(out)
  }

  if (is.function(a)) {
    if (!identical(deparse(a), deparse(b))) note("function", "bodies differ")
    return(out)
  }

  if (is.environment(a)) {
    note("environment", "environments are compared by identity, not by contents")
    return(out)
  }

  if (is.list(a)) {
    if (length(a) != length(b)) {
      note("length", sprintf("%d vs %d", length(a), length(b)))
    }
    nms_a <- names(a); nms_b <- names(b)
    if (!identical(nms_a, nms_b)) {
      # Report the first position where the name differs: an element swap is the
      # common cause and it is invisible in a set comparison.
      n <- max(length(nms_a), length(nms_b))
      nms_a2 <- rep_len(nms_a %||% character(0), n)
      nms_b2 <- rep_len(nms_b %||% character(0), n)
      for (j in seq_len(n)) {
        if (!identical(nms_a2[[j]], nms_b2[[j]])) {
          note("names", sprintf("position %d: %s vs %s", j,
                                nms_a2[[j]] %||% "<unnamed>",
                                nms_b2[[j]] %||% "<unnamed>"))
          break
        }
      }
    }
    for (j in seq_len(min(length(a), length(b)))) {
      if (length(out) >= MAX_REPORTED) break
      nm <- (nms_a %||% character(0))
      lbl <- if (j <= length(nm) && !is.na(nm[[j]]) && nzchar(nm[[j]])) nm[[j]] else j
      out <- diff_objects(a[[j]], b[[j]], path = paste0(path, "$", lbl), out = out)
    }
    return(out)
  }

  if (isS4(a)) {
    if (!identical(a, b)) note("S4", "slot contents differ")
    return(out)
  }
  note("unsupported", sprintf("cannot diff an object of class %s", class_tag(a)))
  out
}

# A diagnostic that prints a 10,000-element vector has replaced one unreadable
# failure with another, so every rendering is capped and says how much it hid.
MAX_SHOWN <- 24L

describe <- function(x) {
  if (is.null(x)) return("NULL")
  cap <- function(txt, n) {
    if (n <= MAX_SHOWN) return(txt)
    sprintf("%s ... (%d of %d shown)", paste(txt[seq_len(MAX_SHOWN)], collapse = " "), MAX_SHOWN, n)
  }
  if (is.character(x)) return(cap(vapply(x, function(v) paste0('"', v, '"'), character(1)), length(x)))
  if (is.logical(x)) return(cap(vapply(x, function(v) if (is.na(v)) "NA" else as.character(v), character(1)), length(x)))
  if (is.integer(x)) return(cap(vapply(x, function(v) if (is.na(v)) "NA" else as.character(v), character(1)), length(x)))
  if (is.double(x)) {
    if (length(x) == 0L) return("numeric(0)")
    tags <- vapply(x, function(v) {
      if (is.nan(v)) "NaN" else if (is.na(v)) "NA_real_" else if (is.infinite(v)) {
        if (v > 0) "Inf" else "-Inf"
      } else if (v == 0) {
        if (1 / v < 0) "-0" else "0"
      } else sprintf("%.17g", v)
    }, character(1))
    return(cap(paste0(tags, " [", hex64(x), "]"), length(x)))
  }
  if (is.list(x)) return(sprintf("list(%d)", length(x)))
  sprintf("<%s>", class_tag(x))
}

# ---------------------------------------------------------------------------
# Capture comparison
# ---------------------------------------------------------------------------

read_kv <- function(path) {
  out <- character(0)
  if (!file.exists(path)) return(out)
  for (line in readLines(path, warn = FALSE)) {
    if (grepl("=", line, fixed = TRUE)) {
      k <- sub("=.*$", "", line); v <- sub("^[^=]*=", "", line)
      out[[k]] <- v
    }
  }
  out
}

capture_payload <- function(dir) {
  rb <- file.path(dir, "result.bytes")
  eb <- file.path(dir, "error.bytes")
  st <- file.path(dir, "state.rds")
  list(
    dir = dir,
    outcome = read_kv(file.path(dir, "meta.tsv"))[["outcome"]] %||% "<none>",
    result = if (file.exists(rb) && file.info(rb)$size > 0) readBin(rb, "raw", file.info(rb)$size) else NULL,
    error = if (file.exists(eb) && file.info(eb)$size > 0) readBin(eb, "raw", file.info(eb)$size) else NULL,
    state = if (file.exists(st)) readRDS(st) else NULL,
    meta = read_kv(file.path(dir, "meta.tsv"))
  )
}

compare_captures <- function(dir_a, dir_b, verbose = TRUE) {
  A <- capture_payload(dir_a)
  B <- capture_payload(dir_b)
  problems <- character(0)
  diffs <- list()

  if (!identical(A$outcome, B$outcome)) {
    problems <- c(problems, sprintf("outcome %s vs %s", A$outcome, B$outcome))
    return(list(ok = FALSE, problems = problems, diffs = diffs))
  }

  if (identical(A$outcome, "success")) {
    pa <- A$result; pb <- B$result
    if (is.null(pa) || is.null(pb)) {
      problems <- c(problems, "one side has no result.bytes")
    } else if (!identical(pa, pb)) {
      problems <- c(problems, sprintf(
        "result bytes differ: %d vs %d bytes", length(pa), length(pb)))
      if (length(pa) == length(pb)) {
        j <- which(as.integer(pa) != as.integer(pb))
        problems <- c(problems, sprintf(
          "  %d differing byte(s), first at offset %d (0x%s vs 0x%s)",
          length(j), j[1], paste(sprintf("%02x", as.integer(pa[j[1]])), collapse = ""),
          paste(sprintf("%02x", as.integer(pb[j[1]])), collapse = "")))
      }
      # Unserialize both and name the first differing field. This never decides
      # acceptance; it only makes the failure readable.
      if (!is.null(pa) && !is.null(pb) && length(pa) == length(pb)) {
        oa <- tryCatch(unserialize(pa), error = function(e) NULL)
        ob <- tryCatch(unserialize(pb), error = function(e) NULL)
        if (!is.null(oa) && !is.null(ob)) diffs <- diff_objects(oa, ob)
      }
    }
  } else {
    if (!identical(A$error, B$error)) {
      problems <- c(problems, "error-condition bytes differ")
      oa <- tryCatch(unserialize(A$error), error = function(e) NULL)
      ob <- tryCatch(unserialize(B$error), error = function(e) NULL)
      if (!is.null(oa) && !is.null(ob)) diffs <- diff_objects(oa, ob)
    }
  }

  # Conditions: class, message, call and order are all part of the contract.
  ca <- A$state$conditions %||% list()
  cb <- B$state$conditions %||% list()
  if (length(ca) != length(cb)) {
    problems <- c(problems, sprintf("condition count %d vs %d", length(ca), length(cb)))
  } else if (length(ca)) {
    for (j in seq_along(ca)) {
      if (!identical(ca[[j]], cb[[j]])) {
        diffs <- c(diffs, list(list(path = sprintf("$.conditions[%d]", j),
                                   kind = "condition",
                                   detail = sprintf("%s %s vs %s %s", ca[[j]]$type,
                                                    describe(ca[[j]]$message),
                                                    cb[[j]]$type,
                                                    describe(cb[[j]]$message)))))
        problems <- c(problems, sprintf("condition %d differs", j))
      }
    }
  }

  # Reference-observable state after the call.
  for (field in c("rng_kind", "dot_random_seed", "options", "globalenv_names",
                  "search_path", "loaded_namespaces")) {
    va <- A$state$after[[field]]
    vb <- B$state$after[[field]]
    if (!identical(va, vb)) {
      problems <- c(problems, sprintf("post-call state differs: %s", field))
      diffs <- c(diffs, diff_objects(va, vb, path = paste0("$.after.", field)))
    }
  }

  # The arms must at least agree about which package they ran.
  for (field in c("pkg_version", "outcome", "fn", "call_style")) {
    va <- A$meta[[field]]; vb <- B$meta[[field]]
    if (!identical(va, vb)) {
      problems <- c(problems, sprintf("run metadata differs: %s (%s vs %s)",
                                      field, va %||% "<none>", vb %||% "<none>"))
    }
  }

  if (verbose) {
    cat(sprintf("%-24s %s\n", "compare", paste(basename(dir_a), "vs", basename(dir_b))))
    cat(sprintf("%-24s %s\n", "outcome", A$outcome))
    if (length(problems) == 0L) {
      cat("VERDICT: identical\n")
    } else {
      cat("VERDICT: DIFFERENT\n")
      for (p in problems) cat("  - ", p, "\n", sep = "")
      if (length(diffs)) {
        cat("  first differences:\n")
        for (d in diffs) {
          cat(sprintf("    %s [%s] %s\n", d$path, d$kind, d$detail))
        }
      }
    }
  }
  list(ok = length(problems) == 0L, problems = problems, diffs = diffs)
}

# ---------------------------------------------------------------------------
# Negative checks. Every one of these MUST be rejected; identical objects MUST
# pass. A comparator that cannot fail is not a comparator.
# ---------------------------------------------------------------------------

flip_one_bit <- function(x) {
  n <- length(x)
  x[[n]] <- bitwXor(x[[n]], 1L)
  x
}

selftest <- function() {
  failures <- 0L
  report <- function(name, ok, why = "") {
    cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", name,
                if (nzchar(why)) paste0(" -- ", why) else ""))
    if (!ok) failures <<- failures + 1L
  }

  # Negative zero is built arithmetically, never written as the literal `-0`.
  # R's parser folds a negated numeric literal in some contexts and not others --
  # measured here, `list(nz = -0)` inside this function folded to +0 while the
  # same literal in a two-line script did not -- so a selftest that used the
  # literal would sometimes be testing +0 against +0 and would pass for the
  # wrong reason. `identical(-0, 0)` is also TRUE in R, which is precisely why
  # the contract compares serialized bytes.
  neg_zero <- function() 0 * -1

  base <- list(
    alpha = 0.05,
    beta = c(1.5, -2.25, 0),
    names = c("a", "b", "c"),
    flag = TRUE,
    na = NA_real_, nan = NaN, inf = Inf, ninf = -Inf,
    nz = neg_zero(), pz = 0,
    i = 42L, d = 42,
    empty = numeric(0),
    nested = list(x = 1:3, y = list(z = "s")),
    df = data.frame(g = factor(c("p", "q")), v = c(1.5, 2.5),
                    row.names = c("r1", "r2")),
    mat = matrix(1:6, nrow = 2, dimnames = list(c("m1", "m2"), c("c1", "c2", "c3"))),
    nullfield = NULL,
    withattr = structure(list(v = 1), an = "kept")
  )
  cat("check_exact.R selftest\n")

  # -- identical objects pass
  report("identical_objects_pass", is_byte_identical(base, base),
         "the comparator must accept a copy of itself")

  # -- one mantissa bit
  t1 <- base; t1$beta <- flip_one_bit(t1$beta)
  report("one_mantissa_bit_rejected", !is_byte_identical(base, t1))

  # -- signed zero, both directions
  if (!identical(1 / base$nz, -Inf) || !identical(1 / base$pz, Inf)) {
    report("signed_zero_fixture_is_sane", FALSE,
           "base must hold one -0 and one +0 for the next two checks to mean anything")
  } else {
    report("signed_zero_fixture_is_sane", TRUE)
  }
  t2 <- base; t2$nz <- 0
  report("negative_zero_to_positive_rejected", !is_byte_identical(base, t2))
  t3 <- base; t3$pz <- neg_zero()
  report("positive_zero_to_negative_rejected", !is_byte_identical(base, t3))

  # -- NA vs NaN, both directions, and inf
  t4 <- base; t4$nan <- NA_real_
  report("nan_to_na_rejected", !is_byte_identical(base, t4))
  t5 <- base; t5$na <- NaN
  report("na_to_nan_rejected", !is_byte_identical(base, t5))
  t6 <- base; t6$inf <- .Machine$double.xmax
  report("inf_to_finite_rejected", !is_byte_identical(base, t6))

  # -- integer vs double with the same value
  t7 <- base; t7$i <- 42
  report("integer_vs_double_rejected", !is_byte_identical(base, t7))

  # -- a name
  t8 <- base; t8$names <- c("a", "b", "d")
  report("changed_name_rejected", !is_byte_identical(base, t8))

  # -- row order in a data.frame and in a matrix
  t9 <- base; t9$df <- t9$df[c(2, 1), , drop = FALSE]
  report("df_row_order_rejected", !is_byte_identical(base, t9))
  t10 <- base; t10$mat <- t10$mat[c(2, 1), , drop = FALSE]
  report("matrix_row_order_rejected", !is_byte_identical(base, t10))

  # -- a NULL field becoming an empty list, and vice versa
  t11 <- base; t11$nullfield <- list()
  report("null_field_rejected", !is_byte_identical(base, t11))
  t12 <- base; t12$nullfield <- NULL
  report("null_field_dropped_rejected", !is_byte_identical(base, t12))

  # -- list order
  t13 <- base; t13$nested <- list(y = base$nested$y, x = base$nested$x)
  report("list_element_order_rejected", !is_byte_identical(base, t13))

  # -- an attribute
  t14 <- base; attr(t14$withattr, "an") <- "changed"
  report("changed_attribute_rejected", !is_byte_identical(base, t14))
  t15 <- base; attr(t15$beta, "units") <- "log"
  report("added_attribute_rejected", !is_byte_identical(base, t15))

  # -- factor levels
  t16 <- base; t16$df$g <- factor(c("q", "p"), levels = c("p", "q"))
  report("factor_value_swap_rejected", !is_byte_identical(base, t16))
  t17 <- base; levels(t17$df$g) <- c("q", "p", "unused")
  report("extra_factor_level_rejected", !is_byte_identical(base, t17))

  # -- a changed warning, in place of the object
  ca <- list(list(class = "simpleWarning", message = "one", call = NA_character_, type = "warning"))
  cb <- list(list(class = "simpleWarning", message = "two", call = NA_character_, type = "warning"))
  report("changed_warning_rejected", !identical(result_bytes(ca), result_bytes(cb)))

  # -- a changed RNG state
  ra <- list(dot_random_seed = c(1L, 2L, 3L), rng_kind = RNGkind())
  rb <- list(dot_random_seed = c(1L, 2L, 4L), rng_kind = RNGkind())
  report("changed_rng_state_rejected", !identical(result_bytes(ra), result_bytes(rb)))

  # -- a reordering of the conditions themselves
  c2 <- list(ca[[1]], list(class = "simpleMessage", message = "note", call = NA_character_, type = "message"))
  c3 <- list(list(class = "simpleMessage", message = "note", call = NA_character_, type = "message"), ca[[1]])
  report("condition_order_rejected", !identical(result_bytes(c2), result_bytes(c3)))

  # -- the diagnostic must actually name things, without deciding anything
  d <- diff_objects(base, t1)
  ok <- length(d) > 0 && any(grepl("double-bits", vapply(d, function(x) x$kind, character(1))))
  report("diagnostic_reports_double_bits", ok,
         if (length(d)) paste("kinds:", paste(vapply(d, function(x) x$kind, character(1)), collapse = ",")) else "no differences reported")

  # Adding a `names` attribute where there was none is an attribute change, and
  # is reported as one rather than as 3 mysterious value diffs.
  t8a <- base; names(t8a$beta) <- c("x", "y", "z")
  d2 <- diff_objects(base, t8a)
  # Adding `names` where there was none changes the *set* of attribute names, so
  # the accurate report is `attribute-set`; a `names` attribute that already
  # existed and changed value is reported as `attribute:names`. Both name the
  # difference, which is all the diagnostic owes.
  ok2 <- length(d2) > 0 &&
    any(grepl("^attribute(:names|-set)$", vapply(d2, function(x) x$kind, character(1))))
  report("diagnostic_reports_added_names_attribute", ok2,
         if (length(d2)) paste("kinds:", paste(vapply(d2, function(x) x$kind, character(1)), collapse = ",")) else "")

  # Changed dimnames on a matrix: also an attribute, and it is the one that makes
  # every element read differently.
  t8d <- base; dimnames(t8d$mat)[[1]] <- c("m2", "m1")
  d2d <- diff_objects(base, t8d)
  ok2d <- length(d2d) > 0 &&
    any(grepl("^attribute:dimnames$", vapply(d2d, function(x) x$kind, character(1))))
  report("diagnostic_reports_dimnames", ok2d)

  # A changed *value* of an existing names attribute.
  t8e <- base; t8e$df <- t8e$df; attr(t8e$df$g, "levels") <- c("p", "q", "r")
  d2e <- diff_objects(base, t8e)
  ok2e <- length(d2e) > 0 &&
    any(grepl("^attribute:levels$", vapply(d2e, function(x) x$kind, character(1))))
  report("diagnostic_reports_changed_factor_levels", ok2e)

  # A renamed list element is reported under `names` with its position.
  t8b <- base; names(t8b)[2] <- "BETA"
  d2b <- diff_objects(base, t8b)
  ok2b <- length(d2b) > 0 &&
    any(grepl("names", vapply(d2b, function(x) x$kind, character(1))))
  report("diagnostic_reports_renamed_element", ok2b)

  d3 <- diff_objects(base, t11)
  ok3 <- length(d3) > 0 &&
    any(grepl("^null$", vapply(d3, function(x) x$kind, character(1))))
  report("diagnostic_reports_null_field", ok3,
         if (length(d3)) paste("kinds:", paste(vapply(d3, function(x) x$kind, character(1)), collapse = ",")) else "")

  # -- the smallest possible differences, which is what a tolerance would absorb
  tiny <- list(v = 1)
  tiny2 <- list(v = 1 + 2^-52)
  report("one_ulp_rejected", !is_byte_identical(tiny, tiny2),
         "adjacent doubles above 1 differ by one mantissa bit")
  near <- list(v = 1)
  near2 <- list(v = 1 + 1e-12)
  report("tiny_relative_difference_rejected", !is_byte_identical(near, near2),
         "isTRUE(all.equal()) accepts this; the exact contract must not")
  report("tiny_relative_difference_is_invisible_to_all_equal",
         isTRUE(all.equal(near$v, near2$v)),
         "this check exists to document what all.equal() would have let through")

  cat(sprintf("  %d check(s) failed\n", failures))
  if (failures > 0L) quit(save = "no", status = 1L)
  quit(save = "no", status = 0L)
}

# ---------------------------------------------------------------------------
# Case manifest
# ---------------------------------------------------------------------------

# A tiny JSON reader for the case manifest. It is deliberately strict: an
# unrecognised key is an error, because a manifest that silently ignores a
# field is a manifest that silently drops a required case.
parse_cases_json <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to read the case manifest", call. = FALSE)
  }
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

check_manifest <- function(path, case_set, observed_path = NULL) {
  m <- parse_cases_json(path)
  known_sets <- names(m$case_sets)
  if (!case_set %in% known_sets) {
    stop("unknown case set ", sQuote(case_set), "; known: ",
         paste(known_sets, collapse = ", "), call. = FALSE)
  }
  # The loop variable is deliberately not named `c`: shadowing base::c inside a
  # function that also uses `c(...)` to build its problem list is a bug that
  # fails at run time, not at read time.
  want <- m$case_sets[[case_set]]$cases
  if (is.null(want) || !length(want)) {
    cat("MANIFEST: INVALID\n")
    cat("  - case set ", case_set, " names no cases\n", sep = "")
    return(list(ok = FALSE, problems = "empty case set"))
  }
  want_ids <- vapply(want, function(cs) cs$id, character(1))
  problems <- character(0)

  if (anyDuplicated(want_ids)) {
    dup <- unique(want_ids[duplicated(want_ids)])
    problems <- c(problems, sprintf("duplicate case id(s) in the manifest: %s",
                                    paste(dup, collapse = ", ")))
  }
  for (cs in want) {
    for (k in c("id", "input", "expect", "tags")) {
      if (is.null(cs[[k]])) {
        problems <- c(problems, sprintf("case %s: missing required key %s",
                                        cs$id %||% "<no id>", k))
      }
    }
    if (!is.null(cs$input) && !file.exists(file.path(dirname(path), cs$input)) &&
        !file.exists(cs$input)) {
      problems <- c(problems, sprintf("case %s: input not found: %s", cs$id, cs$input))
    }
  }

  observed <- character(0)
  if (!is.null(observed_path) && file.exists(observed_path)) {
    obs <- parse_cases_json(observed_path)
    observed <- if (!is.null(obs$cases)) {
      vapply(obs$cases, function(x) if (is.null(x$id)) "<case with no id>" else x$id, character(1))
    } else {
      character(0)
    }
  }
  missing_cases <- setdiff(want_ids, observed)
  if (length(missing_cases)) {
    problems <- c(problems, sprintf("%d required case(s) not observed: %s",
                                    length(missing_cases),
                                    paste(missing_cases, collapse = ", ")))
  }
  # A duplicated observation is not a stronger result, it is two runs of the same
  # case being counted twice -- which is how a suite reports more coverage than
  # it executed. It is rejected rather than collapsed.
  if (anyDuplicated(observed)) {
    dup <- unique(observed[duplicated(observed)])
    problems <- c(problems, sprintf("case(s) observed more than once: %s",
                                    paste(dup, collapse = ", ")))
  }
  extra <- setdiff(observed, want_ids)
  if (length(extra) && !isTRUE(m$case_sets[[case_set]]$allow_extra)) {
    problems <- c(problems, sprintf("%d observed case(s) not in the set: %s",
                                    length(extra), paste(extra, collapse = ", ")))
  }

  cat(sprintf("case set %-10s %d required, %d observed\n", case_set,
              length(want_ids), length(observed)))
  if (length(problems)) {
    cat("MANIFEST: INVALID\n")
    for (p in problems) cat("  - ", p, "\n", sep = "")
    return(list(ok = FALSE, problems = problems))
  }
  cat("MANIFEST: complete\n")
  list(ok = TRUE, problems = character(0))
}

# ---------------------------------------------------------------------------

main <- function() {
  if (isTRUE(opt$selftest)) return(selftest())

  if (!is.null(opt$cases)) {
    if (is.null(opt[["case-set"]])) stop("--cases requires --case-set", call. = FALSE)
    r <- check_manifest(opt$cases, opt[["case-set"]], opt$observed)
    return(if (r$ok) 0L else 1L)
  }

  if (!is.null(opt[["left-obj"]])) {
    if (is.null(opt[["right-obj"]])) stop("--left-obj requires --right-obj", call. = FALSE)
    a <- readRDS(opt[["left-obj"]]); b <- readRDS(opt[["right-obj"]])
    if (is_byte_identical(a, b)) { cat("VERDICT: identical\n"); return(0L) }
    cat("VERDICT: DIFFERENT\n")
    for (d in diff_objects(a, b)) cat(sprintf("  %s [%s] %s\n", d$path, d$kind, d$detail))
    return(1L)
  }

  if (!is.null(opt$left)) {
    if (is.null(opt$right)) stop("--left requires --right", call. = FALSE)
    r <- compare_captures(opt$left, opt$right)
    return(if (r$ok) 0L else 1L)
  }

  cat("usage: check_exact.R --left <dir> --right <dir>\n")
  cat("       check_exact.R --left-obj a.rds --right-obj b.rds\n")
  cat("       check_exact.R --selftest\n")
  cat("       check_exact.R --cases cases.json --case-set <name> [--observed obs.json]\n")
  2L
}

status <- main()
quit(save = "no", status = status)
