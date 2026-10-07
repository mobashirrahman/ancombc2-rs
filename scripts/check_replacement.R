#!/usr/bin/env Rscript
# C2 acceptance: the replacement is a drop-in for the original
# (IMPROVED_PLAN.md S05).
#
#   R_LIBS=<replacement lib> Rscript --vanilla scripts/check_replacement.R \
#       --original-lib <dir> --candidate-lib <dir>
#
# What it proves, and what it deliberately does not:
#
#   * the two installed packages declare the same 32 formals, in the same order,
#     with the same default expressions -- compared as source text, because
#     `identical()` on two promises built at different times is not a statement
#     about defaults;
#   * `pseudo_sens` defaults to TRUE;
#   * every original export is present, and nothing extra is exported;
#   * named and positional calls reach the same code;
#   * the replacement does not resolve, load or read the original.
#
# It does NOT compare results. Byte equality of whole results is S02/S03's job
# with both arms actually installed; this file is the API and the independence
# check, and reporting it as numerical parity would be the exact confusion
# IMPROVED_PLAN.md S07 warns about.

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
  k <- sub("^--", "", a)
  if (i == length(args)) stop("option --", k, " needs a value", call. = FALSE)
  opt[[k]] <- args[[i + 1L]]
  i <- i + 2L
}

orig_lib <- opt[["original-lib"]]
cand_lib <- opt[["candidate-lib"]]
if (is.null(orig_lib) || is.null(cand_lib)) {
  stop("usage: check_replacement.R --original-lib <dir> --candidate-lib <dir>",
       call. = FALSE)
}

failures <- 0L
ok <- function(cond, name, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (cond) "PASS" else "FAIL", name,
              if (nzchar(detail) && !cond) paste0("\n         ", detail) else ""))
  if (!cond) failures <<- failures + 1L
}

# ---------------------------------------------------------------------------
# Load each package from its own library, in separate processes' worth of
# isolation: the original is unlinked from .libPaths() before the candidate is
# loaded, and the candidate's search path is inspected afterwards.
# ---------------------------------------------------------------------------

.libPaths(c(orig_lib, .libPaths()))
ok(length(find.package("ANCOMBC", quiet = TRUE)) == 1L, "original_resolves_uniquely",
   paste(find.package("ANCOMBC", quiet = TRUE), collapse = " "))

orig_path <- normalizePath(find.package("ANCOMBC"))
orig_ns <- asNamespace("ANCOMBC")
orig_formals <- formals(ANCOMBC::ancombc2)
orig_exports <- sort(getNamespaceExports("ANCOMBC"))
orig_all <- paste(deparse(args(ANCOMBC::ancombc2)), collapse = "\n")

# Independence: the original must be gone from the candidate's search path.
unloadNamespace("ANCOMBC")
.libPaths(cand_lib)
cand_ns <- asNamespace("ANCOMBC")
cand_path <- normalizePath(find.package("ANCOMBC"))
cand_formals <- formals(ANCOMBC::ancombc2)
cand_exports <- sort(getNamespaceExports("ANCOMBC"))
cand_all <- paste(deparse(args(ANCOMBC::ancombc2)), collapse = "\n")

cat("package replacement acceptance\n")

# ---- identity of the two installations -----------------------------------
ok(!startsWith(cand_path, normalizePath(orig_lib)),
   "candidate_is_not_the_original_installation",
   sprintf("both resolved to %s", cand_path))
# `R CMD INSTALL` moves an inst/ directory's contents to the installed package
# root, so the source-tree path and the installed path differ. Checking the
# source path here would pass a package whose marker never made it into the
# build, which is the case the guard exists to catch.
ok(file.exists(file.path(cand_path, "REPLACEMENT_PROVENANCE.md")),
   "candidate_declares_itself_a_replacement",
   sprintf("looked for %s", file.path(cand_path, "REPLACEMENT_PROVENANCE.md")))
ok(!startsWith(orig_path, normalizePath(cand_lib)),
   "original_is_not_the_candidate_installation",
   sprintf("both resolved to %s", orig_path))

# No loaded namespace may come from the original's library.
from_orig <- Filter(function(nm) {
  p <- tryCatch(getNamespaceInfo(asNamespace(nm), "path"), error = function(e) NA_character_)
  is.character(p) && length(p) == 1L && !is.na(p) && startsWith(p, normalizePath(orig_lib))
}, loadedNamespaces())
ok(length(from_orig) == 0L, "no_namespace_loaded_from_the_original_library",
   paste(from_orig, collapse = ", "))
ok(!any(grepl(normalizePath(orig_lib), .libPaths(), fixed = TRUE)),
   "original_library_absent_from_candidate_search_path",
   paste(.libPaths(), collapse = " | "))

# ---- formals ---------------------------------------------------------------
ok(identical(names(orig_formals), names(cand_formals)),
   "ancombc2_formal_names_and_order_match",
   sprintf("original: %s\n         candidate: %s",
           paste(names(orig_formals), collapse = ","),
           paste(names(cand_formals), collapse = ",")))
ok(identical(orig_all, cand_all), "ancombc2_default_expressions_match",
   sprintf("original:\n%s\n         candidate:\n%s", orig_all, cand_all))

d_eval <- function(f) {
  # Default expressions as text. Two details matter:
  #
  #  * `identical()` on two promises is a statement about how they were built,
  #    not about what the defaults are, so the comparison is on deparsed source.
  #  * A formal with no default is the empty symbol, and *evaluating* it raises
  #    "argument is missing, with no default". `deparse` does not evaluate it and
  #    returns "", which is exactly the distinction being recorded.
  vapply(seq_along(f), function(i) {
    txt <- paste(deparse(f[[i]]), collapse = " ")
    if (!nzchar(txt)) "<required>" else txt
  }, character(1))
}
do <- d_eval(orig_formals); dc <- d_eval(cand_formals)
diff_defaults <- names(do)[do != dc[names(do)]]
ok(length(diff_defaults) == 0L, "every_default_expression_is_identical",
   paste(sprintf("%s: %s vs %s", diff_defaults, do[diff_defaults],
                 dc[diff_defaults]), collapse = "; "))

ps <- deparse(orig_formals$pseudo_sens)
ok(identical(ps, "TRUE"),
   "pseudo_sens_defaults_to_TRUE",
   sprintf("the pinned original's pseudo_sens default is %s", ps))

ok(!"..." %in% names(orig_formals) && !"..." %in% names(cand_formals),
   "no_dot_dot_dot_added_to_the_entrypoint",
   "adding an accepting `...` would make calls the original rejects succeed here")

# ---- exports ---------------------------------------------------------------
missing_exports <- setdiff(orig_exports, cand_exports)
ok(length(missing_exports) == 0L, "every_original_export_is_present",
   paste(missing_exports, collapse = ", "))
extra_exports <- setdiff(cand_exports, orig_exports)
ok(length(extra_exports) == 0L, "no_export_added_beyond_the_original",
   paste(extra_exports, collapse = ", "))
ok(identical(orig_exports,
             sort(c("ancom", "ancombc", "ancombc2", "data_sanity_check",
                    "secom_dist", "secom_linear", "sim_plnm"))),
   "the_export_set_is_the_original_seven")

# ---- other entry points keep their signatures ------------------------------
sig_bad <- character(0)
for (nm in orig_exports) {
  a <- paste(deparse(args(get(nm, envir = orig_ns))), collapse = "\n")
  b <- paste(deparse(args(get(nm, envir = cand_ns))), collapse = "\n")
  if (!identical(a, b)) sig_bad <- c(sig_bad, nm)
}
ok(length(sig_bad) == 0L, "every_export_keeps_its_signature",
   paste(sig_bad, collapse = ", "))

# ---- DESCRIPTION is honest -------------------------------------------------
cand_desc <- read.dcf(file.path(cand_path, "DESCRIPTION"))
ok(!is.null(cand_desc[1, "X-Replacement-Of"]),
   "candidate_description_records_what_it_replaces",
   "a compatible API is not a claim to be the upstream build")
ok(trimws(cand_desc[1, "Authors@R"]) == trimws(
     read.dcf(file.path(orig_path, "DESCRIPTION"))[1, "Authors@R"]),
   "upstream_attribution_retained")
ok(trimws(cand_desc[1, "License"]) == "Artistic-2.0",
   "upstream_licence_retained", trimws(cand_desc[1, "License"]))

# ---- calls reach the same code -------------------------------------------
# A named call and a positional call with only the required arguments, on a tiny
# deterministic input. The point is that neither is rejected by argument
# handling; the numbers themselves are S02/S03's business.
mk <- function() {
  set.seed(11)
  counts <- matrix(as.integer(rpois(60, lambda = 40)) + 1L, nrow = 6,
                   dimnames = list(paste0("T", 1:6), paste0("S", 1:10)))
  meta <- data.frame(group = factor(rep(c("a", "b", "c"), length.out = 10)),
                     x1 = round(rnorm(10), 6))
  rownames(meta) <- colnames(counts)
  list(data = counts, meta_data = meta, fix_formula = "group + x1")
}
inp <- mk()
named <- tryCatch(
  do.call(ANCOMBC::ancombc2, c(inp, list(p_adj_method = "BH", verbose = FALSE))),
  error = function(e) e)
ok(!inherits(named, "error"), "named_call_succeeds",
   if (inherits(named, "error")) conditionMessage(named) else "")

# Positional: data, meta_data, fix_formula in formal order, everything else
# default. `pseudo_sens` therefore takes the original's own default, which is
# TRUE -- the defect this whole check exists to keep fixed.
pos <- tryCatch(
  do.call(ANCOMBC::ancombc2, list(inp$data, TRUE, "counts", "counts", NULL, NULL,
                                   NULL, inp$meta_data, inp$fix_formula),
          quote = FALSE),
  error = function(e) e)
ok(!inherits(pos, "error"), "positional_call_succeeds",
   if (inherits(pos, "error")) conditionMessage(pos) else "")

if (!inherits(pos, "error")) {
  # A positional call omits pseudo_sens, so it must be the sensitivity analysis.
  # The evidence is the *effect*: `passed_ss` and `diff_robust` are the
  # sensitivity columns, and the same call with pseudo_sens = FALSE does not
  # produce them. Checking a default by asserting the default's value would be
  # circular; this checks what the default changed.
  # The columns are suffixed per coefficient -- passed_ss_(Intercept),
  # diff_robust_groupb and so on -- so the test is on the prefix, not the name.
  res_cols <- names(pos$res)
  ok(any(startsWith(res_cols, "passed_ss")) && any(startsWith(res_cols, "diff_robust")),
     "positional_call_defaults_to_pseudo_sens_TRUE",
     sprintf("res columns: %s", paste(res_cols, collapse = ",")))

  nosens <- tryCatch(
    suppressMessages(do.call(ANCOMBC::ancombc2,
                             c(inp, list(pseudo_sens = FALSE, verbose = FALSE)))),
    error = function(e) e)
  if (inherits(nosens, "error")) {
    ok(FALSE, "the_pseudo_sens_FALSE_control_call_succeeds",
       conditionMessage(nosens))
  } else {
    ok(!any(startsWith(names(nosens$res), "passed_ss")),
       "pseudo_sens_false_really_drops_the_sensitivity_columns",
       "without this control, the previous check could pass on a package that always adds the columns")
  }
}


# ---- an unsupported argument is still an error ----------------------------
bad <- tryCatch(
  do.call(ANCOMBC::ancombc2, c(inp, list(not_an_argument = 1, verbose = FALSE))),
  error = function(e) e)
ok(inherits(bad, "error"), "unknown_argument_is_still_refused",
   "the entrypoint must not have grown a `...`")

cat(sprintf("  %d check(s) failed\n", failures))
if (failures > 0L) quit(save = "no", status = 1L)
quit(save = "no", status = 0L)
