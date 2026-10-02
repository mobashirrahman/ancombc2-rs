#!/usr/bin/env Rscript
# The R arm of the simulation grid: ANCOMBC 2.15.2 over the Rust generator's tables.
#
#   Rscript --vanilla scripts/sim_r.R --data DIR --out results.jsonl
#
# `DIR` is what `ancombc2-sim generate` wrote: one `cellNNNN/repNNNN/` per
# replicate, each holding `counts.tsv`, `meta.tsv`, `truth.tsv` and
# `replicate.json`.
#
# Why the data is written out and read back rather than regenerated in R
# ------------------------------------------------------------------------
# The acceptance rule is "Rust's empirical FDR is within Monte-Carlo error of
# R's". That comparison is only meaningful if both arms analysed the same
# numbers, and the only way to guarantee that is to analyse the same *file*. A
# shared seed is not a shared table: any change to the generator on either side
# silently decouples the arms, and the resulting comparison would be noise
# dressed as a result. So the generator is Rust's, the table is written once, and
# both arms read it.
#
# What this script computes
# -------------------------
# For each replicate: run `ancombc2()`, then for the retained taxa compare
# `beta` (the group term) against the truth in `truth.tsv` and count
# true/false `diff_abn` calls. It writes one JSON object per replicate with the
# *same field names* as the Rust arm's rows, because `ancombc2-sim summarise`
# pools both from the vectors and would silently mis-pair a differently named
# field. The per-taxon vectors are written in the same order as `truth.tsv`'s
# taxa *restricted to the retained ones*, which is the order `ancombc2()`'s
# result rows use.
#
# It writes the numbers whether or not each replicate succeeded, recording the
# failure as `error`, so a replicate that R cannot analyse stays in the
# denominator instead of quietly improving the average.

suppressWarnings(suppressMessages({
  library(utils)
}))

# ---- arguments -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  if (i == length(args)) stop("--", name, " needs a value")
  args[[i + 1L]]
}
data_dir <- arg("--data")
out_path <- arg("--out")
grid_path <- arg("--grid")
if (is.null(data_dir) || is.null(out_path)) {
  stop("usage: sim_r.R --data DIR --out results.jsonl [--grid grid.json]")
}

# The analysis-level factors come from the grid file, which is the same file the
# Rust arm reads. They used to be literals here, and they did not match what the
# Rust arm was doing: `struc_zero`/`neg_lb` were hardcoded TRUE here and left at
# `AncombcConfig::default()`'s FALSE there, so at 90% zero inflation the oracle
# retained 11 of 500 taxa and Rust retained 500. Comparing two analyses that are
# not the same analysis produces a number, and the number means nothing.
#
# So this reads them. `jsonlite` is used here and nowhere else in the oracle path;
# the replicate *rows* are still written by hand, because that writer has no
# dependency at all and should keep not having one.
grid <- if (is.null(grid_path)) {
  list()
} else {
  # `fromJSON` takes a file path, but only for a file: passing a directory fails
  # with a message that does not name the argument, which is a poor way to learn
  # that a `--grid` was pointed at the wrong thing.
  if (!file.exists(grid_path) || dir.exists(grid_path)) {
    stop("--grid must be a JSON file, not a directory: ", grid_path)
  }
  jsonlite::fromJSON(grid_path, simplifyVector = TRUE)
}
gval <- function(key, default) {
  v <- grid[[key]]
  if (is.null(v)) default else v
}
if (!dir.exists(data_dir)) stop("no such directory: ", data_dir)

# `jsonlite` is not a dependency here, and the replicate files carry four or
# five scalar fields, so they are read with a regex rather than a JSON parser.
# The writer is `ancombc2-sim generate`, so the format is fixed; if it ever
# changes, this fails loudly rather than reading zeros.
read_replicate_meta <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = " ")
  one <- function(key) {
    m <- regmatches(txt, regexpr(sprintf('"%s"[[:space:]]*:[[:space:]]*(-?[0-9.eE+]+)', key), txt))
    if (!length(m)) stop("no key ", key, " in ", path)
    as.numeric(sub(sprintf('"%s"[[:space:]]*:[[:space:]]*', key), "", m))
  }
  oneb <- function(key) {
    m <- regmatches(txt, regexpr(sprintf('"%s"[[:space:]]*:[[:space:]]*(true|false)', key), txt))
    if (!length(m)) stop("no key ", key, " in ", path)
    identical(sub(sprintf('"%s"[[:space:]]*:[[:space:]]*', key), "", m), "true")
  }
  list(
    cell = one("cell"),
    rep = one("rep"),
    n_taxa = one("n_taxa"),
    n_samp = one("n_samp"),
    da_proportion = one("da_proportion"),
    log_fc = one("log_fc"),
    zero_inflation = one("zero_inflation"),
    lib_mean = one("lib_mean"),
    lib_cv = one("lib_cv"),
    confound = oneb("confound")
  )
}

# ---- a minimal JSON writer -------------------------------------------------
# Only what the row schema needs: strings, numbers, booleans, and arrays of each.
# `NA` is written as `null` so a missing scalar matches the Rust arm's `null`,
# which is how `ancombc2-sim summarise` recognises an undefined rate rather than
# a zero.
jnum <- function(x) {
  if (length(x) == 0L) return("null")
  if (is.na(x[1])) return("null")
  if (is.infinite(x[1])) return("null")
  sprintf("%.17g", as.numeric(x[1]))
}
jvec <- function(x) {
  if (length(x) == 0L) return("[]")
  paste0("[", paste(vapply(x, jnum, ""), collapse = ","), "]")
}
jbool <- function(x) if (isTRUE(x[1])) "true" else "false"
jstr <- function(x) paste0("\"", gsub("\"", "\\\\\"", x[1]), "\"")

# The row schema. The field names must match the Rust arm's exactly: `summarise`
# deserialises both arms into the same struct.
blank_row <- function(m, arm) {
  list(
    grid = "full", cell = m$cell, rep = m$rep, arm = arm,
    n_taxa = m$n_taxa, n_samp = m$n_samp,
    da_proportion = m$da_proportion, log_fc = m$log_fc,
    zero_inflation = m$zero_inflation, lib_mean = m$lib_mean,
    lib_cv = m$lib_cv, confound = m$confound,
    names = character(0), truth_log_fc = numeric(0), is_da = logical(0),
    has_effect = logical(0), confounded = logical(0),
    beta = numeric(0), se = numeric(0), p = numeric(0), q = numeric(0),
    diff_abn = logical(0),
    n_retained = 0, n_da_retained = 0, n_diff_abn = 0,
    n_confounded_retained = 0, confounded_false_positive_rate = NA_real_,
    empirical_fdr = NA_real_, power = NA_real_,
    sign_concordance = NA_real_, null_sign_rate = NA_real_,
    lfc_bias = NA_real_, lfc_rmse = NA_real_, lfc_mae = NA_real_,
    f1 = NA_real_, jaccard = NA_real_,
    error = NULL
  )
}

render_row <- function(r) {
  paste0(
    "{",
    "\"grid\":", jstr(r$grid),
    ",\"cell\":", jnum(r$cell),
    ",\"rep\":", jnum(r$rep),
    ",\"arm\":", jstr(r$arm),
    ",\"n_taxa\":", jnum(r$n_taxa),
    ",\"n_samp\":", jnum(r$n_samp),
    ",\"da_proportion\":", jnum(r$da_proportion),
    ",\"log_fc\":", jnum(r$log_fc),
    ",\"zero_inflation\":", jnum(r$zero_inflation),
    ",\"lib_mean\":", jnum(r$lib_mean),
    ",\"lib_cv\":", jnum(r$lib_cv),
    ",\"confound\":", jbool(r$confound),
    ",\"names\":[", paste(vapply(r$names, jstr, ""), collapse = ","), "]",
    ",\"truth_log_fc\":", jvec(r$truth_log_fc),
    ",\"is_da\":[", paste(vapply(r$is_da, jbool, ""), collapse = ","), "]",
    ",\"has_effect\":[", paste(vapply(r$has_effect, jbool, ""), collapse = ","), "]",
    ",\"confounded\":[", paste(vapply(r$confounded, jbool, ""), collapse = ","), "]",
    ",\"beta\":", jvec(r$beta),
    ",\"se\":", jvec(r$se),
    ",\"p\":", jvec(r$p),
    ",\"q\":", jvec(r$q),
    ",\"diff_abn\":[", paste(vapply(r$diff_abn, jbool, ""), collapse = ","), "]",
    ",\"n_retained\":", jnum(r$n_retained),
    ",\"n_da_retained\":", jnum(r$n_da_retained),
    ",\"n_diff_abn\":", jnum(r$n_diff_abn),
    ",\"n_confounded_retained\":", jnum(r$n_confounded_retained),
    ",\"confounded_false_positive_rate\":", jnum(r$confounded_false_positive_rate),
    ",\"empirical_fdr\":", jnum(r$empirical_fdr),
    ",\"power\":", jnum(r$power),
    ",\"sign_concordance\":", jnum(r$sign_concordance),
    ",\"null_sign_rate\":", jnum(r$null_sign_rate),
    ",\"lfc_bias\":", jnum(r$lfc_bias),
    ",\"lfc_rmse\":", jnum(r$lfc_rmse),
    ",\"lfc_mae\":", jnum(r$lfc_mae),
    ",\"f1\":", jnum(r$f1),
    ",\"jaccard\":", jnum(r$jaccard),
    ",\"error\":", if (is.null(r$error)) "null" else jstr(r$error),
    "}"
  )
}

# ---- ANCOMBC ----------------------------------------------------------------

# The oracle is *sourced* from the pinned checkout, not installed as a package:
# it declares R >= 4.5.0 and the interpreter here is 4.3.3, and the analysis
# needs the unexported internals. `reference/R/harness.R` is the same loader the
# golden contract and the benchmark arm use, so the simulation's R arm, the
# goldens, and the benchmarks all run identical oracle code. See
# reference/env/ORACLE.md.
oracle_dir <- Sys.getenv("ANCOMBC_ORACLE_DIR", unset = "reference/ANCOMBC")
harness_dir <- file.path("reference", "R")
for (f in c("stubs.R", "harness.R")) {
  if (!file.exists(file.path(harness_dir, f))) {
    stop("cannot find ", file.path(harness_dir, f),
         "; run this script from the repository root")
  }
  source(file.path(harness_dir, f))
}
H <- tryCatch(
  load_harness(oracle_dir),
  error = function(e) {
    stop("cannot load the oracle from ", oracle_dir, ": ", conditionMessage(e),
         "\nSet ANCOMBC_ORACLE_DIR to the pinned ANCOMBC 2.15.2 checkout.")
  }
)

# ---- one replicate ----------------------------------------------------------

analyse_one <- function(dir) {
  meta <- read_replicate_meta(file.path(dir, "replicate.json"))
  r <- blank_row(meta, "r")

  # A matrix, not a data frame: the oracle's internal helpers dispatch on the
  # type, and a data frame makes `is.infinite()` fail with "default method not
  # implemented for type 'list'" deep inside the structural-zero scan. The
  # benchmark arm coerces for the same reason.
  counts <- as.matrix(utils::read.delim(file.path(dir, "counts.tsv"), row.names = 1,
                                        check.names = FALSE))
  storage.mode(counts) <- "double"

  truth <- utils::read.delim(file.path(dir, "truth.tsv"), stringsAsFactors = FALSE)
  if (nrow(truth) != nrow(counts)) {
    stop("truth.tsv has ", nrow(truth), " taxa but counts.tsv has ", nrow(counts),
         ": ", dir)
  }
  # The generator writes the taxa in table order. This is checked rather than
  # assumed, because the metrics join on the retained taxa's order and a
  # mis-ordered truth file would silently shift every truth value by a row.
  if (!identical(as.character(truth$taxon), rownames(counts))) {
    stop("truth.tsv taxa do not match counts.tsv rows, or are out of order: ", dir)
  }
  # The sample column becomes the row names, exactly as the benchmark arm reads
  # its metadata: the oracle matches metadata rows to count-table columns by
  # name, and a data frame with a `sample` column and no row names fails that
  # check even when the samples agree.
  meta_df <- utils::read.delim(file.path(dir, "meta.tsv"), row.names = 1,
                               check.names = FALSE)
  # The generator writes group1/group2; the oracle takes a bare group name and
  # treats the first level as the reference, which is the control group, so the
  # reported `group` coefficient is the treated-minus-control contrast.
  meta_df$group <- factor(sub("^group", "", meta_df$group), levels = c("1", "2"))

  # `ref_run` is the same entry point the golden contract and the benchmark arm
  # use, so all three exercise identical oracle code. The defaults are spelled
  # out rather than inherited, because a default that drifts between the oracle
  # versions would silently change this arm's results.
  res <- try(suppressWarnings(suppressMessages(
    H$ref_run(
      counts,
      meta_df,
      "group",
      p_adj_method = gval("p_adjust", "BH"),
      pseudo = gval("pseudo", 0.5),
      pseudo_sens = FALSE,
      conservative = TRUE,
      prv_cut = gval("prevalence", 0.0),
      lib_cut = gval("lib_size", 0.0),
      s0_perc = gval("s0_perc", 0.05),
      group = "group",
      struc_zero = gval("struc_zero", TRUE),
      neg_lb = gval("neg_lb", TRUE),
      alpha = gval("alpha", 0.05),
      global = FALSE,
      pairwise = FALSE
    )
  )), silent = TRUE)
  if (inherits(res, "try-error")) {
    r$error <- trimws(as.character(res))
    return(r)
  }

  # `ref_run` returns per-taxon by per-coefficient matrices, not a long frame:
  # `beta` and `se` are `n_retained x n_coef` and `diff_abn` is a logical matrix
  # of the same shape. The estimand is the group *column*, which is found by name
  # rather than by position: the intercept is first only by convention, and a
  # silent column shift would compare the intercept to a log fold change.
  taxa <- as.character(res$taxa_retained)
  coef_names <- colnames(res$beta)
  # The generator's group levels are `1` and `2`, so R names the treatment
  # coefficient `group2`.
  group_col <- match("group2", coef_names)
  if (is.na(group_col)) {
    r$error <- paste0("no group2 coefficient among: ", paste(coef_names, collapse = ", "))
    return(r)
  }
  keep <- match(taxa, truth$taxon)
  if (anyNA(keep)) {
    r$error <- "the retained taxa are not all present in truth.tsv"
    return(r)
  }
  tr <- truth[keep, , drop = FALSE]

  r$names <- taxa
  r$truth_log_fc <- as.numeric(tr$log_fc)
  r$is_da <- as.logical(tr$is_da)
  # `has_effect` and `confounded` are separate from the label: in a confounded
  # cell the labelled taxa have no compositional effect, only a group-dependent
  # zero rate. Counting the label would score the confounded taxa as true
  # positives, which is indistinguishable from the sampling-fraction correction
  # doing nothing.
  if (is.null(tr$has_effect)) {
    r$error <- "truth.tsv has no `has_effect` column; regenerate it with `ancombc2-sim generate`"
    return(r)
  }
  r$has_effect <- as.logical(tr$has_effect)
  r$confounded <- if (is.null(tr$confounded)) rep(FALSE, nrow(tr)) else as.logical(tr$confounded)
  r$beta <- as.numeric(res$beta[, group_col])
  r$se <- as.numeric(res$se[, group_col])
  r$p <- as.numeric(res$p[, group_col])
  r$q <- as.numeric(res$q[, group_col])
  # The treatment contrast of the pairwise test is what `diff_abn` decides; the
  # intercept row of that matrix is never a divergence call.
  r$diff_abn <- as.logical(res$diff_abn[, group_col])

  n <- length(r$beta)
  r$n_retained <- n
  r$n_da_retained <- sum(r$has_effect)
  called <- !is.na(r$diff_abn) & r$diff_abn
  r$n_diff_abn <- sum(called)
  r$empirical_fdr <- if (r$n_diff_abn > 0) sum(called & !r$has_effect) / r$n_diff_abn else NA_real_
  r$power <- if (r$n_da_retained > 0) sum(called & r$has_effect) / r$n_da_retained else NA_real_

  # The direct test of the sampling-fraction correction: the confounded taxa have
  # no compositional effect, so any call on one is a false positive.
  n_conf <- sum(r$confounded)
  r$n_confounded_retained <- n_conf
  r$confounded_false_positive_rate <-
    if (n_conf > 0) sum(called & r$confounded) / n_conf else NA_real_

  ok <- is.finite(r$beta)
  da <- r$has_effect & !is.na(r$has_effect)
  r$sign_concordance <- if (any(da)) {
    mean(ifelse(r$truth_log_fc[da] > 0, r$beta[da] > 0, r$beta[da] < 0))
  } else NA_real_
  nul <- !da & ok
  r$null_sign_rate <- if (any(nul)) mean(r$beta[nul] > 0) else NA_real_

  err <- r$beta[ok] - r$truth_log_fc[ok]
  r$lfc_bias <- if (length(err)) mean(err) else NA_real_
  r$lfc_rmse <- if (length(err)) sqrt(mean(err^2)) else NA_real_
  r$lfc_mae <- if (length(err)) mean(abs(err)) else NA_real_

  tp <- sum(called & da); fp <- sum(called & !da); fn <- sum(!called & da)
  r$f1 <- if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_
  r$jaccard <- if (tp + fp + fn > 0) tp / (tp + fp + fn) else NA_real_
  r
}

# ---- the run ----------------------------------------------------------------

# The output directory is created here rather than assumed: a CI job or a `make`
# target that has not written to this path before should not fail on `file()`,
# and failing at the point of first write -- after the analysis is set up -- is
# a worse place to discover a missing directory.
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

dirs <- sort(list.dirs(data_dir, recursive = TRUE, full.names = TRUE))
dirs <- dirs[file.exists(file.path(dirs, "counts.tsv"))]

# `--cells 126,158` restricts the run to those cell indices, and `--reps N` to the
# first N replicates of each. Both exist because the full grid is 252,000
# replicates and about 39 hours of R; answering "does the oracle agree with Rust
# on *this* cell?" should not require waiting for the whole surface. The
# replicate directories carry their own `cell` in `replicate.json`, so this reads
# identity from the data rather than inferring it from the path.
cells_arg <- arg("--cells")
reps_arg <- arg("--reps")
if (!is.null(cells_arg)) {
  keep_cells <- as.integer(strsplit(cells_arg, ",", fixed = TRUE)[[1]])
  if (anyNA(keep_cells)) stop("--cells must be a comma-separated list of integers")
  meta_of <- vapply(dirs, function(d) {
    m <- regmatches(paste(readLines(file.path(d, "replicate.json"), warn = FALSE),
                          collapse = " "),
                    regexpr('"cell"[[:space:]]*:[[:space:]]*[0-9]+',
                            paste(readLines(file.path(d, "replicate.json"),
                                            warn = FALSE), collapse = " ")))
    if (!length(m)) NA_integer_ else as.integer(sub('.*:[[:space:]]*', "", m))
  }, integer(1), USE.NAMES = FALSE)
  dirs <- dirs[!is.na(meta_of) & meta_of %in% keep_cells]
  cat(sprintf("restricted to %d cell(s): %s\n", length(unique(meta_of[!is.na(meta_of) & meta_of %in% keep_cells])),
              cells_arg))
}
if (!is.null(reps_arg)) {
  keep_reps <- as.integer(reps_arg)
  meta_of2 <- vapply(dirs, function(d) {
    txt <- paste(readLines(file.path(d, "replicate.json"), warn = FALSE), collapse = " ")
    m <- regmatches(txt, regexpr('"rep"[[:space:]]*:[[:space:]]*[0-9]+', txt))
    if (!length(m)) NA_integer_ else as.integer(sub('.*:[[:space:]]*', "", m))
  }, integer(1), USE.NAMES = FALSE)
  dirs <- dirs[!is.na(meta_of2) & meta_of2 < keep_reps]
}
if (!length(dirs)) {
  stop("no replicate directories under ", data_dir,
       "; run `ancombc2-sim generate` first")
}
# --shard K --nshards N takes every Nth replicate directory, counting from
# --skip. The replicates are independent by construction -- `simulate(cell, seed,
# rep)` is a pure function of its three arguments, and the oracle is sequential
# per replicate because the reference harness stubs `foreach` -- so sharding
# changes nothing about any row except which process wrote it. The oracle cannot
# be parallelised *internally* without changing its behaviour, which is the point:
# this parallelises around it rather than inside it.
#
# It is worth doing because the full grid is 252,000 replicates at about 1.7 per
# second on one core -- 39 hours. Sharded across the cores this host has, it is
# hours.
shard <- as.integer(arg("--shard", "1"))
nshards <- as.integer(arg("--nshards", "1"))
skip <- as.integer(arg("--skip", "0"))
if (is.na(nshards) || nshards < 1L) stop("--nshards must be at least 1")
if (is.na(shard) || shard < 1L || shard > nshards) {
  stop("--shard must be in 1..--nshards")
}
if (!is.na(skip) && skip > 0L) dirs <- dirs[-seq_len(min(skip, length(dirs)))]
if (nshards > 1L) dirs <- dirs[seq(shard, length(dirs), by = nshards)]
cat(sprintf("r arm: %d replicate(s) under %s (shard %d/%d, skipped %d)\n",
            length(dirs), data_dir, shard, nshards, skip))
con <- file(out_path, open = "wt")
on.exit(close(con), add = TRUE)
t0 <- Sys.time()
for (i in seq_along(dirs)) {
  row <- try(analyse_one(dirs[[i]]), silent = TRUE)
  if (inherits(row, "try-error")) {
    # A failure that is not about the data still has to be recorded against the
    # right cell, so the replicate's own metadata is read for identity.
    meta <- try(read_replicate_meta(file.path(dirs[[i]], "replicate.json")), silent = TRUE)
    row <- if (inherits(meta, "try-error")) {
      blank_row(list(cell = NA, rep = NA, n_taxa = NA, n_samp = NA,
                     da_proportion = NA, log_fc = NA, zero_inflation = NA,
                     lib_mean = NA, lib_cv = NA, confound = FALSE), "r")
    } else blank_row(meta, "r")
    row$error <- trimws(as.character(row))
  }
  writeLines(render_row(row), con)
  if (i %% 200L == 0L || i == length(dirs)) {
    cat(sprintf("  %d/%d (%.0fs)\n", i, length(dirs),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush(stdout())
  }
}
cat(sprintf("wrote %d row(s) to %s in %.0fs\n", length(dirs), out_path,
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
