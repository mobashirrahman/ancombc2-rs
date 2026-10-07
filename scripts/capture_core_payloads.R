#!/usr/bin/env Rscript
# Capture the pinned original's INTERNAL `.ancombc2_core` payloads.
#
#   Rscript --vanilla scripts/capture_core_payloads.R --profile <id> \
#       --case <input.rds> --out validation/exact/payloads/<case>.rds
#
# Why this exists, and why it is a separate generator
# ---------------------------------------------------
# `.ancombc2_core` returns *assembled tables*, not the numerical arrays they were
# assembled from. `beta_hat` survives as the `lfc_*` columns of `res` and `se_hat`
# as `se_*`, but `dof`, `vcov_hat`, `var_delta` and `s02` are in none of the ten
# returned fields and cannot be recovered from them.
#
# So the acceptance test for S07 -- "a transport/assembly test with original
# numerical payloads passes whole-result bytes" -- needs the arrays themselves,
# and the only place they exist is inside a running `.ancombc2_core`.
#
# How, without reimplementing it
# -------------------------------
# The pinned sources are sourced into an environment (as `reference/R/oracle.R`
# already does for the legacy harness), and `.ancombc2_core`'s **last statement**
# -- `return(out)` -- is replaced with a `return()` that also hands back the
# internals. Nothing else about the function is touched, so the assembly code being
# tested downstream is the original's, and the numbers are the original's.
#
# The division of responsibility this preserves
# ---------------------------------------------
# This script may read the pinned source. Nothing the *candidate* runs may.
# `scripts/exact_runner.R` and `r/ANCOMBC/R/*.R` contain no reference to
# `reference/` or to `validation/exact/payloads/`, and the generated payloads live
# in a directory that only the acceptance test reads.

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
`%||%` <- function(a, b) if (is.null(a)) b else a

for (k in c("out", "case")) if (is.null(opt[[k]])) {
  stop("usage: capture_core_payloads.R --out <rds> --case <input.rds> [--tag <name>]",
       call. = FALSE)
}
tag <- opt$tag %||% tools::file_path_sans_ext(basename(opt[["case"]]))

# ---- source the pinned R files --------------------------------------------
# Order matters only in that `ancombc_prep.R` is where `.ancombc2_core` lives;
# the others supply what it calls.
oracle_env <- new.env(parent = globalenv())
pinned <- file.path("reference", "ANCOMBC", "R")
for (f in c("utils.R", "ancombc_prep.R", "ancombc_bias_correct.R", "ancombc_mult.R")) {
  sys.source(file.path(pinned, f), envir = oracle_env)
}

# `.ancombc2_core` uses `%dorng%` and `%dopar%`. Sourcing the pinned files does
# not bring those in -- the upstream NAMESPACE's `importFrom` only applies to an
# *installed* package -- so they have to be supplied here.
#
# The **real** operators are used, not the sequential stubs in
# `reference/R/stubs.R`. `registerDoSEQ()` is what the original itself calls when
# `n_cl = 1` (`ancombc2.R:421`), so a sequential backend is the original's own
# one-worker configuration rather than a substitute for it; the stubs would make a
# stub the reference for the fixed-effects path, which is the one thing
# IMPROVED_PLAN.md S15 names as invalid.
if (!requireNamespace("foreach", quietly = TRUE) ||
    !requireNamespace("doRNG", quietly = TRUE) ||
    !requireNamespace("rngtools", quietly = TRUE)) {
  stop("the pinned `.ancombc2_core` needs foreach, doRNG and rngtools; ",
       "capturing payloads with a stubbed backend would not be the original's result",
       call. = FALSE)
}
assign("%dorng%", doRNG::`%dorng%`, envir = oracle_env)
assign("%dopar%", foreach::`%dopar%`, envir = oracle_env)
assign("%:%", foreach::`%:%`, envir = oracle_env)
assign("foreach", foreach::foreach, envir = oracle_env)
assign("registerDoSEQ", foreach::registerDoSEQ, envir = oracle_env)
assign("registerDoSEQ", NULL, envir = oracle_env)
foreach::registerDoSEQ()
# The original's own bookkeeping around `%dorng%`: it seeds per iteration so a
# sequential run and a parallel run agree. Reproduced so the captured payloads do
# not depend on ambient RNG state.
if (!exists(".Random.seed", envir = globalenv())) set.seed(42L)

# ---- patch only the return statement ---------------------------------------
core <- get(".ancombc2_core", envir = oracle_env)
b <- body(core)
last <- b[[length(b)]]
if (!identical(as.character(last[[1L]]), "return")) {
  stop("`.ancombc2_core` no longer ends in `return(out)`; the capture patch needs ",
       "to be updated rather than applied blindly. Its last statement is: ",
       paste(deparse(last), collapse = " "), call. = FALSE)
}
b[[length(b)]] <- quote(return(list(
  # The ten fields the original returns, exactly as it returns them. This is the
  # acceptance target.
  feature_table          = O2,
  bias_correct_log_table = y_bias_crt,
  samp_frac              = theta_hat,
  delta_em               = delta_em,
  delta_wls              = delta_wls,
  res                    = res,
  res_global             = res_global,
  res_pair               = res_pair,
  res_dunn               = res_dunn,
  res_trend              = res_trend,
  # The internals, for the transport and the assembly. Not part of the original's
  # return value and not part of any candidate result.
  .internals = list(
    O1            = O1,
    O2            = O2,
    x             = x,
    fix_eff       = fix_eff,
    tax_name      = tax_name,
    n_tax         = n_tax,
    beta_hat      = beta_hat,
    var_hat       = var_hat,
    dof           = dof,
    vcov_hat      = vcov_hat,
    y_bias_crt    = y_bias_crt,
    theta_hat     = theta_hat,
    delta_em      = delta_em,
    delta_wls      = delta_wls,
    var_delta     = var_delta,
    s02           = s02,
    meta_data     = meta_data,
    n_samp        = n_samp,
    p_adj_method  = p_adj_method,
    alpha         = alpha,
    global        = global,
    pairwise      = pairwise,
    dunnet        = dunnet,
    trend         = trend,
    group         = group,
    mdfdr_control = mdfdr_control,
    trend_control = trend_control,
    fix_formula   = fix_formula,
    rand_formula  = rand_formula
  )
)))
body(core) <- b
assign(".ancombc2_core", core, envir = oracle_env)

# ---- run it on a saved exact input ----------------------------------------
inp <- readRDS(opt[["case"]])
if (!identical(inp$schema, "ancombc2-exact-input/1")) {
  stop("not an ancombc2-exact-input/1 object", call. = FALSE)
}
# The core is a plain function call, not an `ancombc2()` call: `data_sanity_check`
# and `.data_core` have already run by the time `ancombc2()` reaches the core, and
# re-running them here would mean testing them twice.
# `data_sanity_check` sets `aggregate_data = data` when the caller omitted it
# (data_sanity_check.R:205-207), and `ancombc2()` then passes the *checked* pair.
# An omitted `aggregate_data` here would reach `.ancombc2_core` as NULL and
# `O2 = NULL + pseudo` is `numeric(0)`, so the capture must do what the sanity
# check does.
core_args <- c(list(data = inp$args$data,
                    aggregate_data = inp$args$aggregate_data %||% inp$args$data,
                    meta_data = inp$args$meta_data,
                    fix_formula = inp$args$fix_formula),
               inp$args[names(inp$args) %in%
                          c("rand_formula", "p_adj_method", "pseudo", "s0_perc",
                            "group", "alpha", "verbose", "global", "pairwise",
                            "dunnet", "trend", "iter_control", "em_control",
                            "lme_control", "mdfdr_control", "trend_control")])

t0 <- proc.time()[["elapsed"]]
res <- suppressMessages(do.call(get(".ancombc2_core", envir = oracle_env), core_args))
t1 <- proc.time()[["elapsed"]]

expected <- res[setdiff(names(res), ".internals")]
int <- res$.internals

# ---- record the shapes, so the transport's plan can be checked against them -
describe <- function(x) {
  list(type = typeof(x), length = length(x),
       dim = if (!is.null(dim(x))) as.integer(dim(x)) else NULL,
       class = class(x))
}
inventory <- list(
  beta_hat  = describe(int$beta_hat),
  var_hat   = describe(int$var_hat),
  dof       = describe(int$dof),
  vcov_hat  = list(type = typeof(int$vcov_hat), length = length(int$vcov_hat),
                   per_item = describe(int$vcov_hat[[1L]])),
  y_bias_crt = describe(int$y_bias_crt),
  theta_hat = describe(int$theta_hat),
  delta_em  = describe(int$delta_em),
  delta_wls = describe(int$delta_wls),
  var_delta = describe(int$var_delta),
  s02       = describe(int$s02),
  O2        = describe(int$O2),
  x         = describe(int$x)
)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  schema = "ancombc2-exact-core-payloads/1",
  tag = tag,
  case = opt[["case"]],
  generator = "scripts/capture_core_payloads.R",
  provenance = list(
    oracle = "ANCOMBC 2.15.2 @ dc4febdf59badb3a8dfe0c767ef2186323c2199a",
    source = "reference/ANCOMBC/R/{utils,ancombc_prep,ancombc_bias_correct,ancombc_mult}.R",
    patch = "`.ancombc2_core`'s final `return(out)` replaced by a return that also exposes the internals; nothing else changed",
    warning = "These payloads are test fixtures. They must never be read by a candidate run, and using them is not evidence of Rust numerical parity."
  ),
  inventory = inventory,
  elapsed_s = t1 - t0,
  internals = int,
  expected = expected
), opt$out, version = 3)

cat(sprintf("captured %s in %.1fs: n_tax=%d n_samp=%d p=%d -> %s\n", tag, t1 - t0,
            int$n_tax, int$n_samp, length(int$fix_eff), opt$out))
for (k in names(inventory)) {
  v <- inventory[[k]]
  cat(sprintf("  %-12s %-8s %s\n", k, v$type,
              if (!is.null(v$dim)) paste(v$dim, collapse = "x") else
              if (!is.null(v$per_item)) paste0(v$length, " x ", paste(v$per_item$dim, collapse = "x")) else
              paste0("len ", v$length)))
}
