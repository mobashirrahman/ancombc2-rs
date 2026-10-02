# Sequential stand-ins for the foreach/doRNG machinery used by the oracle.
#
# The oracle evaluates expressions of the form
#
#     foreach(i = seq_len(n), .combine = rbind) %dorng% { ... }
#
# i.e. `foreach(...)` builds a loop specification and `%dorng%` iterates it. With
# the package's own `n_cl = 1` default, `foreach::registerDoSEQ()` is called, so
# the real `%dorng%` runs the loop sequentially; the RNG stream is irrelevant for
# the fixed-effects path because no routine it wraps draws random numbers
# (`.iter_mle` consumes `rnorm` only inside a discarded smoke-test fit).
#
# These stubs therefore reproduce `n_cl = 1` behaviour exactly.

foreach_stub <- function(..., .combine = NULL) {
  args <- list(...)
  if (length(args) != 1L) stop("stub supports a single loop variable only")
  nm <- names(args)
  if (is.null(nm) || !nzchar(nm)) stop("stub requires a named loop variable")
  list(var = nm, value = args[[1L]], combine = .combine)
}

# The body arrives as a promise; substitute it once, then evaluate the captured
# expression once per value with the loop variable bound in the caller's frame.
`%dorng%` <- function(obj, expr) {
  # `expr` is a promise; substitute once to get the unevaluated body, then
  # evaluate that language object once per value. (as.list() here would turn the
  # call into a plain list, which eval() returns as-is.)
  body <- substitute(expr)
  env <- parent.frame()
  vals <- obj$value
  out <- vector("list", length(vals))
  for (k in seq_along(vals)) {
    assign(obj$var, vals[[k]], envir = env)
    out[[k]] <- eval(body, envir = env)
  }
  if (is.null(out)) return(NULL)
  if (!is.null(obj$combine)) {
    # `.combine` arrives either as a symbol (rbind / cbind) or, because
    # `foreach_stub` captures arguments by value, as the function object itself.
    cmb <- obj$combine
    if (is.character(cmb)) cmb <- get(cmb, envir = env)
    if (identical(cmb, base::rbind)) return(do.call(rbind, out))
    if (identical(cmb, base::cbind)) return(do.call(cbind, out))
    return(Reduce(cmb, out))
  }
  out
}

registerDoSEQ <- function(...) invisible(NULL)

# Called by `ancombc2()` when n_cl == 1.
stub_cluster <- function(n_cl) {
  if (n_cl > 1L) stop("golden generation runs with n_cl = 1 only")
  invisible(NULL)
}
