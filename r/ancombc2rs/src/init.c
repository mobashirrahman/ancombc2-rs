#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <string.h>

/* The Rust side allocates with libc's malloc, so the returned pointer is released
 * with free() and not with R_free(). R_free is only correct for memory R itself
 * allocated. */
extern char *ancombc2_rs_run(const char *request);
extern char *ancombc2_rs_version(void);
extern void ancombc2_rs_free_string(char *s);

/* Take ownership of a Rust-allocated string and hand R a CHARSXP. The pointer is
 * released either way: on the success path R has its own copy, and on the error
 * path the memory would otherwise leak. */
static SEXP take_r_string(char *s) {
  if (s == NULL) {
    return Rf_mkString("ancombc2-rs returned a null pointer");
  }
  SEXP out = PROTECT(Rf_mkString(s));
  ancombc2_rs_free_string(s);
  UNPROTECT(1);
  return out;
}

SEXP C_ancombc2_rs_run(SEXP request) {
  if (TYPEOF(request) != STRSXP || Rf_length(request) != 1) {
    Rf_error("the request must be a character vector of length 1");
  }
  const char *text = CHAR(STRING_ELT(request, 0));
  if (text == NULL) {
    Rf_error("the request string is NA");
  }
  return take_r_string(ancombc2_rs_run(text));
}

SEXP C_ancombc2_rs_version(void) {
  return take_r_string(ancombc2_rs_version());
}

/* A panic in the Rust core is caught at the FFI boundary and returned as an
 * `error` key, so a stack-overflow or an out-of-bounds read cannot abort the R
 * session. These entry points exist so `.Call` registration is explicit rather
 * than relying on symbol lookup, which would let a renamed Rust function become
 * a runtime error instead of a build error. */
static const R_CallMethodDef CallEntries[] = {
  {"C_ancombc2_rs_run",    (DL_FUNC) &C_ancombc2_rs_run,    1},
  {"C_ancombc2_rs_version",(DL_FUNC) &C_ancombc2_rs_version,0},
  {NULL, NULL, 0}
};

void R_init_ancombc2rs(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);
}
