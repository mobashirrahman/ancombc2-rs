/* The R side of the typed native bridge (IMPROVED_PLAN.md S06).
 *
 * Everything here is about three things and nothing else:
 *
 *   1. **R owns every buffer.** Each output vector is allocated here with
 *      `Rf_alloc*`, `PROTECT`ed, filled by Rust, and handed to R. Rust never
 *      allocates anything R has to free, and never returns a pointer R has to
 *      release. That is the whole ownership story, and it is why the C ABI has no
 *      free function.
 *
 *   2. **Types and shapes are checked before the bridge is called.** A wrong
 *      REALSXP where an INTSXP belongs, or a length that disagrees with a dim,
 *      is an R error with a message naming the argument -- not a read past the end
 *      of a buffer on the other side of an FFI boundary.
 *
 *   3. **A failure becomes an R condition.** The bridge returns a status code and
 *      leaves a message in its own buffer; this file raises that message with
 *      `Rf_error`. Nothing unwinds through the FFI boundary: a Rust panic is
 *      caught in Rust (`catch_unwind`) and arrives here as
 *      `ANCOMBC2_RB_ERR_PANIC`, so it cannot take the R session with it.
 *
 * R's NA representation
 * ---------------------
 * `NA_integer_` is `INT_MIN` and `NA_real_` is a *specific* NaN payload,
 * 0x7FF00000000007A2. Both are carried through untouched: the integer buffer is
 * memcpy'd as i32 and the real buffer as f64, so the payloads survive verbatim.
 * Nothing here tests `R_IsNA` and reconstructs, which is how a `NaN` becomes an
 * `NA` and `-0.0` becomes `0.0`.
 *
 * Names
 * -----
 * A `STRSXP` is a vector of `CHARSXP`, which means reaching into R's internal
 * pointer table to read one. That is R's business, so this file does it: names are
 * packed into a flat byte buffer with offsets, handed to Rust as bytes, and
 * unpacked back with `Rf_mkChar` on the way out. Rust never sees a SEXP.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Status codes, mirrored in crates/ancombc2-rbridge/src/abi.rs. Kept as an enum
 * rather than bare numbers so a mismatch is a compile error here rather than a
 * silent wrong answer at run time. */
enum {
  RB_OK                 = 0,
  RB_ERR_ARGS           = 1,
  RB_ERR_STAGE          = 2,
  RB_ERR_PANIC          = 3,
  RB_ERR_STAGE_UNKNOWN  = 4
};

/* The probe selector, mirrored in crates/ancombc2-rbridge/src/abi.rs. */
enum {
  ECHO_DATA_INT         = 0,
  ECHO_DATA_REAL        = 1,
  ECHO_AGGREGATE        = 2,
  ECHO_DESIGN           = 3,
  ECHO_DESIGN_COMPLETE  = 4,
  ECHO_GROUP_INDEX      = 5
};

/* Which kind of R vector a probe's payload should be built as, decided by the
 * bridge rather than guessed here. The aggregate probe in particular echoes
 * whichever type it was given, so this file cannot assume. */
enum {
  KIND_NONE  = 0,
  KIND_INT   = 1,
  KIND_REAL  = 2,
  KIND_RAW   = 3
};

extern int32_t ancombc2_rb_probe(
    const void *req, int32_t which,
    int32_t *out_int, double *out_real, uint8_t *out_u8,
    int64_t *out_n, int64_t out_cap, int32_t *out_kind);
extern int64_t ancombc2_rb_last_error_len(void);
extern int64_t ancombc2_rb_copy_last_error(uint8_t *dst, int64_t cap);
struct rb_reducer;
extern int32_t ancombc2_rb_reducer(struct rb_reducer *out);
struct rb_stage;
struct rb_preprocess;
extern int32_t ancombc2_rb_preprocess_probe(
    const struct rb_preprocess *in, double *f64_out, int32_t *i32_out,
    struct rb_stage *out, int64_t slot_count, int64_t *out_f64_len,
    int64_t *out_i32_len, int64_t *out_n_tax, int64_t *out_n_samp,
    int32_t *out_has_group, int8_t *out_names, int64_t *out_name_offsets,
    int64_t *out_names_len);
extern int32_t ancombc2_rb_emit_plan(const void *pay,
                                     int64_t *f64_len, int64_t *i32_len,
                                     int64_t *u8_len, int64_t *plan_len);
extern int32_t ancombc2_rb_emit(const void *pay,
                                double *out_f64, int32_t *out_i32,
                                int64_t *offsets_f64, int64_t *offsets_i32,
                                int64_t *plan_names_offsets, uint8_t *plan_names_bytes,
                                int64_t plan_len);
extern int64_t ancombc2_rb_oracle_sha(uint8_t *dst, int64_t cap);
extern int64_t ancombc2_rb_version(uint8_t *dst, int64_t cap);

/* Raise the bridge's last error as an R condition. `code` goes in the call so a
 * panic is distinguishable in the traceback from a validation failure. */
static void rb_error(int32_t code, const char *what) {
  int64_t n = ancombc2_rb_last_error_len();
  char *buf = (char *) R_alloc((size_t) n + 1, sizeof(char));
  ancombc2_rb_copy_last_error((uint8_t *) buf, n + 1);
  const char *kind =
    code == RB_ERR_PANIC          ? "panicked" :
    code == RB_ERR_STAGE         ? "stage failure" :
    code == RB_ERR_STAGE_UNKNOWN ? "unknown probe" : "invalid argument";
  Rf_error("ancombc2-rbridge: %s in %s: %s", kind, what, buf);
}

/* ------------------------------------------------------------------------ */
/* packed names                                                              */
/* ------------------------------------------------------------------------ */

typedef struct {
  int64_t *offsets;  /* n + 1 */
  uint8_t *bytes;    /* offsets[n] */
  int64_t  n;
  SEXP     keep;    /* the STRSXP, kept alive for the call */
} packed_names;

static void pack_names(packed_names *out, SEXP x, const char *field) {
  out->keep = x;
  if (TYPEOF(x) != STRSXP) {
    Rf_error("`%s` must be a character vector, not %s", field, Rf_type2char(TYPEOF(x)));
  }
  out->n = Rf_xlength(x);
  out->offsets = (int64_t *) R_alloc((size_t) out->n + 1, sizeof(int64_t));
  int64_t total = 0;
  for (int64_t i = 0; i < out->n; i++) {
    out->offsets[i] = total;
    SEXP e = STRING_ELT(x, i);
    if (e == NA_STRING) {
      Rf_error("`%s`[%lld] is NA; every name must be present", field, (long long) (i + 1));
    }
    /* `translateCharUTF8` gives the UTF-8 bytes; `strlen` on them is the byte
     * length, which is what the offsets must count. Not `LENGTH`, which counts
     * characters. */
    const char *s = Rf_translateCharUTF8(e);
    total += (int64_t) strlen(s);
  }
  out->offsets[out->n] = total;
  out->bytes = (uint8_t *) R_alloc((size_t) (total > 0 ? total : 1), sizeof(uint8_t));
  int64_t at = 0;
  for (int64_t i = 0; i < out->n; i++) {
    const char *s = Rf_translateCharUTF8(STRING_ELT(x, i));
    size_t len = strlen(s);
    memcpy(out->bytes + at, s, len);
    at += (int64_t) len;
  }
}

/* Rebuild a STRSXP from a packed buffer. `Rf_mkChar` takes UTF-8 and marks it
 * accordingly, so a name that was UTF-8 on the way in is UTF-8 on the way out
 * with the encoding flag set -- which is what `serialize()` compares. */
static SEXP unpack_names(const packed_names *p, const char *field) {
  SEXP out = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) p->n));
  for (int64_t i = 0; i < p->n; i++) {
    int64_t lo = p->offsets[i], hi = p->offsets[i + 1];
    SET_STRING_ELT(out, i,
                   Rf_mkCharLen((const char *) p->bytes + lo, (int) (hi - lo)));
  }
  UNPROTECT(1);
  (void) field;
  return out;
}

/* ------------------------------------------------------------------------ */
/* argument extraction                                                       */
/* ------------------------------------------------------------------------ */

static double scalar_real(SEXP x, const char *field) {
  if (Rf_length(x) != 1) {
    Rf_error("`%s` must be length 1, not %d", field, Rf_length(x));
  }
  if (TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) {
    Rf_error("`%s` must be a number, not %s", field, Rf_type2char(TYPEOF(x)));
  }
  return (TYPEOF(x) == REALSXP) ? REAL(x)[0] : (double) INTEGER(x)[0];
}

static int32_t scalar_flag(SEXP x, const char *field) {
  if (Rf_length(x) != 1) {
    Rf_error("`%s` must be length 1, not %d", field, Rf_length(x));
  }
  if (TYPEOF(x) != LGLSXP && TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP) {
    Rf_error("`%s` must be TRUE or FALSE, not %s", field, Rf_type2char(TYPEOF(x)));
  }
  int v = (TYPEOF(x) == LGLSXP) ? LOGICAL(x)[0]
          : (TYPEOF(x) == INTSXP) ? INTEGER(x)[0]
                                   : (int) REAL(x)[0];
  if (v == NA_LOGICAL) Rf_error("`%s` is NA", field);
  return v ? 1 : 0;
}

/* A numeric matrix's payload, with its type recorded so the bridge keeps the
 * integer-ness. `is_int` is set to 1 for INTSXP. `rows`/`cols` come from the dim
 * attribute; a vector is treated as one column, which is what the R side sends
 * for a single covariate. */
static void matrix_arg(SEXP x, const char *field,
                       const int32_t **int_out, const double **real_out,
                       int64_t *rows, int64_t *cols, int32_t *is_int) {
  if (TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) {
    Rf_error("`%s` must be a numeric matrix, not %s", field, Rf_type2char(TYPEOF(x)));
  }
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if (dim == R_NilValue) {
    *rows = Rf_xlength(x);
    *cols = 1;
  } else {
    if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
      Rf_error("`%s` must have a two-element dim attribute", field);
    }
    *rows = INTEGER(dim)[0];
    *cols = INTEGER(dim)[1];
  }
  if (*rows < 0 || *cols < 0) {
    Rf_error("`%s` has a negative dimension", field);
  }
  /* The out-parameters are optional: the design is only ever read as REALSXP and
   * its integer-ness is not sent on, so the caller passes NULL for one of them.
   * Writing through a NULL out-parameter segfaults at address nil, which is what
   * the first run did. */
  if (int_out)  *int_out  = (TYPEOF(x) == INTSXP) ? (const int32_t *) INTEGER(x) : NULL;
  if (real_out) *real_out = (TYPEOF(x) == REALSXP) ? (const double *) REAL(x) : NULL;
  if (is_int)   *is_int   = (TYPEOF(x) == INTSXP) ? 1 : 0;
}

/* ------------------------------------------------------------------------ */
/* C_ancombc2_rs_probe                                                       */
/* ------------------------------------------------------------------------ */

/* Arguments:
 *   0 data              numeric matrix, n_tax x n_samp
 *   1 aggregate_data    numeric matrix, same shape and type family
 *   2 design            numeric matrix, n_samp x p
 *   3 design_complete   raw, n_samp bytes, each 0 or 1
 *   4 group_index       integer, n_samp, 0-based-with-minus-one? no: 1-based, 0 = absent
 *   5 group_labels      character, may be character(0)
 *   6 fix_eff           character, length p
 *   7 taxon_names       character, length n_tax
 *   8 sample_names      character, length n_samp
 *   9 pseudo            number
 *  10 s0_perc           number
 *  11 prv_cut           number
 *  12 lib_cut           number
 *  13 alpha             number
 *  14 iter_tol          number
 *  15 iter_max          integer
 *  16 em_tol            number
 *  17 em_max            integer
 *  18 mdfdr_b           integer
 *  19 global            flag
 *  20 pairwise          flag
 *  21 dunnet            flag
 *  22 trend             flag
 *  23 pseudo_sens       flag
 *  24 conservative      flag
 *  25 struc_zero        flag
 *  26 neg_lb            flag
 *  27 which             integer, the probe selector
 *  28 out_n             integer, capacity of the out buffer
 *  29 out_i64           integer(1), receives the true length
 */
SEXP C_ancombc2_rb_probe(SEXP args) {
  /* One PROTECT stack for the whole function, balanced at each exit. The counter
   * is explicit rather than relying on the caller, because `Rf_error` longjmps and
   * an unbalanced stack would leak every protection on the way out. */
  int nprot = 0;
  #define PROTECT_AT(x) do { PROTECT(x); nprot++; } while (0)
  #define UNPROTECT_TO(n) do { while (nprot > (n)) { UNPROTECT(1); nprot--; } } while (0)

  const int nargs = 29;
  if (TYPEOF(args) != VECSXP) {
    Rf_error("internal: the probe argument list must be a VECSXP");
  }
  if (Rf_length(args) != nargs) {
    Rf_error("internal: the probe takes %d arguments, got %d", nargs, Rf_length(args));
  }

  SEXP data   = VECTOR_ELT(args, 0);
  SEXP agg    = VECTOR_ELT(args, 1);
  SEXP design = VECTOR_ELT(args, 2);
  SEXP dcomp  = VECTOR_ELT(args, 3);
  SEXP gidx   = VECTOR_ELT(args, 4);
  SEXP glab   = VECTOR_ELT(args, 5);
  SEXP fixeff = VECTOR_ELT(args, 6);
  SEXP taxnm  = VECTOR_ELT(args, 7);
  SEXP sampnm = VECTOR_ELT(args, 8);

  const int32_t *data_i = NULL, *agg_i = NULL;
  const double  *data_r = NULL, *agg_r = NULL;
  int64_t n_tax = 0, n_samp = 0, a_tax = 0, a_samp = 0, p_rows = 0, p_cols = 0;
  int32_t data_is_int = 0, agg_is_int = 0, design_is_int = 0;
  const double *design_raw = NULL;

  matrix_arg(data, "data", &data_i, &data_r, &n_tax, &n_samp, &data_is_int);
  matrix_arg(agg, "aggregate_data", &agg_i, &agg_r, &a_tax, &a_samp, &agg_is_int);
  /* The design is read only as REALSXP below, so its integer payload is not
   * requested; `design_is_int` still comes back so the widening can happen. */
  matrix_arg(design, "design", NULL, &design_raw, &p_rows, &p_cols, &design_is_int);

  if (TYPEOF(dcomp) != RAWSXP) {
    Rf_error("`design_complete` must be a raw vector, not %s", Rf_type2char(TYPEOF(dcomp)));
  }
  if (TYPEOF(gidx) != INTSXP) {
    Rf_error("`group_index` must be an integer vector, not %s", Rf_type2char(TYPEOF(gidx)));
  }

  /* An INTSXP design is widened to double *here*, on the R main thread, and the
   * widened copy is what crosses. `.ancombc2_core` calls `model.matrix`, which is
   * always double, so widening at the boundary is what the original would have
   * done anyway -- and doing it in R keeps `REAL()` out of Rust. */
  const double *design_r = NULL;
  if (design_is_int) {
    R_xlen_t n = p_rows * p_cols;
    SEXP dbl = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t k = 0; k < n; k++) REAL(dbl)[k] = (double) INTEGER(design)[k];
    nprot++;
    design_r = REAL(dbl);
  } else {
    design_r = design_raw;
  }

  packed_names pn_glab, pn_fix, pn_tax, pn_samp;
  pack_names(&pn_glab, glab, "group_labels");
  pack_names(&pn_fix, fixeff, "fix_eff");
  pack_names(&pn_tax, taxnm, "taxon_names");
  pack_names(&pn_samp, sampnm, "sample_names");

  double pseudo   = scalar_real(VECTOR_ELT(args,  9), "pseudo");
  double s0_perc  = scalar_real(VECTOR_ELT(args, 10), "s0_perc");
  double prv_cut  = scalar_real(VECTOR_ELT(args, 11), "prv_cut");
  double lib_cut  = scalar_real(VECTOR_ELT(args, 12), "lib_cut");
  double alpha    = scalar_real(VECTOR_ELT(args, 13), "alpha");
  double iter_tol = scalar_real(VECTOR_ELT(args, 14), "iter_control$tol");
  int32_t iter_max = (int32_t) scalar_real(VECTOR_ELT(args, 15), "iter_control$max_iter");
  double em_tol   = scalar_real(VECTOR_ELT(args, 16), "em_control$tol");
  int32_t em_max  = (int32_t) scalar_real(VECTOR_ELT(args, 17), "em_control$max_iter");
  int32_t mdfdr_b = (int32_t) scalar_real(VECTOR_ELT(args, 18), "mdfdr_control$B");
  int32_t f_global    = scalar_flag(VECTOR_ELT(args, 19), "global");
  int32_t f_pairwise  = scalar_flag(VECTOR_ELT(args, 20), "pairwise");
  int32_t f_dunnet    = scalar_flag(VECTOR_ELT(args, 21), "dunnet");
  int32_t f_trend     = scalar_flag(VECTOR_ELT(args, 22), "trend");
  int32_t f_pseudo_sens   = scalar_flag(VECTOR_ELT(args, 23), "pseudo_sens");
  int32_t f_conservative = scalar_flag(VECTOR_ELT(args, 24), "conservative");
  int32_t f_struc_zero   = scalar_flag(VECTOR_ELT(args, 25), "struc_zero");
  int32_t f_neg_lb       = scalar_flag(VECTOR_ELT(args, 26), "neg_lb");
  int32_t which      = (int32_t) scalar_real(VECTOR_ELT(args, 27), "which");
  int64_t  out_cap   = (int64_t) scalar_real(VECTOR_ELT(args, 28), "out_n");

  /* The bridge's request struct. Declared here and filled field by field so that
   * a field-order mismatch against `RawRequest` in
   * crates/ancombc2-rbridge/src/abi.rs is visible in this file rather than
   * producing a plausible wrong answer. */
  struct rb_request {
    const int32_t *data_int;
    const double  *data_real;
    int64_t data_rows, data_cols;
    int32_t data_is_int;
    const int32_t *aggregate_int;
    const double  *aggregate_real;
    int64_t aggregate_rows, aggregate_cols;
    int32_t aggregate_is_int;
    const double  *design;
    int64_t design_rows, design_cols;
    const uint8_t *design_complete;
    const int32_t *group_index;
    int64_t n_group_levels;
    const int64_t *group_labels_offsets; const uint8_t *group_labels_bytes; int64_t n_group_labels;
    const int64_t *fix_eff_offsets;      const uint8_t *fix_eff_bytes;      int64_t n_fix_eff;
    const int64_t *taxon_names_offsets;  const uint8_t *taxon_names_bytes;  int64_t n_taxon_names;
    const int64_t *sample_names_offsets; const uint8_t *sample_names_bytes; int64_t n_sample_names;
    double pseudo, s0_perc, prv_cut, lib_cut, alpha, iter_tol;
    int32_t iter_max;
    double em_tol;
    int32_t em_max, mdfdr_b;
    int32_t global, pairwise, dunnet, trend, pseudo_sens, conservative, struc_zero, neg_lb;
  } req;

  memset(&req, 0, sizeof(req));
  req.data_int = data_i; req.data_real = data_r;
  req.data_rows = n_tax; req.data_cols = n_samp; req.data_is_int = data_is_int;
  req.aggregate_int = agg_i; req.aggregate_real = agg_r;
  req.aggregate_rows = a_tax; req.aggregate_cols = a_samp; req.aggregate_is_int = agg_is_int;
  req.design = design_r; req.design_rows = p_rows; req.design_cols = p_cols;
  req.design_complete = (const uint8_t *) RAW(dcomp);
  req.group_index = (const int32_t *) INTEGER(gidx);
  req.n_group_levels = pn_glab.n;
  req.group_labels_offsets = pn_glab.offsets; req.group_labels_bytes = pn_glab.bytes; req.n_group_labels = pn_glab.n;
  req.fix_eff_offsets = pn_fix.offsets; req.fix_eff_bytes = pn_fix.bytes; req.n_fix_eff = pn_fix.n;
  req.taxon_names_offsets = pn_tax.offsets; req.taxon_names_bytes = pn_tax.bytes; req.n_taxon_names = pn_tax.n;
  req.sample_names_offsets = pn_samp.offsets; req.sample_names_bytes = pn_samp.bytes; req.n_sample_names = pn_samp.n;
  req.pseudo = pseudo; req.s0_perc = s0_perc; req.prv_cut = prv_cut; req.lib_cut = lib_cut;
  req.alpha = alpha; req.iter_tol = iter_tol; req.iter_max = iter_max;
  req.em_tol = em_tol; req.em_max = em_max; req.mdfdr_b = mdfdr_b;
  req.global = f_global; req.pairwise = f_pairwise; req.dunnet = f_dunnet; req.trend = f_trend;
  req.pseudo_sens = f_pseudo_sens; req.conservative = f_conservative;
  req.struc_zero = f_struc_zero; req.neg_lb = f_neg_lb;

  /* Scratch for the bridge to write into. `R_alloc` is R's allocator, so it is
   * freed by R's own cleanup -- there is nothing for Rust to release and nothing
   * for this file to free. */
  int64_t want = (which == ECHO_DATA_INT || which == ECHO_DATA_REAL ||
                  which == ECHO_AGGREGATE)          ? n_tax * n_samp
               : (which == ECHO_DESIGN)              ? p_rows * p_cols
               : (which == ECHO_DESIGN_COMPLETE)      ? n_samp
               : (which == ECHO_GROUP_INDEX)          ? n_samp
               : 0;
  int64_t cap = (out_cap > 0) ? out_cap : want;
  if (cap < want) {
    UNPROTECT_TO(0);
    Rf_error("the probe buffer holds %lld element(s) but %lld were requested",
             (long long) cap, (long long) want);
  }

  int32_t *scratch_i32 = (int32_t *) R_alloc((size_t) (cap > 0 ? cap : 1), sizeof(int32_t));
  double  *scratch_f64 = (double *)  R_alloc((size_t) (cap > 0 ? cap : 1), sizeof(double));
  uint8_t *scratch_u8  = (uint8_t *) R_alloc((size_t) (cap > 0 ? cap : 1), sizeof(uint8_t));
  int64_t *out_len     = (int64_t *) R_alloc(1, sizeof(int64_t));
  int32_t  out_kind    = KIND_NONE;

  /* Every scratch buffer is passed for every probe. They are separate
   * allocations, so a probe that writes the wrong one is visible as a wrong
   * payload rather than as an aliased one; passing NULL here instead would make a
   * wrong-kind probe a null dereference. */
  int32_t status = ancombc2_rb_probe(
      &req, which,
      scratch_i32, scratch_f64, scratch_u8,
      out_len, cap, &out_kind);

  if (status != RB_OK) {
    UNPROTECT_TO(0);
    rb_error(status, "probe");
  }

  int64_t true_len = out_len[0];

  SEXP payload;
  switch (out_kind) {
    case KIND_INT:
      payload = PROTECT(Rf_allocVector(INTSXP, true_len)); nprot++;
      if (true_len > 0) memcpy(INTEGER(payload), scratch_i32, (size_t) true_len * sizeof(int32_t));
      break;
    case KIND_RAW:
      payload = PROTECT(Rf_allocVector(RAWSXP, true_len)); nprot++;
      if (true_len > 0) memcpy(RAW(payload), scratch_u8, (size_t) true_len);
      break;
    case KIND_REAL:
      payload = PROTECT(Rf_allocVector(REALSXP, true_len)); nprot++;
      if (true_len > 0) memcpy(REAL(payload), scratch_f64, (size_t) true_len * sizeof(double));
      break;
    default:
      UNPROTECT_TO(0);
      Rf_error("ancombc2-rbridge: probe %d reported no payload kind", which);
  }

  /* Names come back rebuilt from the same packed buffers the request used, so a
   * name that survived the request is proved to have survived intact -- through
   * the pack, through Rust's monotone-partition check, and back out. */
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
  SET_VECTOR_ELT(out, 0, payload);
  SET_VECTOR_ELT(out, 1, unpack_names(&pn_glab, "group_labels"));
  SET_VECTOR_ELT(out, 2, unpack_names(&pn_fix,  "fix_eff"));
  SET_VECTOR_ELT(out, 3, unpack_names(&pn_tax,  "taxon_names"));
  SET_VECTOR_ELT(out, 4, unpack_names(&pn_samp, "sample_names"));
  SEXP dims = PROTECT(Rf_allocVector(INTSXP, 2)); nprot++;
  INTEGER(dims)[0] = (int) n_tax;
  INTEGER(dims)[1] = (int) n_samp;
  SET_VECTOR_ELT(out, 5, dims);

  UNPROTECT_TO(0);
  #undef PROTECT_AT
  #undef UNPROTECT_TO
  return out;
}

/* ------------------------------------------------------------------------ */
/* ancombc2_rb_reducer                                                       */
/* ------------------------------------------------------------------------ */

/* R's `rowSums`/`colSums`/`rowMeans`/`colMeans` with `na.rm = TRUE`, called across
 * the boundary rather than re-derived in Rust.
 *
 * R accumulates these in C `long double`, which is 64 bits of mantissa on x86-64,
 * binary128 on aarch64, and `double` on Windows. A Rust re-derivation would pin
 * one platform's answer as the oracle's; asking `base::rowMeans` is exact by
 * construction on whichever machine R is on, and it is literally the function the
 * oracle called. The measurement that forced this is in
 * `crates/ancombc2-core/src/reduce.rs`: for `c(1e16, 1, -1e16, 1e-17, 1)`,
 * `rowMeans` is 0x1.999999999999ap-2 where an `f64` sum gives 0x1.999999999999ap-3.
 *
 * `base::rowMeans` is a closure wrapping `.Internal(rowMeans(x, m, n, na.rm))`,
 * and `R_xlen_t`'s `long double` accumulator is inside that internal. Evaluating
 * the closure is therefore the way to reach it from C: the internal itself is not
 * exported.
 */

/* Build and evaluate `fn_(x, na.rm = TRUE)` for a base closure.
 *
 * The call is constructed rather than parsed because a cached LANGSXP would be
 * mutated by evaluation, and two reductions in flight must not share one. It is
 * cheap: three CONS cells, and the function cell is interned.
 */
/* Call a reduction wrapper with `x`: `fn(x)`, evaluated.
 *
 * `fn` is handed in from R rather than looked up here. Three reasons, in order of
 * how much they cost:
 *
 *  * `Rf_findNamespace` is not in `Rinternals.h` -- only the internal
 *    `R_FindNamespace` is -- so resolving a name from C means depending on a symbol
 *    the public API does not promise. The first attempt at this did that and
 *    segfaulted on the first call.
 *  * A wrapper with one formal makes this a positional call. Writing
 *    `rowMeans(x, na.rm = TRUE)` by hand needs a named argument, and two distinct
 *    mistakes both surface as "could not find function na.rm" rather than as a name
 *    mismatch.
 *  * It costs nothing. The closures are created once per R session; this is four
 *    calls per run.
 */
static SEXP eval_na_rm_reduction(SEXP fn, const char *fn_name, SEXP x) {
    int nprot = 0;
    if (TYPEOF(fn) != CLOSXP) {
      Rf_error("the `%s` reduction was not passed as a closure, but as %s", fn_name,
               Rf_type2char(TYPEOF(fn)));
    }
    SEXP call = PROTECT(Rf_lang2(fn, x)); nprot++;
    SEXP out = PROTECT(Rf_eval(call, R_BaseEnv)); nprot++;

    if (TYPEOF(out) != REALSXP) {
      UNPROTECT(nprot);
      Rf_error("`%s` returned %s, not a double vector", fn_name,
               Rf_type2char(TYPEOF(out)));
    }
    UNPROTECT(nprot);
    return out;
}

/* Wrap a row-major `rows x cols` buffer as an R matrix and hand it to `fn_`.
 *
 * The bridge's buffers are row-major and R's are column-major, so the copy here
 * is the transpose. It is not optional: passing the buffer straight to R would
 * read element `(i, j)` as `(j, i)`, which for a symmetric matrix is invisible and
 * for anything else is wrong.
 */
static SEXP reduce_rows(SEXP fn, const char *fn_name, const double *x, int64_t rows,
                        int64_t cols, int64_t want_len) {
  int nprot = 0;
  SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int) rows, (int) cols)); nprot++;
  for (int64_t j = 0; j < cols; j++) {
    for (int64_t i = 0; i < rows; i++) {
      REAL(m)[i + j * rows] = x[i * cols + j];
    }
  }
  SEXP out = PROTECT(eval_na_rm_reduction(fn, fn_name, m)); nprot++;
  int64_t got = (int64_t) Rf_xlength(out);
  if (got != want_len) {
    UNPROTECT(nprot);
    Rf_error("`%s` on a %lld x %lld matrix returned %lld values, expected %lld",
             fn_name, (long long) rows, (long long) cols,
             (long long) got, (long long) want_len);
  }
  /* Drop both `m` and `out`, then re-protect `out` alone.
   *
   * The caller reads `REAL(out)` after this returns and unprotects it there, so it
   * must be protected on return -- returning it unprotected is how a `.Call` hands
   * back a collectable object, and the next allocation frees it. And `m` cannot stay
   * on the stack: leaving it there is a PROTECT imbalance of one per reduction, which
   * R reports as "stack imbalance in '<-', 2 then 4" rather than as anything local.
   * Popping both and re-protecting the one that escapes is the balanced way to say
   * "only this leaves".
   */
  UNPROTECT(nprot);
  return PROTECT(out);
}

/* `x` is `rows x cols` row-major. `out` receives `rows` values for a row-wise
 * reduction and `cols` for a column-wise one -- stated per reduction rather than
 * inferred from `min(rows, cols)`, which is right only while one axis happens to be
 * the smaller one. */
#define DEFINE_ROW_REDUCTION(cname, rname, slot)                                  \
  static int32_t cname(const double *x, int64_t rows, int64_t cols, double *out) { \
    if (x == NULL || out == NULL) return RB_ERR_ARGS;                             \
    if (rows <= 0 || cols <= 0) return RB_ERR_ARGS;                              \
    SEXP r = reduce_rows(VECTOR_ELT(rb_reduce_args, slot), rname, x, rows, cols,   \
                         rows);                                                    \
    memcpy(out, REAL(r), (size_t) rows * sizeof(double));                         \
    UNPROTECT(1);                                                                 \
    return RB_OK;                                                                 \
  }

#define DEFINE_COL_REDUCTION(cname, rname, slot)                                  \
  static int32_t cname(const double *x, int64_t rows, int64_t cols, double *out) { \
    if (x == NULL || out == NULL) return RB_ERR_ARGS;                             \
    if (rows <= 0 || cols <= 0) return RB_ERR_ARGS;                              \
    SEXP r = reduce_rows(VECTOR_ELT(rb_reduce_args, slot), rname, x, rows, cols,   \
                         cols);                                                    \
    memcpy(out, REAL(r), (size_t) cols * sizeof(double));                         \
    UNPROTECT(1);                                                                 \
    return RB_OK;                                                                 \
  }

/* The R caller parks the four wrappers here for the duration of one probe call.
 * It is a plain `static SEXP`, set and cleared around the call, because the Rust side
 * reaches these functions through a function-pointer table with a fixed signature
 * and this is how the closures reach it without a lookup. */
static SEXP rb_reduce_args = NULL;

/* Argument indices of the four wrappers, in the order of `struct rb_reducer`, which
 * the bridge fills positionally. These are indices into the *caller's* argument
 * vector, where they sit after the eleven preprocessing arguments -- so they are not
 * the same numbers as the struct's field order, and the `+ 11` is not optional. */
#define RB_SLOT_ROW_MEANS 11
#define RB_SLOT_COL_MEANS 12
#define RB_SLOT_COL_SUMS  13
#define RB_SLOT_ROW_SUMS  14

DEFINE_ROW_REDUCTION(rb_row_means, "rb_row_means_na_rm", RB_SLOT_ROW_MEANS)
DEFINE_COL_REDUCTION(rb_col_means, "rb_col_means_na_rm", RB_SLOT_COL_MEANS)
DEFINE_COL_REDUCTION(rb_col_sums, "rb_col_sums_na_rm", RB_SLOT_COL_SUMS)
DEFINE_ROW_REDUCTION(rb_row_sums, "rb_row_sums_na_rm", RB_SLOT_ROW_SUMS)

/* Field order must match `Reducer` in crates/ancombc2-rbridge/src/output.rs. */
struct rb_reducer {
  int32_t (*row_means)(const double *, int64_t, int64_t, double *);
  int32_t (*col_means)(const double *, int64_t, int64_t, double *);
  int32_t (*col_sums)(const double *, int64_t, int64_t, double *);
  int32_t (*row_sums)(const double *, int64_t, int64_t, double *);
};

int32_t ancombc2_rb_reducer(struct rb_reducer *out) {
  if (out == NULL) return RB_ERR_ARGS;
  out->row_means = rb_row_means;
  out->col_means = rb_col_means;
  out->col_sums = rb_col_sums;
  out->row_sums = rb_row_sums;
  return RB_OK;
}

/* ------------------------------------------------------------------------ */
/* C_ancombc2_rb_fit_probe                                                   */
/* ------------------------------------------------------------------------ */

struct rb_fit {
  const double *x;      /* n_samp x p, column-major */
  const double *y;      /* n_taxa x n_samp, row-major */
  int64_t n_samp;
  int64_t p;
  int64_t n_taxa;
  const uint8_t *observed;
};

extern int32_t ancombc2_rb_fit_probe(const struct rb_fit *in, double *beta,
                                     double *fitted, double *dof, int *groups);

/* Arguments:
 *   0 x  1 y  2 observed  3 n_samp  4 p  5 n_taxa
 *
 * `y` arrives as an `n_taxa x n_samp` matrix -- R's column-major -- and the bridge
 * wants it row-major, because `.lm_fit_all` indexes `Ymat[i, ]` per taxon and the
 * core's buffers are taxon-major. `observed` arrives as a logical matrix of the same
 * shape and is flattened the same way.
 *
 * Returns a list of `beta`, `fitted`, `dof`, in `.lm_fit_all`'s own shapes.
 */
SEXP C_ancombc2_rb_fit_probe(SEXP args) {
  int nprot = 0;
  #define UNPROTECT_TO(n) do { while (nprot > (n)) { UNPROTECT(1); nprot--; } } while (0)

  const int nargs = 6;
  if (TYPEOF(args) != VECSXP) Rf_error("internal: fit probe takes a VECSXP");
  if (Rf_length(args) != nargs) {
    Rf_error("internal: fit probe takes %d arguments, got %d", nargs, Rf_length(args));
  }

  SEXP x = VECTOR_ELT(args, 0);
  SEXP y = VECTOR_ELT(args, 1);
  SEXP obs = VECTOR_ELT(args, 2);
  int64_t n_samp = (int64_t) scalar_real(VECTOR_ELT(args, 3), "n_samp");
  int64_t p = (int64_t) scalar_real(VECTOR_ELT(args, 4), "p");
  int64_t n_taxa = (int64_t) scalar_real(VECTOR_ELT(args, 5), "n_taxa");
  if (n_samp <= 0 || p <= 0 || n_taxa <= 0) {
    UNPROTECT_TO(0);
    Rf_error("n_samp, p and n_taxa must all be positive, got %lld, %lld, %lld",
             (long long) n_samp, (long long) p, (long long) n_taxa);
  }

  /* `x` is widened to `REALSXP` here rather than pointed at in place, because
   * `matrix_arg` hands back `NULL` for a `real_out` when the matrix is
   * `INTSXP` -- and then `dp` is NULL and every design value reads as garbage.
   * `rb_fit_probe` coerces in R, so in practice `x` is already `REALSXP` and this
   * is a no-op copy; the guard is here because a NULL design pointer is a silent
   * wrong answer rather than a crash, and because the probe should not depend on
   * which layer happened to widen it. */
  SEXP xmat = PROTECT(Rf_coerceVector(x, REALSXP)); nprot++;
  if (TYPEOF(xmat) != REALSXP) {
    UNPROTECT_TO(nprot);
    Rf_error("`x` must be a numeric matrix, not %s", Rf_type2char(TYPEOF(x)));
  }
  const double *dp = NULL, *yp = NULL;
  int64_t xr = 0, xc = 0;
  matrix_arg(xmat, "x", NULL, &dp, &xr, &xc, NULL);
  if (dp == NULL) {
    UNPROTECT_TO(nprot);
    Rf_error("`x` could not be read as doubles");
  }
  if (xr != n_samp || xc != p) {
    UNPROTECT_TO(0);
    Rf_error("`x` is %lld x %lld, not %lld x %lld", (long long) xr, (long long) xc,
             (long long) n_samp, (long long) p);
  }
  /* `as.vector(t(y))` is the row-major flattening, because `t()` then `as.vector()`
   * reads the transposed matrix column-major, i.e. the original rows in order. */
  {
    /* `as.matrix` at the C level: an integer matrix is widened, which is what R's
     * `as.matrix` does to an integer one, and a data.frame is refused rather than
     * guessed at. */
    SEXP ymat = PROTECT(Rf_coerceVector(y, REALSXP)); nprot++;
    if (TYPEOF(ymat) != REALSXP) {
      UNPROTECT_TO(nprot);
      Rf_error("`y` must be a numeric matrix, not %s", Rf_type2char(TYPEOF(y)));
    }
    SEXP dim = Rf_getAttrib(ymat, R_DimSymbol);
    if (dim == R_NilValue || TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
      UNPROTECT_TO(nprot);
      Rf_error("`y` must be a matrix");
    }
    if ((int64_t) INTEGER(dim)[0] != n_taxa || (int64_t) INTEGER(dim)[1] != n_samp) {
      UNPROTECT_TO(nprot);
      Rf_error("`y` is %d x %d, not %lld x %lld", INTEGER(dim)[0], INTEGER(dim)[1],
               (long long) n_taxa, (long long) n_samp);
    }
    SEXP flat = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (n_taxa * n_samp))); nprot++;
    for (int64_t t = 0; t < n_taxa; t++) {
      for (int64_t j = 0; j < n_samp; j++) {
        /* R's column-major index of `(t, j)` is `t + j * n_taxa` -- the stride is
         * the number of *rows*. Writing `j + t * n_samp` instead leaves the buffer
         * in R's own column-major order, which is not a permutation that looks
         * wrong: it hands the core taxon 1's response interleaved with taxon 2's,
         * 3's and so on, so every taxon gets a plausible-looking response and every
         * fit is quietly the wrong one. The sentinel `y[t, j] = 100 t + j` is what
         * made it visible -- the core reported a first-half mean of 351, which is
         * `mean(101, 201, 301, 401, 501, 601)`. */
        REAL(flat)[t * n_samp + j] = REAL(ymat)[t + j * n_taxa];
      }
    }
    yp = REAL(flat);
  }
  {
    SEXP omat = PROTECT(Rf_coerceVector(obs, LGLSXP)); nprot++;
    SEXP odim = Rf_getAttrib(omat, R_DimSymbol);
    if (odim == R_NilValue || TYPEOF(odim) != INTSXP || Rf_length(odim) != 2 ||
        (int64_t) INTEGER(odim)[0] != n_taxa || (int64_t) INTEGER(odim)[1] != n_samp) {
      UNPROTECT_TO(nprot);
      Rf_error("`observed` must be an %lld x %lld logical matrix", (long long) n_taxa,
               (long long) n_samp);
    }
    /* R's `LGLSXP` is an `int32_t` array -- four bytes per cell, *not* one, despite
     * the name and despite `TRUE`/`FALSE` being 0/1. Casting `LOGICAL()` to
     * `uint8_t *` and handing it over as the core's byte-per-cell mask therefore
     * reads only the first *quarter* of the buffer, in a stride of one: an all-`TRUE`
     * 72-cell mask arrived with 18 cells set, and the four bytes of the first cell
     * (`01 00 00 00`) read as `TRUE`, `FALSE`, `FALSE`, `FALSE`. It was not a
     * transposition and not a stride mistake in the loop; it was the element size.
     *
     * So the mask is narrowed here, into a real byte buffer, rather than cast. */
    uint8_t *oflat = (uint8_t *) R_alloc((size_t) (n_taxa * n_samp), sizeof(uint8_t));
    if (oflat == NULL) {
      UNPROTECT_TO(nprot);
      Rf_error("could not allocate the %lld-cell observed mask",
               (long long) (n_taxa * n_samp));
    }
    for (int64_t t = 0; t < n_taxa; t++) {
      for (int64_t j = 0; j < n_samp; j++) {
        /* The same column-major index as `y` above: `(t, j)` is at `t + j * n_taxa`.
         * Getting this wrong is invisible on an all-`TRUE` mask and only shows up
         * once a cell is missing, which is exactly when it matters. */
        /* R's `TRUE` is 1 and `FALSE` is 0, but `NA` is `INT_MIN`, so this
         * normalises rather than casts: `INT_MIN` would set every bit. */
        oflat[t * n_samp + j] = (LOGICAL(omat)[t + j * n_taxa] == TRUE) ? 1u : 0u;
      }
    }
    struct rb_fit in;
    memset(&in, 0, sizeof(in));
    in.x = dp; in.y = yp;
    in.n_samp = n_samp; in.p = p; in.n_taxa = n_taxa;
    in.observed = oflat;

    SEXP beta = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (n_taxa * p))); nprot++;
    SEXP fitted = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (n_taxa * n_samp))); nprot++;
    SEXP dof = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) n_taxa)); nprot++;
    /* 3 + 1 + 2 per group: n_samp, p, n_taxa, the group count, then (rows, taxa). */
    SEXP groups = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (7 + 2 * n_taxa))); nprot++;
    int32_t st = ancombc2_rb_fit_probe(&in, REAL(beta), REAL(fitted), REAL(dof),
                                       INTEGER(groups));
    if (st != RB_OK) { UNPROTECT_TO(0); rb_error(st, "fit probe"); }

    /* `beta` and `fitted` come back **row-major** (taxon-major), because that is the
     * core's own layout and `lm_fit_all` writes a taxon's row at a time. R reads a
     * matrix **column-major**, so giving these a `dim` without transposing first
     * makes R read down the taxon axis instead of along the sample axis. That is
     * not a rounding difference: R's `fitted[1, ]` was the interleaved responses of
     * taxa 1, 3, 5, ... , and it looked like a plausible fit of the wrong thing.
     * The fix is to transpose into R's order here rather than to make the core lay
     * its output out for R, which would put a presentation concern in the numerics. */
    /* R's column-major index for element `(t, j)` of an `n_taxa x n_samp` matrix is
     * `j * n_taxa + t`: the stride is the number of *rows*, which is the detail that
     * makes this easy to get wrong in a way that still looks like data. Writing
     * `j + t * n_samp` instead -- the core's row-major stride -- produces a buffer
     * whose every row is a plausible-looking mixture of the other rows. */
    SEXP b_tr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (n_taxa * p))); nprot++;
    for (int64_t t = 0; t < n_taxa; t++) {
      for (int64_t a = 0; a < p; a++) {
        REAL(b_tr)[a * n_taxa + t] = REAL(beta)[t * p + a];
      }
    }
    SEXP f_tr = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (n_taxa * n_samp))); nprot++;
    for (int64_t t = 0; t < n_taxa; t++) {
      for (int64_t j = 0; j < n_samp; j++) {
        REAL(f_tr)[j * n_taxa + t] = REAL(fitted)[t * n_samp + j];
      }
    }
    SEXP bdim = PROTECT(Rf_allocVector(INTSXP, 2)); nprot++;
    INTEGER(bdim)[0] = (int) n_taxa; INTEGER(bdim)[1] = (int) p;
    Rf_setAttrib(b_tr, R_DimSymbol, bdim);
    SEXP fdim = PROTECT(Rf_allocVector(INTSXP, 2)); nprot++;
    INTEGER(fdim)[0] = (int) n_taxa; INTEGER(fdim)[1] = (int) n_samp;
    Rf_setAttrib(f_tr, R_DimSymbol, fdim);

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_VECTOR_ELT(out, 0, b_tr); SET_STRING_ELT(names_sexp, 0, Rf_mkChar("beta"));
    SET_VECTOR_ELT(out, 1, f_tr); SET_STRING_ELT(names_sexp, 1, Rf_mkChar("fitted"));
    SET_VECTOR_ELT(out, 2, dof); SET_STRING_ELT(names_sexp, 2, Rf_mkChar("dof"));
    SET_VECTOR_ELT(out, 3, groups); SET_STRING_ELT(names_sexp, 3, Rf_mkChar("groups"));
    Rf_setAttrib(out, R_NamesSymbol, names_sexp);
    UNPROTECT_TO(0);
    #undef UNPROTECT_TO
    return out;
  }
}

/* ------------------------------------------------------------------------ */
/* C_ancombc2_rb_preprocess_probe                                            */
/* ------------------------------------------------------------------------ */

/* The slot the bridge fills per stage. Must match `RawStage` in
 * crates/ancombc2-rbridge/src/abi.rs: two i64 offsets and three i64s, 40 bytes.
 *
 * Offsets, not pointers. The first version of this held `const double *` and
 * `const int64_t *` pointing into the bridge's own `Vec`s, and the shim read them
 * after `ancombc2_rb_preprocess_probe` returned -- by which time Rust had dropped the
 * locals they pointed at. Nothing crashed and every shape was right; the first two
 * elements of each stage were simply whatever the allocator had left there. The
 * bridge now copies the values into buffers allocated here. */
struct rb_stage {
  int64_t f64_at;
  int64_t i32_at;
  int64_t n_rows;
  int64_t n_cols;
  int64_t n;
};

/* The stage names, in the bridge's order. `abi::preprocess_probe_layout_tests`
 * pins that order; this list is the same specification on the R side, and the
 * length check below makes a disagreement an error rather than a shift. */
static const char *preprocess_stage_names[] = {
  "prevalence1", "tax_keep1", "lib_size1", "samp_keep1",
  "O1", "log1", "means1", "y1",
  "tax_keep2", "tax_keep2_absolute", "samp_keep2",
  "O2", "log2", "means2", "y2",
  "zero_keep",
  "group_prevalence", "group_size", "group_lower", "zero_ind"
};
#define N_PREPROCESS_STAGES \
  ((int64_t) (sizeof(preprocess_stage_names) / sizeof(preprocess_stage_names[0])))

struct rb_preprocess {
  const double *data;
  const double *aggregate;
  const double *group_indicator;
  int64_t n_tax;
  int64_t n_samp;
  int64_t n_groups;
  int32_t struc_zero;
  double prv_cut;
  double lib_cut;
  double pseudo;
  int32_t neg_lb;
};

/* Turn one slot into an R object: a matrix when the slot has two dimensions, an
 * integer vector when the payload is in the `i32` buffer, a double vector otherwise.
 *
 * A zero-length stage is `integer(0)` when it is an index stage and `numeric(0)`
 * otherwise, because an empty `which()` is `integer(0)` in the reference and `NULL` is
 * not the same value. An uncomputed stage is `integer(0)` with a class, so "the screen
 * did not run" is distinguishable from "the screen kept nothing".
 */
static SEXP stage_to_sexp(const struct rb_stage *s, const double *f64_buf,
                          const int32_t *i32_buf) {
  int nprot = 0;
  SEXP v;
  int is_int = s->i32_at >= 0;
  int computed = s->f64_at >= 0 || is_int;
  if (s->n == 0) {
    v = PROTECT(Rf_allocVector(is_int ? INTSXP : REALSXP, 0)); nprot++;
    if (!computed) Rf_setAttrib(v, R_ClassSymbol, Rf_mkString("uncomputed-stage"));
  } else if (is_int) {
    v = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) s->n)); nprot++;
    for (int64_t k = 0; k < s->n; k++) INTEGER(v)[k] = (int) i32_buf[s->i32_at + k];
  } else {
    v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) s->n)); nprot++;
    memcpy(REAL(v), f64_buf + s->f64_at, (size_t) s->n * sizeof(double));
  }
  /* `n_cols > 0`, not `> 1`: the reference's `feature_table[tax_keep, , drop = FALSE]`
   * is `n x 1` and keeps its `dim`, and a stage that declares a shape must get one.
   * `n_cols == 0` is what marks a vector, and that is the only thing the `0` means. */
  if (s->n_rows > 0 && s->n_cols > 0) {
    if (s->n_rows * s->n_cols != s->n) {
      UNPROTECT(nprot);
      Rf_error("a stage claims %lld x %lld but carries %lld values",
               (long long) s->n_rows, (long long) s->n_cols, (long long) s->n);
    }
    SEXP d = PROTECT(Rf_allocVector(INTSXP, 2)); nprot++;
    INTEGER(d)[0] = (int) s->n_rows;
    INTEGER(d)[1] = (int) s->n_cols;
    Rf_setAttrib(v, R_DimSymbol, d);
    UNPROTECT(1); nprot--;
  }
  /* `v` is returned *still protected*, and the caller unprotects it once it has
   * stored `v` in a protected VECSXP. Unprotecting here instead is a PROTECT
   * imbalance of one per stage -- nineteen of them -- which R reports as "stack
   * imbalance in '<-', 2 then 21" and which points at nothing in particular. */
  return v;
}

/* Arguments, one VECSXP:
 *   0 feature_table   1 feature_table_aggregate   2 group_indicator
 *   3 n_tax  4 n_samp  5 n_groups
 *   6 struc_zero  7 prv_cut  8 lib_cut  9 pseudo  10 neg_lb
 *  11 rb_row_means_na_rm  12 rb_col_means_na_rm  13 rb_col_sums_na_rm
 *  14 rb_row_sums_na_rm
 *
 * Returns a named list of the nineteen stages, plus `n_tax`, `n_samp`, `n_groups`
 * and `has_group` so the R side can rebuild the shapes it expects.
 *
 * This is a diagnostic. It computes the same arithmetic `ancombc2()` will, with the
 * same R-backed reduction, and hands each stage back separately so a difference can
 * be attributed to a stage rather than to the result.
 */
SEXP C_ancombc2_rb_preprocess_probe(SEXP args) {
  int nprot = 0;
  #define UNPROTECT_TO(n) do { while (nprot > (n)) { UNPROTECT(1); nprot--; } } while (0)

  const int nargs = 15;
  if (TYPEOF(args) != VECSXP) Rf_error("internal: preprocess probe takes a VECSXP");
  if (Rf_length(args) != nargs) {
    Rf_error("internal: preprocess probe takes %d arguments, got %d", nargs,
             Rf_length(args));
  }

  SEXP data = VECTOR_ELT(args, 0);
  SEXP agg = VECTOR_ELT(args, 1);
  SEXP gind = VECTOR_ELT(args, 2);
  int64_t n_tax = (int64_t) scalar_real(VECTOR_ELT(args, 3), "n_tax");
  int64_t n_samp = (int64_t) scalar_real(VECTOR_ELT(args, 4), "n_samp");
  int64_t n_groups = (int64_t) scalar_real(VECTOR_ELT(args, 5), "n_groups");
  int32_t struc_zero = scalar_flag(VECTOR_ELT(args, 6), "struc_zero");
  double prv_cut = scalar_real(VECTOR_ELT(args, 7), "prv_cut");
  double lib_cut = scalar_real(VECTOR_ELT(args, 8), "lib_cut");
  double pseudo = scalar_real(VECTOR_ELT(args, 9), "pseudo");
  int32_t neg_lb = scalar_flag(VECTOR_ELT(args, 10), "neg_lb");

  if (n_tax <= 0 || n_samp <= 0 || n_groups < 0) {
    UNPROTECT_TO(0);
    Rf_error("n_tax and n_samp must be positive and n_groups non-negative, got %lld, %lld, %lld",
             (long long) n_tax, (long long) n_samp, (long long) n_groups);
  }
  if (struc_zero && n_groups == 0) {
    UNPROTECT_TO(0);
    Rf_error("struc_zero needs a group: the reference stops with \"group must be "
             "specified\" before reaching .get_struc_zero");
  }

  const double *dp = NULL, *ap = NULL, *gp = NULL;
  matrix_arg(data, "feature_table", NULL, &dp, &n_tax, &n_samp, NULL);
  if (n_tax != (int64_t) scalar_real(VECTOR_ELT(args, 3), "n_tax")) {
    UNPROTECT_TO(0);
    Rf_error("`feature_table` has %lld rows, not the %lld the caller declared",
             (long long) n_tax, (long long) scalar_real(VECTOR_ELT(args, 3), "n_tax"));
  }
  {
    int64_t ar = 0, ac = 0;
    matrix_arg(agg, "feature_table_aggregate", NULL, &ap, &ar, &ac, NULL);
    if (ar != n_tax || ac != n_samp) {
      UNPROTECT_TO(0);
      Rf_error("`feature_table_aggregate` is %lld x %lld but `feature_table` is %lld x %lld",
               (long long) ar, (long long) ac, (long long) n_tax, (long long) n_samp);
    }
  }
  if (n_groups > 0) {
    if (TYPEOF(gind) != REALSXP && TYPEOF(gind) != INTSXP) {
      UNPROTECT_TO(0);
      Rf_error("`group_indicator` must be a numeric matrix");
    }
    SEXP gd = Rf_getAttrib(gind, R_DimSymbol);
    if (gd == R_NilValue || TYPEOF(gd) != INTSXP || Rf_length(gd) != 2) {
      UNPROTECT_TO(0);
      Rf_error("`group_indicator` must have a two-element dim attribute");
    }
    if ((int64_t) INTEGER(gd)[0] != n_samp || (int64_t) INTEGER(gd)[1] != n_groups) {
      UNPROTECT_TO(0);
      Rf_error("`group_indicator` is %d x %d but should be %lld x %lld",
               INTEGER(gd)[0], INTEGER(gd)[1], (long long) n_samp, (long long) n_groups);
    }
    gp = TYPEOF(gind) == REALSXP ? REAL(gind) : NULL;
    if (gp == NULL) {
      UNPROTECT_TO(0);
      Rf_error("`group_indicator` must be double; it is an indicator, not a count");
    }
  }

  struct rb_preprocess in;
  memset(&in, 0, sizeof(in));
  in.data = dp; in.aggregate = ap; in.group_indicator = gp;
  in.n_tax = n_tax; in.n_samp = n_samp; in.n_groups = n_groups;
  in.struc_zero = struc_zero; in.prv_cut = prv_cut; in.lib_cut = lib_cut;
  in.pseudo = pseudo; in.neg_lb = neg_lb;

  /* Arguments 11-14 are the four reduction wrappers, in `struct rb_reducer` order.
   * They are parked for the bridge's call and cleared immediately after, so a later
   * call that somehow reached a reduction without going through this shim would find
   * NULL and say so rather than dereference it. */
  if (Rf_length(args) != 15) {
    UNPROTECT_TO(0);
    Rf_error("internal: preprocess probe takes 15 arguments, got %d", Rf_length(args));
  }
  for (int k = 11; k < 15; k++) {
    if (TYPEOF(VECTOR_ELT(args, k)) != CLOSXP) {
      UNPROTECT_TO(0);
      Rf_error("argument %d must be a reduction closure", k + 1);
    }
  }

  /* Two passes over the stages: the first asks how many doubles and integers the
   * bridge needs, the second allocates exactly that and fills it. The alternative --
   * handing the bridge buffers it sized itself -- would mean two calls and a second
   * entry point, and this way the allocation and the write cannot disagree. */
  int64_t f64_len = 0, i32_len = 0, names_len = 0;
  int64_t out_tax = 0, out_samp = 0;
  int32_t has_group = 0;
  SEXP prev_reduce_args = rb_reduce_args;
  rb_reduce_args = args;

  struct rb_stage *slots = (struct rb_stage *) R_alloc((size_t) N_PREPROCESS_STAGES,
                                                      sizeof(struct rb_stage));
  int64_t *name_off = (int64_t *) R_alloc((size_t) N_PREPROCESS_STAGES + 1,
                                          sizeof(int64_t));
  int8_t *name_bytes = (int8_t *) R_alloc((size_t) N_PREPROCESS_STAGES * 32 + 1,
                                          sizeof(int8_t));

  int32_t st = ancombc2_rb_preprocess_probe(
      &in, NULL, NULL, slots, N_PREPROCESS_STAGES, &f64_len, &i32_len, &out_tax,
      &out_samp, &has_group, name_bytes, name_off, &names_len);
  if (st != RB_OK) {
    rb_reduce_args = prev_reduce_args;
    UNPROTECT_TO(0);
    rb_error(st, "preprocess probe (sizing)");
  }
  if (f64_len < 0 || i32_len < 0) {
    rb_reduce_args = prev_reduce_args;
    UNPROTECT_TO(0);
    Rf_error("the bridge reported %lld doubles and %lld integers",
             (long long) f64_len, (long long) i32_len);
  }
  /* Real R vectors, not `R_alloc`.
   *
   * The bridge calls back into R for the row means, and those calls allocate. Holding
   * an `R_alloc` buffer across them is what produced the bug this comment replaces: the
   * `lib_size1` stage came back as uninitialised bytes while every other stage was
   * correct, because the buffer the bridge wrote into had been reclaimed by the time
   * the shim read it back. A `PROTECT`ed `REALSXP` is GC-managed, so it cannot be.
   *
   * The alternative -- having the bridge allocate -- would mean a second allocation
   * path in Rust and a second ownership story for no gain. */
  SEXP f64_sexp = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) (f64_len > 0 ? f64_len : 1)));
  nprot++;
  SEXP i32_sexp = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) (i32_len > 0 ? i32_len : 1)));
  nprot++;
  double *f64_buf = REAL(f64_sexp);
  int32_t *i32_buf = INTEGER(i32_sexp);
  st = ancombc2_rb_preprocess_probe(
      &in, f64_buf, i32_buf, slots, N_PREPROCESS_STAGES, &f64_len, &i32_len, &out_tax,
      &out_samp, &has_group, name_bytes, name_off, &names_len);
  rb_reduce_args = prev_reduce_args;
  if (st != RB_OK) { UNPROTECT_TO(0); rb_error(st, "preprocess probe"); }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, N_PREPROCESS_STAGES + 4)); nprot++;
  SEXP names_sexp = PROTECT(Rf_allocVector(STRSXP, N_PREPROCESS_STAGES + 4)); nprot++;
  static const char *tail_names[] = {"n_tax", "n_samp", "n_groups", "has_group"};
  /* The bridge's own names, checked against this file's table. Two lists of stage
   * names is two things to keep in step; comparing them means a change on either side
   * is an error naming both, rather than a stage silently reported under the wrong
   * name -- which would look like a numerical difference. */
  for (int64_t i = 0; i < N_PREPROCESS_STAGES; i++) {
    int64_t lo = name_off[i], hi = name_off[i + 1];
    if (hi < lo || hi > names_len) {
      UNPROTECT_TO(nprot);
      Rf_error("the bridge's stage names are not a monotone partition (%lld..%lld of %lld)",
               (long long) lo, (long long) hi, (long long) names_len);
    }
    char *nm = (char *) R_alloc((size_t) (hi - lo) + 1, sizeof(char));
    memcpy(nm, name_bytes + lo, (size_t) (hi - lo));
    nm[hi - lo] = '\0';
    if (strcmp(nm, preprocess_stage_names[i]) != 0) {
      UNPROTECT_TO(nprot);
      Rf_error("the bridge's stage %lld is \"%s\" and this shim expected \"%s\"",
               (long long) (i + 1), nm, preprocess_stage_names[i]);
    }
    SEXP stage = stage_to_sexp(&slots[i], f64_buf, i32_buf);
    SET_VECTOR_ELT(out, i, stage);
    /* Release the protect `stage_to_sexp` left behind, but do *not* touch `nprot`:
     * that protect is the callee's, so the counter never saw it go up and must not
     * see it come down. Counting it here drove `nprot` to -15 by the end of the loop,
     * which made the final `UNPROTECT_TO(0)` a no-op and leaked four real PROTECTs --
     * reported by R as "stack imbalance in '<-', 2 then 6". */
    UNPROTECT(1);
    SET_STRING_ELT(names_sexp, i, Rf_mkChar(nm));
  }
  for (int k = 0; k < 4; k++) {
    SET_VECTOR_ELT(out, N_PREPROCESS_STAGES + k, ScalarInteger(0));
    SET_STRING_ELT(names_sexp, N_PREPROCESS_STAGES + k, Rf_mkChar(tail_names[k]));
  }
  SET_VECTOR_ELT(out, N_PREPROCESS_STAGES + 0, ScalarInteger((int) out_tax));
  SET_VECTOR_ELT(out, N_PREPROCESS_STAGES + 1, ScalarInteger((int) out_samp));
  SET_VECTOR_ELT(out, N_PREPROCESS_STAGES + 2, ScalarInteger((int) n_groups));
  SET_VECTOR_ELT(out, N_PREPROCESS_STAGES + 3, ScalarLogical(has_group));
  Rf_setAttrib(out, R_NamesSymbol, names_sexp);

  UNPROTECT_TO(0);
  #undef UNPROTECT_TO
  return out;
}

/* ------------------------------------------------------------------------ */
/* C_ancombc2_rb_emit                                                        */
/* ------------------------------------------------------------------------ */

/* A double buffer of exactly `n` values, widening an INTSXP if that is what
 * arrived. The widening happens here, on the R main thread, so the wire only ever
 * has one numeric representation; the *original's* type is recorded separately
 * and put back by the plan.
 *
 * `store` is R_NilValue unless this call allocated. The allocation is returned
 * rather than protected inside the helper because a helper cannot reach the
 * caller's PROTECT counter, and an unbalanced PROTECT on an error path is worse
 * than a plain store. */
static const double *wire_doubles(SEXP x, const char *field, int64_t n, SEXP *store) {
  *store = R_NilValue;
  if (TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) {
    Rf_error("`%s` must be an integer or double vector, not %s", field,
             Rf_type2char(TYPEOF(x)));
  }
  if (Rf_xlength(x) != (R_xlen_t) n) {
    Rf_error("`%s` has %lld values, expected %lld", field,
             (long long) Rf_xlength(x), (long long) n);
  }
  if (TYPEOF(x) == REALSXP) return REAL(x);
  SEXP d = Rf_allocVector(REALSXP, (R_xlen_t) n);
  for (int64_t k = 0; k < n; k++) REAL(d)[k] = (double) INTEGER(x)[k];
  *store = d;
  return REAL(d);
}

/* Arguments, all in one VECSXP so the registration arity stays 1:
 *   0 n_tax  1 n_samp  2 p  3 n_taxa_bias
 *   4 beta_hat   5 var_hat   6 dof   7 vcov_hat (flat, row-major per taxon)
 *   8 y_bias_crt (as a matrix)   9 theta_hat
 *  10 delta_em  11 delta_wls  12 var_delta  13 s02
 *  14 O2
 *  15 taxa  16 taxa_bias
 *  17 taxon_names  18 sample_names  19 fix_eff
 *
 * Returns a list: `payload` (named by the plan), `plan_names`, `types`, `rtype`,
 * `dims`, `taxa`, `taxa_bias`, the three name vectors, and `dof_is_int`.
 *
 * Nothing here computes anything. It is the same discipline as the probe: show
 * the path loses nothing before trusting what travels over it.
 */
SEXP C_ancombc2_rb_emit(SEXP args) {
  int nprot = 0;
  #define UNPROTECT_TO(n) do { while (nprot > (n)) { UNPROTECT(1); nprot--; } } while (0)

  const int nargs = 20;
  if (TYPEOF(args) != VECSXP) Rf_error("internal: emit takes a VECSXP");
  if (Rf_length(args) != nargs) {
    Rf_error("internal: emit takes %d arguments, got %d", nargs, Rf_length(args));
  }

  int64_t n_tax  = (int64_t) scalar_real(VECTOR_ELT(args, 0), "n_tax");
  int64_t n_samp = (int64_t) scalar_real(VECTOR_ELT(args, 1), "n_samp");
  int64_t p      = (int64_t) scalar_real(VECTOR_ELT(args, 2), "p");
  int64_t n_tb   = (int64_t) scalar_real(VECTOR_ELT(args, 3), "n_taxa_bias");
  if (n_tax <= 0 || n_samp <= 0 || p <= 0) {
    Rf_error("n_tax, n_samp and p must all be positive, got %lld, %lld, %lld",
             (long long) n_tax, (long long) n_samp, (long long) p);
  }
  if (n_tb < n_tax) {
    UNPROTECT_TO(0);
    Rf_error("taxa_bias has %lld entries and taxa has %lld: O1 is a superset of O2, "
             "not the reverse", (long long) n_tb, (long long) n_tax);
  }
  const int64_t nx = n_tax * p;
  const int64_t ns = n_tax * n_samp;
  const int64_t nv = n_tax * p * p;

  /* `dof`'s original type is what decides whether it comes back as INTSXP, so it
   * is read before anything is widened. */
  SEXP dof_in = VECTOR_ELT(args, 6);
  const int dof_is_int = (TYPEOF(dof_in) == INTSXP);

  /* Eleven inputs, eleven slots. They are *not* pooled: every one of these
   * pointers stays live until `ancombc2_rb_emit` has read them, so sharing a slot
   * would leave the earlier buffer unprotected and still in use. */
  SEXP keep[11];
  memset(keep, 0, sizeof(keep));
  const double *beta   = wire_doubles(VECTOR_ELT(args, 4),  "beta_hat",   nx,     &keep[0]);
  const double *varhat = wire_doubles(VECTOR_ELT(args, 5),  "var_hat",    nx,     &keep[1]);
  const double *dof    = wire_doubles(dof_in,                 "dof",        nx,     &keep[2]);
  const double *vcov   = wire_doubles(VECTOR_ELT(args, 7),  "vcov_hat",   nv,     &keep[3]);
  const double *ybc    = wire_doubles(VECTOR_ELT(args, 8),  "y_bias_crt", ns,     &keep[4]);
  const double *theta  = wire_doubles(VECTOR_ELT(args, 9),  "theta_hat",  n_samp, &keep[5]);
  const double *dem    = wire_doubles(VECTOR_ELT(args, 10), "delta_em",   p,      &keep[6]);
  const double *dwl    = wire_doubles(VECTOR_ELT(args, 11), "delta_wls",  p,      &keep[7]);
  const double *vd     = wire_doubles(VECTOR_ELT(args, 12), "var_delta",  p,      &keep[8]);
  const double *s02v   = wire_doubles(VECTOR_ELT(args, 13), "s02",        p,      &keep[9]);
  const double *o2     = wire_doubles(VECTOR_ELT(args, 14), "O2",         ns,     &keep[10]);

  for (int k = 0; k < 11; k++) {
    if (keep[k] != R_NilValue) { PROTECT(keep[k]); nprot++; }
  }

  SEXP taxa = VECTOR_ELT(args, 15);
  SEXP taxa_bias = VECTOR_ELT(args, 16);
  if (TYPEOF(taxa) != INTSXP || Rf_xlength(taxa) != (R_xlen_t) n_tax) {
    UNPROTECT_TO(0);
    Rf_error("`taxa` must be an integer vector of length n_tax (%lld)", (long long) n_tax);
  }
  if (TYPEOF(taxa_bias) != INTSXP || Rf_xlength(taxa_bias) != (R_xlen_t) n_tb) {
    UNPROTECT_TO(0);
    Rf_error("`taxa_bias` must be an integer vector of length %lld", (long long) n_tb);
  }

  packed_names pn_tax, pn_samp, pn_fix;
  pack_names(&pn_tax,  VECTOR_ELT(args, 17), "taxon_names");
  pack_names(&pn_samp, VECTOR_ELT(args, 18), "sample_names");
  pack_names(&pn_fix,  VECTOR_ELT(args, 19), "fix_eff");
  if (pn_tax.n != n_tax) {
    UNPROTECT_TO(0);
    Rf_error("`taxon_names` has %lld entries but n_tax is %lld",
             (long long) pn_tax.n, (long long) n_tax);
  }
  if (pn_samp.n != n_samp) {
    UNPROTECT_TO(0);
    Rf_error("`sample_names` has %lld entries but n_samp is %lld",
             (long long) pn_samp.n, (long long) n_samp);
  }
  if (pn_fix.n != p) {
    UNPROTECT_TO(0);
    Rf_error("`fix_eff` has %lld entries but p is %lld", (long long) pn_fix.n, (long long) p);
  }

  /* The shapes, written out here rather than taken from the bridge. Two sources
   * of truth is the point: if the bridge's plan and this table ever disagree, the
   * loop below fails loudly instead of quietly reshaping the result. */
  struct shape { const char *name; int64_t rows, cols; int is_int; int is_list; };
  struct shape shapes[] = {
    {"theta_hat",  n_samp, 1,      0, 0},
    {"beta_hat",   n_tax,  p,      0, 0},
    {"var_hat",    n_tax,  p,      0, 0},
    {"dof",        n_tax,  p, dof_is_int, 0},
    {"vcov_hat",   n_tax,  p,      0, 1},
    {"delta_em",   p,      1,      0, 0},
    {"delta_wls",  p,      1,      0, 0},
    {"var_delta",  p,      1,      0, 0},
    {"s02",        p,      1,      0, 0},
    {"y_bias_crt", n_tax,  n_samp, 0, 0},
    {"O2",         n_tax,  n_samp, 0, 0},
  };
  const int64_t n_shapes = (int64_t) (sizeof(shapes) / sizeof(shapes[0]));

  /* The wire struct. Its field order matches `RawPayloads` in abi.rs, and the
   * compiler checks that here rather than leaving it to a reader. */
  struct rb_payloads {
    const double *beta_hat, *var_hat, *dof;
    int32_t dof_is_int;
    const double *vcov, *y_bias_crt, *theta_hat;
    const double *delta_em, *delta_wls, *var_delta, *s02, *o2;
    const int64_t *taxa, *taxa_bias;
    const int64_t *taxon_names_offsets; const uint8_t *taxon_names_bytes; int64_t n_taxon_names;
    const int64_t *sample_names_offsets; const uint8_t *sample_names_bytes; int64_t n_sample_names;
    const int64_t *fix_eff_offsets;      const uint8_t *fix_eff_bytes;      int64_t n_fix_eff;
    int64_t n_tax, n_samp, p;
    int64_t n_taxa_bias;
  } pay;

  memset(&pay, 0, sizeof(pay));
  pay.beta_hat = beta; pay.var_hat = varhat; pay.dof = dof; pay.dof_is_int = dof_is_int;
  pay.vcov = vcov; pay.y_bias_crt = ybc; pay.theta_hat = theta;
  pay.delta_em = dem; pay.delta_wls = dwl; pay.var_delta = vd; pay.s02 = s02v; pay.o2 = o2;
  pay.taxa = (const int64_t *) INTEGER(taxa);
  pay.taxa_bias = (const int64_t *) INTEGER(taxa_bias);
  pay.taxon_names_offsets = pn_tax.offsets;    pay.taxon_names_bytes = pn_tax.bytes;    pay.n_taxon_names = pn_tax.n;
  pay.sample_names_offsets = pn_samp.offsets; pay.sample_names_bytes = pn_samp.bytes; pay.n_sample_names = pn_samp.n;
  pay.fix_eff_offsets = pn_fix.offsets;       pay.fix_eff_bytes = pn_fix.bytes;       pay.n_fix_eff = pn_fix.n;
  pay.n_tax = n_tax; pay.n_samp = n_samp; pay.p = p; pay.n_taxa_bias = n_tb;

  /* Ask for the plan, allocate from it, then fill. The allocation and the plan
   * cannot disagree about a length, because the allocation came from the plan. */
  int64_t f64_len = 0, i32_len = 0, u8_len = 0, plan_len = 0;
  int32_t st = ancombc2_rb_emit_plan(&pay, &f64_len, &i32_len, &u8_len, &plan_len);
  if (st != RB_OK) { UNPROTECT_TO(0); rb_error(st, "emit plan"); }
  if (plan_len != n_shapes) {
    UNPROTECT_TO(0);
    Rf_error("the bridge's plan has %lld entries and this shim knows %lld",
             (long long) plan_len, (long long) n_shapes);
  }
  /* Only the f64 region must be non-empty. `i32_len == 0` is a legitimate plan --
   * it is what a double `dof` produces, since `dof` is the only integer entry --
   * and refusing it here would have made a fractional `df.residual` unrunnable.
   * `max(1, len)` below is for the same reason: `R_alloc(0, ...)` is not something
   * to rely on, and the buffer is never read when the region is empty. */
  if (f64_len <= 0 || i32_len < 0 || u8_len < 0 || plan_len <= 0) {
    UNPROTECT_TO(0);
    Rf_error("the plan reserved %lld doubles, %lld integers and %lld bytes across "
             "%lld entries, which cannot hold a result", (long long) f64_len,
             (long long) i32_len, (long long) u8_len, (long long) plan_len);
  }

  double  *out_f64 = (double *)  R_alloc((size_t) (f64_len > 0 ? f64_len : 1), sizeof(double));
  int32_t *out_i32 = (int32_t *) R_alloc((size_t) (i32_len > 0 ? i32_len : 1), sizeof(int32_t));
  int64_t *off64  = (int64_t *) R_alloc((size_t) plan_len, sizeof(int64_t));
  int64_t *off32  = (int64_t *) R_alloc((size_t) plan_len, sizeof(int64_t));
  /* The plan's names are short identifiers. The reservation is not a guess about
   * their content: what the bridge actually wrote is checked against it below. */
  const int64_t names_cap = plan_len * 64 + 1;
  int64_t *plan_names_off = (int64_t *) R_alloc((size_t) plan_len + 1, sizeof(int64_t));
  uint8_t *plan_names_bts = (uint8_t *) R_alloc((size_t) names_cap, sizeof(uint8_t));

  st = ancombc2_rb_emit(&pay, out_f64, out_i32, off64, off32,
                        plan_names_off, plan_names_bts, plan_len);
  if (st != RB_OK) { UNPROTECT_TO(0); rb_error(st, "emit"); }
  if (plan_names_off[plan_len] > names_cap) {
    UNPROTECT_TO(0);
    Rf_error("the bridge wrote %lld bytes of plan names into %lld reserved",
             (long long) plan_names_off[plan_len], (long long) names_cap);
  }

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 11)); nprot++;
  SEXP payload = PROTECT(Rf_allocVector(VECSXP, plan_len)); nprot++;
  SEXP plan_names = PROTECT(Rf_allocVector(STRSXP, plan_len)); nprot++;
  SEXP types = PROTECT(Rf_allocVector(STRSXP, plan_len)); nprot++;
  SEXP rtype = PROTECT(Rf_allocVector(INTSXP, plan_len)); nprot++;
  SEXP dims  = PROTECT(Rf_allocVector(INTSXP, plan_len * 2)); nprot++;

  for (int64_t i = 0; i < plan_len; i++) {
    const int64_t lo = plan_names_off[i], hi = plan_names_off[i + 1];
    SEXP nm = Rf_mkCharLen((const char *) plan_names_bts + lo, (int) (hi - lo));
    SET_STRING_ELT(plan_names, i, nm);
    const char *got = CHAR(nm);
    if (strcmp(got, shapes[i].name) != 0) {
      UNPROTECT_TO(nprot);
      Rf_error("the bridge's plan entry %lld is \"%s\" and this shim expected \"%s\"",
               (long long) (i + 1), got, shapes[i].name);
    }
    /* `Rf_mkString` returns a length-1 STRSXP, which `SET_STRING_ELT` refuses;
     * it wants the CHARSXP itself. `Rf_mkChar` is the one that returns it. */
    const char *tname = shapes[i].is_list ? "list-of-matrix"
                     : (shapes[i].is_int  ? "integer" : "double");
    SET_STRING_ELT(types, i, Rf_mkChar(tname));
    INTEGER(rtype)[i] = shapes[i].is_list ? 3L : (shapes[i].is_int ? 1L : 2L);
    INTEGER(dims)[2 * i]     = (int) shapes[i].rows;
    INTEGER(dims)[2 * i + 1] = (int) shapes[i].cols;

    const int64_t len = shapes[i].rows * shapes[i].cols;
    SEXP v;
    if (shapes[i].is_list) {
      /* A list of `rows` square matrices. The wire is row-major within a taxon
       * and R is column-major, so the transpose happens exactly here, in one
       * visible place. For a symmetric matrix it is invisible; for `p >= 3` it is
       * not, which is why it is written out rather than left implicit. */
      v = PROTECT(Rf_allocVector(VECSXP, shapes[i].rows)); nprot++;
      const double *base = out_f64 + off64[i];
      for (int64_t t = 0; t < shapes[i].rows; t++) {
        SEXP m = PROTECT(Rf_allocMatrix(REALSXP, (int) shapes[i].cols, (int) shapes[i].cols));
        nprot++;
        const double *src = base + t * p * p;
        /* R's element (i, j) lives at `i + j*p`. The wire is row-major, so its
         * element (i, j) is at `i*p + j`. Writing one into the other position --
         * `b + a*p` from `a*p + b` -- transposes every matrix, and a symmetric
         * covariance hides it. */
        for (int64_t i = 0; i < p; i++) {
          for (int64_t j = 0; j < p; j++) {
            REAL(m)[i + j * p] = src[i * p + j];
          }
        }
        SET_VECTOR_ELT(v, t, m);
        UNPROTECT(1); nprot--;
      }
    } else if (shapes[i].is_int) {
      v = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t) len)); nprot++;
      memcpy(INTEGER(v), out_i32 + off32[i], (size_t) len * sizeof(int32_t));
    } else {
      v = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) len)); nprot++;
      memcpy(REAL(v), out_f64 + off64[i], (size_t) len * sizeof(double));
    }
    /* No transpose for the matrices: the wire order already *is* R's
     * column-major order, because that is how the payload was captured. */
    if (!shapes[i].is_list && shapes[i].cols > 1) {
      SEXP d = PROTECT(Rf_allocVector(INTSXP, 2)); nprot++;
      INTEGER(d)[0] = (int) shapes[i].rows;
      INTEGER(d)[1] = (int) shapes[i].cols;
      Rf_setAttrib(v, R_DimSymbol, d);
      UNPROTECT(1); nprot--;
    }
    SET_VECTOR_ELT(payload, i, v);
    UNPROTECT(1); nprot--;
  }
  Rf_setAttrib(payload, R_NamesSymbol, plan_names);

  SET_VECTOR_ELT(out, 0, payload);
  SET_VECTOR_ELT(out, 1, plan_names);
  SET_VECTOR_ELT(out, 2, types);
  SET_VECTOR_ELT(out, 3, rtype);
  SET_VECTOR_ELT(out, 4, dims);
  SET_VECTOR_ELT(out, 5, taxa);
  SET_VECTOR_ELT(out, 6, taxa_bias);
  SET_VECTOR_ELT(out, 7, unpack_names(&pn_tax,  "taxon_names"));
  SET_VECTOR_ELT(out, 8, unpack_names(&pn_samp, "sample_names"));
  SET_VECTOR_ELT(out, 9, unpack_names(&pn_fix,  "fix_eff"));
  SET_VECTOR_ELT(out, 10, ScalarInteger(dof_is_int ? 1 : 0));

  static const char *out_names[] = {
    "payload", "plan_names", "types", "rtype", "dims",
    "taxa", "taxa_bias", "taxon_names", "sample_names", "fix_eff", "dof_is_int"
  };
  SEXP out_names_sexp = PROTECT(Rf_allocVector(STRSXP, 11)); nprot++;
  for (int k = 0; k < 11; k++) SET_STRING_ELT(out_names_sexp, k, Rf_mkChar(out_names[k]));
  Rf_setAttrib(out, R_NamesSymbol, out_names_sexp);

  UNPROTECT_TO(0);
  #undef UNPROTECT_TO
  return out;
}

/* ------------------------------------------------------------------------ */

/* Read one NUL-terminated string out of the bridge.
 *
 * `Rf_ScalarString(Rf_mkCharLen(...))` rather than the bare `Rf_mkCharLen`: R
 * coerces a returned CHARSXP into a length-1 STRSXP on the way back, but the
 * coercion is implicit and R rejects a CHARSXP that has picked up an attribute
 * along the way ("cannot have attributes on a CHARSXP"). Returning a STRSXP makes
 * the type a fact of this file rather than a property of the caller. */
static SEXP copy_bridge_string(int64_t (*fn)(uint8_t *, int64_t), const char *what) {
  int64_t n = fn(NULL, 0);
  uint8_t *buf = (uint8_t *) R_alloc((size_t) n + 1, sizeof(uint8_t));
  int64_t got = fn(buf, n + 1);
  (void) what;
  return Rf_ScalarString(Rf_mkCharLen((const char *) buf, (int) got));
}

SEXP C_ancombc2_rb_oracle_sha(void) {
  return copy_bridge_string(ancombc2_rb_oracle_sha, "oracle sha");
}

SEXP C_ancombc2_rb_version(void) {
  return copy_bridge_string(ancombc2_rb_version, "version");
}

static const R_CallMethodDef CallEntries[] = {
  /* Arity 1: the probe takes a single VECSXP, not 30 separate arguments.
   * Registering it as 30 would make `.Call` expect 30 arguments and reject the
   * list -- which is exactly what happened on the first run. */
  {"C_ancombc2_rb_probe",      (DL_FUNC) &C_ancombc2_rb_probe,       1},
  {"C_ancombc2_rb_emit",       (DL_FUNC) &C_ancombc2_rb_emit,        1},
  {"C_ancombc2_rb_preprocess_probe", (DL_FUNC) &C_ancombc2_rb_preprocess_probe, 1},
  {"C_ancombc2_rb_fit_probe",        (DL_FUNC) &C_ancombc2_rb_fit_probe,        1},
  {"C_ancombc2_rb_oracle_sha", (DL_FUNC) &C_ancombc2_rb_oracle_sha,  0},
  {"C_ancombc2_rb_version",    (DL_FUNC) &C_ancombc2_rb_version,     0},
  {NULL, NULL, 0}
};

void R_init_ANCOMBC(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  /* Symbols are looked up by name from `useDynLib(ANCOMBC, .registration = TRUE)`,
   * so a renamed Rust function is a load-time error rather than a run-time one. */
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);
}
