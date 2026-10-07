#!/usr/bin/env Rscript
# Dump the runtime facts a certified profile must match, as a single
# KEY=VALUE stream on stdout. scripts/verify_profile.py parses it.
#
# Run through the profile's own interpreter, e.g.
#   R_LIBS=<isolated original lib> \
#     micromamba run -p <env> Rscript scripts/probe_runtime.R
#
# Nothing here consults the oracle or a golden. It reports what the process
# actually is. No non-base package is used, so a missing dependency cannot
# make the probe itself unrunnable.

q <- function(x) paste0("'", gsub("'", "'\\''", x), "'")
kv <- function(key, value) cat(key, "=", q(value), "\n", sep = "")

kv("r_version_string", R.version.string)
kv("r_version", paste(R.version$major, R.version$minor, sep = "."))
kv("r_platform", R.version$platform)
kv("r_arch", R.version$arch)
kv("r_build", R.version$os)

kv("libpaths", paste(.libPaths(), collapse = "|"))

ext <- extSoftVersion()
kv("blas_vendor", if ("BLAS" %in% names(ext)) ext[["BLAS"]] else "<absent>")
kv("laplack_vendor", if ("LAPACK" %in% names(ext)) ext[["LAPACK"]] else "<absent>")

kv("locale", Sys.getlocale())
kv("rng_kind", paste(RNGkind(), collapse = "|"))
for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "RCPP_PARALLEL_NUM_THREADS",
            "ANCOMBC_ORACLE_DIR")) {
  kv(paste0("env_", v), Sys.getenv(v, unset = "<unset>"))
}

ip <- as.data.frame(utils::installed.packages(), stringsAsFactors = FALSE)
kv("installed_packages", paste(sort(paste0(ip$Package, "@", ip$Version)),
                               collapse = "|"))

found <- tryCatch(find.package("ANCOMBC", quiet = TRUE),
                  error = function(e) character())
kv("ancombc_found", if (length(found)) found[1] else "<not-found>")
if (length(found)) {
  d <- read.dcf(file.path(found[1], "DESCRIPTION"))
  kv("ancombc_version", as.character(utils::packageVersion("ANCOMBC")))
  kv("ancombc_depends", d[1, "Depends"])
  kv("ancombc_imports", d[1, "Imports"])
}
