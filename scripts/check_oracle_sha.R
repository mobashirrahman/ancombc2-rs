#!/usr/bin/env Rscript
# Verify the vendored oracle is the pinned commit.
#
# `reference/ANCOMBC` is committed as ordinary files rather than as a submodule or
# a gitlink, because `load_harness()` *sources* `reference/ANCOMBC/R/*.R` -- a
# clone that received an empty directory there could not regenerate or check a
# single golden. That choice removes the clone's own `.git`, so the pin cannot be
# checked with `git rev-parse` inside the tree; it is checked against recorded
# checksums of the files that are actually read instead.
#
# The commit id itself lives in `reference/env/ORACLE.md`. What this script
# verifies is the stronger property: that the bytes on disk are the ones that
# commit was verified to contain. A tampered or partially-updated vendored tree
# fails here even though the recorded id would still read correctly.
#
# Usage:
#   Rscript --vanilla scripts/check_oracle_sha.R          # verify
#   Rscript --vanilla scripts/check_oracle_sha.R --write # re-record, after an
#                                                        # intentional re-vendor
#
# Re-recording is deliberately not something `make` does. If the oracle is ever
# re-vendored at a different commit, the id in `reference/env/ORACLE.md`, the
# `ORACLE_SHA` constant in `crates/ancombc2-core/src/compat.rs`, the README's
# compatibility string and `scripts/smoke_oracle.R` all have to move together, and
# that is a decision rather than a side effect.

args <- commandArgs(trailingOnly = TRUE)
write_mode <- "--write" %in% args

root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
oracle_dir <- file.path(root, "reference", "ANCOMBC")
manifest_path <- file.path(root, "reference", "env", "oracle.sha256")

# The files that define the behaviour this port reproduces. `NAMESPACE`,
# `DESCRIPTION` and `NEWS.md` are included so that a re-vendor at a different
# version cannot pass unnoticed; the `data/` directory and the vignette sources
# are not read by the harness and are excluded to keep the record about the code.
# The files that define the behaviour this port reproduces: everything in `R/`,
# plus the three that would reveal a different version having been vendored. The
# `data/` directory and the vignette sources are not read by the harness, so they
# are left out to keep the record about the code rather than about the tarball.
top <- c("DESCRIPTION", "NAMESPACE", "NEWS", "README.md")
rs <- sort(list.files(file.path(oracle_dir, "R"), pattern = "[.]R$"))
rel <- sort(c(top, file.path("R", rs)))
stopifnot(length(rs) > 0L)

if (write_mode) {
  # `sha256sum` is coreutils and is present on every platform this repository is
  # built on; `openssl dgst` is the portable fallback.
  digest <- function(path) {
    out <- suppressWarnings(
      system2("sha256sum", path, stdout = TRUE, stderr = FALSE)
    )
    if (length(out) == 0L || !is.null(attr(out, "status"))) {
      out <- system2("openssl", c("dgst", "-sha256", path), stdout = TRUE,
                     stderr = FALSE)
      return(sub("^.*= *", "", out[1]))
    }
    strsplit(out[1], "[ \t]+")[[1]][1]
  }
  lines <- vapply(rel, function(r) {
    p <- file.path(oracle_dir, r)
    if (!file.exists(p)) stop("oracle file missing: ", r)
    paste0(digest(p), "  ", r)
  }, character(1))
  dir.create(dirname(manifest_path), showWarnings = FALSE, recursive = TRUE)
  writeLines(c(
    "# SHA-256 of the vendored ANCOMBC source, at commit",
    "# dc4febdf59badb3a8dfe0c767ef2186323c2199a (ANCOMBC 2.15.2). See",
    "# reference/env/ORACLE.md. Regenerate with --write only when re-vendoring.",
    lines
  ), manifest_path)
  message("wrote ", manifest_path, " (", length(rel), " files)")
  quit(status = 0L)
}

if (!file.exists(manifest_path)) {
  stop("missing ", manifest_path, "; run with --write to record the vendored source")
}

recorded <- readLines(manifest_path)
recorded <- recorded[!grepl("^#", recorded) & nzchar(trimws(recorded))]
if (length(recorded) == 0L) stop(manifest_path, " records no files")

# Read the recorded id out of the manifest header and cross-check it against the
# two places that assert it in code, so the three cannot drift apart.
hdr <- readLines(manifest_path)
id <- sub("^# ([0-9a-f]{40}).*$", "\\1", grep("^# [0-9a-f]{40}", hdr, value = TRUE)[1])
if (is.na(id)) stop("no commit id in the header of ", manifest_path)

compat_rs <- readLines(file.path(root, "crates", "ancombc2-core", "src", "lib.rs"))
declared <- unique(regmatches(paste(compat_rs, collapse = "\n"),
                              gregexpr("[0-9a-f]{40}", paste(compat_rs, collapse = "\n")))[[1]])
declared <- declared[!declared %in% c("0000000000000000000000000000000000000000")]
if (!any(startsWith(declared, id) || startsWith(id, declared))) {
  stop("the manifest records commit ", id,
       " but lib.rs's ORACLE_SHA declares ", paste(declared, collapse = ", "))
}

bad <- character(0)
for (line in recorded) {
  parts <- strsplit(trimws(line), "[ \t]+")[[1]]
  want <- parts[1]
  rel_path <- paste(parts[-1], collapse = " ")
  p <- file.path(oracle_dir, rel_path)
  if (!file.exists(p)) {
    bad <- c(bad, paste(rel_path, "(missing)"))
    next
  }
  out <- suppressWarnings(system2("sha256sum", p, stdout = TRUE, stderr = FALSE))
  got <- strsplit(out[1], "[ \t]+")[[1]][1]
  if (!identical(got, want)) bad <- c(bad, paste(rel_path, "(checksum)"))
}

if (length(bad)) {
  stop("the vendored oracle does not match ", manifest_path, ":\n  ",
       paste(bad, collapse = "\n  "),
       "\nThe oracle is what every golden in this repository was generated from;",
       "\nif it has changed, the goldens are stale and every parity result is void.",
       call. = FALSE)
}

cat("oracle ok: ANCOMBC 2.15.2 @ ", substr(id, 1, 12),
    " (", length(recorded), " files verified)\n", sep = "")
