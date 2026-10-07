# Draw least-squares problems and record what R's own `Cdqrls` returns for them,
# as raw little-endian doubles, for the Rust transcription to be compared against
# bit for bit (`crates/ancombc2-core/tests/dqrls_vs_r.rs`).
#
#   Rscript scripts/make_dqrls_cases.R OUT.bin [SEED]
#
# Per case, in order: n, p, x (n*p), y (n), then R's coefficients (p), residuals
# (n), effects (n), qraux (p), pivot (p, 1-based), rank.
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
set.seed(if (length(args) >= 2) as.integer(args[[2]]) else 20261007L)

cases <- list()
add <- function(x, y) cases[[length(cases) + 1L]] <<- list(x = x, y = y)

for (n in c(3L, 5L, 8L, 12L, 17L, 31L, 32L, 33L, 64L, 100L, 257L, 1000L, 4097L, 20000L)) {
  for (p in unique(pmin(c(1L, 2L, 3L, 5L, 8L, 18L), n))) {
    x <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))[, seq_len(p), drop = FALSE]
    add(x, rnorm(n))                                    # full rank
    if (p >= 3L) {
      x2 <- x; x2[, 3L] <- x2[, 2L] * 2 + 1e-9 * rnorm(n)   # near collinear
      add(x2, rnorm(n))
      x3 <- x; x3[, 3L] <- 2 * x3[, 2L]                      # exactly aliased
      add(x3, rnorm(n))
      x4 <- x; x4[, 2L] <- x4[, 2L] * 1e12                   # badly scaled
      add(x4, rnorm(n))
      x5 <- x; x5[, 2L] <- 0                                 # zero column
      add(x5, rnorm(n))
    }
    add(x, rnorm(n) * 1e8)                              # big response
    add(x, round(rnorm(n), 1))                          # short decimals
  }
}
# indicator designs, the shape ANCOM-BC's own fits take
for (n in c(12L, 40L, 200L)) {
  g <- factor(rep_len(c("a", "b", "c"), n))
  add(model.matrix(~ g), rnorm(n))
  add(model.matrix(~ g + rnorm(n)), log(rpois(n, 20) + 0.5))
}

con <- file(out, "wb")
for (cs in cases) {
  x <- cs$x; storage.mode(x) <- "double"; y <- as.double(cs$y)
  z <- .Call(stats:::C_Cdqrls, x, y, 1e-7, FALSE)
  n <- nrow(x); p <- ncol(x)
  writeBin(as.double(c(n, p, x, y, z$coefficients, z$residuals, z$effects,
                       z$qraux, z$pivot, z$rank)), con)
}
close(con)
cat(length(cases), "cases ->", out, "\n")
