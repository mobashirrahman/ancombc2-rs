#!/usr/bin/env Rscript
# Install a minimal stand-in for `phyloseq`: the S4 class definitions only.
#
#   Rscript --vanilla scripts/install_phyloseq_shim.R
#   Rscript --vanilla scripts/install_phyloseq_shim.R --lib /path/to/rlib
#
# # Why this exists
# ----------------
# `dietswap` and `atlas1006`, two of the four Layer 4 datasets, ship as
# serialized `phyloseq` S4 objects in Bioconductor's `microbiome` package. Reading
# one needs the *class definitions*, because every base generic on the object
# dispatches on its class and looks for the package that defines it -- so
# `dim()`, `as.vector()` and `length()` all fail with "unable to find required
# package 'phyloseq'" on an object whose payload is a perfectly ordinary matrix.
#
# The real `phyloseq` cannot be installed on every machine: it depends on
# `RCurl`, which needs libcurl development headers, and installing those needs
# root. So this writes a package containing the class definitions and nothing
# else.
#
# # What it is not
# ---------------
# * It contains **no code from `phyloseq`**. Six `setClass` calls, written here.
# * It implements **no** `phyloseq` method. `as(otu_table(x), "matrix")` does not
#   work with it, and nothing in this repository calls it:
#   `scripts/prepare_realdata.R` reads the slots through `attributes()` instead,
#   because the payload is already a base matrix and data frame with its shape
#   and names in their attributes.
# * It is **shadowed by the real package** wherever that is installed, so this is
#   a fallback and not an override. The name is deliberately the real one: the
#   class attributes in a serialized object record `package = "phyloseq"`, and R
#   resolves them by that name.
#
# The alternative -- not having the datasets -- means the real-data layer covers
# two datasets instead of four, and loses the three-group case and the
# sensitivity analysis. That is a real loss of coverage, so the shim is worth its
# thirty lines.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L) return(default)
  args[[i + 1L]]
}
lib <- arg("--lib", Sys.getenv("R_LIBS_USER", unset = .libPaths()[1]))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("phyloseq", quietly = TRUE)) {
  cat("phyloseq is already installed; nothing to do.\n")
  cat("  (If that is the shim itself, delete ", file.path(lib, "phyloseq"),
      " and rerun to reinstall.)\n", sep = "")
  quit(status = 0)
}

src <- file.path(tempdir(), "phyloseq-shim")
unlink(src, recursive = TRUE)
dir.create(file.path(src, "R"), recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  "Package: phyloseq",
  "Type: Package",
  "Title: Class definitions only, for reading local .rda files",
  "Version: 0.0.1",
  "Description: A stand-in for the real phyloseq package, providing only the S4",
  "    class definitions needed to read a serialized phyloseq object's slots. It",
  "    contains no code from phyloseq and implements none of its methods. It",
  "    exists so that scripts/prepare_realdata.R can extract the count and",
  "    sample-data matrices from a .rda file on a machine where phyloseq itself",
  "    cannot be installed, because its dependency RCurl needs libcurl",
  "    development headers and installing those needs root. The representations",
  "    are deliberately loose and the reader goes through attributes() rather than",
  "    any phyloseq method, so no behaviour of the real package is relied upon.",
  "    If the real phyloseq is installed it takes precedence and nothing here is",
  "    used. See scripts/install_phyloseq_shim.R.",
  "License: Artistic-2.0",
  "Author: ancombc2-rs contributors",
  "Maintainer: ancombc2-rs contributors <noreply@example.invalid>"
), file.path(src, "DESCRIPTION"))

writeLines(
  "exportClasses(otu_table, sample_data, refseq, tree, taxonomyTable, physeq)",
  file.path(src, "NAMESPACE"))

writeLines(c(
  "# Class definitions only; see DESCRIPTION.",
  "#",
  "# The slot names are the union of the ones phyloseq has used. The serialized",
  "# objects in microbiome 1.24 carry `sam_data` and `phy_tree` (the older names)",
  "# while the current package calls them `sample_data` and `tree`, so a",
  "# definition matching only one set fails to find the other. Declaring the union",
  "# means both serializations read, and reading a slot an object does not carry",
  "# yields NULL -- which the caller sees rather than an error.",
  "setClass(\"otu_table\",",
  "         representation(data = \"ANY\", taxa_are_rows = \"logical\"),",
  "         prototype(data = matrix(numeric(0), 0, 0), taxa_are_rows = TRUE))",
  "setClass(\"sample_data\", representation(data = \"ANY\"))",
  "setClass(\"taxonomyTable\", representation(data = \"ANY\"))",
  "setClass(\"refseq\", representation(sequences = \"ANY\", names = \"ANY\"))",
  "setClass(\"tree\", representation(phylo = \"ANY\"))",
  "setClass(\"physeq\",",
  "         representation(otu_table = \"otu_table\",",
  "                        sample_data = \"sample_data\",",
  "                        sam_data = \"sample_data\",",
  "                        refseq = \"refseq\",",
  "                        tree = \"tree\",",
  "                        phy_tree = \"tree\",",
  "                        tax_table = \"taxonomyTable\"))"
), file.path(src, "R", "classes.R"))

cat("building the shim...\n")
build <- file.path(tempdir(), "phyloseq_0.0.1.tar.gz")
st <- system2("R", c("CMD", "build", "--no-build-vignettes", shQuote(src)),
              stdout = FALSE, stderr = FALSE)
if (!identical(st, 0L) || !file.exists(build)) {
  # `R CMD build` writes into the working directory, so look there too.
  build <- "phyloseq_0.0.1.tar.gz"
}
if (!file.exists(build)) stop("R CMD build did not produce a tarball")

cat("installing into ", lib, "\n", sep = "")
st <- system2("R", c("CMD", "INSTALL", "-l", shQuote(lib), shQuote(build)),
              stdout = FALSE, stderr = FALSE)
if (!identical(st, 0L)) stop("R CMD INSTALL failed")

.libPaths(c(lib, .libPaths()))
ok <- suppressWarnings(suppressMessages(requireNamespace("phyloseq", quietly = TRUE)))
cat("phyloseq class definitions loadable:", ok, "\n")
if (!ok) stop("the shim did not install cleanly")
cat("done\n")
