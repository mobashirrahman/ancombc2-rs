# Private package library (nloptr, ...) required by the ANCOMBC oracle.
# See reference/env/ORACLE.md.
# Set ANCOMBC_RLIB to the directory holding them; unset, nothing is prepended.
local({
  lib <- Sys.getenv("ANCOMBC_RLIB", "")
  if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))
})
