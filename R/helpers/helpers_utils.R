# =============================================================================
# R/helpers/helpers_utils.R
# Utility helpers extracted from 04_R_management_zones_combined.R
# Phase 1: Pure functions — zero side effects, zero global access
# =============================================================================

.timer <- function(label, expr) {
  t0 <- proc.time()
  result <- force(expr)
  dt <- (proc.time() - t0)[3]
  cat(sprintf("  [%s] %.1f s\n", label, dt))
  invisible(result)
}
