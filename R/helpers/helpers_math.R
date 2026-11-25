# =============================================================================
# R/helpers/helpers_math.R
# Math helpers extracted from 04_R_management_zones_combined.R
# Phase 1: Pure functions — zero side effects, zero global access
# =============================================================================

normalise <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

safe_scale_matrix <- function(x) {
  center <- colMeans(x, na.rm = TRUE)
  spread <- apply(x, 2, sd, na.rm = TRUE)
  spread[is.na(spread) | spread < .Machine$double.eps] <- 1
  scale(x, center = center, scale = spread)
}
