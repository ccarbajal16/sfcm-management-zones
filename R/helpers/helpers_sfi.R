# =============================================================================
# R/helpers/helpers_sfi.R
# Soil Fertility Index helpers extracted from 04_R_management_zones_combined.R
# Phase 1: Pure functions — zero side effects, zero global access
# =============================================================================

mf_s_shape <- function(x, a, b) {
  mid <- (a + b) / 2
  u   <- numeric(length(x))
  u[x <= a] <- 0
  u[x >= b] <- 1
  lo <- x > a & x <= mid
  hi <- x > mid & x < b
  u[lo] <- 2 * ((x[lo] - a) / (b - a))^2
  u[hi] <- 1 - 2 * ((x[hi] - b) / (b - a))^2
  u
}

mf_s_reverse <- function(x, a, b) 1 - mf_s_shape(x, a, b)

mf_bell <- function(x, a, b) {
  mid   <- (a + b) / 2
  sigma <- (b - a)
  exp(-((x - mid)^2) / (2 * sigma^2))
}

score_indicator <- function(x, params) {
  switch(params$role,
         more_better = mf_s_shape(x,   params$a, params$b),
         less_better = mf_s_reverse(x, params$a, params$b),
         optimum     = mf_bell(x,      params$a, params$b),
         stop("Unknown role: ", params$role))
}

classify_thresh <- function(val, low, high, low_label, high_label,
                            mid_label = "Marginal") {
  case_when(val < low ~ low_label, val > high ~ high_label, TRUE ~ mid_label)
}
