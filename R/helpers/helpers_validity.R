# =============================================================================
# R/helpers/helpers_validity.R
# Fuzzy clustering validity indices extracted from 04_R_management_zones_combined.R
# Phase 1: Pure functions — zero side effects, zero global access
# =============================================================================

compute_fpi <- function(U, C) {
  n <- nrow(U)
  (C / (C - 1)) * (1 - (1 / n) * sum(U^2))
}

compute_nce <- function(U, C) {
  n <- nrow(U)
  U_safe <- pmax(U, .Machine$double.eps)
  -(1 / (n * log(C))) * sum(U_safe * log(U_safe))
}

compute_xie_beni <- function(U, V, X, m) {
  n <- nrow(X)
  C <- nrow(V)
  Jm <- 0
  for (k in seq_len(C)) {
    diff <- sweep(X, 2, V[k, ])
    Jm <- Jm + sum(U[, k]^m * rowSums(diff^2))
  }
  if (C < 2) return(Inf)
  min_sep <- Inf
  for (j in 1:(C - 1)) {
    for (l in (j + 1):C) {
      d2 <- sum((V[j, ] - V[l, ])^2)
      if (d2 < min_sep) min_sep <- d2
    }
  }
  if (min_sep < .Machine$double.eps) return(Inf)
  Jm / (n * min_sep)
}

compute_davies_bouldin <- function(X, clusters, centroids) {
  C <- nrow(centroids)
  if (C < 2) return(NA_real_)
  S <- numeric(C)
  for (j in seq_len(C)) {
    idx_j <- which(clusters == j)
    if (length(idx_j) < 2) { S[j] <- 0; next }
    S[j] <- mean(sqrt(rowSums(sweep(X[idx_j, , drop = FALSE],
                                     2, centroids[j, ])^2)))
  }
  R_mat <- matrix(0, C, C)
  for (j in seq_len(C)) {
    for (l in seq_len(C)) {
      if (j == l) next
      d_jl <- sqrt(sum((centroids[j, ] - centroids[l, ])^2))
      if (d_jl < .Machine$double.eps) d_jl <- 1e-12
      R_mat[j, l] <- (S[j] + S[l]) / d_jl
    }
  }
  mean(apply(R_mat, 1, max))
}

compute_calinski_harabasz <- function(X, clusters, centroids) {
  n <- nrow(X)
  C <- nrow(centroids)
  if (C < 2 || C >= n) return(NA_real_)
  grand_mean <- colMeans(X)
  Bk <- sum(sapply(seq_len(C), function(j) {
    nj <- sum(clusters == j)
    nj * sum((centroids[j, ] - grand_mean)^2)
  }))
  Wk <- sum(sapply(seq_len(C), function(j) {
    idx_j <- which(clusters == j)
    if (length(idx_j) < 1) return(0)
    sum(sweep(X[idx_j, , drop = FALSE], 2, centroids[j, ])^2)
  }))
  if (Wk < .Machine$double.eps) return(NA_real_)
  (Bk / (C - 1)) / (Wk / (n - C))
}
