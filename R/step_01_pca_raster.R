# =============================================================================
# R/step_01_pca_raster.R
# Pipeline Step 01 — Load and prepare PCA raster bands for clustering
# Extracted from 04_R_management_zones_combined.R lines 101-118
# =============================================================================

run_pca_raster <- function(cfg) {
  cat("Loading PC raster...\n")
  pca_rast_all  <- rast(cfg$pca_raster)
  n_pcs_cluster <- min(cfg$n_pcs_cluster, nlyr(pca_rast_all))
  pca_rast      <- pca_rast_all[[1:n_pcs_cluster]]
  names(pca_rast) <- paste0("PC", 1:n_pcs_cluster)

  df_all  <- .timer("Extract raster to data.frame", {
    as.data.frame(pca_rast, xy = TRUE, na.rm = TRUE)
  })
  xy_all  <- df_all[, c("x", "y")]
  mat_all <- as.matrix(df_all[, paste0("PC", 1:n_pcs_cluster)])

  center_vals <- colMeans(mat_all)
  scale_vals  <- apply(mat_all, 2, sd)
  mat_scl     <- scale(mat_all, center = center_vals, scale = scale_vals)

  n_pixels <- nrow(mat_all)
  cat("Pixels available:", n_pixels, "\n")

  template <- pca_rast[[1]]

  list(
    pca_rast   = pca_rast,
    mat_scl    = mat_scl,
    df_all     = df_all,
    xy_all     = xy_all,
    mat_all    = mat_all,
    template   = template,
    n_pixels   = n_pixels,
    n_pcs_cluster = n_pcs_cluster
  )
}
