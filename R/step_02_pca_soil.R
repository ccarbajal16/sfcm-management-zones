# =============================================================================
# R/step_02_pca_soil.R
# Pipeline Step 02 — Load soil CSV observations and transform to UTM
# Extracted from 04_R_management_zones_combined.R lines 120-134
# =============================================================================

run_pca_soil <- function(cfg) {
  cat("Loading soil data...\n")
  soil_raw    <- read_csv(cfg$soil_csv, show_col_types = FALSE)
  target_vars <- intersect(cfg$target_vars, names(soil_raw))
  cat("Target variables:", paste(target_vars, collapse = ", "), "\n")

  pts_sf <- st_as_sf(soil_raw,
                     coords = c(cfg$lon_col, cfg$lat_col),
                     crs    = 4326) |>
    st_transform(cfg$crs_code)

  soil_raw$x <- st_coordinates(pts_sf)[, 1]
  soil_raw$y <- st_coordinates(pts_sf)[, 2]

  pts_sf_utm <- st_as_sf(soil_raw, coords = c("x", "y"),
                         crs = cfg$crs_code, remove = FALSE)

  list(
    soil_raw    = soil_raw,
    pts_sf_utm  = pts_sf_utm,
    target_vars = target_vars
  )
}
