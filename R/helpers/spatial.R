# =============================================================================
# R/helpers/spatial.R
# Shared spatial operations extracted from 04_R_management_zones_combined.R
# Phase 1: Pure functions — zero side effects, zero global access
# =============================================================================

make_raster_fast <- function(values, template_rast, cells, name) {
  r <- template_rast
  r[] <- NA
  r[cells] <- values
  names(r) <- name
  r
}

resample_to_template <- function(r, template, method = "bilinear") {
  r_resampled <- terra::resample(r, template, method = method)
  terra::mask(r_resampled, template)
}

crs_transform <- function(sf_obj, crs_code) {
  sf::st_transform(sf_obj, crs = crs_code)
}

extract_raster_values <- function(r, pts) {
  terra::extract(r, terra::vect(pts), ID = FALSE)[, 1, drop = TRUE]
}

stack_bands <- function(raster_list) {
  terra::rast(raster_list)
}
