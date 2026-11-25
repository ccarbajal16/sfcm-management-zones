# =============================================================================
# R/step_05_kriging.R
# Pipeline Step 05 — Soil Fertility Index (SFI) + AHP-weighted kriging
# Extracted from 04_R_management_zones_combined.R lines 778-924
# =============================================================================

run_kriging <- function(cfg) {

  # ===========================================================================
  # Load soil data internally (needed for SFI computation)
  # ===========================================================================

  soil_raw    <- read_csv(cfg$soil_csv, show_col_types = FALSE)
  target_vars <- intersect(cfg$target_vars, names(soil_raw))

  pts_sf <- st_as_sf(soil_raw,
                     coords = c(cfg$lon_col, cfg$lat_col),
                     crs    = 4326) |>
    st_transform(cfg$crs_code)

  soil_raw$x <- st_coordinates(pts_sf)[, 1]
  soil_raw$y <- st_coordinates(pts_sf)[, 2]

  pts_sf_utm <- st_as_sf(soil_raw, coords = c("x", "y"),
                         crs = cfg$crs_code, remove = FALSE)

  # Load template from pca raster (needed for kriging grid)
  pca_rast_all <- rast(cfg$pca_raster)
  template <- pca_rast_all[[1]]

  # ===========================================================================
  # STEP 5 — PATHWAY B: SOIL FERTILITY INDEX
  # ===========================================================================

  cat("\n=== STEP 5: Pathway B — Soil Fertility Index ===\n")

  source("R/helpers/helpers_sfi.R")

  sfi_vars_avail <- intersect(names(cfg$sfi_vars), names(soil_raw))
  if (length(sfi_vars_avail) == 0)
    stop("No SFI variables found in soil CSV.")

  cat("  SFI variables:", paste(sfi_vars_avail, collapse = ", "), "\n")

  score_mat <- matrix(NA, nrow = nrow(soil_raw), ncol = length(sfi_vars_avail),
                      dimnames = list(NULL, sfi_vars_avail))
  for (v in sfi_vars_avail) {
    score_mat[, v] <- score_indicator(soil_raw[[v]], cfg$sfi_vars[[v]])
  }

  ahp_w <- cfg$sfi_ahp_weights[sfi_vars_avail]
  if (any(is.na(ahp_w))) {
    cat("  WARNING: missing AHP weights — using equal weights.\n")
    ahp_w[is.na(ahp_w)] <- 1
  }
  ahp_w <- ahp_w / sum(ahp_w)
  cat("  AHP weights:\n")
  print(round(ahp_w, 3))

  soil_raw$SFI <- as.numeric(score_mat %*% ahp_w)
  soil_raw$SFI_class <- cut(soil_raw$SFI, breaks = cfg$sfi_breaks,
                            labels = cfg$sfi_labels, include.lowest = TRUE)

  cat(sprintf("  SFI: %.3f - %.3f (mean = %.3f)\n",
              min(soil_raw$SFI, na.rm = TRUE),
              max(soil_raw$SFI, na.rm = TRUE),
              mean(soil_raw$SFI, na.rm = TRUE)))

  # -- Kriging SFI surface -------------------------------------------------------

  pts_sfi_utm <- st_as_sf(soil_raw |> filter(!is.na(SFI)),
                          coords = c("x", "y"), crs = cfg$crs_code)
  pts_sp <- as(pts_sfi_utm["SFI"], "Spatial")

  coarse_fact <- max(1, ceiling(sqrt(ncell(template) / 500000)))
  template_coarse <- aggregate(template, fact = coarse_fact)
  pred_grid_df <- as.data.frame(template_coarse, xy = TRUE, na.rm = TRUE)
  pred_grid_sf <- st_as_sf(pred_grid_df, coords = c("x", "y"),
                           crs = cfg$crs_code)
  pred_grid_sp <- as(pred_grid_sf, "Spatial")

  cat("  Kriging SFI (grid:", nrow(pred_grid_sp), "cells)...\n")
  krige_result <- tryCatch(
    .timer("autoKrige SFI", {
      automap::autoKrige(SFI ~ 1, input_data = pts_sp,
                         new_data = pred_grid_sp)
    }),
    error = function(e) {
      cat("  autoKrige failed:", conditionMessage(e), "— IDW fallback.\n")
      NULL
    }
  )

  sfi_cols <- c("Very Low" = "#d73027", "Low" = "#fc8d59",
                "Moderate" = "#fee08b", "High" = "#91cf60",
                "Very High" = "#1a9850")

  sfi_krige_rast <- NULL
  sfi_var_rast   <- NULL
  sfi_class_rast <- NULL

  if (!is.null(krige_result)) {
    krige_sf <- st_as_sf(krige_result$krige_output)
    krige_coords <- st_coordinates(krige_sf)
    krige_df <- st_drop_geometry(krige_sf)
    krige_df$x <- krige_coords[, 1]
    krige_df$y <- krige_coords[, 2]

    sfi_krige_rast <- rast(krige_df[, c("x", "y", "var1.pred")], type = "xyz")
    sfi_var_rast   <- rast(krige_df[, c("x", "y", "var1.var")],  type = "xyz")
    crs(sfi_krige_rast) <- paste0("EPSG:", cfg$crs_code)
    crs(sfi_var_rast)   <- paste0("EPSG:", cfg$crs_code)
    sfi_krige_rast <- resample(sfi_krige_rast, template, method = "bilinear")
    sfi_var_rast   <- resample(sfi_var_rast,   template, method = "bilinear")
    sfi_krige_rast <- clamp(sfi_krige_rast, 0, 1)

    rcl_mat <- matrix(
      c(cfg$sfi_breaks[-length(cfg$sfi_breaks)],
        cfg$sfi_breaks[-1],
        seq_along(cfg$sfi_labels)),
      ncol = 3
    )
    sfi_class_rast <- terra::classify(sfi_krige_rast, rcl_mat)
    levels(sfi_class_rast) <- data.frame(
      value = seq_along(cfg$sfi_labels), label = cfg$sfi_labels
    )

    writeRaster(sfi_krige_rast,  file.path(cfg$out_data, "zones_sfi.tif"),
                overwrite = TRUE)
    writeRaster(sfi_class_rast,  file.path(cfg$out_data, "zones_sfi_class.tif"),
                overwrite = TRUE)
    writeRaster(sfi_var_rast,    file.path(cfg$out_data, "sfi_kriging_variance.tif"),
                overwrite = TRUE)
  }

  # -- Figures 33-35: Pathway B maps --------------------------------------------

  if (!is.null(sfi_krige_rast)) {
    source("R/helpers/visualization.R")

    sfi_plot_rast <- downsample_continuous(sfi_krige_rast, cfg$plot_max_cells)

    p_sfi <- ggplot() +
      geom_spatraster(data = sfi_plot_rast) +
      scale_fill_distiller(palette = "RdYlGn", direction = 1,
                           name = "SFI", na.value = "grey90", limits = c(0, 1)) +
      geom_sf(data = pts_sfi_utm, aes(colour = SFI), size = 2, shape = 16) +
      scale_colour_distiller(palette = "RdYlGn", direction = 1,
                             name = "SFI (obs)", limits = c(0, 1)) +
      labs(title = "Soil Fertility Index — kriged surface (Pathway B)",
           subtitle = paste0("AHP: ",
                             paste(names(ahp_w), round(ahp_w, 2),
                                   sep = "=", collapse = ", "))) +
      theme_map()
    save_figure(p_sfi, file.path(cfg$out_figures, "33_sfi_map.png"),
                width = 8, height = 6.5, dpi = 300)

    sfi_class_plot <- downsample_for_plot(sfi_class_rast, cfg$plot_max_cells)
    p_sfi_class <- ggplot() +
      geom_spatraster(data = sfi_class_plot) +
      scale_fill_manual(values = sfi_cols, name = "Fertility\nclass",
                        na.value = "grey90") +
      labs(title = "SFI fertility classification",
           subtitle = "Breaks: 0 | 0.25 | 0.50 | 0.75 | 0.90 | 1.0") +
      theme_map()
    save_figure(p_sfi_class, file.path(cfg$out_figures, "34_sfi_class_map.png"),
                width = 8, height = 6.5, dpi = 300)

    sfi_var_plot <- downsample_continuous(sfi_var_rast, cfg$plot_max_cells)
    p_krig_var <- ggplot() +
      geom_spatraster(data = sfi_var_plot) +
      scale_fill_viridis_c(option = "B", name = "Kriging\nvariance",
                           na.value = "grey90") +
      labs(title = "SFI kriging prediction variance",
           subtitle = "Higher = greater interpolation uncertainty") +
      theme_map()
    save_figure(p_krig_var, file.path(cfg$out_figures, "35_sfi_kriging_variance.png"),
                width = 7, height = 6, dpi = 300)
  }

  list(
    sfi_krige_rast  = sfi_krige_rast,
    sfi_var_rast    = sfi_var_rast,
    sfi_class_rast  = sfi_class_rast,
    soil_raw        = soil_raw,
    pts_sf_utm      = pts_sf_utm,
    ahp_w           = ahp_w,
    pts_sfi_utm     = pts_sfi_utm,
    sfi_cols        = sfi_cols,
    target_vars     = target_vars
  )
}
