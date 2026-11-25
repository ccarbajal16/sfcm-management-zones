# =============================================================================
# R/step_08_final_figure.R
# Pipeline Step 08 — Final combined figure + summary
# Updated for the SFCM (geocmeans) pipeline. Panels: SFCM zones, SFCM
# uncertainty, SFI surface. No more grid-ClustGeo or ensemble cross-tab.
# =============================================================================

build_final_figure <- function(cfg, sfcm_result, sfcm_alpha_result,
                                kriging_result, soil_result, validation_result) {

  k_opt        <- sfcm_alpha_result$k_opt
  zone_rast    <- sfcm_alpha_result$zone_rast
  uncert_rast  <- sfcm_alpha_result$uncert_rast
  pts_sf_utm   <- sfcm_result$pts_sf_utm
  sfi_krige    <- kriging_result$sfi_krige_rast
  alpha_opt    <- sfcm_alpha_result$alpha_opt
  window_opt   <- sfcm_alpha_result$window_opt
  soil_raw     <- soil_result$soil_raw
  moran_results <- validation_result$moran_results
  kw_results   <- validation_result$kw_results
  zone_cols    <- hcl.colors(k_opt, palette = "Set2")

  cat("\n=== STEP 8: Final combined figure ===\n")

  source("R/helpers/visualization.R")

  # Panel 1: SFCM hard zones
  p_map_a <- ggplot() +
    geom_spatraster(data = as.factor(zone_rast)) +
    scale_fill_manual(values = setNames(zone_cols, as.character(1:k_opt)),
                      name = "Zone", na.value = "grey90") +
    geom_sf(data = pts_sf_utm, colour = "black", fill = "white",
            size = 1.5, shape = 21, stroke = 0.5) +
    labs(title = sprintf("SFCM zones (k = %d, alpha = %.2f, window = %d)",
                         k_opt, alpha_opt, window_opt),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  # Panel 2: SFCM uncertainty (1 - max membership)
  p_map_uncert <- ggplot() +
    geom_spatraster(data = uncert_rast) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         name = "Uncertainty", na.value = "grey90") +
    geom_sf(data = pts_sf_utm, colour = "black", fill = "white",
            size = 1.5, shape = 21, stroke = 0.5) +
    labs(title = "SFCM uncertainty", x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  # Panel 3: SFI surface (Pathway B)
  p_map_b <- if (!is.null(sfi_krige)) {
    sfi_plot2 <- downsample_continuous(sfi_krige, cfg$plot_max_cells)
    ggplot() +
      geom_spatraster(data = sfi_plot2) +
      scale_fill_distiller(palette = "RdYlGn", direction = 1,
                           name = "SFI", na.value = "grey90",
                           limits = c(0, 1)) +
      geom_sf(data = pts_sf_utm, colour = "black", fill = "white",
              size = 1.5, shape = 21, stroke = 0.5) +
      labs(title = "Pathway B: SFI surface", x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else NULL

  plot_list <- list(p_map_a, p_map_uncert)
  if (!is.null(p_map_b)) plot_list <- c(plot_list, list(p_map_b))

  wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = "SFCM Management Zone Delineation (geocmeans)",
      subtitle = paste0(
        "k = ", k_opt, " (composite validity) | ",
        "alpha = ", sprintf("%.2f", alpha_opt),
        " | window = ", window_opt, "x", window_opt,
        " | n = ", nrow(soil_raw), " SoilGrids 250m points")
    )

  save_figure(last_plot(), file.path(cfg$out_figures, "40_final_figure.png"),
              width = 6 * length(plot_list), height = 6, dpi = 300)

  # ===========================================================================
  # == SUMMARY ==================================================================
  # ===========================================================================

  cat("\n=== SFCM pipeline complete ===\n")
  cat(sprintf("  k selected         : %d (composite: FPI+NCE+Sil+XB+DB+CH)\n", k_opt))
  cat(sprintf("  Pathway A (SFCM)   : alpha = %.2f | window = %dx%d\n",
              alpha_opt, window_opt, window_opt))
  if (!is.null(sfcm_alpha_result$full_silhouette) &&
      is.finite(sfcm_alpha_result$full_silhouette)) {
    cat(sprintf("  SFCM silhouette    : %.3f\n", sfcm_alpha_result$full_silhouette))
  }
  if (!is.null(sfcm_alpha_result$full_inertia) &&
      is.finite(sfcm_alpha_result$full_inertia)) {
    cat(sprintf("  SFCM inertia       : %.3f\n", sfcm_alpha_result$full_inertia))
  }
  if (!is.null(sfcm_alpha_result$full_spconsist) &&
      is.finite(sfcm_alpha_result$full_spconsist)) {
    cat(sprintf("  SFCM spConsistency : %.3f\n", sfcm_alpha_result$full_spconsist))
  }
  cat("  Pathway B          : AHP-SFI + kriging\n")
  cat("  Outputs:\n")
  cat("    Rasters  :", file.path(cfg$out_data, "zones_*.tif"), "\n")
  cat("    Tables   :", file.path(cfg$out_results, "*.csv"), "\n")
  cat("    Figures  :", file.path(cfg$out_figures, "29-40_*.png"), "\n")
  if (!is.null(moran_results) && nrow(moran_results) > 0) {
    cat("\n  Moran's I (significant):\n")
    print(moran_results |> filter(significant))
  }
  if (!is.null(kw_results) && nrow(kw_results) > 0) {
    cat("\n  Kruskal-Wallis:\n")
    print(kw_results)
  }

  invisible(NULL)
}
