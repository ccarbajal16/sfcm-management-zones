# =============================================================================
# R/step_05_profiling.R
# Pipeline Step 05 — Zone profiling: statistical characterization per zone
# Consumes sfcm_alpha_result (from step_04), not ensemble_result.
# =============================================================================

run_zone_profiling <- function(cfg, sfcm_alpha_result) {

  zone_rast   <- sfcm_alpha_result$zone_rast
  pts_sf_utm  <- sfcm_alpha_result$pts_sf_utm
  target_vars <- intersect(cfg$target_vars, names(pts_sf_utm))
  k_opt       <- sfcm_alpha_result$k_opt

  # ===========================================================================
  # STEP 5 — ZONE PROFILING
  # ===========================================================================

  cat("\n=== STEP 5: Zone profiling ===\n")

  obs_zone <- terra::extract(zone_rast, vect(pts_sf_utm), ID = FALSE)[, "zone"]

  soil_zones <- pts_sf_utm |>
    st_drop_geometry() |>
    select(all_of(target_vars), x, y) |>
    mutate(zone = paste0("Zone ", obs_zone)) |>
    drop_na(zone) |>
    filter(!is.na(obs_zone))

  zone_summary <- soil_zones |>
    group_by(zone) |>
    summarise(
      n = n(),
      across(all_of(target_vars),
             list(mean   = ~mean(.x,   na.rm = TRUE),
                  sd     = ~sd(.x,     na.rm = TRUE),
                  median = ~median(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )

  print(zone_summary)
  write_csv(zone_summary, file.path(cfg$out_results, "zone_profiles.csv"))

  # -- Figure 33: Box plots ------------------------------------------------------

  zone_cols <- hcl.colors(k_opt, palette = "Set2")

  soil_long <- soil_zones |>
    pivot_longer(all_of(target_vars), names_to = "property", values_to = "value") |>
    mutate(property = factor(property, levels = target_vars))

  ggplot(soil_long, aes(zone, value, fill = zone)) +
    geom_boxplot(outlier.shape = 21, linewidth = 0.6, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", shape = 18,
                 size = 3, colour = "white") +
    facet_wrap(~property, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
    labs(title = "Soil property distribution by management zone",
         subtitle = "Diamond = mean", x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_angle_30())

  save_figure(last_plot(), file.path(cfg$out_figures, "36_zone_boxplots.png"),
              width = 11, height = 7, dpi = 300)

  # -- Figure 34: Radar chart ----------------------------------------------------

  radar_df <- zone_summary |>
    select(zone, ends_with("_mean")) |>
    rename_with(~str_remove(.x, "_mean"), ends_with("_mean")) |>
    pivot_longer(-zone, names_to = "property", values_to = "value") |>
    group_by(property) |>
    mutate(norm = (value - min(value)) / (max(value) - min(value) + 1e-9)) |>
    ungroup()

  ggplot(radar_df, aes(x = property, y = norm, colour = zone, group = zone)) +
    geom_polygon(aes(fill = zone), alpha = 0.15, linewidth = 1) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
    scale_fill_manual(values   = setNames(zone_cols, paste0("Zone ", 1:k_opt))) +
    coord_polar() +
    labs(title = "Zone profile radar chart (normalised)", x = NULL, y = NULL,
         colour = NULL, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom")

  save_figure(last_plot(), file.path(cfg$out_figures, "37_zone_radar.png"),
              width = 7, height = 7, dpi = 300)

  list(
    zone_summary = zone_summary,
    soil_zones   = soil_zones,
    k_opt        = k_opt,
    target_vars  = target_vars,
    pts_sf_utm   = pts_sf_utm,
    zone_cols    = zone_cols
  )
}

element_angle_30 <- function() {
  element_text(angle = 30, hjust = 1, size = 9)
}
