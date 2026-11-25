# =============================================================================
# R/step_04_sfcm_alpha.R
# Pipeline Step 04 — Sweep alpha (spatial weight) and window size for SFCM,
# pick the combo that maximises a composite score, then re-fit the full grid
# with the chosen (alpha, window) so the final zone raster reflects the tuned
# spatial context. Produces figures 30, 31, 32.
# =============================================================================

run_sfcm_alpha <- function(cfg, sfcm_result) {

  cat("\n=== STEP 4: SFCM alpha + window selection ===\n")

  suppressPackageStartupMessages(library(geocmeans))

  k_opt    <- sfcm_result$k_opt
  pca_rast <- sfcm_result$pca_rast
  pca_list <- lapply(seq_len(nlyr(pca_rast)), function(j) pca_rast[[j]])
  names(pca_list) <- paste0("PC", seq_len(nlyr(pca_rast)))

  # -- Sub-sample the PCA raster for the alpha/window sweep ---------------------
  # We use a coarser aggregation to keep the sweep fast. The chosen (alpha,
  # window) is then re-fitted on the full grid at the end of this step.
  pca_sub <- .sfcm_subsample(pca_rast, n_target = cfg$sfcm_alpha_sel_size, seed = 7)
  pca_sub_list <- lapply(seq_len(nlyr(pca_sub)), function(j) pca_sub[[j]])
  names(pca_sub_list) <- paste0("PC", seq_len(nlyr(pca_sub)))
  n_sub <- sum(!is.na(terra::values(pca_sub[[1]])))
  cat("  Subsample size:", n_sub, "cells\n")

  alphas   <- cfg$sfcm_alpha_grid
  windows  <- cfg$sfcm_window_sizes
  n_runs   <- length(alphas) * length(windows)
  cat("  Sweep grid :", n_runs, "runs (", length(alphas), "alphas x",
      length(windows), "windows)\n")

  eval_df <- vector("list", n_runs)
  idx <- 1L
  for (a in alphas) {
    for (w in windows) {
      cat(sprintf("    alpha = %.2f  window = %dx%d ... ", a, w, w))
      window_mat <- matrix(1, nrow = w, ncol = w)

      t0 <- Sys.time()
      fit <- tryCatch(
        geocmeans::SFCMeans(
          data       = pca_sub_list,
          k          = k_opt,
          m          = cfg$fuzziness,
          alpha      = a,
          window     = window_mat,
          lag_method = cfg$sfcm_lag_method,
          maxiter    = cfg$sfcm_maxiter,
          tol        = cfg$sfcm_tol,
          standardize = FALSE,
          init       = cfg$sfcm_init,
          seed       = 42,
          verbose    = FALSE
        ),
        error = function(e) {
          cat("FAILED:", conditionMessage(e), "\n")
          NULL
        }
      )
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(fit)) {
        eval_df[[idx]] <- tibble(alpha = a, window = w, status = "error",
                                 silhouette = NA_real_, inertia = NA_real_,
                                 spconsist = NA_real_, composite = NA_real_,
                                 seconds = elapsed)
        idx <- idx + 1L
        next
      }

      sil  <- tryCatch(geocmeans::calcSilhouetteIdx(fit$Data, fit$Belongings),
                       error = function(e) NA_real_)
      iner <- tryCatch(geocmeans::calcexplainedInertia(fit$Data, fit$Belongings),
                       error = function(e) NA_real_)
      spc  <- tryCatch({
        w3 <- matrix(1, nrow = 3, ncol = 3)  # fixed for the spatial index calc
        geocmeans::spConsistency(fit, window = w3, nrep = 5)$Mean
      }, error = function(e) NA_real_)

      cat(sprintf("sil=%.3f  inertia=%.3f  spC=%.3f  (%.1fs)\n",
                  sil, iner, spc, elapsed))

      eval_df[[idx]] <- tibble(alpha = a, window = w, status = "ok",
                               silhouette = sil, inertia = iner,
                               spconsist = spc, seconds = elapsed)
      idx <- idx + 1L
    }
  }

  eval_df <- bind_rows(eval_df)
  write_csv(eval_df, file.path(cfg$out_results, "sfcm_alpha_sweep.csv"))

  # Composite score: rank-based + averaged. Higher = better.
  ok <- eval_df |> filter(status == "ok", is.finite(silhouette),
                          is.finite(inertia), is.finite(spconsist))
  if (nrow(ok) == 0) {
    stop("SFCM alpha sweep produced no valid runs.")
  }
  ok <- ok |>
    mutate(r_sil = rank(silhouette),
           r_in  = rank(inertia),
           r_spc = rank(-spconsist)) |>
    mutate(composite = 0.25 * r_sil + 0.25 * r_in + 0.50 * r_spc)

  spatial_ok <- ok |> filter(alpha > 0)
  best <- if (nrow(spatial_ok) > 0) {
    spatial_ok |> slice_max(composite, n = 1, with_ties = FALSE)
  } else {
    ok |> slice_max(composite, n = 1, with_ties = FALSE)
  }
  alpha_opt  <- best$alpha
  window_opt <- best$window
  cat(sprintf("\n  Selected: alpha = %.2f, window = %dx%d\n",
              alpha_opt, window_opt, window_opt))

  # -- Figure 30, 31, 32 are produced AFTER the final re-fit (so they show the
  #    chosen parameter combination, not the subsampled version). We stash the
  #    chosen (alpha, window) and continue.

  # -- Re-fit the FULL grid with the chosen (alpha, window) ---------------------
  cat(sprintf("\n  Re-fitting full grid at alpha=%.2f, window=%dx%d ...\n",
              alpha_opt, window_opt, window_opt))

  window_full <- matrix(1, nrow = window_opt, ncol = window_opt)
  sfcm_full <- .timer("Full-grid SFCMeans (tuned)", {
    geocmeans::SFCMeans(
      data       = pca_list,
      k          = k_opt,
      m          = cfg$fuzziness,
      alpha      = alpha_opt,
      window     = window_full,
      lag_method = cfg$sfcm_lag_method,
      maxiter    = cfg$sfcm_maxiter,
      tol        = cfg$sfcm_tol,
      standardize = FALSE,
      init       = cfg$sfcm_init,
      seed       = 42,
      verbose    = FALSE
    )
  })

  # Re-rasterize the tuned full-grid SFCM
  member_mat    <- sfcm_full$Belongings
  hard_zones_v  <- as.integer(sub("^V", "", sfcm_full$Groups))
  centroids     <- sfcm_full$Centers

  valid_mask    <- complete.cases(sfcm_result$mat_all)
  template      <- sfcm_result$template
  xy_valid      <- sfcm_result$xy_all[valid_mask, , drop = FALSE]
  cells_valid   <- cellFromXY(template, as.matrix(xy_valid))

  ord           <- t(apply(-member_mat, 1, order))
  max_member    <- member_mat[cbind(seq_len(nrow(member_mat)), ord[, 1])]
  second_max    <- member_mat[cbind(seq_len(nrow(member_mat)), ord[, 2])]
  uncertainty   <- 1 - max_member
  confusion_idx <- 1 - (max_member - second_max)

  source("R/helpers/spatial.R")
  zone_rast   <- make_raster_fast(hard_zones_v,  template, cells_valid, "zone")
  uncert_rast <- make_raster_fast(uncertainty,   template, cells_valid, "uncertainty")
  ci_rast     <- make_raster_fast(confusion_idx, template, cells_valid, "confusion_index")
  member_rasts <- lapply(seq_len(k_opt), function(j) {
    make_raster_fast(member_mat[, j], template, cells_valid, paste0("zone", j))
  })
  member_stack <- rast(member_rasts)

  writeRaster(zone_rast,    file.path(cfg$out_data, "zones_sfcm.tif"),
              overwrite = TRUE)
  writeRaster(uncert_rast,  file.path(cfg$out_data, "zones_uncertainty.tif"),
              overwrite = TRUE)
  writeRaster(ci_rast,      file.path(cfg$out_data, "zones_confusion.tif"),
              overwrite = TRUE)
  writeRaster(member_stack, file.path(cfg$out_data, "zones_membership.tif"),
              overwrite = TRUE)

  # -- Quality indices on the tuned full grid ----------------------------------
  full_sil  <- tryCatch(geocmeans::calcSilhouetteIdx(sfcm_full$Data,
                                                     sfcm_full$Belongings),
                        error = function(e) NA_real_)
  full_iner <- tryCatch(geocmeans::calcexplainedInertia(sfcm_full$Data,
                                                       sfcm_full$Belongings),
                        error = function(e) NA_real_)
  full_spc  <- tryCatch({
    w3 <- matrix(1, nrow = 3, ncol = 3)
    geocmeans::spConsistency(sfcm_full, window = w3, nrep = 5)$Mean
  }, error = function(e) NA_real_)
  cat(sprintf("  Full-grid:  silhouette = %.3f  inertia = %.3f  spConsist = %.3f\n",
              full_sil, full_iner, full_spc))

  # -- Figure 30: SFCM hard zone map (tuned) ------------------------------------

  zone_cols <- hcl.colors(k_opt, palette = "Set2")
  pts_sf_utm <- sfcm_result$pts_sf_utm

  p_zones <- ggplot() +
    geom_spatraster(data = as.factor(zone_rast)) +
    scale_fill_manual(values = setNames(zone_cols, as.character(1:k_opt)),
                      name = "Zone", na.value = "grey90") +
    geom_sf(data = pts_sf_utm, colour = "black", fill = "white",
            size = 1.5, shape = 21, stroke = 0.5) +
    labs(title = sprintf("SFCM zones (alpha = %.2f, window = %d)",
                         alpha_opt, window_opt),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  # -- Figure 31: SFCM confusion index ------------------------------------------

  p_ci <- ggplot() +
    geom_spatraster(data = ci_rast) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Confusion",
                         na.value = "grey90") +
    geom_sf(data = pts_sf_utm, colour = "black", fill = "white",
            size = 1.5, shape = 21, stroke = 0.5) +
    labs(title = "SFCM confusion index", x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  # -- Figure 32: alpha x window composite-score heatmap ------------------------

  p_alpha <- ggplot(ok, aes(factor(alpha), factor(window), fill = composite)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_tile(data = best, aes(x = factor(alpha), y = factor(window)),
              fill = NA, colour = "#d73027", linewidth = 1.4) +
    geom_text(aes(label = sprintf("%.1f", composite),
                  colour = ifelse(composite > 12, "black", "white")),
              fontface = "bold", size = 3.3, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_viridis_c(name = "Composite\nscore") +
    labs(title = "SFCM alpha x window composite score",
         subtitle = sprintf("best: alpha = %.2f, window = %dx%d",
                            alpha_opt, window_opt, window_opt),
         x = "alpha (spatial weight)", y = "window (cells)") +
    theme_bw(base_size = 11)

  source("R/helpers/visualization.R")
  save_figure(p_zones, file.path(cfg$out_figures, "30_zone_maps_pathwayA.png"),
              width = 7, height = 6, dpi = 300)
  save_figure(p_ci, file.path(cfg$out_figures, "31_confusion_index.png"),
              width = 7, height = 6, dpi = 300)
  save_figure(p_alpha, file.path(cfg$out_figures, "32_sfcm_alpha.png"),
              width = 8.5, height = 5, dpi = 300)

  list(
    sfcm_res       = sfcm_full,
    alpha_opt      = alpha_opt,
    window_opt     = window_opt,
    eval_df        = eval_df,
    best           = best,
    zone_rast      = zone_rast,
    uncert_rast    = uncert_rast,
    ci_rast        = ci_rast,
    member_stack   = member_stack,
    membership_mat = member_mat,
    centroids      = centroids,
    k_opt          = k_opt,
    pts_sf_utm     = pts_sf_utm,
    template       = sfcm_result$template,
    n_pixels       = sfcm_result$n_pixels,
    pca_rast       = sfcm_result$pca_rast,
    n_pcs_cluster  = sfcm_result$n_pcs_cluster,
    full_silhouette = full_sil,
    full_inertia   = full_iner,
    full_spconsist = full_spc
  )
}

# ---- helper: random subsample of a SpatRaster ------------------------------
# Uses terra::spatSample() for an unbiased random pick of n cells, then
# returns a SpatRaster with only those cells (rest masked to NA).
.sfcm_subsample <- function(rast, n_target, seed) {
  set.seed(seed)
  n_total <- sum(!is.na(terra::values(rast[[1]])))
  if (n_total <= n_target) return(rast)
  size <- min(n_target, n_total)
  s <- terra::spatSample(rast, size = size, method = "random",
                         na.rm = TRUE, as.raster = TRUE)
  s
}
