# =============================================================================
# R/step_03_sfcm.R
# Pipeline Step 03 — k-selection (plain CMeans from geocmeans).
# The SFCM final fit is done in step_04 after the alpha+window sweep.
# =============================================================================

run_sfcm <- function(cfg, pca_result, soil_result) {

  # ===========================================================================
  # STEP 3 — k-SELECTION (plain FCM, no spatial context)
  # ===========================================================================

  cat("\n=== STEP 3: k-selection (geocmeans::CMeans, no spatial) ===\n")

  source("R/helpers/helpers_validity.R")
  source("R/helpers/helpers_math.R")

  mat_all  <- pca_result$mat_all
  n_pixels <- pca_result$n_pixels

  set.seed(99)
  idx_sel <- sample(n_pixels, min(cfg$sample_size_ksel, n_pixels))
  sub_mat <- mat_all[idx_sel, , drop = FALSE]

  idx_sil <- sample(length(idx_sel), min(cfg$sample_size_sil, length(idx_sel)))
  sub_sil <- sub_mat[idx_sil, , drop = FALSE]

  k_range <- cfg$k_min:cfg$k_max
  n_k <- length(k_range)

  fpi_vec <- numeric(n_k); nce_vec <- numeric(n_k); sil_vec <- numeric(n_k)
  xb_vec  <- numeric(n_k); db_vec  <- numeric(n_k); ch_vec  <- numeric(n_k)

  cat("Computing 6 validity indices for k =", cfg$k_min, "to", cfg$k_max, "...\n")

  suppressPackageStartupMessages(library(geocmeans))

  .timer("k-selection loop (6 indices)", {
    for (i in seq_along(k_range)) {
      k <- k_range[i]
      cat("  k =", k, "...")

      res_fcm <- geocmeans::CMeans(sub_mat, k = k, m = cfg$fuzziness,
                                   standardize = FALSE, maxiter = 200,
                                   tol = 0.01, seed = 99,
                                   verbose = FALSE, init = "kpp")
      U      <- res_fcm$Belongings
      V      <- res_fcm$Centers
      cl_sub <- max.col(U, ties.method = "first")

      fpi_vec[i] <- compute_fpi(U, k)
      nce_vec[i] <- compute_nce(U, k)
      xb_vec[i]  <- compute_xie_beni(U, V, sub_mat, cfg$fuzziness)
      db_vec[i]  <- compute_davies_bouldin(sub_mat, cl_sub, V)
      ch_vec[i]  <- compute_calinski_harabasz(sub_mat, cl_sub, V)

      # Silhouette on the sil subsample (extract the matching rows from U)
      sil_vec[i] <- tryCatch(
        mean(geocmeans::calcSilhouetteIdx(sub_sil, U[idx_sil, , drop = FALSE]),
             na.rm = TRUE),
        error = function(e) NA_real_
      )

      cat(" done\n")
    }
  })

  fpi_n <- normalise(fpi_vec); nce_n <- normalise(nce_vec); sil_n <- normalise(sil_vec)
  xb_n  <- normalise(xb_vec);  db_n  <- normalise(db_vec);  ch_n  <- normalise(ch_vec)

  composite <- (1 - fpi_n) + (1 - nce_n) + sil_n + (1 - xb_n) + (1 - db_n) + ch_n

  val_df <- tibble(
    k         = k_range,
    FPI       = round(fpi_vec, 4),
    NCE       = round(nce_vec, 4),
    Silhouette = round(sil_vec, 4),
    XieBeni   = round(xb_vec, 4),
    DB        = round(db_vec, 4),
    CH        = round(ch_vec, 4),
    Composite = round(composite, 4)
  )
  print(val_df)
  write_csv(val_df, file.path(cfg$out_results, "validity_indices.csv"))

  k_opt <- if (!is.null(cfg$k_final)) cfg$k_final else k_range[which.max(composite)]
  cat("Selected k =", k_opt, "(composite validity score)\n")

  # -- Figure 29: Validity indices ------------------------------------------------

  p_fpi <- ggplot(val_df, aes(k, FPI)) +
    geom_line(colour = "#2166ac", linewidth = 1) +
    geom_point(size = 2.5, colour = "#2166ac") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "FPI (minimise)", y = "FPI") + theme_bw(base_size = 10)
  p_nce <- ggplot(val_df, aes(k, NCE)) +
    geom_line(colour = "#4393c3", linewidth = 1) +
    geom_point(size = 2.5, colour = "#4393c3") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "NCE (minimise)", y = "NCE") + theme_bw(base_size = 10)
  p_sil <- ggplot(val_df, aes(k, Silhouette)) +
    geom_line(colour = "#e36b2d", linewidth = 1) +
    geom_point(size = 2.5, colour = "#e36b2d") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Silhouette (maximise)", y = "Silhouette") + theme_bw(base_size = 10)
  p_xb <- ggplot(val_df, aes(k, XieBeni)) +
    geom_line(colour = "#b2182b", linewidth = 1) +
    geom_point(size = 2.5, colour = "#b2182b") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Xie-Beni (minimise)", y = "Xie-Beni") + theme_bw(base_size = 10)
  p_db <- ggplot(val_df, aes(k, DB)) +
    geom_line(colour = "#1a9850", linewidth = 1) +
    geom_point(size = 2.5, colour = "#1a9850") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Davies-Bouldin (minimise)", y = "DB") + theme_bw(base_size = 10)
  p_ch <- ggplot(val_df, aes(k, CH)) +
    geom_line(colour = "#762a83", linewidth = 1) +
    geom_point(size = 2.5, colour = "#762a83") +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Calinski-Harabasz (maximise)", y = "CH") + theme_bw(base_size = 10)
  p_composite <- ggplot(val_df, aes(k, Composite)) +
    geom_line(colour = "black", linewidth = 1.2) +
    geom_point(size = 3, colour = "black") +
    geom_vline(xintercept = k_opt, linetype = "dashed", colour = "red") +
    annotate("text", x = k_opt + 0.18, y = max(composite) * 0.96,
             label = paste("k =", k_opt), hjust = 0,
             colour = "red", fontface = "bold", size = 4) +
    scale_x_continuous(breaks = k_range) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.14))) +
    labs(title = "Composite validity (maximise)", y = "Score") +
    theme_bw(base_size = 10)

  p_validity <- (p_fpi | p_nce | p_sil) / (p_xb | p_db | p_ch) / p_composite +
    plot_annotation(title = "Extended cluster validity indices for zone number selection",
                    subtitle = "Composite = normalised sum of (1-FPI)+(1-NCE)+Sil+(1-XB)+(1-DB)+CH")
  save_figure(p_validity, file.path(cfg$out_figures, "29_validity_indices.png"),
              width = 14, height = 12, dpi = 300)

  list(
    k_opt        = k_opt,
    val_df       = val_df,
    composite    = composite,
    mat_scl      = pca_result$mat_scl,
    template     = pca_result$template,
    n_pixels     = n_pixels,
    pca_rast     = pca_result$pca_rast,
    n_pcs_cluster = pca_result$n_pcs_cluster,
    df_all       = pca_result$df_all,
    xy_all       = pca_result$xy_all,
    mat_all      = pca_result$mat_all,
    pts_sf_utm   = soil_result$pts_sf_utm
  )
}
