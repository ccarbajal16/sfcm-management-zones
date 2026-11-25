# =============================================================================
# R/step_06_validation.R
# Pipeline Step 06 — Spatial validation: Kruskal-Wallis, Dunn tests, Moran's I
# (renamed from step_08; also fixes the Moran's I figure bug by switching to
# scale() with a fallback for length-1 outputs and tighter NA handling)
# =============================================================================

run_spatial_validation <- function(cfg, profiling_result) {

  soil_zones  <- profiling_result$soil_zones
  target_vars <- profiling_result$target_vars

  # ===========================================================================
  # STEP 6 — SPATIAL VALIDATION
  # ===========================================================================

  cat("\n=== STEP 6: Spatial validation ===\n")

  if (length(target_vars) == 0 || nrow(soil_zones) == 0) {
    cat("  No target variables or empty soil_zones — skipping validation.\n")
    return(list(kw_results = NULL, moran_results = NULL, dunn_results = NULL))
  }

  kw_results <- purrr::map_dfr(target_vars, function(v) {
    df_v <- tryCatch(
      soil_zones |> select(zone, value = all_of(v)) |> drop_na(),
      error = function(e) return(NULL)
    )
    if (is.null(df_v) || n_distinct(df_v$zone) < 2 || nrow(df_v) < 5)
      return(NULL)
    kt <- kruskal.test(value ~ zone, data = df_v)
    tibble(property = v, KW_chi2 = round(kt$statistic, 3),
           KW_df = kt$parameter, KW_p = round(kt$p.value, 4),
           significant = kt$p.value < 0.05)
  })

  cat("\nKruskal-Wallis:\n")
  print(kw_results)

  dunn_vars <- if (nrow(kw_results) > 0 && "significant" %in% names(kw_results)) {
    kw_results |> filter(significant) |> pull(property)
  } else character(0)

  dunn_results <- if (length(dunn_vars) > 0) {
    purrr::map_dfr(dunn_vars, function(v) {
      df_v <- soil_zones |> select(zone, value = all_of(v)) |> drop_na()
      dt <- tryCatch(
        dunn.test::dunn.test(df_v$value, df_v$zone,
                             method = "bonferroni", altp = TRUE),
        error = function(e) {
          cat(sprintf("  dunn.test failed for %s: %s\n", v, conditionMessage(e)))
          return(NULL)
        }
      )
      if (is.null(dt)) return(NULL)
      tibble(property = v, comparison = dt$comparisons,
             Z = round(dt$Z, 3), p_adj = round(dt$altP.adjusted, 4))
    })
  } else {
    tibble(property = character(), comparison = character(),
           Z = numeric(), p_adj = numeric())
  }

  write_csv(bind_rows(kw_results),
            file.path(cfg$out_results, "kruskal_dunn_tests.csv"))
  write_csv(dunn_results,
            file.path(cfg$out_results, "dunn_posthoc_tests.csv"))

  # -- Moran's I -----------------------------------------------------------------

  cat("\nMoran's I...\n")

  coords_mat <- as.matrix(soil_zones[, c("x", "y")])
  n_pts <- nrow(coords_mat)
  k_nn  <- min(8, n_pts - 1)

  nb_knn <- tryCatch(
    spdep::knearneigh(coords_mat, k = k_nn) |> spdep::knn2nb(),
    error = function(e) {
      cat(sprintf("  knearneigh failed: %s\n", conditionMessage(e)))
      NULL
    }
  )

  lw_knn <- if (!is.null(nb_knn)) {
    tryCatch(
      spdep::nb2listw(nb_knn, style = "W"),
      error = function(e) {
        cat(sprintf("  nb2listw failed: %s\n", conditionMessage(e)))
        NULL
      }
    )
  } else NULL

  moran_results <- if (!is.null(lw_knn)) {
    purrr::map_dfr(target_vars, function(v) {
      x_v <- soil_zones[[v]]
      if (sum(!is.na(x_v)) < 10) return(NULL)
      mt <- tryCatch(
        spdep::moran.test(x_v, lw_knn, na.action = na.exclude),
        error = function(e) {
          cat(sprintf("  moran.test failed for %s: %s\n",
                      v, conditionMessage(e)))
          return(NULL)
        }
      )
      if (is.null(mt)) return(NULL)
      # Residual Moran's I after removing zone means
      zone_means <- ave(x_v, soil_zones$zone, FUN = function(z) mean(z, na.rm = TRUE))
      resid_v <- x_v - zone_means
      mt_res <- tryCatch(
        spdep::moran.test(resid_v, lw_knn, na.action = na.exclude),
        error = function(e) NULL
      )
      tibble(property = v,
             MoranI = round(mt$estimate["Moran I statistic"], 4),
             Expectation = round(mt$estimate["Expectation"], 6),
             Variance = round(mt$estimate["Variance"], 8),
             p_value = round(mt$p.value, 4),
             significant = mt$p.value < 0.05,
             MoranI_resid = if (!is.null(mt_res)) round(mt_res$estimate["Moran I statistic"], 4) else NA_real_,
             p_resid = if (!is.null(mt_res)) round(mt_res$p.value, 4) else NA_real_)
    })
  } else {
    tibble(property = character(), MoranI = numeric(),
           Expectation = numeric(), Variance = numeric(),
           p_value = numeric(), significant = logical())
  }

  cat("\nMoran's I results:\n")
  print(moran_results)
  write_csv(moran_results, file.path(cfg$out_results, "moranI_results.csv"))

  # -- Figure 35: Moran's I ------------------------------------------------------
  # Bug fix vs step_08_validation.R: we use a vector-safe scaler that always
  # returns a length-N numeric vector, and we build the scatter plots with
  # explicit `[ , drop = FALSE]` to avoid matrix-degenerate cases.

  safe_scale_vec <- function(x) {
    x <- as.numeric(x)
    finite_x <- x[is.finite(x)]
    if (length(finite_x) < 2L || sd(finite_x) < .Machine$double.eps) {
      return(rep(0, length(x)))
    }
    mu <- mean(finite_x)
    sdv <- stats::sd(finite_x)
    (x - mu) / sdv
  }

  tryCatch({

    if (nrow(moran_results) > 0) {
      p_moran_bar <- ggplot(moran_results,
                            aes(x = reorder(property, MoranI), y = MoranI,
                                fill = significant)) +
        geom_col() +
        geom_hline(yintercept = 0, linewidth = 0.6) +
        scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#b2b2b2"),
                          labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05"),
                          name = NULL) +
        labs(title = "Moran's I by soil property", x = NULL, y = "Moran's I") +
        coord_flip() + theme_bw(base_size = 11)

      sig_vars <- if ("significant" %in% names(moran_results)) {
        vv <- moran_results |> filter(significant) |> arrange(desc(MoranI)) |> pull(property)
        if (length(vv) == 0)
          moran_results |> slice_max(MoranI, n = 4) |> pull(property)
        else vv
      } else {
        moran_results |> slice_max(MoranI, n = 4) |> pull(property)
      }
      sig_vars <- head(sig_vars, 4)

      moran_scatter_list <- if (!is.null(lw_knn) && length(sig_vars) > 0) {
        Filter(Negate(is.null), purrr::map(sig_vars, function(v) {
          x_raw <- soil_zones[[v]]
          if (length(x_raw) == 0) return(NULL)
          x_fill <- ifelse(is.na(x_raw), mean(x_raw, na.rm = TRUE), x_raw)
          x_v    <- safe_scale_vec(x_fill)
          if (length(x_v) < 5L) return(NULL)

          wz_v <- tryCatch(
            spdep::lag.listw(lw_knn, x_v),
            error = function(e) rep(NA_real_, length(x_v))
          )
          wz_v <- as.numeric(wz_v)
          keep <- is.finite(x_v) & is.finite(wz_v)
          if (sum(keep) < 5) return(NULL)
          df_sc <- tibble(x_std = x_v[keep], wx = wz_v[keep])

          mi <- moran_results |> filter(property == v)
          I_val <- if (nrow(mi) > 0) mi$MoranI[1] else NA_real_
          p_val <- if (nrow(mi) > 0) mi$p_value[1] else NA_real_

          ggplot(df_sc, aes(x_std, wx)) +
            geom_point(alpha = 0.4, size = 1.2, colour = "#4393c3") +
            geom_smooth(method = "lm", se = FALSE, colour = "#d6604d",
                        linewidth = 0.8) +
            geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
            geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
            annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3,
                     label = sprintf("I = %.3f\np = %.4f", I_val, p_val),
                     size = 3, fontface = "italic") +
            labs(title = v, x = "z", y = "Wz") +
            theme_bw(base_size = 10)
        }))
      } else list()

      # Render bar + scatter separately, then combine via cowplot-free
      # patchwork `plot_layout` (which is more robust than the `/` operator
      # for heterogeneous plot lists).
      if (length(moran_scatter_list) > 0) {
        n_sc <- length(moran_scatter_list)
        # build the scatter grid in a 2-column layout
        p_scatter_combined <- wrap_plots(moran_scatter_list, ncol = 2)
        p_moran <- (p_moran_bar / p_scatter_combined) +
          plot_layout(heights = c(1, 2)) +
          plot_annotation(title = "Moran's I spatial autocorrelation",
                          theme = theme(plot.title = element_text(size = 13, face = "bold")))
      } else {
        p_moran <- p_moran_bar +
          plot_annotation(title = "Moran's I spatial autocorrelation")
      }

      save_figure(p_moran,
                  file.path(cfg$out_figures, "38_moran_correlogram.png"),
                  width = 9, height = 10, dpi = 300)

    } else {
      cat("  No Moran's I results to plot.\n")
    }

  }, error = function(e) {
    cat(sprintf("  WARNING: Moran's I figure generation failed: %s\n",
                conditionMessage(e)))
    cat("  Call:", deparse(conditionCall(e)), "\n")
    cat("  (Numerical results saved — skipping figure only.)\n")
  })

  list(
    kw_results    = kw_results,
    moran_results = moran_results,
    dunn_results  = dunn_results
  )
}
