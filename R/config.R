# =============================================================================
# R/config.R
# Configuration loader — reads YAML, validates, returns immutable list
# Validates the SFCM (geocmeans) config schema. The old ClustGeo keys
# (alpha_grid, geo_unit, sample_size_geo, grid_clustgeo) are gone.
# =============================================================================

if (!"yaml" %in% rownames(installed.packages())) install.packages("yaml")

load_config <- function(yaml_path = "config/default_config.yaml",
                        overrides = list()) {
  if (!file.exists(yaml_path)) {
    stop("Config file not found: ", yaml_path, call. = FALSE)
  }

  cfg <- tryCatch(
    yaml::read_yaml(yaml_path),
    error = function(e) stop(
      "Failed to parse config file: ", conditionMessage(e), call. = FALSE
    )
  )

  if (length(overrides) > 0) {
    cfg <- modifyList(cfg, overrides)
  }

  if (is.null(cfg[["plot_max_cells"]])) {
    cfg$plot_max_cells <- 200000L
  }

  cfg$sfi_ahp_weights <- unlist(cfg$sfi_ahp_weights)

  validate_config(cfg)

  cfg
}

validate_config <- function(cfg) {
  required <- c(
    "pca_raster", "soil_csv",
    "crs_code", "lon_col", "lat_col", "target_vars",
    "n_pcs_cluster", "k_min", "k_max",
    "fuzziness", "fcm_nstart",
    "sample_size_ksel", "sample_size_fcm", "sample_size_sil",
    "sfcm_alpha_grid", "sfcm_window_sizes", "sfcm_lag_method",
    "sfcm_nstart", "sfcm_maxiter", "sfcm_tol", "sfcm_chunk_size",
    "sfcm_alpha_sel_size", "sfcm_init",
    "sfi_vars", "sfi_ahp_weights", "sfi_breaks", "sfi_labels",
    "thresh",
    "out_data", "out_figures", "out_results"
  )

  missing <- setdiff(required, names(cfg))
  if (length(missing) > 0) {
    stop("Missing required config keys: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  stopifnot(
    is.character(cfg$pca_raster) && length(cfg$pca_raster) == 1,
    is.character(cfg$soil_csv)   && length(cfg$soil_csv) == 1,
    is.numeric(cfg$crs_code)     && length(cfg$crs_code) == 1,
    is.character(cfg$lon_col)    && length(cfg$lon_col) == 1,
    is.character(cfg$lat_col)    && length(cfg$lat_col) == 1,
    is.character(cfg$target_vars) && length(cfg$target_vars) >= 1,
    is.numeric(cfg$n_pcs_cluster) && cfg$n_pcs_cluster >= 1,
    is.numeric(cfg$k_min)         && cfg$k_min >= 2,
    is.numeric(cfg$k_max)         && cfg$k_max > cfg$k_min,
    is.numeric(cfg$fuzziness)     && cfg$fuzziness > 1,
    is.numeric(cfg$fcm_nstart)    && cfg$fcm_nstart >= 1,
    is.numeric(cfg$sample_size_ksel) && cfg$sample_size_ksel > 0,
    is.numeric(cfg$sample_size_fcm)  && cfg$sample_size_fcm > 0,
    is.numeric(cfg$sample_size_sil)  && cfg$sample_size_sil > 0,
    is.numeric(cfg$sfcm_alpha_grid)   && all(cfg$sfcm_alpha_grid >= 0),
    is.numeric(cfg$sfcm_window_sizes) && all(cfg$sfcm_window_sizes >= 3) &&
      all(cfg$sfcm_window_sizes %% 2 == 1),  # must be odd
    is.character(cfg$sfcm_lag_method) && length(cfg$sfcm_lag_method) == 1,
    is.numeric(cfg$sfcm_nstart)    && cfg$sfcm_nstart >= 1,
    is.numeric(cfg$sfcm_maxiter)   && cfg$sfcm_maxiter >= 1,
    is.numeric(cfg$sfcm_tol)       && cfg$sfcm_tol > 0,
    is.numeric(cfg$sfcm_chunk_size) && cfg$sfcm_chunk_size >= 1,
    is.numeric(cfg$sfcm_alpha_sel_size) && cfg$sfcm_alpha_sel_size > 0,
    is.character(cfg$sfcm_init),
    is.list(cfg$sfi_vars),
    is.numeric(cfg$sfi_ahp_weights),
    is.numeric(cfg$sfi_breaks) && length(cfg$sfi_breaks) >= 2,
    is.character(cfg$sfi_labels) &&
      length(cfg$sfi_labels) == length(cfg$sfi_breaks) - 1,
    is.list(cfg$thresh),
    is.character(cfg$out_data)    && length(cfg$out_data) == 1,
    is.character(cfg$out_figures) && length(cfg$out_figures) == 1,
    is.character(cfg$out_results) && length(cfg$out_results) == 1
  )

  if (!is.null(cfg[["k_final"]])) {
    stopifnot(is.numeric(cfg$k_final), cfg$k_final >= cfg$k_min)
  }

  if (!is.null(cfg[["plot_max_cells"]])) {
    stopifnot(is.numeric(cfg$plot_max_cells), cfg$plot_max_cells > 0)
  }

  known <- c(required, "k_final", "plot_max_cells")
  unknown <- setdiff(names(cfg), known)
  if (length(unknown) > 0) {
    warning("Unknown config keys: ", paste(unknown, collapse = ", "),
            call. = FALSE)
  }

  invisible(cfg)
}
