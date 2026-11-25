# =============================================================================
# main.R — Orchestrator: SFCM Management Zone Delineation Pipeline
#
# Entry point for the modular pipeline. Sources all modules and runs the
# 8-step workflow sequentially.
#
# Usage: Rscript main.R
# =============================================================================

# -- 0. Helpers ----------------------------------------------------------------
source("R/helpers/helpers_utils.R")
source("R/helpers/visualization.R")
source("R/helpers/spatial.R")
source("R/helpers/helpers_math.R")
source("R/helpers/helpers_sfi.R")
source("R/helpers/helpers_validity.R")

# -- 1. Configuration ---------------------------------------------------------
source("R/config.R")
CFG <- load_config("config/default_config.yaml")

# -- 2. Environment -----------------------------------------------------------
source("R/setup.R")
load_required_packages(CFG)
setup_environment(CFG)

# -- 3. Pipeline step modules -------------------------------------------------
source("R/step_01_pca_raster.R")
source("R/step_02_pca_soil.R")
source("R/step_03_sfcm.R")
source("R/step_04_sfcm_alpha.R")
source("R/step_05_profiling.R")
source("R/step_06_kriging.R")
source("R/step_06_validation.R")
source("R/step_07_prescriptions.R")
source("R/step_08_final_figure.R")

# -- 4. Run pipeline ----------------------------------------------------------

main <- function(cfg_override = list()) {
  cfg <- if (length(cfg_override) > 0L) {
    load_config("config/default_config.yaml", overrides = cfg_override)
  } else {
    CFG
  }

  pca_result <- tryCatch(
    run_pca_raster(cfg),
    error = function(e) stop("Step 01 (PCA Raster) failed: ", e$message, call. = FALSE)
  )

  soil_result <- tryCatch(
    run_pca_soil(cfg),
    error = function(e) stop("Step 02 (PCA Soil) failed: ", e$message, call. = FALSE)
  )

  sfcm_result <- tryCatch(
    run_sfcm(cfg, pca_result, soil_result),
    error = function(e) stop("Step 03 (SFCM) failed: ", e$message, call. = FALSE)
  )

  sfcm_alpha_result <- tryCatch(
    run_sfcm_alpha(cfg, sfcm_result),
    error = function(e) stop("Step 04 (SFCM alpha sweep) failed: ",
                             e$message, call. = FALSE)
  )

  profiling_result <- tryCatch(
    run_zone_profiling(cfg, sfcm_alpha_result),
    error = function(e) stop("Step 05 (Zone Profiling) failed: ",
                             e$message, call. = FALSE)
  )

  kriging_result <- tryCatch(
    run_kriging(cfg),
    error = function(e) stop("Step 06 (Kriging / SFI) failed: ",
                             e$message, call. = FALSE)
  )

  validation_result <- tryCatch(
    run_spatial_validation(cfg, profiling_result),
    error = function(e) stop("Step 07 (Spatial Validation) failed: ",
                             e$message, call. = FALSE)
  )

  prescriptions_result <- tryCatch(
    run_prescriptions(cfg, profiling_result),
    error = function(e) stop("Step 07 (Prescriptions) failed: ",
                             e$message, call. = FALSE)
  )

  build_final_figure(cfg, sfcm_result, sfcm_alpha_result,
                     kriging_result, soil_result, validation_result)

  cat("\nPipeline complete.\n")

  invisible(list(
    pca            = pca_result,
    soil           = soil_result,
    sfcm           = sfcm_result,
    sfcm_alpha     = sfcm_alpha_result,
    profiling      = profiling_result,
    kriging        = kriging_result,
    validation     = validation_result,
    prescriptions  = prescriptions_result
  ))
}

main()
