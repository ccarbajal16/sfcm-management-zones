# =============================================================================
# R/setup.R
# Environment setup — packages and output directories
# Phase 2: replaces package loading + dir creation in 04_R_management_zones_combined.R
# =============================================================================

load_required_packages <- function(cfg) {
  required <- c(
    "yaml",
    "terra",
    "sf",
    "tidyverse",
    "geocmeans",
    "ppclust",
    "cluster",
    "viridis",
    "patchwork",
    "tidyterra",
    "ClusterR",
    "matrixStats",
    "parallelDist",
    "spdep",
    "automap",
    "dunn.test"
  )

  new_pkgs <- required[!(required %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) install.packages(new_pkgs)
  # Load purrr LAST so its `map` / `map_dfr` are not shadowed by other
  # packages (notably mclust, which was previously in this list and breaks
  # downstream map() calls with "invalid 'length' argument").
  invisible(lapply(required, library, character.only = TRUE))
  suppressPackageStartupMessages(library(purrr))
}

setup_environment <- function(cfg) {
  dirs <- c(cfg$out_data, cfg$out_figures, cfg$out_results)
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  invisible(cfg)
}
