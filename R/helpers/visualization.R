# =============================================================================
# R/helpers/visualization.R
# Shared ggplot2 helpers extracted from 04_R_management_zones_combined.R
# Phase 1: Pure functions — zero side effects, zero global access
# =============================================================================

downsample_for_plot <- function(r, max_cells) {
  if (ncell(r) > max_cells) {
    fact <- ceiling(sqrt(ncell(r) / max_cells))
    aggregate(r, fact = fact, fun = "modal", na.rm = TRUE)
  } else r
}

downsample_continuous <- function(r, max_cells) {
  if (ncell(r) > max_cells) {
    fact <- ceiling(sqrt(ncell(r) / max_cells))
    aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
  } else r
}

zone_map_plot <- function(r, k_opt, title = NULL, point_sf = NULL,
                          show_legend = TRUE) {
  zone_cols <- hcl.colors(k_opt, palette = "Set2")
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = as.factor(r)) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(zone_cols, as.character(seq_len(k_opt))),
      name = "Zone", na.value = "grey90"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  if (!is.null(point_sf)) {
    p <- p + ggplot2::geom_sf(data = point_sf, colour = "white",
                              size = 1.2, shape = 16)
  }

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title, x = NULL, y = NULL)
  } else {
    p <- p + ggplot2::labs(x = NULL, y = NULL)
  }

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}

continuous_plot <- function(r, palette = "viridis", limits = NULL,
                            title = NULL) {
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = r)

  if (palette == "distiller") {
    p <- p + ggplot2::scale_fill_distiller(
      palette = "RdYlGn", direction = 1,
      name = NULL, na.value = "grey90", limits = limits
    )
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(
      option = if (palette == "viridis") "D" else palette,
      name = NULL, na.value = "grey90", limits = limits
    )
  }

  p <- p + ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title, x = NULL, y = NULL)
  } else {
    p <- p + ggplot2::labs(x = NULL, y = NULL)
  }

  p
}

theme_map <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
}

save_figure <- function(fig, path, width, height, dpi = 300) {
  ggplot2::ggsave(path, fig, width = width, height = height, dpi = dpi)
}
