# =============================================================================
# R/step_09_prescriptions.R
# Pipeline Step 09 — Fertility prescriptions per management zone
# Extracted from 04_R_management_zones_combined.R lines 1297-1368
# =============================================================================

run_prescriptions <- function(cfg, profiling_result) {

  zone_summary <- profiling_result$zone_summary
  target_vars  <- profiling_result$target_vars
  k_opt        <- profiling_result$k_opt

  if (is.null(target_vars)) target_vars <- cfg$target_vars
  if (is.null(k_opt)) k_opt <- nrow(zone_summary)

  # ===========================================================================
  # STEP 9 — FERTILITY PRESCRIPTIONS
  # ===========================================================================

  cat("\n=== STEP 9: Fertility prescriptions ===\n")

  source("R/helpers/helpers_sfi.R")

  prescriptions <- zone_summary |>
    select(zone, n, ends_with("_mean")) |>
    rename_with(~str_remove(.x, "_mean"), ends_with("_mean")) |>
    mutate(
      SOC_status = if ("SOC" %in% names(pick(everything())))
        classify_thresh(SOC, cfg$thresh$SOC_low, cfg$thresh$SOC_high,
                        "Deficient — add organic amendments", "Adequate")
      else NA_character_,

      pH_status = if ("pH" %in% names(pick(everything())))
        case_when(
          pH < cfg$thresh$pH_low  ~ "Acidic — lime application recommended",
          pH > cfg$thresh$pH_high ~ "Alkaline — gypsum or S may be needed",
          TRUE                    ~ "Optimal pH range")
      else NA_character_,

      CEC_status = if ("CEC" %in% names(pick(everything())))
        classify_thresh(CEC, cfg$thresh$CEC_low, cfg$thresh$CEC_high,
                        "Low CEC — low nutrient retention", "High CEC")
      else NA_character_,

      VR_priority = case_when(
        str_detect(coalesce(SOC_status, ""), "Deficient") |
          str_detect(coalesce(pH_status,  ""), "Acidic|Alkaline") |
          str_detect(coalesce(CEC_status, ""), "Low")              ~ "HIGH",
        str_detect(coalesce(SOC_status, ""), "Marginal") |
          str_detect(coalesce(CEC_status, ""), "Marginal")         ~ "MEDIUM",
        TRUE                                                        ~ "LOW"
      )
    ) |>
    select(zone, n, any_of(target_vars),
           any_of(c("SOC_status", "pH_status", "CEC_status")), VR_priority)

  print(prescriptions)
  write_csv(prescriptions,
            file.path(cfg$out_results, "fertility_prescriptions.csv"))

  # -- Figure 41: Prescription summary -------------------------------------------

  zone_cols <- hcl.colors(k_opt, palette = "Set2")

  prescriptions |>
    mutate(zone = fct_reorder(zone, as.integer(str_extract(zone, "\\d+")))) |>
    pivot_longer(any_of(c("SOC_status", "pH_status", "CEC_status")),
                 names_to = "indicator", values_to = "status") |>
    mutate(
      indicator = recode(indicator,
                         "SOC_status" = "Soil Organic Carbon",
                         "pH_status"  = "Soil pH",
                         "CEC_status" = "CEC"),
      severity = case_when(
        str_detect(status, "Deficient|Acidic|Alkaline|Low") ~ "Action needed",
        str_detect(status, "Marginal")                      ~ "Monitor",
        TRUE                                                ~ "Adequate"
      ),
      severity = factor(severity, levels = c("Action needed", "Monitor", "Adequate"))
    ) |>
    ggplot(aes(zone, fill = severity)) +
    geom_bar(position = "stack") +
    facet_wrap(~indicator, nrow = 1) +
    scale_fill_manual(values = c("Action needed" = "#d73027",
                                 "Monitor"       = "#fee090",
                                 "Adequate"      = "#1a9850"), name = NULL) +
    labs(title = "Soil fertility status by management zone",
         x = NULL, y = "Count of indicators") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top")

  save_figure(last_plot(), file.path(cfg$out_figures, "39_prescription_summary.png"),
              width = 10, height = 5, dpi = 300)

  list(
    prescriptions = prescriptions
  )
}
