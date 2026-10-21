######################################################################
# 03_presentation_figures.R
# UNSW ASOC x Taylor Fry Case Competition 2026 | Task 1
#
# Run AFTER 02_modelling_updated.R from the project root.
# Uses the existing model outputs; it does NOT fit or select models.
# Exports 16:9, 300-dpi, white-background PNGs to:
#   temp_output/modelling/presentation_figures/
#
# Every chart uses an accessible light-blue -> navy palette.
# FY2027 forecasts are estimates of PEOPLE with each outcome, not
# distinct incident events. Risk rankings/driver coefficients are
# predictive/associational and are NOT evidence of causal effects.
######################################################################

required_packages <- c("dplyr", "tidyr", "readr", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install these packages: install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)

model_dir <- "temp_output/modelling"
figure_dir <- file.path(model_dir, "presentation_figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_output <- function(name, required = TRUE) {
  path <- file.path(model_dir, name)
  if (!file.exists(path)) {
    if (required) {
      stop("Missing modelling output: ", path, "\n",
           "Run source('02_modelling_updated.R') before this script.")
    }
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

# Colours: intentionally within one blue family (not the EDA's red/orange).
BLUE <- list(
  navy = "#103C68", dark = "#155E9A", mid = "#318AC5",
  sky = "#71BCE6", light = "#A8DBF4", ice = "#DDF2FC",
  pale = "#F2F9FD", ink = "#213E59", muted = "#64819A",
  grid = "#E1EEF5", white = "#FFFFFF"
)
sev_colours <- c(Mild = BLUE$light, CloseCall = BLUE$mid, Severe = BLUE$navy)
sev_levels <- c("Mild", "CloseCall", "Severe")
month_order <- month.abb[c(9:12, 1:8)]

# Large, high-contrast typography for projected presentation slides.
theme_slide <- function() {
  theme_minimal(base_size = 15, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = BLUE$white, colour = NA),
      panel.background = element_rect(fill = BLUE$white, colour = NA),
      plot.title = element_text(colour = BLUE$navy, face = "bold", size = 23,
                                margin = margin(b = 7)),
      plot.subtitle = element_text(colour = BLUE$muted, size = 13.5,
                                   margin = margin(b = 17)),
      plot.caption = element_text(colour = BLUE$muted, size = 10, hjust = 0,
                                  margin = margin(t = 16)),
      axis.title = element_text(colour = BLUE$ink, size = 13),
      axis.text = element_text(colour = BLUE$ink, size = 12),
      axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = BLUE$grid, linewidth = 0.45),
      panel.grid.major.y = element_blank(),
      strip.background = element_rect(fill = BLUE$pale, colour = NA),
      strip.text = element_text(colour = BLUE$navy, face = "bold", size = 13,
                                margin = margin(7, 7, 7, 7)),
      legend.title = element_text(colour = BLUE$ink, face = "bold", size = 12),
      legend.text = element_text(colour = BLUE$ink, size = 11),
      legend.position = "bottom",
      plot.margin = margin(18, 28, 18, 22)
    )
}

save_slide <- function(plot, filename, width = 13.333, height = 7.5) {
  ggplot2::ggsave(
    filename = file.path(figure_dir, filename), plot = plot,
    width = width, height = height, units = "in", dpi = 300,
    bg = BLUE$white, limitsize = FALSE
  )
  cat("Saved: ", file.path(figure_dir, filename), "\n", sep = "")
}

# -------------------------------------------------------------------
# FIGURE 1 | FY2025 model-selection comparison
# Relative excess deviance allows outcomes with very different
# absolute incident frequencies to appear on one comparable chart.
# -------------------------------------------------------------------
validation <- read_output("01_validation_comparison.csv")
stopifnot(all(c("outcome", "model", "poisson_deviance") %in% names(validation)))

model_lab <- c(
  Poisson = "Poisson GLM", NegBin = "Negative binomial",
  NegBin_GLMM = "NB mixed model", TwoStage = "Two-stage"
)
validation_plot_data <- validation %>%
  filter(is.finite(poisson_deviance), poisson_deviance > 0) %>%
  group_by(outcome) %>%
  mutate(
    best_deviance = min(poisson_deviance),
    excess_pct = 100 * (poisson_deviance / best_deviance - 1),
    is_winner = abs(poisson_deviance - best_deviance) < 1e-8
  ) %>%
  ungroup() %>%
  mutate(
    outcome = factor(outcome, levels = c(sev_levels, "Total")),
    model_label = factor(
      if_else(model %in% names(model_lab), unname(model_lab[model]), model),
      levels = rev(unname(model_lab))
    ),
    performance = if_else(is_winner, "Selected", "Alternative"),
    label = if_else(is_winner, "BEST", sprintf("+%.1f%%", excess_pct))
  )

p1 <- ggplot(validation_plot_data,
             aes(x = excess_pct, y = model_label, fill = performance)) +
  geom_col(width = 0.58, show.legend = FALSE) +
  geom_text(aes(label = label), hjust = -0.12,
            colour = BLUE$ink, size = 4.3, fontface = "bold") +
  facet_wrap(~ outcome, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = c(Selected = BLUE$dark, Alternative = BLUE$light)) +
  scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1),
                     expand = expansion(mult = c(0, 0.22), add = c(0, 0.5))) +
  labs(
    title = "Which model predicts incident counts best?",
    subtitle = "FY2025 validation | Excess Poisson deviance relative to the best model in each category",
    x = "Excess deviance vs best model (lower is better)", y = NULL,
    caption = "Model selection uses FY2025 only; FY2026 is reserved for an independent diagnostic test."
  ) + theme_slide()
save_slide(p1, "01_validation_model_comparison.png")

# -------------------------------------------------------------------
# FIGURE 2 | Independent FY2026 holdout: predicted vs observed totals
# Panels use separate y scales so rare Severe outcomes remain legible.
# -------------------------------------------------------------------
test <- read_output("03_final_test_metrics.csv")
stopifnot(all(c("outcome", "model", "actual_total", "predicted_total") %in% names(test)))

test_long <- test %>%
  select(outcome, model, actual_total, predicted_total) %>%
  pivot_longer(c(actual_total, predicted_total),
               names_to = "series", values_to = "people") %>%
  mutate(
    outcome = factor(outcome, levels = c(sev_levels, "Total")),
    series = factor(series, levels = c("actual_total", "predicted_total"),
                    labels = c("Observed", "Predicted")),
    value_label = if_else(people >= 100, comma(round(people)),
                          number(people, accuracy = 0.1))
  )

p2 <- ggplot(test_long, aes(x = series, y = people, fill = series)) +
  geom_col(width = 0.60, show.legend = FALSE) +
  geom_text(aes(label = value_label), vjust = -0.45, size = 4.5,
            colour = BLUE$ink, fontface = "bold") +
  facet_wrap(~ outcome, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(Observed = BLUE$light, Predicted = BLUE$dark)) +
  scale_y_continuous(labels = label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Independent test: predicted versus actual",
    subtitle = "FY2026 holdout | Only models chosen on FY2025 validation",
    x = NULL, y = "People involved in incidents",
    caption = "Panels have different y-axis scales. 'Total' is predicted independently of the three severity models."
  ) + theme_slide() +
  theme(axis.text.x = element_text(face = "bold", size = 11))
save_slide(p2, "02_FY2026_holdout_calibration.png")

# -------------------------------------------------------------------
# FIGURE 3 | FY2027 expected incidence by month and severity
# -------------------------------------------------------------------
monthly <- read_output("10_FY2027_by_month.csv")
stopifnot(all(c("month", "visits", sev_levels) %in% names(monthly)))

monthly_long <- monthly %>%
  mutate(month = factor(as.character(month), levels = month_order)) %>%
  pivot_longer(all_of(sev_levels), names_to = "severity", values_to = "expected") %>%
  mutate(severity = factor(severity, levels = sev_levels))

p3 <- ggplot(monthly_long, aes(x = month, y = expected, fill = severity)) +
  geom_col(width = 0.76, show.legend = FALSE) +
  facet_wrap(~ severity, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = sev_colours) +
  scale_y_continuous(labels = label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "When will beach safety incidents occur?",
    subtitle = "FY2027 baseline forecast | September 2026 – August 2027",
    x = "Month", y = "Expected people affected",
    caption = "Predictions assume FY2024–2026 average monthly visits/weather and continuing FY2026 safety coverage."
  ) + theme_slide()
save_slide(p3, "03_FY2027_incidents_by_month.png")

# -------------------------------------------------------------------
# FIGURE 4 | Beaches with the most expected deaths (absolute burden)
# This ranking is by expected NUMBER, not per-visitor risk.
# -------------------------------------------------------------------
beaches <- read_output("11_FY2027_beach_ranking.csv")
stopifnot(all(c("beach_nm", "visits", "Severe", "severe_per_million_visits") %in%
              names(beaches)))

top_beaches <- beaches %>%
  arrange(desc(Severe)) %>%
  slice_head(n = 12) %>%
  mutate(beach_nm = reorder(as.character(beach_nm), Severe))

p4 <- ggplot(top_beaches, aes(x = Severe, y = beach_nm, fill = Severe)) +
  geom_col(width = 0.69) +
  geom_text(aes(label = sprintf("%.2f", Severe)),
            hjust = -0.18, size = 4.2, fontface = "bold", colour = BLUE$ink) +
  scale_fill_gradient(low = BLUE$sky, high = BLUE$navy, guide = "none") +
  scale_x_continuous(labels = label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.20))) +
  labs(
    title = "Where is the greatest expected loss of life?",
    subtitle = "Top 12 beaches by forecast drowning deaths | FY2027",
    x = "Expected deaths in FY2027", y = NULL,
    caption = "Ranked by absolute expected deaths, not deaths per million visits."
  ) + theme_slide() +
  theme(axis.text.y = element_text(size = 12.5), panel.grid.major.x = element_line(colour = BLUE$grid))
save_slide(p4, "04_FY2027_top12_beaches_deaths.png")

# -------------------------------------------------------------------
# FIGURE 5 | Beach comparison: death burden vs visitors and risk rate
# Helps explain why absolute demand and per-visit risk differ.
# -------------------------------------------------------------------
beach_scatter <- beaches %>%
  mutate(visits_m = visits / 1e6) %>%
  arrange(desc(Severe))
label_beaches <- beach_scatter %>% slice_head(n = 6)

p5 <- ggplot(beach_scatter,
             aes(x = visits_m, y = Severe, colour = severe_per_million_visits)) +
  geom_point(size = 3.5, alpha = 0.88) +
  scale_colour_gradient(low = BLUE$light, high = BLUE$navy,
                        name = "Expected deaths\nper million visits") +
  scale_x_continuous(labels = label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_continuous(labels = label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0.02, 0.16))) +
  labs(
    title = "High incident burden is not always high risk per visit",
    subtitle = "Each point represents one beach | FY2027 baseline forecast",
    x = "Forecast annual visits (millions)", y = "Expected drowning deaths",
    caption = "Colour shows the expected death rate per million visits. Labels identify the top six by absolute deaths."
  ) + theme_slide() +
  theme(legend.position = "right")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p5 <- p5 + ggrepel::geom_text_repel(
    data = label_beaches, aes(label = beach_nm),
    colour = BLUE$ink, size = 3.3, fontface = "bold",
    box.padding = 0.35, point.padding = 0.25,
    min.segment.length = 0, max.overlaps = Inf,
    show.legend = FALSE, seed = 2026
  )
} else {
  p5 <- p5 + geom_text(
    data = label_beaches, aes(label = beach_nm),
    hjust = 0, vjust = -0.6, size = 3.1,
    colour = BLUE$ink, check_overlap = TRUE, show.legend = FALSE
  )
}
save_slide(p5, "05_FY2027_beach_risk_vs_footfall.png")

# -------------------------------------------------------------------
# FIGURE 6 | Drivers in the selected Severe model
# If two-stage wins: show both stages separately with CORRECT units:
#   Stage 1: NB total incident-rate ratios
#   Stage 2: conditional odds ratios for Severe among incidents
# If direct model wins: Severe incident-rate ratios only.
# Scaling continuous covariates avoids treating a 1-point change in
# first-aid coverage as comparable to a full campaign-index change.
# -------------------------------------------------------------------
coeff <- read_output("08_driver_associations.csv")
selected <- read_output("02_selected_models.csv")
stopifnot(all(c("outcome", "model", "effect_type", "term", "estimate_log", "std_error") %in%
              names(coeff)))
sev_winner <- selected %>% filter(outcome == "Severe") %>% pull(model)
if (length(sev_winner) != 1) stop("Expected exactly one selected Severe model.")

if (sev_winner == "TwoStage") {
  sev_coeff <- coeff %>%
    filter(outcome %in% c("Severe / Total (frequency stage)",
                          "Severe | incident (share stage)")) %>%
    mutate(stage = if_else(grepl("frequency stage", outcome, fixed = TRUE),
                           "Stage 1 · Total incident frequency (rate ratio)",
                           "Stage 2 · Severe conditional odds (odds ratio)"))
} else {
  sev_coeff <- coeff %>%
    filter(outcome == "Severe", model == sev_winner) %>%
    mutate(stage = "Direct severe incident frequency (rate ratio)")
}

# Each label describes the X change corresponding to the plotted ratio.
driver_steps <- tibble::tibble(
  term = c("swell_max_m", "n_rips", "tourist_ratio", "log_lg_hrs",
           "log_drone_hrs", "campaign_idx", "first_aid_cap", "pa_system", "flags"),
  driver = c("Maximum swell (+1 m)", "Average rips (+1)",
             "Non-local visitor share (+10 pp)",
             "Log(1 + lifeguard hours) (+1)",
             "Log(1 + drone hours) (+1)",
             "Campaign index (+0.25)",
             "First-aid coverage (+10 pts)",
             "PA installed (0 → 1)", "Warning flags installed (0 → 1)"),
  step = c(1, 1, 0.10, 1, 1, 0.25, 10, 1, 1)
)

driver_data <- sev_coeff %>%
  inner_join(driver_steps, by = "term") %>%
  mutate(
    driver = factor(driver, levels = rev(driver_steps$driver)),
    ratio_scaled = exp(estimate_log * step),
    lower_scaled = exp((estimate_log - 1.96 * std_error) * step),
    upper_scaled = exp((estimate_log + 1.96 * std_error) * step)
  ) %>%
  filter(is.finite(ratio_scaled), is.finite(lower_scaled),
         is.finite(upper_scaled), lower_scaled > 0)

if (nrow(driver_data) == 0) {
  warning("No matching Severe coefficients were found; skipping drivers plot.")
} else {
  stages <- unique(driver_data$stage)
  driver_data$stage <- factor(driver_data$stage, levels = stages)
  stage_colors <- setNames(c(BLUE$mid, BLUE$navy)[seq_along(stages)], stages)
  p6 <- ggplot(driver_data, aes(x = ratio_scaled, y = driver, colour = stage)) +
    geom_vline(xintercept = 1, colour = BLUE$muted, linetype = "dashed",
               linewidth = 0.7) +
    geom_segment(aes(x = lower_scaled, xend = upper_scaled, yend = driver),
                 linewidth = 1.1) +
    geom_point(size = 3.2) +
    facet_wrap(~ stage, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = stage_colors, guide = "none") +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    labs(
      title = "What factors are associated with severe incidents?",
      subtitle = paste0("Selected model: ", sev_winner,
                        " | Points: relative ratio; whiskers: approximate 95% CI"),
      x = "Rate ratio / odds ratio (log scale; dashed line = 1)", y = NULL,
      caption = "Units of change are stated on the y-axis. Stage 1 and Stage 2 ratios are NOT interchangeable.\nAssociations are adjusted, not causal intervention effects."
    ) + theme_slide() +
    theme(axis.text.y = element_text(size = 11),
          plot.caption = element_text(size = 9.3, colour = BLUE$muted))
  save_slide(p6, "06_severe_model_driver_associations.png")
}

# -------------------------------------------------------------------
# FIGURE 7 | Severe incidents: expected count vs expected risk per visit
# Helps answer both 'when most deaths' and 'when highest risk per visit'.
# -------------------------------------------------------------------
monthly_severe <- monthly %>%
  mutate(
    month = factor(as.character(month), levels = month_order),
    deaths_per_million = Severe / visits * 1e6
  ) %>%
  select(month, Severe, deaths_per_million) %>%
  pivot_longer(c(Severe, deaths_per_million),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("Severe", "deaths_per_million"),
                         labels = c("Expected deaths", "Deaths per million visits")))

p7 <- ggplot(monthly_severe, aes(x = month, y = value, group = 1)) +
  geom_line(colour = BLUE$dark, linewidth = 1.5) +
  geom_point(colour = BLUE$navy, size = 3.4) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0.07, 0.16))) +
  labs(
    title = "Seasonal demand and per-visit risk are different",
    subtitle = "Predicted deaths and exposure-adjusted rate | FY2027",
    x = "Month", y = NULL,
    caption = "Death rate = expected deaths / forecast visits × 1 million. Month-specific forecasts depend on future-input assumptions."
  ) + theme_slide()
save_slide(p7, "07_FY2027_severe_seasonality_and_risk.png")

# -------------------------------------------------------------------
# FIGURE 8 | Forecast sensitivity to future beach visitor volume
# Assumes +/- 10% visitor scenarios and all other predictors held fixed.
# -------------------------------------------------------------------
sensitivity <- read_output("14_visitor_sensitivity.csv", required = FALSE)
if (!is.null(sensitivity) &&
    all(c("visitor_multiplier", "Mild", "CloseCall", "Severe") %in% names(sensitivity))) {
  sensitivity_long <- sensitivity %>%
    pivot_longer(all_of(sev_levels), names_to = "severity", values_to = "expected") %>%
    group_by(severity) %>%
    mutate(baseline = expected[which.min(abs(visitor_multiplier - 1))],
           change_pct = 100 * (expected / baseline - 1)) %>%
    ungroup() %>%
    mutate(severity = factor(severity, levels = sev_levels),
           visitor_change_pct = 100 * (visitor_multiplier - 1))

  p8 <- ggplot(sensitivity_long,
               aes(x = visitor_change_pct, y = change_pct, colour = severity, group = severity)) +
    geom_hline(yintercept = 0, colour = BLUE$grid, linewidth = 0.9) +
    geom_vline(xintercept = 0, colour = BLUE$grid, linewidth = 0.9) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = sev_colours) +
    scale_x_continuous(breaks = c(-10, 0, 10), labels = label_number(suffix = "%")) +
    scale_y_continuous(labels = label_number(suffix = "%", accuracy = 0.1)) +
    labs(
      title = "How sensitive is the forecast to visitor volumes?",
      subtitle = "FY2027 scenario test | Forecast visits change by −10%, 0%, +10%",
      x = "Change in forecast beach visits", y = "Change in expected people affected",
      colour = "Outcome",
      caption = "Other explanatory variables and safety resources are held constant. This is a sensitivity analysis, not a causal claim."
    ) + theme_slide() +
    theme(legend.position = "bottom")
  save_slide(p8, "08_FY2027_visitor_sensitivity.png")
}

# -------------------------------------------------------------------
# OPTIONAL FIGURE 9 | FY2026 all-candidate diagnostic comparison
# Does not retroactively change the winner selected on FY2025.
# -------------------------------------------------------------------
all_test <- read_output("03b_all_candidate_test_metrics.csv", required = FALSE)
if (!is.null(all_test) &&
    all(c("outcome", "model", "poisson_deviance") %in% names(all_test))) {
  holdout_models <- all_test %>%
    left_join(selected %>% select(outcome, selected_model = model), by = "outcome") %>%
    group_by(outcome) %>%
    mutate(excess_pct = 100 * (poisson_deviance / min(poisson_deviance) - 1)) %>%
    ungroup() %>%
    mutate(
      outcome = factor(outcome, levels = c(sev_levels, "Total")),
      model_label = factor(
        if_else(model %in% names(model_lab), unname(model_lab[model]), model),
        levels = rev(unname(model_lab))
      ),
      selection = if_else(model == selected_model,
                          "FY2025-selected model", "Other candidate"),
      value_label = sprintf("+%.1f%%", excess_pct)
    )

  p9 <- ggplot(holdout_models,
               aes(x = excess_pct, y = model_label, fill = selection)) +
    geom_col(width = 0.57) +
    geom_text(aes(label = value_label), hjust = -0.12,
              colour = BLUE$ink, size = 4) +
    facet_wrap(~ outcome, ncol = 2, scales = "free_x") +
    scale_fill_manual(values = c("FY2025-selected model" = BLUE$dark,
                                 "Other candidate" = BLUE$light)) +
    scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1),
                       expand = expansion(mult = c(0, 0.22), add = c(0, 0.5))) +
    labs(
      title = "Do the selected models hold up on unseen data?",
      subtitle = "FY2026 diagnostic only | Excess Poisson deviance versus best FY2026 candidate",
      x = "Excess test deviance (lower is better)", y = NULL,
      fill = NULL,
      caption = "Dark-blue bars indicate FY2025-selected models. FY2026 results are diagnostic, not used for model selection."
    ) + theme_slide()
  save_slide(p9, "09_FY2026_all_candidate_diagnostic.png")
}

cat("\nDone. Presentation figures: ", normalizePath(figure_dir), "\n", sep = "")
