######################################################################
# 02_modelling.R
# 2026 UNSW ASOC x Taylor Fry Case Competition — Task 1
#
# Input:  temp_output/clean/model_df.rds (from 01_setup_eda_cleaning.R)
#         OR model_df.rds in your working directory
# Output: temp_output/modelling/*.csv, *.png, *.rds
#
# One outcome represents a PERSON involved, not a distinct incident_id.
# Historical fiscal years run Sep-Aug (FY2026 = Sep 2025 to Aug 2026).
# Forecast FY2027 = Sep 2026 to Aug 2027.
#
# Workflow:
#   (1) Train FY2017-FY2024; select model using FY2025
#   (2) Refit winner using FY2017-FY2025; test on FY2026
#   (3) Refit winner on all data; forecast FY2027
#
# TwoStage is now compared for ALL three severity outcomes: Mild,
# CloseCall, and Severe (but NOT Total, because Total / Total = 1).
# Its first stage predicts Total with a negative binomial model;
# its second stage models the probability of the target severity
# among people involved in incidents, using grouped binomial GLM.
#
# FY2026 comparisons for all candidates are exported for diagnostic
# transparency; model selection STILL uses only FY2025 validation.
#
# Optional: install glmmTMB and set RUN_GLMM <- TRUE to add a
# beach-random-intercept negative binomial GLMM as a candidate.
######################################################################

set.seed(2026)

required <- c("dplyr", "tidyr", "readr", "ggplot2", "MASS")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing packages first: install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))")
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
# MASS::glm.nb is explicitly namespaced to avoid masking dplyr::select.

RUN_GLMM <- FALSE  # optional, may take considerably longer
if (RUN_GLMM && !requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("Install glmmTMB first, or set RUN_GLMM <- FALSE")
}

source_candidates <- c("temp_output/clean/model_df.rds", "model_df.rds",
                       "/mnt/data/model_df.rds")
existing <- source_candidates[file.exists(source_candidates)]
if (!length(existing)) stop("Cannot find model_df.rds. Run 01_setup_eda_cleaning.R first.")
model_df <- readRDS(existing[1])

out_dir <- "temp_output/modelling"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_cols <- c("beach_nm", "date", "fy", "cal_m", "season",
                   "visitors", "Mild", "CloseCall", "Severe", "Total",
                   "lg_hrs", "drone_hrs", "log_lg_hrs", "log_drone_hrs",
                   "campaign_idx", "first_aid_cap", "pa_system", "flags",
                   "swell_max_m", "n_rips", "tourist_ratio", "surfer_prop",
                   "rocky", "fishing")
missing_cols <- setdiff(required_cols, names(model_df))
if (length(missing_cols)) {
  stop("model_df.rds is missing these columns: ", paste(missing_cols, collapse = ", "))
}

model_df <- model_df %>%
  mutate(
    date = as.Date(date),
    fy = as.integer(fy),
    beach_nm = factor(beach_nm),
    season = factor(as.character(season),
                    levels = c("Spring", "Summer", "Autumn", "Winter"))
  ) %>%
  arrange(beach_nm, date)

stopifnot(nrow(model_df) == 9600,
          all(model_df$visitors > 0),
          !anyNA(model_df[, required_cols]),
          all(model_df$Total == model_df$Mild + model_df$CloseCall + model_df$Severe),
          all(table(model_df$beach_nm) == 120))

# ---- Predictors -------------------------------------------------------
# Safety levers retained for Task 2; associations are NOT automatically causal.
# `patrolled` omitted: it is exactly lg_hrs > 0.
# `swell_avg_m` omitted: high collinearity with swell_max_m.
# Incident-only variables (response time, swimmer age, etc.) omitted to avoid leakage.

frequent_vars <- c(
  "season", "swell_max_m", "n_rips", "tourist_ratio", "surfer_prop",
  "rocky", "fishing", "log_lg_hrs", "log_drone_hrs",
  "campaign_idx", "first_aid_cap", "pa_system", "flags"
)

# Leaner formula because there are only 220 deaths in the entire dataset.
severe_vars <- c(
  "season", "swell_max_m", "n_rips", "tourist_ratio",
  "log_lg_hrs", "log_drone_hrs", "campaign_idx", "first_aid_cap",
  "pa_system", "flags"
)

# Stage 2 of the frequency-severity challenger: chance of the
# specified severity given that a person was involved in an incident.
# Models Mild, CloseCall and Severe separately as one-vs-rest shares.
share_vars <- c(
  "season", "swell_max_m", "n_rips", "tourist_ratio",
  "log_lg_hrs", "log_drone_hrs", "campaign_idx", "first_aid_cap"
)

count_formula <- function(outcome) {
  vars <- if (outcome == "Severe") severe_vars else frequent_vars
  stats::reformulate(c(vars, "offset(log(visitors))"), response = outcome)
}

share_formula <- function(outcome) {
  if (!outcome %in% c("Mild", "CloseCall", "Severe")) {
    stop("TwoStage share is only meaningful for severity outcomes, not Total")
  }
  stats::reformulate(
    share_vars,
    response = sprintf("cbind(%s, Total - %s)", outcome, outcome)
  )
}

# ---- Model fitting / prediction -------------------------------------
fit_candidate <- function(outcome, kind, data) {
  if (kind == "Poisson") {
    return(stats::glm(count_formula(outcome), family = poisson(link = "log"),
                      data = data, control = glm.control(maxit = 100)))
  }
  if (kind == "NegBin") {
    return(MASS::glm.nb(count_formula(outcome), data = data,
                        control = glm.control(maxit = 100)))
  }
  if (kind == "NegBin_GLMM") {
    f <- update(count_formula(outcome), . ~ . + (1 | beach_nm))
    return(glmmTMB::glmmTMB(f, family = glmmTMB::nbinom2(link = "log"),
                            data = data))
  }
  if (kind == "TwoStage" && outcome %in% c("Mild", "CloseCall", "Severe")) {
    # Stage 1: expected TOTAL person-level incidents per beach-month.
    frequency <- MASS::glm.nb(count_formula("Total"), data = data,
                              control = glm.control(maxit = 100))
    # Stage 2: conditional share of the requested severity. Zero-incident
    # months cannot inform the conditional share, but remain in Stage 1.
    severity_share <- stats::glm(share_formula(outcome),
                                 family = binomial(link = "logit"),
                                 data = data[data$Total > 0, ],
                                 control = glm.control(maxit = 100))
    return(list(frequency = frequency, share = severity_share))
  }
  stop("Unsupported model: ", outcome, "/", kind)
}

predict_candidate <- function(model, outcome, kind, newdata) {
  if (kind == "TwoStage") {
    expected_total <- as.numeric(predict(model$frequency, newdata = newdata,
                                         type = "response"))
    severity_probability <- as.numeric(predict(model$share, newdata = newdata,
                                               type = "response"))
    pred <- expected_total * severity_probability
  } else {
    pred <- as.numeric(predict(model, newdata = newdata, type = "response"))
  }
  if (any(!is.finite(pred)) || any(pred < 0)) {
    stop("Nonfinite/negative predictions for ", outcome, "/", kind)
  }
  pred
}

kinds_for <- function(outcome) {
  options <- c("Poisson", "NegBin")
  if (RUN_GLMM && outcome != "Severe") options <- c(options, "NegBin_GLMM")
  if (outcome %in% c("Mild", "CloseCall", "Severe")) {
    options <- c(options, "TwoStage")
  }
  options
}

# ---- Out-of-time evaluation -----------------------------------------
# Compare means on the SAME scale using Poisson deviance for all candidates.
# This is a predictive comparison, even when a candidate uses a NB distribution.
poisson_deviance <- function(observed, expected) {
  expected <- pmax(expected, 1e-10)
  positive <- observed > 0
  2 * (sum(expected[!positive]) +
         sum(observed[positive] * log(observed[positive] / expected[positive]) -
               observed[positive] + expected[positive]))
}

evaluate <- function(observed, predicted, outcome, kind, split) {
  tibble(
    outcome = outcome, model = kind, split = split,
    poisson_deviance = poisson_deviance(observed, predicted),
    MAE = mean(abs(observed - predicted)),
    RMSE = sqrt(mean((observed - predicted)^2)),
    actual_total = sum(observed), predicted_total = sum(predicted),
    actual_to_expected = ifelse(sum(predicted) > 0,
                                sum(observed) / sum(predicted), NA_real_)
  )
}

train <- model_df %>% filter(fy <= 2024)
valid <- model_df %>% filter(fy == 2025)
train_plus_valid <- model_df %>% filter(fy <= 2025)
test <- model_df %>% filter(fy == 2026)
stopifnot(nrow(train) > 0, nrow(valid) == 960, nrow(test) == 960)

outcomes <- c("Mild", "CloseCall", "Severe", "Total")

cat("\n=== 1. FITTING CANDIDATES; FY2025 VALIDATION ===\n")
validation_rows <- list()
for (y in outcomes) {
  for (kind in kinds_for(y)) {
    cat("Fit", y, kind, "...\n")
    candidate <- tryCatch(fit_candidate(y, kind, train),
                          error = function(e) {
                            warning("Failed: ", y, "/", kind, ": ", conditionMessage(e))
                            NULL
                          })
    if (is.null(candidate)) next
    preds <- tryCatch(predict_candidate(candidate, y, kind, valid),
                      error = function(e) {
                        warning("Predict failed: ", y, "/", kind, ": ",
                                conditionMessage(e))
                        NULL
                      })
    if (is.null(preds)) next
    validation_rows[[length(validation_rows) + 1L]] <-
      evaluate(valid[[y]], preds, y, kind, "FY2025_validation")
  }
}
validation <- dplyr::bind_rows(validation_rows)
if (nrow(validation) == 0) stop("All model candidates failed.")

# Winning model per outcome is chosen ONLY from validation year FY2025.
selected <- validation %>%
  group_by(outcome) %>%
  arrange(poisson_deviance, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(outcome, model, validation_deviance = poisson_deviance)
if (nrow(selected) != length(outcomes)) {
  stop("At least one outcome has no valid model. See warnings above.")
}

cat("\nSelected models (based on FY2025 deviance):\n")
print(selected)
write_csv(validation, file.path(out_dir, "01_validation_comparison.csv"))
write_csv(selected, file.path(out_dir, "02_selected_models.csv"))

cat("\n=== 2. FINAL HOLDOUT TEST: FY2026 ===\n")
test_rows <- list()
all_test_rows <- list()
test_predictions <- test %>%
  select(beach_nm, date, fy, visitors, Mild, CloseCall, Severe, Total)

# Refit EVERY candidate on train + validation to report FY2026 diagnostics.
# Importantly, FY2026 is never used to choose the winner: selected was
# already locked using FY2025 validation deviance.
for (y in outcomes) {
  winner <- selected$model[selected$outcome == y]
  for (kind in kinds_for(y)) {
    cat("Test", y, kind, "...\n")
    fit <- tryCatch(fit_candidate(y, kind, train_plus_valid),
                    error = function(e) {
                      warning("Holdout fit failed: ", y, "/", kind, ": ",
                              conditionMessage(e))
                      NULL
                    })
    if (is.null(fit)) {
      if (kind == winner) stop("Chosen model failed on FY2026: ", y, "/", kind)
      next
    }
    preds <- tryCatch(predict_candidate(fit, y, kind, test),
                      error = function(e) {
                        warning("Holdout predict failed: ", y, "/", kind, ": ",
                                conditionMessage(e))
                        NULL
                      })
    if (is.null(preds)) {
      if (kind == winner) stop("Chosen model could not predict FY2026: ", y, "/", kind)
      next
    }
    diagnostic <- evaluate(test[[y]], preds, y, kind, "FY2026_diagnostic")
    all_test_rows[[length(all_test_rows) + 1L]] <- diagnostic
    if (kind == winner) {
      test_predictions[[paste0("pred_", y)]] <- preds
      test_rows[[length(test_rows) + 1L]] <-
        evaluate(test[[y]], preds, y, kind, "FY2026_final_test")
    }
  }
}

all_test_metrics <- bind_rows(all_test_rows)
test_metrics <- bind_rows(test_rows)
stopifnot(nrow(test_metrics) == length(outcomes))
print(test_metrics)
write_csv(test_metrics, file.path(out_dir, "03_final_test_metrics.csv"))
write_csv(all_test_metrics,
          file.path(out_dir, "03b_all_candidate_test_metrics.csv"))
write_csv(test_predictions, file.path(out_dir, "04_test_predictions.csv"))

# Holdout-year calibration check by month (not used for model selection).
test_by_month <- test_predictions %>%
  mutate(cal_m = as.integer(format(date, "%m"))) %>%
  group_by(cal_m) %>%
  summarise(actual_severe = sum(Severe),
            predicted_severe = sum(pred_Severe), .groups = "drop") %>%
  arrange(match(cal_m, c(9:12, 1:8)))
write_csv(test_by_month, file.path(out_dir, "04b_FY2026_severe_month_calibration.csv"))

# ---- Build 2027 baseline inputs -------------------------------------
# IMPORTANT: 2027's visits/weather are not observed yet.
# Assumptions:
#   - FY2027 monthly VISITS and WEATHER = average of FY2024-FY2026,
#     calculated separately for EACH BEACH AND CALENDAR MONTH.
#   - Safety provision and mostly-static characteristics = FY2026 values
#     for the corresponding beach/month (no extra Task 2 funding).
#   - No new unmodelled shocks or structural changes.
# Replace these assumptions if the team develops a better forecast.

cat("\n=== 3. CONSTRUCTING FY2027 BASELINE ===\n")
climate_cols <- c("swell_max_m", "tourist_ratio", "surfer_prop")
# These extra weather features are retained in the data for scenarios/EDA.
extra_climate <- intersect(c("air_temp_hi", "swell_avg_m", "wind_days",
                             "water_temp", "rain_mm", "turbidity"), names(model_df))
climate_cols <- c(climate_cols, extra_climate)

climatology <- model_df %>%
  filter(fy %in% 2024:2026) %>%
  group_by(beach_nm, cal_m) %>%
  summarise(
    visitors_forecast = mean(visitors),
    across(all_of(climate_cols), mean, .names = "forecast_{.col}"),
    .groups = "drop"
  )

future <- model_df %>%
  filter(fy == 2026) %>%
  select(-any_of(c("Mild", "CloseCall", "Severe", "Total", "n_events")),
         -starts_with("cause_")) %>%
  left_join(climatology, by = c("beach_nm", "cal_m")) %>%
  mutate(visitors = pmax(1, round(visitors_forecast)))

for (v in climate_cols) {
  future[[v]] <- future[[paste0("forecast_", v)]]
}
future <- future %>%
  select(-visitors_forecast, -starts_with("forecast_")) %>%
  mutate(
    # Input dates are the 1st of each month, so adding 1 to calendar year is safe.
    date = as.Date(sprintf("%04d-%02d-01",
                           as.integer(format(date, "%Y")) + 1L, cal_m)),
    fy = 2027L,
    month = format(date, "%m/%Y"),
    log_visitors = log(visitors),
    log_lg_hrs = log1p(lg_hrs),
    log_drone_hrs = log1p(drone_hrs)
  ) %>%
  arrange(beach_nm, date)

stopifnot(nrow(future) == 960,
          all(future$visitors > 0),
          all(!is.na(future[, c(frequent_vars, "visitors")])) )

write_csv(future, file.path(out_dir, "05_FY2027_assumed_inputs.csv"))

# ---- Refit chosen models using ALL observed years --------------------
cat("\n=== 4. REFIT WINNERS; FORECAST FY2027 ===\n")
final_models <- list()
forecast <- future %>%
  select(beach_nm, date, month, cal_m, season, fy, postcode, visitors,
         lg_hrs, drone_hrs, campaign_idx, first_aid_cap, pa_system, flags)

for (i in seq_len(nrow(selected))) {
  y <- selected$outcome[i]
  kind <- selected$model[i]
  fitted_model <- fit_candidate(y, kind, model_df)
  final_models[[y]] <- list(kind = kind, fit = fitted_model)
  forecast[[paste0("pred_", y)]] <- predict_candidate(fitted_model, y, kind, future)
}

# The sum of independently selected severity forecasts is the
# accounting total. It need not equal pred_Total from the independently
# fitted Total model. If different outcomes select two-stage models,
# their one-vs-rest probabilities are estimated separately and may not
# sum exactly to one (even when using the same frequency stage).
# pred_Total is retained as an independent model-check forecast.
forecast <- forecast %>%
  mutate(pred_Total_from_severities = pred_Mild + pred_CloseCall + pred_Severe)

saveRDS(final_models, file.path(out_dir, "06_fitted_models.rds"))
write_csv(forecast, file.path(out_dir, "07_FY2027_beach_month_forecast.csv"))

# ---- Factor estimates (ASSOCIATIONS, not proven treatment effects) ---
coefficient_rows <- list()
append_coefficients <- function(fit, label, kind, interpretation) {
  if (kind == "NegBin_GLMM") {
    est_table <- as.matrix(summary(fit)$coefficients$cond)
  } else {
    est_table <- as.matrix(summary(fit)$coefficients)
  }
  tibble(
    outcome = label, model = kind, effect_type = interpretation,
    term = rownames(est_table),
    estimate_log = as.numeric(est_table[, 1]),
    std_error = as.numeric(est_table[, 2]),
    ratio = exp(estimate_log),
    lower_95 = exp(estimate_log - 1.96 * std_error),
    upper_95 = exp(estimate_log + 1.96 * std_error)
  )
}
for (y in outcomes) {
  item <- final_models[[y]]
  if (item$kind == "TwoStage") {
    coefficient_rows[[length(coefficient_rows) + 1L]] <-
      append_coefficients(item$fit$frequency,
                          paste0(y, " / Total (frequency stage)"),
                          "TwoStage_NB", "Incident-rate ratio")
    coefficient_rows[[length(coefficient_rows) + 1L]] <-
      append_coefficients(item$fit$share,
                          paste0(y, " | incident (share stage)"),
                          "TwoStage_logistic", "Conditional odds ratio")
  } else {
    coefficient_rows[[length(coefficient_rows) + 1L]] <-
      append_coefficients(item$fit, y, item$kind, "Incident-rate ratio")
  }
}
coefficients <- bind_rows(coefficient_rows)
write_csv(coefficients, file.path(out_dir, "08_driver_associations.csv"))

# ---- Tables for the presentation ------------------------------------
statewide <- forecast %>%
  summarise(
    Mild = sum(pred_Mild),
    CloseCall = sum(pred_CloseCall),
    Severe = sum(pred_Severe),
    Total_from_severities = sum(pred_Total_from_severities),
    Total_independent_model = sum(pred_Total),
    Forecast_visits = sum(visitors)
  )
write_csv(statewide, file.path(out_dir, "09_FY2027_statewide_summary.csv"))

by_month <- forecast %>%
  group_by(cal_m) %>%
  summarise(
    visits = sum(visitors), Mild = sum(pred_Mild),
    CloseCall = sum(pred_CloseCall), Severe = sum(pred_Severe),
    .groups = "drop"
  ) %>%
  arrange(match(cal_m, c(9:12, 1:8))) %>%
  mutate(month = factor(month.abb[cal_m], levels = month.abb[c(9:12, 1:8)]))
write_csv(by_month, file.path(out_dir, "10_FY2027_by_month.csv"))

by_beach <- forecast %>%
  group_by(beach_nm, postcode) %>%
  summarise(
    visits = sum(visitors), Mild = sum(pred_Mild),
    CloseCall = sum(pred_CloseCall), Severe = sum(pred_Severe),
    .groups = "drop"
  ) %>%
  mutate(
    Total = Mild + CloseCall + Severe,
    severe_per_million_visits = Severe / visits * 1e6
  ) %>%
  arrange(desc(Severe))
write_csv(by_beach, file.path(out_dir, "11_FY2027_beach_ranking.csv"))

# Monthly predicted outcomes by severity.
plot_month <- by_month %>%
  pivot_longer(c(Mild, CloseCall, Severe), names_to = "severity", values_to = "expected") %>%
  ggplot(aes(x = month, y = expected)) +
  geom_col(fill = "#2471A3") +
  facet_wrap(~ severity, scales = "free_y", ncol = 1) +
  labs(title = "FY2027 expected beach incidents by month",
       subtitle = "Person-level outcomes; Sep 2026 to Aug 2027",
       x = "Calendar month", y = "Expected people") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "12_monthly_forecast.png"), plot_month,
       width = 10, height = 8, dpi = 200)

plot_beaches <- by_beach %>%
  slice_max(Severe, n = 15, with_ties = FALSE) %>%
  mutate(beach_nm = reorder(as.character(beach_nm), Severe)) %>%
  ggplot(aes(x = Severe, y = beach_nm)) +
  geom_col(fill = "#C0392B") +
  labs(title = "FY2027: beaches with the highest expected deaths",
       x = "Expected deaths", y = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "13_top_beaches_expected_deaths.png"), plot_beaches,
       width = 10, height = 6, dpi = 200)

# ---- Sensitivity to unknown visitor volumes -------------------------
# Illustrative +/-10% visitor scenarios, holding other drivers fixed.
# This is NOT a causal estimate of changing safety interventions.
visitor_scenarios <- bind_rows(lapply(c(0.9, 1, 1.1), function(multiplier) {
  scenario <- future
  scenario$visitors <- pmax(1, round(scenario$visitors * multiplier))
  scenario$log_visitors <- log(scenario$visitors)
  out <- tibble(visitor_multiplier = multiplier)
  for (y in c("Mild", "CloseCall", "Severe")) {
    obj <- final_models[[y]]
    out[[y]] <- sum(predict_candidate(obj$fit, y, obj$kind, scenario))
  }
  out
}))
write_csv(visitor_scenarios, file.path(out_dir, "14_visitor_sensitivity.csv"))

cat("\n=== FY2027 STATEWIDE EXPECTED PEOPLE ===\n")
print(statewide)
cat("\n=== TOP 10 BEACHES BY EXPECTED DEATHS ===\n")
print(by_beach %>% select(beach_nm, Severe, severe_per_million_visits) %>% head(10))
cat("\nSaved all modelling results to: ", normalizePath(out_dir), "\n", sep = "")
cat("\nIMPORTANT: Coefficients describe conditional associations, not\n",
    "causal safety-intervention effects. FY2027 inputs are assumed.\n", sep = "")
