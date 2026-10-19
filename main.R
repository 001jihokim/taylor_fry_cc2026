# Main code file for EDA, Modelling, Model diagnostics and Visualisations
# UNSW x TAYLOR FRY T3 2026 Case Competition
# Team Name - Omega Consulting
# Members - Jiho Kim, Becky Thai, Christina Lou, Timothy Bae

# =============================================================================
# 2026 ASOC x Taylor Fry Case Competition
# Part A: data cleaning   |   Part B: Task 1 (forecast + drivers of incidents)
# Fiscal year (FY) convention used throughout: Sep-Aug, labelled by the year it
# ENDS in. FY2017 = Sep 2016-Aug 2017 ... FY2026 = Sep 2025-Aug 2026.
# "Next year" = FY2027 (Sep 2026-Aug 2027).
# =============================================================================

library(MASS)        # glm.nb, stepAIC  (load BEFORE tidyverse so dplyr::select wins)
library(tidyverse)
library(lubridate)
library(broom)
library(scales)
library(gbm)         # gradient boosting (B9)

set.seed(2026)
data_dir <- "."      # <- folder containing the three csv files
out_dir  <- "outputs"; dir.create(out_dir, showWarnings = FALSE)

# =============================================================================
# PART A: DATA CLEANING
# =============================================================================

# ---- A1. Load ---------------------------------------------------------------
beach <- read_csv(file.path(data_dir, "beach_monthly.csv"), show_col_types = FALSE)
inc   <- read_csv(file.path(data_dir, "incidents.csv"),     show_col_types = FALSE)
# pricing.csv has 3 blank trailing rows - drop them (used in Task 2)
pricing <- read_csv(file.path(data_dir, "pricing.csv"), show_col_types = FALSE) %>%
  filter(!is.na(item))

glimpse(beach); glimpse(inc)

# ---- A2. beach_monthly: types, dates, derived time fields -------------------
beach <- beach %>%
  mutate(
    month_date    = as.Date(paste0("01/", month), format = "%d/%m/%Y"),
    cal_month     = month(month_date),
    fy            = year(month_date) + (cal_month >= 9),
    summer        = cal_month %in% c(12, 1, 2),
    marine_health = factor(marine_health, levels = c("L", "M", "H"))
  )

# ---- A3. beach_monthly: integrity checks (should all be 0 / TRUE) -----------
stopifnot(
  !anyDuplicated(beach[, c("beach_nm", "month")]),                  # unique grain
  nrow(beach) == 80 * 120,                                          # balanced panel
  all(beach$swell_max_m >= beach$swell_avg_m),                      # per dictionary
  all(beach$patrolled == as.integer(beach$lg_hrs > 0))              # per dictionary
)
colSums(is.na(beach))                                               # expect all 0
beach %>% summarise(across(where(is.numeric), list(min = min, max = max))) %>%
  pivot_longer(everything(), names_to = c("var", ".value"), names_pattern = "(.*)_(min|max)")

# Static vs time-varying sanity check: static columns should have 1 value/beach
static_cols <- c("postcode", "shore_len_m", "n_rips", "rocky", "fishing", "rockpool")
beach %>% group_by(beach_nm) %>%
  summarise(across(all_of(static_cols), n_distinct), .groups = "drop") %>%
  summarise(across(all_of(static_cols), max))                       # expect all 1

# drone_hrs capacity check: hardware supports <=120 hrs/month per drone
beach %>% filter(drone_hrs > 0) %>% count(drone_hrs > 120)

# ---- A4. incidents: types, dates, factors -----------------------------------
inc <- inc %>%
  mutate(
    dte_occurred = dmy(dte_occurred),
    dte_reported = dmy(dte_reported),
    month        = format(dte_occurred, "%m/%Y"),            # join key to beach_monthly
    fy           = year(dte_occurred) + (month(dte_occurred) >= 9),
    cal_month    = month(dte_occurred),
    report_lag   = as.numeric(dte_reported - dte_occurred),
    inc_type     = factor(inc_type, levels = c("Mild", "CloseCall", "Severe"), ordered = TRUE),
    severe       = as.integer(inc_type == "Severe"),
    # NA cause (child rows only) and NA surfer -> explicit "Unknown" level
    inc_cause    = factor(replace_na(inc_cause, "Unknown")),
    surfer_f     = factor(case_when(surfer == 1 ~ "surfer", surfer == 0 ~ "non_surfer",
                                    TRUE ~ "unknown")),
    swim_level   = factor(swim_level, levels = c("beginner", "intermediate", "experienced")),
    sex          = factor(sex),
    local_ind    = factor(local_ind, levels = c("local", "not_local")),
    drone_present = as.integer(drone_present == "Y"),
    child        = as.integer(age < 18),
    time_band    = cut(inc_hour, c(-1, 5, 11, 17, 23),
                       labels = c("night", "morning", "afternoon", "evening")),
    # response_min includes reporting delay and has extreme values (max ~5,800).
    # Keep original, add a log version winsorised at the 99.5th percentile.
    log_resp     = log(pmin(response_min, quantile(response_min, 0.995)))
  )

# ---- A5. incidents: integrity checks ----------------------------------------
stopifnot(
  !anyDuplicated(inc[, c("incident_id", "person_no")]),
  all(inc$dte_reported >= inc$dte_occurred),
  min(inc$dte_occurred) >= as.Date("2016-09-01"),
  max(inc$dte_occurred) <= as.Date("2026-08-31")
)
nrow(anti_join(inc, beach, by = c("beach_nm", "month")))            # expect 0
colSums(is.na(inc))
summary(inc$response_min); quantile(inc$response_min, c(.5, .9, .99, .999))
# alone = 1 only valid for single-person incidents
inc %>% group_by(incident_id) %>% mutate(n_ppl = n()) %>% ungroup() %>%
  filter(alone == 1, n_ppl > 1) %>% nrow()                           # expect 0
# is cause constant within an incident (except NA children)?
inc %>% filter(inc_cause != "Unknown") %>% group_by(incident_id) %>%
  summarise(k = n_distinct(inc_cause)) %>% count(k)

# ---- A6. Build the beach-month modelling panel ------------------------------
# Events = distinct incident_id; people counted by severity (inc_type is per person)
inc_bm <- inc %>%
  group_by(beach_nm, month) %>%
  summarise(n_events = n_distinct(incident_id),
            n_people = n(),
            n_mild   = sum(inc_type == "Mild"),
            n_close  = sum(inc_type == "CloseCall"),
            n_severe = sum(inc_type == "Severe"),
            .groups = "drop")

# beach-months with no incidents become 0
panel <- beach %>%
  left_join(inc_bm, by = c("beach_nm", "month")) %>%
  mutate(across(c(n_events, n_people, n_mild, n_close, n_severe), ~ replace_na(.x, 0L)))
stopifnot(sum(panel$n_people) == nrow(inc),
          sum(panel$n_events) == n_distinct(inc$incident_id))

# ---- A7. Feature engineering (same function reused for the forecast data) ---
add_feats <- function(d) {
  d %>% mutate(
    patrolled     = as.integer(lg_hrs > 0),
    lg_per_1k     = lg_hrs / (visitors / 1000),     # lifeguard intensity per 1,000 visits
    lg_log        = log1p(lg_per_1k),
    drone_10      = drone_hrs / 10,
    fa_10         = first_aid_cap / 10,
    equip_10      = equip_rating / 10,
    rain_100      = rain_mm / 100,
    turb_10       = turbidity / 10,
    len_log       = log(shore_len_m),
    marine_health = factor(marine_health, levels = c("L", "M", "H")),
    trend         = fy - 2021                        # linear time trend
  )
}
panel <- add_feats(panel)

# =============================================================================
# PART B: TASK 1
# =============================================================================

# ---- B1. Exploratory analysis ----------------------------------------------
# (i) Totals by severity and year
by_fy <- panel %>% group_by(fy) %>%
  summarise(visitors = sum(visitors), events = sum(n_events), people = sum(n_people),
            mild = sum(n_mild), close = sum(n_close), severe = sum(n_severe))
print(by_fy)

p_trend <- by_fy %>%
  pivot_longer(c(mild, close, severe), names_to = "severity", values_to = "n") %>%
  ggplot(aes(fy, n, colour = severity)) + geom_line() + geom_point() +
  facet_wrap(~severity, scales = "free_y", ncol = 1) +
  labs(title = "Incidents (people) by severity and fiscal year", x = "FY", y = "Count") +
  theme_minimal()
ggsave(file.path(out_dir, "trend_by_severity.png"), p_trend, width = 7, height = 7)

# (ii) Seasonality: incidents and incident RATE per million visits by calendar month
seas <- panel %>% group_by(cal_month) %>%
  summarise(people = sum(n_people), severe = sum(n_severe), visits = sum(visitors),
            rate_per_m = people / visits * 1e6, severe_per_m = severe / visits * 1e6)
print(seas)

# (iii) Beach league table
beach_tbl <- panel %>% group_by(beach_nm, postcode) %>%
  summarise(people = sum(n_people), severe = sum(n_severe), visits = sum(visitors),
            rate_per_m = people / visits * 1e6, severe_per_m = severe / visits * 1e6,
            avg_lg_hrs = mean(lg_hrs), .groups = "drop") %>%
  arrange(desc(severe), desc(people))
print(head(beach_tbl, 15))

# (iv) Patrol effect, raw: severe & total rate by patrolled status
panel %>% group_by(patrolled) %>%
  summarise(visits = sum(visitors), people = sum(n_people), severe = sum(n_severe),
            rate_per_m = people / visits * 1e6, severe_per_m = severe / visits * 1e6)

# (v) Incident-level descriptives
inc %>% count(inc_type) %>% mutate(share = n / sum(n))
inc %>% count(inc_cause, inc_type) %>% group_by(inc_cause) %>% mutate(share = n / sum(n))
inc %>% group_by(time_band) %>% summarise(n = n(), severe_rate = mean(severe))
inc %>% group_by(drone_present) %>% summarise(n = n(), severe_rate = mean(severe),
                                              med_resp = median(response_min))
inc %>% group_by(n_lifeguards) %>% summarise(n = n(), severe_rate = mean(severe),
                                             med_resp = median(response_min))

# (vi) Over-dispersion check: variance >> mean => Poisson is too restrictive
with(panel, c(mean = mean(n_events), var = var(n_events)))
with(panel, c(mean = mean(n_severe), var = var(n_severe)))
mean(panel$n_events == 0); mean(panel$n_severe == 0)

# (vii) Correlation screen of candidate drivers (look for spurious proxies,
#       e.g. icecream_visits / n_bins tend to proxy for heat and crowd size)
panel %>% mutate(rate = n_events / visitors * 1e5) %>%
  dplyr::select(rate, air_temp_hi, swell_avg_m, swell_max_m, wind_days, water_temp, rain_mm,
                turbidity, tourist_ratio, surfer_prop, n_rips, shore_len_m, n_bins,
                icecream_visits, lg_per_1k, drone_hrs, campaign_idx, first_aid_cap,
                equip_rating) %>%
  cor(use = "complete.obs") %>% round(2) %>% .[, "rate"]

# ---- B2. Frequency models: Negative Binomial GLM with log(visitors) offset ---
# Exposure = visitors, so coefficients are effects on the incident RATE per visit.
rhs_full <- paste(
  "trend + air_temp_hi + swell_avg_m + swell_max_m + wind_days + water_temp + rain_100 +",
  "turb_10 + tourist_ratio + surfer_prop + marine_health + n_rips + len_log + rocky +",
  "fishing + rockpool + n_bins + icecream_visits +",
  "patrolled + lg_log + drone_10 + campaign_idx + fa_10 + pa_system + flags + equip_10")

# Severe is sparse (~220 people over 9,600 beach-months), so start from a smaller set
rhs_severe <- paste(
  "trend + swell_max_m + water_temp + tourist_ratio + surfer_prop + marine_health +",
  "n_rips + len_log + rocky + patrolled + lg_log + drone_10 + campaign_idx + fa_10 +",
  "pa_system + flags + equip_10")

outcomes <- list(events = "n_events", mild = "n_mild", close = "n_close", severe = "n_severe")

# Fit with NB; fall back to Poisson if theta is not identifiable (common for sparse counts)
fit_count <- function(y, rhs, data) {
  f   <- as.formula(paste(y, "~", rhs, "+ offset(log(visitors))"))
  fit <- tryCatch(suppressWarnings(glm.nb(f, data = data, control = glm.control(maxit = 200))),
                  error = function(e) NULL)
  if (is.null(fit) || fit$theta > 1e4)
    fit <- glm(f, data = data, family = poisson())
  fit
}

# Backward selection by AIC on the full data to choose the final specification
sel <- function(y, rhs) {
  # NB: fit with data = panel (a global object) rather than via fit_count(), because
  # stepAIC re-evaluates the model call and a function-argument called `data` is not found
  f    <- as.formula(paste(y, "~", rhs, "+ offset(log(visitors))"))
  full <- tryCatch(suppressWarnings(glm.nb(f, data = panel, control = glm.control(maxit = 200))),
                   error = function(e) NULL)
  if (is.null(full) || full$theta > 1e4) full <- glm(f, data = panel, family = poisson())
  st   <- suppressWarnings(stepAIC(full, direction = "backward", trace = FALSE))
  paste(attr(terms(st), "term.labels"), collapse = " + ")
}
rhs_final <- list(
  events = sel("n_events", rhs_full),
  mild   = sel("n_mild",   rhs_full),
  close  = sel("n_close",  rhs_full),
  severe = sel("n_severe", rhs_severe))
print(rhs_final)

models <- map2(outcomes, rhs_final, ~ fit_count(.x, .y, panel))
map(models, summary)

# ---- B3. Drivers: rate ratios (per unit and per 1 SD) -----------------------
sds <- panel %>% summarise(across(where(is.numeric), sd)) %>% unlist()

tidy_rr <- function(fit, label) {
  tidy(fit) %>%                                   # Wald CIs (profile CIs are slow for NB)
    filter(term != "(Intercept)") %>%
    mutate(conf.low = exp(estimate - 1.96 * std.error),
           conf.high = exp(estimate + 1.96 * std.error),
           estimate = exp(estimate)) %>%
    mutate(sd_x = ifelse(term %in% names(sds), sds[term], 1),
           rr_1sd = exp(log(estimate) * sd_x),
           sig = p.value < 0.05, model = label)
}
rr_all <- imap_dfr(models, ~ tidy_rr(.x, .y))
write_csv(rr_all, file.path(out_dir, "rate_ratios.csv"))
rr_all %>% filter(model %in% c("events", "severe")) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, rr_1sd, p.value) %>%
  arrange(model, p.value) %>% print(n = 60)

p_rr <- rr_all %>% filter(model %in% c("events", "severe")) %>%
  ggplot(aes(rr_1sd, reorder(term, rr_1sd), colour = sig)) +
  geom_point() + geom_vline(xintercept = 1, linetype = 2) +
  facet_wrap(~model, scales = "free") + scale_x_log10() +
  labs(title = "Rate ratio per 1 SD change in each driver", x = "Rate ratio (log scale)", y = NULL) +
  theme_minimal()
ggsave(file.path(out_dir, "rate_ratios.png"), p_rr, width = 10, height = 6)

# Model diagnostics: dispersion, actual vs fitted by decile
sapply(models, function(m) c(theta = if (inherits(m, "negbin")) m$theta else NA, AIC = AIC(m)))
panel %>% mutate(fit = fitted(models$events)) %>%
  mutate(decile = ntile(fit, 10)) %>% group_by(decile) %>%
  summarise(actual = sum(n_events), fitted = sum(fit))

# ---- B4. Incident-level severity model: P(Severe | incident, person) --------
# Complements the frequency model by using person/incident-level information
sev_fit <- glm(severe ~ I(age / 10) + sex + surfer_f + swim_level + local_ind + alone +
                 inc_cause + time_band + log_resp + n_lifeguards + drone_present,
               data = inc, family = binomial())
sev_or <- tidy(sev_fit, conf.int = FALSE) %>%
  mutate(OR = exp(estimate), lo = exp(estimate - 1.96 * std.error),
         hi = exp(estimate + 1.96 * std.error)) %>% filter(term != "(Intercept)")
print(sev_or %>% arrange(p.value), n = 30)
write_csv(sev_or, file.path(out_dir, "severity_odds_ratios.csv"))
# Caveat: log_resp includes reporting delay (late-found drownings => long response),
# so its association with Severe is partly reverse causation. Read with care.

# ---- B5. Forward validation by fiscal year ---------------------------------
# Train on FY < k, predict FY k (k = 2022..2026). Specification held fixed.
pois_dev <- function(y, mu) 2 * sum(ifelse(y > 0, y * log(y / mu), 0) - (y - mu))

validate <- function(y, rhs, folds = 2022:2026) {
  map_dfr(folds, function(k) {
    tr <- panel %>% filter(fy < k); te <- panel %>% filter(fy == k)
    fit <- fit_count(y, rhs, tr)
    mu  <- predict(fit, newdata = te, type = "response")
    base <- tr %>% summarise(r = sum(.data[[y]]) / sum(visitors)) %>% pull(r) * te$visitors
    tibble(outcome = y, test_fy = k,
           actual = sum(te[[y]]), predicted = sum(mu), pct_err = predicted / actual - 1,
           mae_bm = mean(abs(te[[y]] - mu)), mae_bm_naive = mean(abs(te[[y]] - base)),
           dev_model = pois_dev(te[[y]], mu), dev_naive = pois_dev(te[[y]], pmax(base, 1e-6)))
  })
}
val <- imap_dfr(outcomes, ~ validate(.x, rhs_final[[.y]]))
print(val, n = 40)
write_csv(val, file.path(out_dir, "forward_validation.csv"))
# "naive" = overall incident rate per visit x visitors (no covariates), as a benchmark

# ---- B6. Build FY2027 (Sep 2026 - Aug 2027) covariates ----------------------
# Assumptions: weather/mix = same-month average of last 3 FYs; visitors = last-3-FY
# same-month average grown at recent annual growth; safety resourcing held at FY2026
# levels (business-as-usual, no new $2M); static beach features unchanged.
last1 <- panel %>% filter(fy == max(fy))
last3 <- panel %>% filter(fy >= max(fy) - 2)

g <- by_fy %>% filter(fy %in% c(max(fy), max(fy) - 2)) %>% arrange(fy) %>%
  { (.$visitors[2] / .$visitors[1])^(1 / 2) - 1 }        # annual visitor growth, last 2 yrs
message("Assumed annual visitor growth: ", percent(g, 0.1))

clim_vars <- c("air_temp_hi", "swell_avg_m", "swell_max_m", "wind_days", "water_temp",
               "rain_mm", "turbidity", "tourist_ratio", "icecream_visits", "surfer_prop")

clim <- last3 %>% group_by(beach_nm, cal_month) %>%
  summarise(across(all_of(clim_vars), mean), visitors = mean(visitors) * (1 + g)^2,
            .groups = "drop")
bau  <- last1 %>% dplyr::select(beach_nm, cal_month, lg_hrs, drone_hrs, campaign_idx,
                                first_aid_cap, marine_health)
stat <- last1 %>% group_by(beach_nm) %>% slice_tail(n = 1) %>% ungroup() %>%
  dplyr::select(beach_nm, postcode, shore_len_m, n_rips, pa_system, flags, equip_rating,
                rocky, fishing, n_bins, rockpool)

future <- expand_grid(beach_nm = unique(beach$beach_nm),
                      month_date = seq(as.Date("2026-09-01"), by = "month", length.out = 12)) %>%
  mutate(cal_month = month(month_date), fy = 2027L,
         month = format(month_date, "%m/%Y")) %>%
  left_join(clim, by = c("beach_nm", "cal_month")) %>%
  left_join(bau,  by = c("beach_nm", "cal_month")) %>%
  left_join(stat, by = "beach_nm") %>%
  add_feats()
stopifnot(nrow(future) == 960, !anyNA(future[, c("visitors", "lg_hrs", "air_temp_hi")]))

# ---- B7. Forecast with parameter + process uncertainty ----------------------
# Simulate: draw coefficients ~ MVN(beta, vcov), then counts ~ NB(mu, theta)
simulate_counts <- function(fit, nd, nsim = 2000) {
  X    <- model.matrix(reformulate(attr(terms(fit), "term.labels")), data = nd)
  beta <- mvrnorm(nsim, coef(fit), vcov(fit))
  mu   <- exp(X %*% t(beta) + log(nd$visitors))               # rows = beach-months
  if (inherits(fit, "negbin")) {
    matrix(rnbinom(length(mu), size = fit$theta, mu = mu), nrow(mu))
  } else {
    matrix(rpois(length(mu), mu), nrow(mu))
  }
}

forecast_all <- function(nd, nsim = 2000) {
  imap(models, ~ list(point = predict(.x, newdata = nd, type = "response"),
                      sims  = simulate_counts(.x, nd, nsim)))
}
fc <- forecast_all(future)

# (a) Annual totals by severity with 90% prediction intervals
fc_total <- imap_dfr(fc, function(f, nm) {
  tot <- colSums(f$sims)
  tibble(outcome = nm, expected = sum(f$point),
         lo90 = quantile(tot, 0.05), hi90 = quantile(tot, 0.95))
})
print(fc_total)
write_csv(fc_total, file.path(out_dir, "forecast_FY2027_totals.csv"))
# Compare with history
by_fy %>% filter(fy >= 2024)

# (b) WHERE: by beach (expected counts + rate)
future_out <- future %>%
  mutate(exp_events = fc$events$point, exp_mild = fc$mild$point,
         exp_close = fc$close$point, exp_severe = fc$severe$point)

fc_beach <- future_out %>% group_by(beach_nm, postcode) %>%
  summarise(visits = sum(visitors), events = sum(exp_events), mild = sum(exp_mild),
            close = sum(exp_close), severe = sum(exp_severe),
            severe_per_m = severe / visits * 1e6, .groups = "drop") %>%
  arrange(desc(severe))
print(head(fc_beach, 15))
write_csv(fc_beach, file.path(out_dir, "forecast_FY2027_by_beach.csv"))

# (c) WHEN: by month
fc_month <- future_out %>% group_by(month_date) %>%
  summarise(events = sum(exp_events), severe = sum(exp_severe), visits = sum(visitors))
print(fc_month)
write_csv(fc_month, file.path(out_dir, "forecast_FY2027_by_month.csv"))

# (d) WHERE x WHEN heatmap of expected severe incidents (top 20 beaches)
top20 <- fc_beach$beach_nm[1:20]
p_heat <- future_out %>% filter(beach_nm %in% top20) %>%
  mutate(beach_nm = factor(beach_nm, levels = rev(top20))) %>%
  ggplot(aes(format(month_date, "%b"), beach_nm, fill = exp_severe)) + geom_tile() +
  scale_x_discrete(limits = format(seq(as.Date("2026-09-01"), by = "month", length.out = 12), "%b")) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(title = "Expected severe incidents, FY2027 - top 20 beaches", x = NULL, y = NULL,
       fill = "Expected") + theme_minimal()
ggsave(file.path(out_dir, "heatmap_severe.png"), p_heat, width = 9, height = 6)

# ---- B8. Scenario / sensitivity analysis -----------------------------------
scen <- list(
  base            = future,
  visitors_plus10 = future %>% mutate(visitors = visitors * 1.10),
  visitors_minus10 = future %>% mutate(visitors = visitors * 0.90),
  rougher_surf    = future %>% mutate(swell_avg_m = swell_avg_m * 1.15,
                                      swell_max_m = swell_max_m * 1.15),
  hotter          = future %>% mutate(air_temp_hi = air_temp_hi + 1.5, water_temp = water_temp + 1),
  lifeguards_minus10 = future %>% mutate(lg_hrs = lg_hrs * 0.90) %>% add_feats()
)
scen_out <- imap_dfr(scen, function(nd, nm) {
  nd <- add_feats(nd %>% dplyr::select(-any_of(c("patrolled", "lg_per_1k", "lg_log", "drone_10",
                                                 "fa_10", "equip_10", "rain_100", "turb_10",
                                                 "len_log", "trend"))))
  tibble(scenario = nm,
         events = sum(predict(models$events, nd, type = "response")),
         severe = sum(predict(models$severe, nd, type = "response")))
})
print(scen_out)
write_csv(scen_out, file.path(out_dir, "scenarios.csv"))

# =============================================================================
# B9. GRADIENT BOOSTING (GBM) - non-linear benchmark and ensemble
# Poisson-loss GBM with log(visitors) offset (same exposure logic as the GLM).
# Uses: (i) check the GLM isn't missing non-linearities/interactions,
#       (ii) alternative forecast + simple GLM/GBM average, (iii) importance & PDPs.
# NB: predict.gbm() ignores the offset, so log(visitors) is added back manually.
# =============================================================================
gbm_vars <- c("trend", "air_temp_hi", "swell_avg_m", "swell_max_m", "wind_days", "water_temp",
              "rain_100", "turb_10", "tourist_ratio", "surfer_prop", "marine_health", "n_rips",
              "len_log", "rocky", "fishing", "rockpool", "n_bins", "patrolled", "lg_log",
              "drone_10", "campaign_idx", "fa_10", "pa_system", "flags", "equip_10")
# Severe is sparse -> smaller, more regularised model
gbm_vars_severe <- setdiff(gbm_vars, c("fishing", "rockpool", "n_bins", "wind_days", "rain_100",
                                       "turb_10", "air_temp_hi", "swell_avg_m"))

gbm_params <- list(
  events = list(vars = gbm_vars, depth = 3, minobs = 20, shrink = 0.05, ntrees = 1500),
  mild   = list(vars = gbm_vars, depth = 3, minobs = 20, shrink = 0.05, ntrees = 1500),
  close  = list(vars = gbm_vars, depth = 3, minobs = 20, shrink = 0.03, ntrees = 1500),
  severe = list(vars = gbm_vars_severe, depth = 2, minobs = 30, shrink = 0.01, ntrees = 1500))

fit_gbm <- function(y, data, prm) {
  f <- as.formula(paste(y, "~", paste(prm$vars, collapse = " + "), "+ offset(log(visitors))"))
  m <- gbm(f, data = as.data.frame(data), distribution = "poisson",
           n.trees = prm$ntrees, interaction.depth = prm$depth, shrinkage = prm$shrink,
           n.minobsinnode = prm$minobs, bag.fraction = 0.7, cv.folds = 5,
           keep.data = FALSE, verbose = FALSE, n.cores = 1)
  m$best <- gbm.perf(m, method = "cv", plot.it = FALSE)   # CV-chosen number of trees
  m
}
pred_gbm <- function(m, nd) {
  exp(predict(m, as.data.frame(nd), n.trees = m$best, type = "link") + log(nd$visitors))
}

# ---- B9.1 Forward validation (same folds as the GLM) -------------------------
validate_gbm <- function(y, prm, folds = 2022:2026) {
  map_dfr(folds, function(k) {
    tr <- panel %>% filter(fy < k); te <- panel %>% filter(fy == k)
    m  <- fit_gbm(y, tr, prm); mu <- pred_gbm(m, te)
    tibble(outcome = y, test_fy = k, pred_gbm = sum(mu),
           mae_bm_gbm = mean(abs(te[[y]] - mu)), dev_gbm = pois_dev(te[[y]], pmax(mu, 1e-6)))
  })
}
y_map <- unlist(outcomes)                                     # name -> column
val_gbm <- imap_dfr(gbm_params, ~ validate_gbm(y_map[[.y]], .x) %>% mutate(outcome = .y))

# Side-by-side with the GLM from B5 (lower deviance / MAE = better)
val_cmp <- val %>%
  mutate(outcome = names(outcomes)[match(outcome, unlist(outcomes))]) %>%
  left_join(val_gbm, by = c("outcome", "test_fy")) %>%
  mutate(pct_err_gbm = pred_gbm / actual - 1) %>%
  dplyr::select(outcome, test_fy, actual, predicted_glm = predicted, pred_gbm,
                pct_err_glm = pct_err, pct_err_gbm,
                mae_glm = mae_bm, mae_gbm = mae_bm_gbm, dev_glm = dev_model, dev_gbm, dev_naive)
print(val_cmp, n = 40)
val_cmp %>% group_by(outcome) %>%
  summarise(across(c(mae_glm, mae_gbm, dev_glm, dev_gbm, dev_naive), mean))
write_csv(val_cmp, file.path(out_dir, "validation_glm_vs_gbm.csv"))

# ---- B9.2 Final GBMs on all data, importance and partial dependence ----------
gbm_models <- imap(gbm_params, ~ fit_gbm(y_map[[.y]], panel, .x))

imp <- imap_dfr(gbm_models, ~ summary(.x, n.trees = .x$best, plotit = FALSE) %>%
                  as_tibble() %>% mutate(model = .y))
write_csv(imp, file.path(out_dir, "gbm_importance.csv"))
p_imp <- imp %>% filter(model %in% c("events", "severe")) %>%
  group_by(model) %>% slice_max(rel.inf, n = 12) %>% ungroup() %>%
  ggplot(aes(rel.inf, reorder(var, rel.inf))) + geom_col() + facet_wrap(~model, scales = "free") +
  labs(title = "GBM relative influence", x = "Relative influence (%)", y = NULL) + theme_minimal()
ggsave(file.path(out_dir, "gbm_importance.png"), p_imp, width = 10, height = 5)

# Partial dependence for the safety levers (shape of the effect, not just a slope)
pdp <- function(m, var) {
  g <- plot(m, i.var = var, n.trees = m$best, return.grid = TRUE, continuous.resolution = 40)
  tibble(var = var, x = as.numeric(g[[1]]), y = g$y)
}
safety_vars <- c("lg_log", "drone_10", "campaign_idx", "fa_10", "equip_10")
pd <- map_dfr(safety_vars, ~ pdp(gbm_models$events, .x) %>% mutate(model = "events")) %>%
  bind_rows(map_dfr(intersect(safety_vars, gbm_vars_severe),
                    ~ pdp(gbm_models$severe, .x) %>% mutate(model = "severe")))
p_pd <- ggplot(pd, aes(x, y)) + geom_line() + facet_grid(model ~ var, scales = "free") +
  labs(title = "GBM partial dependence - safety levers (log-rate scale)", x = NULL, y = NULL) +
  theme_minimal()
ggsave(file.path(out_dir, "gbm_pdp_safety.png"), p_pd, width = 12, height = 5)

# ---- B9.3 FY2027 forecast: GLM vs GBM vs 50/50 ensemble ----------------------
fc_cmp <- imap_dfr(gbm_models, function(m, nm) {
  glm_pt <- sum(fc[[nm]]$point); gbm_pt <- sum(pred_gbm(m, future))
  tibble(outcome = nm, glm = glm_pt, gbm = gbm_pt, ensemble = (glm_pt + gbm_pt) / 2)
})
print(fc_cmp)
write_csv(fc_cmp, file.path(out_dir, "forecast_FY2027_glm_gbm.csv"))

# Beach-level agreement between the two approaches (severe)
agree <- future %>% mutate(glm = fc$severe$point, gbm = pred_gbm(gbm_models$severe, future)) %>%
  group_by(beach_nm) %>% summarise(glm = sum(glm), gbm = sum(gbm), .groups = "drop")
cor(agree$glm, agree$gbm, method = "spearman")             # rank agreement across beaches
# Interpretation guide: if GBM validates materially better, report the ensemble point
# forecast and keep the GLM prediction intervals; if similar, stay with the GLM for
# simplicity and explainability. Prediction intervals still come from the GLM simulation.

# =============================================================================
# KEY ASSUMPTIONS / LIMITATIONS TO STATE ON THE SLIDES
# - FY2027 covariates: climatology (last 3 FY), visitors grown at recent trend,
#   safety resourcing held at FY2026 (no new funding).
# - Associations, not causal effects: lifeguards/drones are placed at riskier beaches
#   (reverse causality), so raw patrol coefficients may be biased towards zero/positive.
# - Variable selection used the full history; validation reuses that specification
#   (mild optimism). Severe counts are sparse (220 deaths) -> wide intervals.
# - response_min includes reporting delay; n_lifeguards / drone_present are measured
#   at incident time (post-treatment) so they are used only in the severity model.
# - GBM tuning is light (fixed depth/shrinkage, CV tree count with random folds, which
#   ignores beach/time structure so the CV-chosen tree count is a little optimistic).
# - Spurious proxies (icecream_visits, n_bins) can survive selection - sanity-check
#   signs/sizes before interpreting them.
# =============================================================================