# =====================================================================
# 01_setup_eda_cleaning.R
# UNSW x Taylor Fry T3 2026 Case Competition - Omega Consulting
# Members - Jiho Kim, Becky Thai, Christina Lou, Timothy Bae
#
# Part 1 of 2: setup, import, validation, cleaning, then EDA.
# Reads:  data/*.csv (never modified)
# Writes: temp_output/clean/*.rds  (model-ready tables, read by 02_modelling.R)
#         temp_output/eda/         (figures f01-f09, tables t01-t10)
# temp_output/ is git-ignored; run this script, then 02_modelling.R, to rebuild everything.
# Open taylor_fry_cc2026.Rproj first so the relative paths below resolve from the project root.
# =====================================================================

# =====================================================================
# PART A. SETUP
# =====================================================================
# ---- Paths ----------------------------------------------------------
data_dir  <- "data"                 # raw csvs: beach_monthly.csv, incidents.csv, pricing.csv (read-only)
out_root  <- "temp_output"          # EVERYTHING the code generates goes here (git-ignored)
clean_dir <- file.path(out_root, "clean")    # cleaned .rds tables (Part B)
eda_dir   <- file.path(out_root, "eda")      # EDA figures + tables (Part C)
model_dir <- file.path(out_root, "models")   # model outputs + forecasts (02_modelling.R)
for (d in c(clean_dir, eda_dir, model_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Packages -------------------------------------------------------
# Wrangling / plotting
pkgs_core  <- c("dplyr", "tidyr", "readr", "lubridate", "ggplot2", "scales", "patchwork", "corrplot")
# Modelling: MASS (negative binomial), statmod + tweedie (Tweedie GLM and its power p),
#            gbm + xgboost (Poisson GBM with exposure), ranger (random forest)
pkgs_model <- c("MASS", "statmod", "tweedie", "gbm", "xgboost", "ranger")

missing_pkgs <- setdiff(c(pkgs_core, pkgs_model), rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)
# Only the core packages are attached. Modelling packages are called with pkg:: or attached in
# 02_modelling.R, because MASS masks dplyr::select().
invisible(lapply(pkgs_core, library, character.only = TRUE))

set.seed(2026)
fy_of <- function(d) year(d) + (month(d) >= 9)   # fiscal year Sep-Aug; FY2026 = Sep-2025 to Aug-2026


# =====================================================================
# PART B. IMPORT, VALIDATION AND CLEANING
# =====================================================================
# ---------------------------------------------------------------------
# 1. IMPORT (explicit column types so nothing is silently mis-parsed)
# ---------------------------------------------------------------------
beach_raw <- read_csv(file.path(data_dir, "beach_monthly.csv"), show_col_types = FALSE,
                      col_types = cols(beach_nm = "c", month = "c", marine_health = "c", .default = "d"))
inc_raw   <- read_csv(file.path(data_dir, "incidents.csv"), show_col_types = FALSE,
                      col_types = cols(dte_reported = "c", dte_occurred = "c", inc_type = "c", inc_cause = "c",
                                       sex = "c", swim_level = "c", local_ind = "c", beach_nm = "c",
                                       drone_present = "c", .default = "d"))
price_raw <- read_csv(file.path(data_dir, "pricing.csv"), show_col_types = FALSE)

# ---------------------------------------------------------------------
# 2. VALIDATION - check the data against the data dictionary
#    Each check prints PASS / FLAG; FLAGs are handled in section 3.
# ---------------------------------------------------------------------
chk <- function(label, ok, detail = "") cat(sprintf("[%s] %-62s %s\n", ifelse(ok, "PASS", "FLAG"), label, detail))
cat("\n== Data validation ==\n")
chk("beach_monthly: 80 beaches x 120 months, no duplicate keys",
    nrow(beach_raw) == 9600 && n_distinct(beach_raw$beach_nm) == 80 && !anyDuplicated(beach_raw[c("beach_nm", "month")]))
chk("incidents: (incident_id, person_no) unique", !anyDuplicated(inc_raw[c("incident_id", "person_no")]))
chk("beach_monthly: no missing values", sum(is.na(beach_raw)) == 0)
na_inc <- colSums(is.na(inc_raw)); na_inc <- na_inc[na_inc > 0]
chk("incidents: no missing values", length(na_inc) == 0, paste(names(na_inc), na_inc, sep = "=", collapse = ", "))
chk("every incident beach exists in beach_monthly", all(inc_raw$beach_nm %in% beach_raw$beach_nm))
chk("patrolled == (lg_hrs > 0)", all(beach_raw$patrolled == as.integer(beach_raw$lg_hrs > 0)))
chk("swell_max_m >= swell_avg_m", all(beach_raw$swell_max_m >= beach_raw$swell_avg_m))
chk("proportions within [0,1]", all(between(beach_raw$tourist_ratio, 0, 1), between(beach_raw$surfer_prop, 0, 1)))
chk("campaign_idx constant within postcode-month",
    beach_raw %>% group_by(postcode, month) %>% summarise(k = n_distinct(campaign_idx), .groups = "drop") %>% with(all(k == 1)))
chk("first_aid_cap constant within postcode-month",
    beach_raw %>% group_by(postcode, month) %>% summarise(k = n_distinct(first_aid_cap), .groups = "drop") %>% with(all(k == 1)))
static_vars <- c("postcode", "shore_len_m", "n_rips", "rocky", "fishing", "rockpool")
chk("static beach attributes never change",
    beach_raw %>% group_by(beach_nm) %>% summarise(across(all_of(static_vars), n_distinct)) %>%
      dplyr::select(-beach_nm) %>% as.matrix() %>% { all(. == 1) })
n_alone_bad <- inc_raw %>% add_count(incident_id, name = "k") %>% filter(alone == 1, k > 1) %>% nrow()
chk("alone = 1 only for single-person incidents", n_alone_bad == 0)
chk("dte_reported >= dte_occurred", all(dmy(inc_raw$dte_reported) >= dmy(inc_raw$dte_occurred)))
chk("pricing: 10 items, no blank rows", sum(!is.na(price_raw$item)) == 10 && nrow(price_raw) == 10,
    sprintf("(%d rows read, %d blank)", nrow(price_raw), sum(is.na(price_raw$item))))
chk("response_min plausible (< 6 hours)", max(inc_raw$response_min) < 360,
    sprintf("(max = %.0f min)", max(inc_raw$response_min)))

# Values sitting exactly on a bound suggest censoring/capping in the source data
bounds <- tibble(var = c("rain_mm", "turbidity", "turbidity", "swell_max_m", "drone_hrs", "first_aid_cap", "equip_rating"),
                 val = c(700, 0, 100, 5, 150, 100, 100))
bounds$n_at_bound <- mapply(function(v, x) sum(beach_raw[[v]] == x), bounds$var, bounds$val)
cat("\nBeach-months sitting exactly on a bound (possible capping):\n"); print(bounds)

# ---------------------------------------------------------------------
# 3. CLEANING
# ---------------------------------------------------------------------
## 3a. Pricing - drop trailing blank rows, tidy text
pricing <- price_raw %>% filter(!is.na(item)) %>%
  mutate(item = trimws(item), effect_per_unit = gsub("\\s+", " ", effect_per_unit))

## 3b. Beach-month table
beach <- beach_raw %>%
  mutate(date      = dmy(paste0("01/", month)),
         cal_m     = month(date),
         month_lab = factor(month.abb[cal_m], levels = month.abb[c(9:12, 1:8)]),   # Sep..Aug ordering
         fy        = fy_of(date),
         season    = case_when(cal_m %in% c(12, 1, 2) ~ "Summer", cal_m %in% 3:5 ~ "Autumn",
                               cal_m %in% 6:8 ~ "Winter", TRUE ~ "Spring"),
         season    = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter")),
         marine_health = factor(marine_health, levels = c("L", "M", "H"), ordered = TRUE),
         postcode  = factor(postcode),
         # derived exposure / intensity measures
         log_visitors   = log(visitors),
         lg_hrs_per_10k = 1e4 * lg_hrs / visitors,         # patrol intensity relative to crowd size
         drone_any      = as.integer(drone_hrs > 0),
         campaign_any   = as.integer(campaign_idx > 0),
         # flag values on a bound (kept as-is; flag lets models/sensitivity tests treat them separately)
         cap_flag = as.integer(rain_mm >= 700 | turbidity %in% c(0, 100) | swell_max_m >= 5),
         # FY2020-21: visitors and incidents both dip then recover. Flag so models can test a level shift.
         period_dip = as.integer(fy %in% c(2020, 2021)))

## 3c. Incident (person-level) table
inc <- inc_raw %>%
  mutate(d_occ = dmy(dte_occurred), d_rep = dmy(dte_reported),
         month = format(d_occ, "%m/%Y"), cal_m = month(d_occ), fy = fy_of(d_occ),
         report_lag_days = as.integer(d_rep - d_occ),
         # response_min includes reporting delay (dictionary). Incidents reported on a later day carry
         # ~1440 min per day of lag. Keep the original and add an on-scene measure.
         response_onscene = response_min - 1440 * report_lag_days,
         late_report = as.integer(report_lag_days > 0),
         inc_type   = factor(inc_type, levels = c("Mild", "CloseCall", "Severe"), ordered = TRUE),
         swim_level = factor(swim_level, levels = c("beginner", "intermediate", "experienced"), ordered = TRUE),
         child      = as.integer(age < 18),
         age_band   = cut(age, c(0, 12, 17, 24, 39, 59, 74, Inf),
                          labels = c("0-12", "13-17", "18-24", "25-39", "40-59", "60-74", "75+")),
         surfer     = factor(ifelse(is.na(surfer), "Unknown", ifelse(surfer == 1, "Yes", "No")),
                             levels = c("No", "Yes", "Unknown")),
         local      = as.integer(local_ind == "local"),
         drone_present = as.integer(drone_present == "Y")) %>%
  # inc_cause is incident-level but missing for some children. Where another person in the same
  # incident has a recorded cause, inherit it; otherwise mark as Unspecified.
  group_by(incident_id) %>%
  mutate(cause_filled = first(na.omit(inc_cause)),
         inc_cause = coalesce(inc_cause, cause_filled, "Unspecified"),
         n_people  = n()) %>%
  ungroup() %>% dplyr::select(-cause_filled)

cat(sprintf("\ninc_cause: %d missing child rows recovered from same-incident adults, %d left 'Unspecified'\n",
            sum(is.na(inc_raw$inc_cause)) - sum(inc$inc_cause == "Unspecified"), sum(inc$inc_cause == "Unspecified")))
cat(sprintf("response_min: %d late-reported people; on-scene response median %.1f min (raw max %.0f -> cleaned max %.0f)\n",
            sum(inc$late_report), median(inc$response_onscene), max(inc$response_min), max(inc$response_onscene)))

## 3d. Beach-month modelling table: counts of PEOPLE by severity (+ events and cause), zero-filled
#  People (not incidents) are the unit for severity because inc_type is person-level and the client
#  objective counts severe outcomes (deaths). Incident (event) counts are kept for reference.
people_counts <- inc %>% count(beach_nm, month, inc_type) %>%
  pivot_wider(names_from = inc_type, values_from = n, values_fill = 0)
event_counts  <- inc %>% distinct(beach_nm, month, incident_id) %>% count(beach_nm, month, name = "n_events")
cause_counts  <- inc %>% count(beach_nm, month, inc_cause) %>%
  pivot_wider(names_from = inc_cause, values_from = n, values_fill = 0, names_prefix = "cause_")

model_df <- beach %>%
  left_join(people_counts, by = c("beach_nm", "month")) %>%
  left_join(event_counts,  by = c("beach_nm", "month")) %>%
  left_join(cause_counts,  by = c("beach_nm", "month")) %>%
  mutate(across(c(Mild, CloseCall, Severe, n_events, starts_with("cause_")), ~ replace_na(.x, 0L)),
         Total = Mild + CloseCall + Severe)

stopifnot(nrow(model_df) == 9600,
          sum(model_df$Total) == nrow(inc),
          sum(model_df$Severe) == sum(inc$inc_type == "Severe"),
          sum(model_df$n_events) == n_distinct(inc$incident_id))

## 3e. Static beach attribute table (one row per beach) - used for 'where' analysis and Task 2 costing
beach_static <- beach %>% group_by(beach_nm, postcode) %>%
  summarise(shore_len_m = first(shore_len_m), n_rips = first(n_rips), rocky = first(rocky),
            fishing = first(fishing), rockpool = first(rockpool),
            pa_system_now = last(pa_system, order_by = date), flags_now = last(flags, order_by = date),
            months_patrolled = sum(patrolled), ever_patrolled = as.integer(any(patrolled == 1)),
            .groups = "drop")

## 3f. Cross-validation fold IDs (defined once here so every model uses identical splits)
#  Plain LOOCV (leave one beach-month out) is not used: rows are correlated within a beach and within
#  a season, so it leaks information and overstates accuracy. Instead:
#    - fold_year : leave-one-fiscal-year-out  -> tests temporal generalisation (what Task 1 needs)
#    - fold_beach: leave-one-beach-out         -> tests generalisation to beaches' own risk levels
#  The forward (expanding window) check - train FY<t, test FY=t for t = 2024-26 - is run in 02_modelling.R.
model_df <- model_df %>%
  mutate(fold_year  = as.integer(factor(fy)),
         fold_beach = as.integer(factor(beach_nm)))

# ---------------------------------------------------------------------
# 4. SAVE CLEANED DATA (all generated files go to temp_output/, never to data/)
# ---------------------------------------------------------------------
saveRDS(model_df,     file.path(clean_dir, "model_df.rds"))        # beach-month: counts + drivers + folds
saveRDS(inc,          file.path(clean_dir, "incidents_clean.rds"))  # person-level, cleaned
saveRDS(beach_static, file.path(clean_dir, "beach_static.rds"))     # one row per beach
saveRDS(pricing,      file.path(clean_dir, "pricing.rds"))
cat("\nSaved cleaned data to", clean_dir, "\n")

# =====================================================================
# PART C. EXPLORATORY DATA ANALYSIS (uses the cleaned objects above)
# =====================================================================
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                                plot.title = element_text(face = "bold")))
sev_cols  <- c(Mild = "#7FB3D5", CloseCall = "#F5B041", Severe = "#C0392B")
save_plot <- function(p, name, w = 8, h = 4.5) ggsave(file.path(eda_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
save_tab  <- function(d, name) write_csv(d, file.path(eda_dir, paste0(name, ".csv")))

rate_tab <- function(d, ...) d %>% group_by(...) %>%
  summarise(beach_months = n(), visitors_m = sum(visitors) / 1e6,
            Mild = sum(Mild), CloseCall = sum(CloseCall), Severe = sum(Severe),
            mild_per_m = Mild / visitors_m, close_per_m = CloseCall / visitors_m, severe_per_m = Severe / visitors_m,
            .groups = "drop")

## 1. Headline volumes and trend by fiscal year
by_fy <- rate_tab(model_df, fy); cat("\n== By fiscal year ==\n"); print(by_fy, width = Inf); save_tab(by_fy, "t01_by_fiscal_year")

p_fy_counts <- by_fy %>% pivot_longer(c(Mild, CloseCall, Severe), names_to = "sev", values_to = "n") %>%
  mutate(sev = factor(sev, names(sev_cols))) %>%
  ggplot(aes(fy, n, colour = sev)) + geom_line(linewidth = 1) + geom_point() +
  facet_wrap(~ sev, scales = "free_y") + scale_colour_manual(values = sev_cols, guide = "none") +
  scale_x_continuous(breaks = 2017:2026) + labs(title = "People involved in incidents by fiscal year", x = "Fiscal year (Sep-Aug)", y = "People")
p_fy_rate <- by_fy %>% ggplot(aes(fy, severe_per_m)) + geom_col(fill = sev_cols["Severe"]) +
  geom_line(aes(y = visitors_m / 30), colour = "grey30", linewidth = 0.8) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 30, name = "Visits (millions)")) + scale_x_continuous(breaks = 2017:2026) +
  labs(title = "Deaths per million visits (bars) vs total visits (line)", x = "Fiscal year", y = "Deaths per million visits")
save_plot(p_fy_counts / p_fy_rate, "f01_trend_by_fy", h = 7)

## 2. Seasonality - when do incidents happen?
by_month <- rate_tab(model_df, month_lab); save_tab(by_month, "t02_by_calendar_month")
p_month <- by_month %>% pivot_longer(c(mild_per_m, close_per_m, severe_per_m), names_to = "sev", values_to = "rate") %>%
  mutate(sev = factor(sev, c("mild_per_m", "close_per_m", "severe_per_m"), names(sev_cols))) %>%
  ggplot(aes(month_lab, rate, fill = sev)) + geom_col() + facet_wrap(~ sev, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = sev_cols, guide = "none") +
  labs(title = "Incident rate per million visits by month", x = NULL, y = "People per million visits")
p_month_vol <- by_month %>% ggplot(aes(month_lab, visitors_m)) + geom_col(fill = "grey60") +
  labs(title = "Visits by month (exposure)", x = NULL, y = "Visits (m), 10 years")
save_plot(p_month_vol / p_month + plot_layout(heights = c(1, 3)), "f02_seasonality", h = 8)

by_hour <- inc %>% count(inc_hour, inc_type) %>% group_by(inc_type) %>% mutate(share = n / sum(n))
save_plot(ggplot(by_hour, aes(inc_hour, share, colour = inc_type)) + geom_line(linewidth = 1) +
            scale_colour_manual(values = sev_cols) + scale_y_continuous(labels = percent) +
            labs(title = "Time of day incidents occur (share within severity)", x = "Hour", y = NULL, colour = NULL),
          "f03_hour_of_day")

## 3. Where - beach-level risk (counts AND rates; small beaches have noisy rates)
by_beach <- rate_tab(model_df, beach_nm) %>% left_join(beach_static, by = "beach_nm") %>% arrange(desc(Severe))
save_tab(by_beach, "t03_by_beach")
cat("\n== Top 10 beaches by deaths (10 yrs) ==\n")
print(head(by_beach %>% dplyr::select(beach_nm, visitors_m, Severe, severe_per_m, ever_patrolled), 10))
p_beach <- by_beach %>% ggplot(aes(visitors_m, Severe, colour = factor(ever_patrolled))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(data = head(by_beach, 8), aes(label = beach_nm), size = 2.8, vjust = -0.8, show.legend = FALSE) +
  scale_colour_manual(values = c(`0` = "#C0392B", `1` = "#2471A3"), labels = c("Never patrolled", "Patrolled (some months)")) +
  labs(title = "Deaths vs visits by beach, 10 years", x = "Visits (millions)", y = "Deaths", colour = NULL)
save_plot(p_beach, "f04_beach_deaths_vs_visits")

## 4. One-way analysis of categorical drivers (crude rates - NOT causal; confounded with season/crowds)
cat_vars <- c("patrolled", "flags", "pa_system", "drone_any", "campaign_any", "fishing", "rocky", "rockpool", "marine_health", "season")
oneway_cat <- bind_rows(lapply(cat_vars, function(v) rate_tab(model_df, level = as.character(.data[[v]])) %>% mutate(var = v, .before = 1)))
save_tab(oneway_cat, "t04_oneway_categorical")
p_cat <- oneway_cat %>% ggplot(aes(level, severe_per_m)) + geom_col(fill = sev_cols["Severe"]) +
  facet_wrap(~ var, scales = "free_x", nrow = 2) + labs(title = "Crude deaths per million visits by beach feature", x = NULL, y = NULL)
save_plot(p_cat, "f05_oneway_categorical", w = 10, h = 5)

## 5. One-way analysis of continuous drivers (rate by quintile)
cont_vars <- c("swell_avg_m", "swell_max_m", "n_rips", "tourist_ratio", "surfer_prop", "water_temp", "air_temp_hi",
               "wind_days", "rain_mm", "turbidity", "lg_hrs_per_10k", "drone_hrs", "campaign_idx", "first_aid_cap",
               "equip_rating", "icecream_visits", "n_bins", "shore_len_m")
oneway_cont <- bind_rows(lapply(cont_vars, function(v) model_df %>%
                                  mutate(bin = ntile(.data[[v]], 5)) %>% group_by(bin) %>%
                                  summarise(mid = median(.data[[v]]), visitors_m = sum(visitors) / 1e6,
                                            total_per_m = sum(Total) / visitors_m, severe_per_m = sum(Severe) / visitors_m, .groups = "drop") %>%
                                  mutate(var = v, .before = 1)))
save_tab(oneway_cont, "t05_oneway_continuous")
p_cont <- oneway_cont %>% pivot_longer(c(total_per_m, severe_per_m)) %>%
  group_by(var, name) %>% mutate(rel = value / mean(value)) %>% ungroup() %>%
  ggplot(aes(bin, rel, colour = name)) + geom_hline(yintercept = 1, colour = "grey70") + geom_line(linewidth = 0.9) + geom_point(size = 1.2) +
  facet_wrap(~ var, ncol = 6) + scale_colour_manual(values = c(severe_per_m = unname(sev_cols["Severe"]), total_per_m = "grey40"),
                                                    labels = c("Deaths", "All incidents")) +
  labs(title = "Relative incident rate by quintile of each variable (1 = average)", x = "Quintile (low -> high)", y = "Relative rate", colour = NULL)
save_plot(p_cont, "f06_oneway_continuous", w = 12, h = 6.5)

## 6. Multicollinearity among candidate predictors
pred_vars <- c("visitors", cont_vars[cont_vars != "lg_hrs_per_10k"], "lg_hrs", "patrolled", "pa_system", "flags", "fishing", "rocky", "rockpool")
cm <- cor(model_df[pred_vars], method = "spearman")
png(file.path(eda_dir, "f07_correlation.png"), width = 2000, height = 2000, res = 220)
corrplot(cm, method = "color", type = "upper", order = "hclust", tl.col = "black", tl.cex = 0.7, diag = FALSE)
invisible(dev.off())
vif <- sapply(pred_vars, function(v) 1 / (1 - summary(lm(reformulate(setdiff(pred_vars, v), v), data = model_df))$r.squared))
high_cor <- as.data.frame(as.table(cm)) %>% filter(as.character(Var1) < as.character(Var2), abs(Freq) > 0.6) %>% arrange(desc(abs(Freq)))
cat("\n== Predictor pairs with |Spearman rho| > 0.6 ==\n"); print(high_cor)
cat("\n== Variance inflation factors (> 5 worth watching) ==\n"); print(round(sort(vif, decreasing = TRUE), 1))
save_tab(tibble(var = names(vif), vif = vif), "t06_vif")

## 7. Distribution of the response - justifies Poisson vs NB vs Tweedie / zero-inflation
disp <- model_df %>% summarise(across(c(Mild, CloseCall, Severe, Total),
                                      list(mean = mean, var = var, zero_share = ~ mean(.x == 0), max = max), .names = "{.col}.{.fn}")) %>%
  pivot_longer(everything(), names_to = c("sev", "stat"), names_sep = "\\.") %>%
  pivot_wider(names_from = stat, values_from = value) %>% mutate(var_to_mean = var / mean)
cat("\n== Response distribution per beach-month (var/mean >> 1 => overdispersion) ==\n"); print(disp)
save_tab(disp, "t07_response_dispersion")
sev_pois <- glm(Severe ~ 1 + offset(log(visitors)), poisson, model_df)
cat(sprintf("Severe: observed zero beach-months %d vs Poisson-expected %.0f (intercept-only rate model)\n",
            sum(model_df$Severe == 0), sum(dpois(0, fitted(sev_pois)))))
p_dist <- model_df %>% pivot_longer(c(Mild, CloseCall, Severe), names_to = "sev", values_to = "n") %>%
  mutate(sev = factor(sev, names(sev_cols))) %>%
  ggplot(aes(n, fill = sev)) + geom_histogram(binwidth = 1) + facet_wrap(~ sev, scales = "free") +
  scale_fill_manual(values = sev_cols, guide = "none") + labs(title = "People per beach-month", x = "Count", y = "Beach-months")
save_plot(p_dist, "f08_response_distribution")

## 8. Person-level profile of incidents (who) - descriptive only, NOT used as model predictors
person_profile <- function(v) inc %>% count(level = as.character(.data[[v]]), inc_type) %>%
  group_by(level) %>% mutate(n_level = sum(n), share = n / n_level) %>% filter(inc_type == "Severe") %>%
  transmute(var = v, level, people = n_level, severe = n, severe_share = share)
person_tab <- bind_rows(lapply(c("swim_level", "local_ind", "age_band", "sex", "surfer", "alone", "inc_cause",
                                 "drone_present", "n_lifeguards"), person_profile))
save_tab(person_tab, "t08_person_profile_severity")
p_person <- person_tab %>% ggplot(aes(level, severe_share)) + geom_col(fill = sev_cols["Severe"]) +
  facet_wrap(~ var, scales = "free_x", ncol = 3) + scale_y_continuous(labels = percent) +
  labs(title = "Share of people in an incident who died, by characteristic", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_person, "f09_person_profile", w = 10, h = 8)

## 9. Response time vs safety resources (on-scene response, late reports excluded)
resp <- inc %>% filter(late_report == 0) %>% distinct(incident_id, .keep_all = TRUE) %>%
  group_by(n_lifeguards = pmin(n_lifeguards, 3), drone_present) %>%
  summarise(incidents = n(), median_resp = median(response_onscene), p90_resp = quantile(response_onscene, .9), .groups = "drop")
cat("\n== On-scene response time (min) by lifeguards on duty (3 = 3+) and drone ==\n"); print(resp)
save_tab(resp, "t09_response_time")

## 10. Natural experiments - beaches where flags / PA were installed mid-period
switchers <- model_df %>% group_by(beach_nm) %>% filter(n_distinct(pa_system) > 1 | n_distinct(flags) > 1) %>%
  group_by(beach_nm, pa_system, flags) %>%
  summarise(months = n(), visitors_m = sum(visitors) / 1e6, total_per_m = sum(Total) / visitors_m,
            severe = sum(Severe), .groups = "drop")
cat("\n== Beaches that installed PA / flags during the period (before vs after) ==\n"); print(switchers, n = 40)
save_tab(switchers, "t10_install_before_after")

cat("\nSaved EDA outputs to", eda_dir, "\n")