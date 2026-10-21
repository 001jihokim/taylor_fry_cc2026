######################################################################
# 01_setup_eda_cleaning.R
# UNSW ASOC x Taylor Fry Case Competition 2026 - Omega Consulting
# Jiho Kim, Becky Thai, Christina Lou, Timothy Bae
#
# Part 1 of 2: setup, import, cleaning, EDA and pre-modelling checks
#   Reads:  data/*.csv               (raw data, never modified)
#   Writes: temp_output/clean/*.rds  (model-ready tables, read by 02_modelling.R)
#           temp_output/eda/         (Figure 3.x / 4.x as .png, tables as .csv)
#
# Assumption: the data provided are accurate to the data dictionary, so no separate
# validation step is run. Cleaning only handles what the dictionary itself describes
# (e.g. NA causes for children, reporting delay in response_min).
#
# Open taylor_fry_cc2026.Rproj first so the relative paths below resolve.
######################################################################


######################################################################
# PART 0 SETUP
######################################################################

pkgs <- c("dplyr", "tidyr", "readr", "lubridate", "ggplot2", "scales",
          "patchwork", "corrplot", "car", "MASS")
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(ggplot2)
library(scales)
library(patchwork)
library(corrplot)
# car (VIF) and MASS (negative binomial) are called as car:: and MASS::
# rather than attached, because they mask dplyr functions (select, recode)

set.seed(2026)

data_dir  <- "data"
clean_dir <- file.path("temp_output", "clean")
eda_dir   <- file.path("temp_output", "eda")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(eda_dir,   recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 11))
sev_cols <- c(Mild = "#7FB3D5", CloseCall = "#F5B041", Severe = "#C0392B")

# helpers
save_fig <- function(plot, name, width = 8, height = 4.5) {
  ggsave(file.path(eda_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 200)
}

save_tab <- function(table, name) {
  write_csv(table, file.path(eda_dir, paste0(name, ".csv")))
}

# crude rates per million visits for any grouping of the beach-month table
rate_tab <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      beach_months = n(),
      visitors_m   = sum(visitors) / 1e6,
      Mild         = sum(Mild),
      CloseCall    = sum(CloseCall),
      Severe       = sum(Severe),
      mild_per_m   = Mild / visitors_m,
      close_per_m  = CloseCall / visitors_m,
      severe_per_m = Severe / visitors_m,
      .groups      = "drop"
    )
}

# fiscal year runs Sep-Aug to match the data period (FY2026 = Sep 2025 - Aug 2026)
fy_of <- function(date) year(date) + (month(date) >= 9)


######################################################################
# PART 1 IMPORT
######################################################################
# explicit column types so nothing is silently mis-parsed

beach_raw <- read_csv(
  file.path(data_dir, "beach_monthly.csv"),
  col_types = cols(beach_nm = "c", month = "c", marine_health = "c", .default = "d")
)

inc_raw <- read_csv(
  file.path(data_dir, "incidents.csv"),
  col_types = cols(
    dte_reported  = "c", dte_occurred = "c", inc_type   = "c", inc_cause = "c",
    sex           = "c", swim_level   = "c", local_ind  = "c", beach_nm  = "c",
    drone_present = "c", .default     = "d"
  )
)

price_raw <- read_csv(file.path(data_dir, "pricing.csv"), show_col_types = FALSE)


######################################################################
# PART 2 CLEANING
######################################################################
# Summary of fixes (detail and reasoning at each step below)
#   1. inc_cause missing (1,764 rows, all children): inherit the cause recorded for an
#      adult in the same incident (140 rows); label the rest "Unspecified" (1,624)  -> 2.3
#   2. surfer missing (1,734 rows): labelled "Unknown", rows kept                 -> 2.3
#   3. response_min outliers (max 5,801 min): remove the reporting delay to give an
#      on-scene response time (max 320 min); original column kept                -> 2.3
#   Smaller fixes:
#   - three blank rows dropped from pricing.csv                                   -> 2.1
#   - FY2020-21 flagged as a dip period                                           -> 2.2
#   - values sitting exactly on a limit (e.g. rain at 700 mm) are kept as-is, with
#     no flag: a sensitivity analysis on them adds complexity for little benefit
#   Variable choice: patrolled is dropped as a model predictor because it is exactly
#   lg_hrs > 0, so it is collinear with lifeguard hours by construction           -> 2.2
#   icecream_visits and n_bins are dropped as predictors: crowd / heat proxies    -> 4.3
#   Modelling table: people counted by beach-month and severity, zero-filled, with
#   checks that all 33,941 people and all 220 deaths survive the join             -> 2.4
#   Excluded from the models: person-level fields and response time               -> 2.3

######################################################################
# 2.1 pricing - drop the three trailing blank rows, tidy text
pricing <- price_raw %>%
  filter(!is.na(item)) %>%
  mutate(
    item            = trimws(item),
    effect_per_unit = gsub("\\s+", " ", effect_per_unit)
  )

######################################################################
# 2.2 beach-month table
beach <- beach_raw %>%
  mutate(
    date          = dmy(paste0("01/", month)),
    cal_m         = month(date),
    month_lab     = factor(month.abb[cal_m], levels = month.abb[c(9:12, 1:8)]),  # Sep..Aug
    fy            = fy_of(date),
    season        = case_when(
      cal_m %in% c(12, 1, 2) ~ "Summer",
      cal_m %in% 3:5         ~ "Autumn",
      cal_m %in% 6:8         ~ "Winter",
      TRUE                   ~ "Spring"
    ),
    season        = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter")),
    # Sep-Apr is the patrol season: ~80% of beach-months patrolled vs ~12% in May-Aug
    patrol_season = as.integer(cal_m %in% c(9:12, 1:4)),
    marine_health = factor(marine_health, levels = c("L", "M", "H")),
    postcode      = factor(postcode),
    
    # exposure and intensity measures
    log_visitors   = log(visitors),
    # lifeguards enter the models as log_lg_hrs only. patrolled (= lg_hrs > 0) stays in
    # the data for description but is not a predictor: log_lg_hrs is 0 exactly when
    # patrolled is 0, so the two are collinear by construction
    log_lg_hrs     = log1p(lg_hrs),       # lg_hrs is 0 at unpatrolled beach-months
    log_drone_hrs  = log1p(drone_hrs),    # drone_hrs is 0 in 77% of beach-months
    lg_hrs_per_10k = 1e4 * lg_hrs / visitors,
    drone_any      = as.integer(drone_hrs > 0),
    campaign_any   = as.integer(campaign_idx > 0),
    
    # FY2020-21: visits and incidents both dip, then recover - flag to test a level shift
    dip_fy2020_21 = as.integer(fy %in% c(2020, 2021))
  )

######################################################################
# 2.3 incident (person-level) table
#
# Fix 1 - missing inc_cause (1,764 rows, all children)
#   What: 140 rows take the cause recorded for an adult in the same incident; the
#         other 1,624 are labelled "Unspecified".
#   Why:  cause is an incident-level field, so a child in a family incident shares the
#         adult's cause. Imputing anything beyond that would be guessing.
#
# Fix 2 - missing surfer (1,734 rows)
#   What: labelled "Unknown" rather than dropping the rows.
#   Why:  the ages and severities of these rows look like the rest of the data, which
#         suggests they are missing at random. Dropping rows would undercount incidents
#         at beach-month level, which feeds the forecast directly.
#
# Fix 3 - response_min outliers (max 5,801 min, about 4 days)
#   Cause: the dictionary says the measure includes reporting delay. All 52 extreme
#          values are incidents reported 1-4 days late, each adding ~1,440 min per day.
#   What:  subtracting that delay gives an on-scene response time (response_onscene)
#          with a clean maximum of 320 min; the original column is kept.
#   Why:   these 52 rows would otherwise dominate any average response time.
#
# Fields deliberately excluded from the models (kept here for description only):
#   person-level fields such as swim_level, age and alone describe who had an incident,
#   not the beach population, and are unknown when forecasting. Response time is
#   excluded too, because it is only observed after an incident has happened.
inc <- inc_raw %>%
  mutate(
    d_occ = dmy(dte_occurred),
    d_rep = dmy(dte_reported),
    month = format(d_occ, "%m/%Y"),   # join key: month OCCURRED, not reported (dictionary s6)
    cal_m = month(d_occ),
    fy    = fy_of(d_occ),
    
    # Fix 3: remove the reporting delay from response_min
    report_lag_days  = as.integer(d_rep - d_occ),
    late_report      = as.integer(report_lag_days > 0),
    response_onscene = response_min - 1440 * report_lag_days,
    
    inc_type   = factor(inc_type, levels = c("Mild", "CloseCall", "Severe"), ordered = TRUE),
    swim_level = factor(swim_level, levels = c("beginner", "intermediate", "experienced"), ordered = TRUE),
    # Fix 2: missing surfer -> "Unknown"
    surfer     = factor(
      case_when(is.na(surfer) ~ "Unknown", surfer == 1 ~ "Yes", TRUE ~ "No"),
      levels = c("No", "Yes", "Unknown")
    ),
    child      = as.integer(age < 18),
    age_band   = cut(age, breaks = c(0, 12, 17, 24, 39, 59, 74, Inf),
                     labels = c("0-12", "13-17", "18-24", "25-39", "40-59", "60-74", "75+")),
    local      = as.integer(local_ind == "local"),
    drone_present   = as.integer(drone_present == "Y"),
    in_patrol_hours = as.integer(n_lifeguards > 0)
  ) %>%
  # Fix 1: inherit inc_cause from the same incident, otherwise "Unspecified"
  group_by(incident_id) %>%
  mutate(
    n_people  = n(),
    inc_cause = coalesce(inc_cause, first(na.omit(inc_cause)), "Unspecified")
  ) %>%
  ungroup()

stopifnot(all(inc$response_onscene > 0))

cat(sprintf("\ninc_cause: %d child rows recovered from the same incident, %d left Unspecified\n",
            sum(is.na(inc_raw$inc_cause)) - sum(inc$inc_cause == "Unspecified"),
            sum(inc$inc_cause == "Unspecified")))
cat(sprintf("surfer: %d people recorded as Unknown (%.1f%%)\n",
            sum(inc$surfer == "Unknown"), 100 * mean(inc$surfer == "Unknown")))
cat(sprintf("response_min: %d late reports; on-scene max %.0f min (raw max %.0f)\n",
            sum(inc$late_report), max(inc$response_onscene), max(inc$response_min)))

######################################################################
# 2.4 beach-month modelling table (dictionary s6)
# Incident rows are aggregated to PEOPLE by beach-month and severity, then joined to
# the monthly data. People are the unit because inc_type is person-level and the client
# objective counts severe outcomes. EVENTS (distinct incident_id) are kept for reference.
# Months with no incidents are set to ZERO, not missing: 98% of beach-months have zero
# deaths, and those zeros are real information. The checks below confirm all 33,941
# people and all 220 deaths survive the join.

people_counts <- inc %>%
  count(beach_nm, month, inc_type) %>%
  pivot_wider(names_from = inc_type, values_from = n, values_fill = 0)

event_counts <- inc %>%
  distinct(beach_nm, month, incident_id) %>%
  count(beach_nm, month, name = "n_events")

cause_counts <- inc %>%
  count(beach_nm, month, inc_cause) %>%
  pivot_wider(names_from = inc_cause, values_from = n, values_fill = 0, names_prefix = "cause_")

model_df <- beach %>%
  left_join(people_counts, by = c("beach_nm", "month")) %>%
  left_join(event_counts,  by = c("beach_nm", "month")) %>%
  left_join(cause_counts,  by = c("beach_nm", "month")) %>%
  mutate(
    across(c(Mild, CloseCall, Severe, n_events, starts_with("cause_")), ~ replace_na(.x, 0L)),
    Total = Mild + CloseCall + Severe
  )

stopifnot(
  nrow(model_df) == 9600,
  sum(model_df$Total)    == nrow(inc),
  sum(model_df$Severe)   == sum(inc$inc_type == "Severe"),
  sum(model_df$n_events) == n_distinct(inc$incident_id)
)

cat(sprintf("\nmodel_df: %d beach-months, %d people and %d deaths after the join; %.1f%% of beach-months have zero deaths\n",
            nrow(model_df), sum(model_df$Total), sum(model_df$Severe), 100 * mean(model_df$Severe == 0)))

######################################################################
# 2.5 postcode-month table (dictionary s6)
# campaigns and community first aid are bought per postcode, so Task 2 costing
# and any postcode-level effects need this view
postcode_month <- model_df %>%
  group_by(postcode, month, date, fy, cal_m) %>%
  summarise(
    n_beaches     = n(),
    visitors      = sum(visitors),
    campaign_idx  = first(campaign_idx),
    first_aid_cap = first(first_aid_cap),
    Mild          = sum(Mild),
    CloseCall     = sum(CloseCall),
    Severe        = sum(Severe),
    Total         = sum(Total),
    .groups       = "drop"
  )

######################################################################
# 2.6 static beach table (one row per beach) - 'where' analysis and Task 2 costing
beach_static <- beach %>%
  group_by(beach_nm, postcode) %>%
  summarise(
    shore_len_m      = first(shore_len_m),
    n_rips           = first(n_rips),
    rocky            = first(rocky),
    fishing          = first(fishing),
    rockpool         = first(rockpool),
    pa_system_now    = last(pa_system, order_by = date),
    flags_now        = last(flags, order_by = date),
    months_patrolled = sum(patrolled),
    ever_patrolled   = as.integer(any(patrolled == 1)),
    .groups          = "drop"
  ) %>%
  add_count(postcode, name = "beaches_in_postcode")

######################################################################
# 2.7 save cleaned data (all generated files go to temp_output/, never to data/)
saveRDS(model_df,       file.path(clean_dir, "model_df.rds"))         # beach-month counts + drivers
saveRDS(inc,            file.path(clean_dir, "incidents_clean.rds"))  # person-level, cleaned
saveRDS(beach_static,   file.path(clean_dir, "beach_static.rds"))     # one row per beach
saveRDS(postcode_month, file.path(clean_dir, "postcode_month.rds"))   # postcode-level view
saveRDS(pricing,        file.path(clean_dir, "pricing.rds"))
cat("\nSaved cleaned data to", clean_dir, "\n")


######################################################################
# PART 3 EXPLORATORY DATA ANALYSIS
######################################################################
# crude rates are per million visits. They are NOT causal - see Part 4 for confounding.

######################################################################
# Figure 3.1: trend by fiscal year
by_fy <- rate_tab(model_df, fy)
save_tab(by_fy, "t01_by_fiscal_year")
cat("\n== By fiscal year ==\n")
print(by_fy, width = Inf)

plot_fy_counts <- by_fy %>%
  pivot_longer(c(Mild, CloseCall, Severe), names_to = "sev", values_to = "people") %>%
  mutate(sev = factor(sev, levels = names(sev_cols))) %>%
  ggplot(aes(x = fy, y = people, colour = sev)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~ sev, scales = "free_y") +
  scale_colour_manual(values = sev_cols, guide = "none") +
  scale_x_continuous(breaks = seq(2017, 2026, by = 3)) +
  labs(title = "Figure 3.1: People involved in incidents by fiscal year",
       x = "Fiscal year (Sep-Aug)", y = "People")

plot_fy_rate <- ggplot(by_fy, aes(x = fy)) +
  geom_col(aes(y = severe_per_m), fill = sev_cols["Severe"]) +
  geom_line(aes(y = visitors_m / 30), colour = "grey30", linewidth = 0.8) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 30, name = "Visits (millions)")) +
  scale_x_continuous(breaks = 2017:2026) +
  labs(title = "Deaths per million visits (bars) vs total visits (line)",
       x = "Fiscal year", y = "Deaths per million visits")

save_fig(plot_fy_counts / plot_fy_rate, "f3_1_trend_by_fy", height = 7)

######################################################################
# Figure 3.2: seasonality - when do incidents happen?
by_month <- rate_tab(model_df, month_lab)
save_tab(by_month, "t02_by_calendar_month")

plot_month_vol <- ggplot(by_month, aes(x = month_lab, y = visitors_m)) +
  geom_col(fill = "grey60") +
  labs(title = "Figure 3.2: Visits by month (exposure)", x = NULL, y = "Visits (m), 10 years")

plot_month_rate <- by_month %>%
  pivot_longer(c(mild_per_m, close_per_m, severe_per_m), names_to = "sev", values_to = "rate") %>%
  mutate(sev = factor(sev, levels = c("mild_per_m", "close_per_m", "severe_per_m"),
                      labels = names(sev_cols))) %>%
  ggplot(aes(x = month_lab, y = rate, fill = sev)) +
  geom_col() +
  facet_wrap(~ sev, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = sev_cols, guide = "none") +
  labs(title = "Incident rate per million visits by month", x = NULL, y = "People per million visits")

save_fig(plot_month_vol / plot_month_rate + plot_layout(heights = c(1, 3)),
         "f3_2_seasonality", height = 8)

######################################################################
# Figure 3.3: time of day, and whether lifeguards were on duty
by_hour <- inc %>%
  group_by(inc_hour, inc_type) %>%
  summarise(people = n(), on_duty = mean(in_patrol_hours), .groups = "drop") %>%
  group_by(inc_type) %>%
  mutate(share = people / sum(people)) %>%
  ungroup()

plot_hour_share <- ggplot(by_hour, aes(x = inc_hour, y = share, colour = inc_type)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = sev_cols) +
  scale_y_continuous(labels = percent) +
  labs(title = "Figure 3.3: Time of day incidents occur (share within severity)",
       x = NULL, y = NULL, colour = NULL)

plot_hour_duty <- inc %>%
  group_by(inc_hour) %>%
  summarise(on_duty = mean(in_patrol_hours)) %>%
  ggplot(aes(x = inc_hour, y = on_duty)) +
  geom_col(fill = "#2471A3") +
  scale_y_continuous(labels = percent) +
  labs(title = "Share of incidents with a lifeguard on duty", x = "Hour", y = NULL)

save_fig(plot_hour_share / plot_hour_duty + plot_layout(heights = c(2, 1)),
         "f3_3_hour_of_day", height = 6)

######################################################################
# Figure 3.4: where - beach-level deaths vs visits (small beaches have noisy rates)
by_beach <- rate_tab(model_df, beach_nm) %>%
  left_join(beach_static, by = "beach_nm") %>%
  arrange(desc(Severe))
save_tab(by_beach, "t03_by_beach")

cat("\n== Top 10 beaches by deaths (10 years) ==\n")
print(by_beach %>% select(beach_nm, visitors_m, Severe, severe_per_m, ever_patrolled) %>% head(10))

plot_beach <- ggplot(by_beach, aes(x = visitors_m, y = Severe, colour = factor(ever_patrolled))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(data = head(by_beach, 8), aes(label = beach_nm),
            size = 2.8, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = c(`0` = "#C0392B", `1` = "#2471A3"),
                      labels = c("Never patrolled", "Patrolled (some months)")) +
  labs(title = "Figure 3.4: Deaths vs visits by beach, 10 years",
       x = "Visits (millions)", y = "Deaths", colour = NULL)

save_fig(plot_beach, "f3_4_beach_deaths_vs_visits")

######################################################################
# Figure 3.5: one-way analysis of categorical drivers
cat_vars <- c("patrolled", "flags", "pa_system", "drone_any", "campaign_any",
              "fishing", "rocky", "rockpool", "marine_health", "season")

oneway_cat <- bind_rows(lapply(cat_vars, function(v) {
  rate_tab(model_df, level = as.character(.data[[v]])) %>%
    mutate(var = v, .before = 1)
}))
save_tab(oneway_cat, "t04_oneway_categorical")

plot_cat <- ggplot(oneway_cat, aes(x = level, y = severe_per_m)) +
  geom_col(fill = sev_cols["Severe"]) +
  facet_wrap(~ var, scales = "free_x", nrow = 2) +
  labs(title = "Figure 3.5: Crude deaths per million visits by beach feature", x = NULL, y = NULL)

save_fig(plot_cat, "f3_5_oneway_categorical", width = 10, height = 5)

######################################################################
# Figure 3.6: one-way analysis of continuous drivers (rate by quintile)
cont_vars <- c("swell_avg_m", "swell_max_m", "n_rips", "tourist_ratio", "surfer_prop",
               "water_temp", "air_temp_hi", "wind_days", "rain_mm", "turbidity",
               "lg_hrs_per_10k", "drone_hrs", "campaign_idx", "first_aid_cap",
               "equip_rating", "shore_len_m")   # crowd proxies excluded (see Part 4)

oneway_cont <- bind_rows(lapply(cont_vars, function(v) {
  model_df %>%
    mutate(bin = ntile(.data[[v]], 5)) %>%
    group_by(bin) %>%
    summarise(
      mid          = median(.data[[v]]),
      visitors_m   = sum(visitors) / 1e6,
      total_per_m  = sum(Total) / visitors_m,
      severe_per_m = sum(Severe) / visitors_m,
      .groups      = "drop"
    ) %>%
    mutate(var = v, .before = 1)
}))
save_tab(oneway_cont, "t05_oneway_continuous")

plot_cont <- oneway_cont %>%
  pivot_longer(c(total_per_m, severe_per_m), names_to = "measure", values_to = "rate") %>%
  group_by(var, measure) %>%
  mutate(rel = rate / mean(rate)) %>%
  ungroup() %>%
  ggplot(aes(x = bin, y = rel, colour = measure)) +
  geom_hline(yintercept = 1, colour = "grey70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  facet_wrap(~ var, ncol = 6) +
  scale_colour_manual(values = c(severe_per_m = unname(sev_cols["Severe"]), total_per_m = "grey40"),
                      labels = c(severe_per_m = "Deaths", total_per_m = "All incidents")) +
  labs(title = "Figure 3.6: Relative incident rate by quintile of each variable (1 = average)",
       x = "Quintile (low -> high)", y = "Relative rate", colour = NULL)

save_fig(plot_cont, "f3_6_oneway_continuous", width = 12, height = 6.5)

######################################################################
# Figure 3.7: who - share of people in an incident who died
# descriptive only: person-level fields are not known in advance, so NOT model predictors
person_profile <- function(v) {
  inc %>%
    count(level = as.character(.data[[v]]), inc_type) %>%
    group_by(level) %>%
    mutate(people = sum(n), severe_share = n / people) %>%
    ungroup() %>%
    filter(inc_type == "Severe") %>%
    transmute(var = v, level, people, severe = n, severe_share)
}

person_vars <- c("swim_level", "local_ind", "age_band", "sex", "surfer", "alone",
                 "inc_cause", "drone_present", "n_lifeguards")
person_tab  <- bind_rows(lapply(person_vars, person_profile))
save_tab(person_tab, "t06_person_profile_severity")

plot_person <- ggplot(person_tab, aes(x = level, y = severe_share)) +
  geom_col(fill = sev_cols["Severe"]) +
  facet_wrap(~ var, scales = "free_x", ncol = 3) +
  scale_y_continuous(labels = percent) +
  labs(title = "Figure 3.7: Share of people in an incident who died, by characteristic",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(plot_person, "f3_7_person_profile", width = 10, height = 8)

######################################################################
# Figure 3.8: response time vs safety resources (on-scene response, late reports excluded)
resp <- inc %>%
  filter(late_report == 0) %>%
  distinct(incident_id, .keep_all = TRUE) %>%
  group_by(n_lifeguards = pmin(n_lifeguards, 3), drone_present) %>%
  summarise(
    incidents   = n(),
    median_resp = median(response_onscene),
    p90_resp    = quantile(response_onscene, 0.9),
    .groups     = "drop"
  )
save_tab(resp, "t07_response_time")

plot_resp <- resp %>%
  mutate(
    lifeguards = factor(n_lifeguards, levels = 0:3, labels = c("0", "1", "2", "3+")),
    drone      = ifelse(drone_present == 1, "Drone flying", "No drone")
  ) %>%
  ggplot(aes(x = lifeguards, y = median_resp, fill = drone)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = median_resp, ymax = p90_resp),
                position = position_dodge(width = 0.8), width = 0.2, colour = "grey30") +
  scale_fill_manual(values = c("No drone" = "grey60", "Drone flying" = "#2471A3")) +
  labs(title = "Figure 3.8: On-scene response time by lifeguards on duty and drone",
       subtitle = "Bar = median, line = up to 90th percentile; late reports excluded",
       x = "Lifeguards on duty", y = "Minutes", fill = NULL)

save_fig(plot_resp, "f3_8_response_time")

######################################################################
# Figure 3.9: natural experiments - beaches that installed a PA system or flags mid-period
installs <- model_df %>%
  group_by(beach_nm) %>%
  filter(n_distinct(pa_system) > 1 | n_distinct(flags) > 1) %>%
  mutate(installed = case_when(
    n_distinct(pa_system) > 1 & n_distinct(flags) > 1 ~ "PA system and flags",
    n_distinct(pa_system) > 1                         ~ "PA system",
    TRUE                                              ~ "Flags"
  )) %>%
  group_by(beach_nm, installed, pa_system, flags) %>%
  summarise(
    start       = min(date),
    months      = n(),
    visitors_m  = sum(visitors) / 1e6,
    total_per_m = sum(Total) / visitors_m,
    severe      = sum(Severe),
    .groups     = "drop"
  ) %>%
  group_by(beach_nm) %>%
  mutate(period = factor(ifelse(start == min(start), "Before", "After"),
                         levels = c("Before", "After"))) %>%
  ungroup()
save_tab(installs, "t08_install_before_after")

plot_installs <- ggplot(installs, aes(x = total_per_m, y = beach_nm)) +
  geom_line(aes(group = beach_nm), colour = "grey70") +
  geom_point(aes(colour = period, size = months)) +
  geom_text(aes(label = paste0(severe, " d")), size = 2.6, vjust = -1.1, colour = "grey30") +
  facet_wrap(~ installed, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(Before = "grey50", After = "#2471A3")) +
  labs(title = "Figure 3.9: Incident rate before vs after PA or flag installs",
       subtitle = "Label = deaths in that period. Short periods (small points) are noisy",
       x = "All incidents per million visits", y = NULL, colour = NULL, size = "Months")

save_fig(plot_installs, "f3_9_install_before_after", height = 6)


######################################################################
# PART 4 PRE-MODELLING CHECKS
######################################################################
# what has to be accounted for before 02_modelling.R:
#   4.1 response distribution (Poisson vs NB, zero inflation)
#   4.2 exposure (is a log(visitors) offset justified?)
#   4.3 correlation and 4.4 VIF (multicollinearity)
#   4.5 confounding with season (patrols and campaigns are seasonal by design)
#   4.6 within- vs between-beach variation (what the data can identify)
#   4.7 time and beach effects not captured by the drivers
#   4.8 sample size for the severe model

# candidate predictors, on the scale they would enter a rate model
# (visitors is the exposure, so it is an offset rather than a predictor)
# patrolled is not a candidate: it is exactly lg_hrs > 0 (see 2.2)
# icecream_visits and n_bins are not candidates: they are crowd / heat proxies that
# track visitors and temperature rather than drive risk (evidence printed in 4.3)
levers   <- c("log_lg_hrs", "log_drone_hrs", "campaign_idx",
              "first_aid_cap", "pa_system", "flags")
env_vars <- c("equip_rating", "swell_avg_m", "swell_max_m", "wind_days", "water_temp",
              "air_temp_hi", "rain_mm", "turbidity", "tourist_ratio", "surfer_prop",
              "n_rips", "shore_len_m", "rocky", "fishing", "rockpool", "marine_health")
crowd_proxies <- c("icecream_visits", "n_bins")
cand_vars <- c(levers, env_vars)

rate_formula <- function(y, vars) {
  reformulate(c(vars, "offset(log(visitors))"), response = y)
}

######################################################################
# 4.1 response distribution: over-dispersion and excess zeros
disp <- model_df %>%
  summarise(across(c(Mild, CloseCall, Severe, Total),
                   list(mean = mean, var = var, zero_share = ~ mean(.x == 0), max = max),
                   .names = "{.col}.{.fn}")) %>%
  pivot_longer(everything(), names_to = c("sev", "stat"), names_sep = "\\.") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(var_to_mean = var / mean)

# crude var/mean overstates dispersion because rates differ by beach and month, so
# also check Pearson dispersion and expected zeros AFTER fitting the drivers + month
disp_fit <- bind_rows(lapply(c("Mild", "CloseCall", "Severe", "Total"), function(y) {
  fit_pois <- glm(rate_formula(y, c(cand_vars, "month_lab")), family = poisson, data = model_df)
  fit_nb   <- suppressWarnings(MASS::glm.nb(rate_formula(y, c(cand_vars, "month_lab")), data = model_df))
  tibble(
    sev               = y,
    pearson_disp      = sum(residuals(fit_pois, type = "pearson")^2) / df.residual(fit_pois),
    zeros_observed    = sum(model_df[[y]] == 0),
    zeros_poisson_exp = sum(dpois(0, fitted(fit_pois))),
    nb_theta          = fit_nb$theta,
    aic_poisson       = AIC(fit_pois),
    aic_nb            = AIC(fit_nb)
  )
}))

disp <- left_join(disp, disp_fit, by = "sev")
save_tab(disp, "t09_response_dispersion")
cat("\n== 4.1 Response distribution per beach-month ==\n")
cat("(Pearson dispersion >> 1 => over-dispersion; observed zeros >> expected => zero inflation)\n")
print(disp, width = Inf)

plot_dist <- model_df %>%
  pivot_longer(c(Mild, CloseCall, Severe), names_to = "sev", values_to = "people") %>%
  mutate(sev = factor(sev, levels = names(sev_cols))) %>%
  ggplot(aes(x = people, fill = sev)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ sev, scales = "free") +
  scale_fill_manual(values = sev_cols, guide = "none") +
  labs(title = "Figure 4.1a: People per beach-month", x = "Count", y = "Beach-months")

save_fig(plot_dist, "f4_1a_response_distribution")

disp_long <- disp %>%
  mutate(sev = factor(sev, levels = c("Mild", "CloseCall", "Severe", "Total")))

plot_disp <- disp_long %>%
  select(sev, `Crude var / mean` = var_to_mean, `Pearson, after drivers` = pearson_disp) %>%
  pivot_longer(-sev, names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = sev, y = value, fill = measure)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("grey60", "#2471A3")) +
  labs(title = "Figure 4.1b: Over-dispersion (1 = Poisson)", x = NULL, y = "Dispersion", fill = NULL)

plot_zeros <- disp_long %>%
  transmute(sev, Observed = zeros_observed / nrow(model_df),
            `Poisson expected` = zeros_poisson_exp / nrow(model_df)) %>%
  pivot_longer(-sev, names_to = "measure", values_to = "share") %>%
  ggplot(aes(x = sev, y = share, fill = measure)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c("grey30", "#F5B041")) +
  labs(title = "Zero beach-months: observed vs Poisson-expected", x = NULL, y = NULL, fill = NULL)

plot_aic <- disp_long %>%
  ggplot(aes(x = sev, y = aic_poisson - aic_nb)) +
  geom_col(fill = "#2471A3") +
  labs(title = "AIC improvement from negative binomial over Poisson", x = NULL, y = "AIC drop")

save_fig(plot_disp / plot_zeros / plot_aic, "f4_1b_dispersion_and_zeros", height = 8)

######################################################################
# 4.2 exposure: an offset assumes incidents are proportional to visits (slope = 1)
# fit log(visitors) as a free covariate and test whether its slope is 1, first with
# month only, then with all candidate drivers (busy beaches also get more lifeguards,
# so the month-only slope mixes the crowd effect with the safety-resource effect)
exposure_slope <- function(y, controls, label) {
  fit <- glm(reformulate(c("log_visitors", controls), response = y),
             family = poisson, data = model_df)
  est <- coef(summary(fit))["log_visitors", ]
  tibble(
    sev      = y,
    controls = label,
    slope    = est[["Estimate"]],
    lower_95 = est[["Estimate"]] - 1.96 * est[["Std. Error"]],
    upper_95 = est[["Estimate"]] + 1.96 * est[["Std. Error"]]
  )
}

exposure <- bind_rows(lapply(c("Mild", "CloseCall", "Severe", "Total"), function(y) {
  bind_rows(exposure_slope(y, "month_lab", "month only"),
            exposure_slope(y, c(cand_vars, "month_lab"), "month + all drivers"))
}))
save_tab(exposure, "t10_exposure_slope")
cat("\n== 4.2 Exposure check: slope on log(visitors) (offset assumes 1) ==\n")
print(exposure)

plot_exposure <- exposure %>%
  mutate(sev = factor(sev, levels = c("Mild", "CloseCall", "Severe", "Total"))) %>%
  ggplot(aes(x = sev, y = slope, colour = controls)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(ymin = lower_95, ymax = upper_95), position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("grey50", "#2471A3")) +
  labs(title = "Figure 4.2: Slope on log(visits) with 95% CI (offset assumes 1)",
       x = NULL, y = "Slope", colour = "Controls")

save_fig(plot_exposure, "f4_2_exposure_slope")

######################################################################
# 4.3 correlation among candidate predictors (Spearman: robust to skew and 0/1 vars)
corr_df  <- model_df %>%
  mutate(marine_health = as.numeric(marine_health)) %>%
  select(all_of(cand_vars), log_visitors)
corr_mat <- cor(corr_df, method = "spearman")

png(file.path(eda_dir, "f4_3_correlation.png"), width = 2000, height = 2000, res = 220)
corrplot(corr_mat, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.cex = 0.7, diag = FALSE,
         title = "Figure 4.3: Spearman correlation of candidate predictors", mar = c(0, 0, 2, 0))
invisible(dev.off())

high_cor <- as.data.frame(as.table(corr_mat)) %>%
  filter(as.character(Var1) < as.character(Var2), abs(Freq) > 0.6) %>%
  arrange(desc(abs(Freq))) %>%
  rename(var_1 = Var1, var_2 = Var2, rho = Freq)
cat("\n== 4.3 Predictor pairs with |Spearman rho| > 0.6 ==\n")
print(high_cor)

# evidence for dropping the crowd proxies: how closely they track crowds and heat
proxy_cor <- cor(model_df[crowd_proxies],
                 model_df[c("log_visitors", "air_temp_hi", "log_lg_hrs")],
                 method = "spearman")
cat("\n== 4.3b Crowd proxies vs visits, heat and lifeguard hours (Spearman rho) ==\n")
print(round(proxy_cor, 2))

######################################################################
# 4.4 variance inflation factors
# VIF_j = 1 / (1 - R^2_j), R^2_j from regressing predictor j on all the others.
# vif_design: plain design-matrix VIF.  vif_glm: weighted by a Poisson rate GLM
# (what the model actually sees). marine_health is a 3-level factor, so car::vif gives a
# GVIF; it is reported as (GVIF^(1/(2 df)))^2 to be comparable with ordinary VIFs.
vif_of <- function(fit) {
  v <- car::vif(fit)
  if (is.matrix(v)) v <- v[, "GVIF^(1/(2*Df))"]^2
  v
}

vif_table <- function(vars) {
  fit_glm <- glm(rate_formula("Total", vars), family = poisson, data = model_df)
  fit_lm  <- lm(reformulate(vars, response = "log_visitors"), data = model_df)  # response irrelevant
  tibble(var = vars) %>%
    mutate(vif_design = vif_of(fit_lm)[var],
           vif_glm    = vif_of(fit_glm)[var]) %>%
    arrange(desc(vif_glm))
}

vif_full <- vif_table(cand_vars)
save_tab(vif_full, "t11_vif_full")
cat("\n== 4.4a VIF, all candidate predictors (> 5 worth watching, > 10 serious) ==\n")
print(vif_full, n = Inf)

# drop the highest-VIF variable one at a time until every VIF <= 5. The safety levers
# are protected: Task 2 needs their effects, so they are never dropped.
keep_vars <- cand_vars
vif_steps <- tibble()
repeat {
  v     <- vif_table(keep_vars)
  worst <- v %>% filter(!var %in% levers) %>% slice_max(vif_glm, n = 1)
  if (nrow(worst) == 0 || worst$vif_glm <= 5) break
  vif_steps <- bind_rows(vif_steps, worst %>% mutate(step = nrow(vif_steps) + 1, .before = 1))
  keep_vars <- setdiff(keep_vars, worst$var)
}
vif_reduced <- vif_table(keep_vars)

save_tab(vif_steps,   "t12_vif_dropped")
save_tab(vif_reduced, "t13_vif_reduced")
cat("\n== 4.4b Variables dropped for VIF > 5, in order ==\n")
print(vif_steps)
cat("\n== 4.4c VIF of the reduced candidate set ==\n")
print(vif_reduced, n = Inf)

plot_vif <- bind_rows(
  vif_full    %>% mutate(set = "All candidates"),
  vif_reduced %>% mutate(set = "After pruning")
) %>%
  mutate(set = factor(set, levels = c("All candidates", "After pruning")),
         var = factor(var, levels = rev(vif_full$var))) %>%
  ggplot(aes(x = vif_glm, y = var, fill = set)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.8) +
  geom_vline(xintercept = 5,  linetype = "dashed") +
  geom_vline(xintercept = 10, linetype = "dotted") +
  scale_fill_manual(values = c("grey60", "#2471A3")) +
  labs(title = "Figure 4.4: Variance inflation factors (Poisson-weighted)",
       subtitle = "Dashed = 5 (worth watching), dotted = 10 (serious). Missing blue bar = dropped",
       x = "VIF", y = NULL, fill = NULL)

save_fig(plot_vif, "f4_4_vif", height = 6.5)

######################################################################
# 4.5 confounding with season
# patrols run mainly Sep-Apr and campaigns only Nov-Feb, so the crude 'patrolled' and
# 'campaign' effects in Figure 3.5 partly measure the season. Compare within season.
coverage <- model_df %>%
  group_by(month_lab) %>%
  summarise(
    patrolled_share = mean(patrolled),
    campaign_share  = mean(campaign_any),
    severe_per_m    = 1e6 * sum(Severe) / sum(visitors),
    .groups         = "drop"
  )

plot_coverage <- coverage %>%
  pivot_longer(c(patrolled_share, campaign_share), names_to = "measure", values_to = "share") %>%
  ggplot(aes(x = month_lab, y = share, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(patrolled_share = "#2471A3", campaign_share = "#F5B041"),
                    labels = c(campaign_share = "Campaign active", patrolled_share = "Patrolled")) +
  scale_y_continuous(labels = percent) +
  labs(title = "Figure 4.5a: Safety coverage by month vs death rate",
       x = NULL, y = "Share of beach-months", fill = NULL)

plot_coverage_rate <- ggplot(coverage, aes(x = month_lab, y = severe_per_m)) +
  geom_col(fill = sev_cols["Severe"]) +
  labs(x = NULL, y = "Deaths per million visits")

save_fig(plot_coverage / plot_coverage_rate, "f4_5a_seasonal_coverage", height = 6)

strat_patrol <- model_df %>%
  mutate(season_grp = ifelse(patrol_season == 1, "Sep-Apr", "May-Aug")) %>%
  group_by(season_grp, patrolled) %>%
  summarise(
    beach_months = n(),
    visitors_m   = sum(visitors) / 1e6,
    Severe       = sum(Severe),
    total_per_m  = sum(Total) / visitors_m,
    severe_per_m = Severe / visitors_m,
    .groups      = "drop"
  )

strat_campaign <- model_df %>%
  filter(cal_m %in% c(11, 12, 1, 2)) %>%
  group_by(campaign_any) %>%
  summarise(
    beach_months = n(),
    visitors_m   = sum(visitors) / 1e6,
    Severe       = sum(Severe),
    total_per_m  = sum(Total) / visitors_m,
    severe_per_m = Severe / visitors_m,
    .groups      = "drop"
  )

save_tab(strat_patrol,   "t14_patrolled_within_season")
save_tab(strat_campaign, "t15_campaign_within_nov_feb")
cat("\n== 4.5a Patrolled vs not, within season ==\n")
print(strat_patrol)
cat("\n== 4.5b Campaign vs none, Nov-Feb only ==\n")
print(strat_campaign)

# crude all-month comparisons, to show how much of each gap is season
crude_rates <- function(df, flag) {
  df %>%
    group_by(with_flag = .data[[flag]]) %>%
    summarise(severe_per_m = 1e6 * sum(Severe) / sum(visitors),
              total_per_m  = 1e6 * sum(Total) / sum(visitors), .groups = "drop")
}

strat_plot_df <- bind_rows(
  crude_rates(model_df, "patrolled") %>%
    mutate(comparison = "Patrolled\nall months (crude)"),
  strat_patrol %>%
    transmute(comparison = paste0("Patrolled\n", season_grp, " only"), with_flag = patrolled,
              severe_per_m, total_per_m),
  crude_rates(model_df, "campaign_any") %>%
    mutate(comparison = "Campaign\nall months (crude)"),
  strat_campaign %>%
    transmute(comparison = "Campaign\nNov-Feb only", with_flag = campaign_any,
              severe_per_m, total_per_m)
) %>%
  mutate(group      = ifelse(with_flag == 1, "With", "Without"),
         comparison = factor(comparison, levels = unique(comparison))) %>%
  pivot_longer(c(severe_per_m, total_per_m), names_to = "measure", values_to = "rate") %>%
  mutate(measure = ifelse(measure == "severe_per_m", "Deaths per million visits",
                          "All incidents per million visits"))

plot_strat <- ggplot(strat_plot_df, aes(x = comparison, y = rate, fill = group)) +
  geom_col(position = "dodge") +
  facet_wrap(~ measure, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(With = "#2471A3", Without = "grey60")) +
  labs(title = "Figure 4.5b: Patrol and campaign gaps, crude vs within season",
       subtitle = "If the gap shrinks within season, part of the crude effect is the season",
       x = NULL, y = NULL, fill = NULL)

save_fig(plot_strat, "f4_5b_within_season", width = 9, height = 6)

######################################################################
# 4.6 within- vs between-beach variation
# a driver that barely changes within a beach (e.g. pa_system, flags) is identified
# almost entirely by comparing DIFFERENT beaches, so its estimate can absorb any
# unmeasured beach-level risk. Read those coefficients with caution.
variation <- bind_rows(lapply(cand_vars, function(v) {
  x          <- as.numeric(model_df[[v]])
  beach_mean <- ave(x, model_df$beach_nm)
  tibble(
    var              = v,
    within_share     = 1 - var(beach_mean) / var(x),
    beaches_changing = model_df %>%
      group_by(beach_nm) %>%
      summarise(k = n_distinct(.data[[v]])) %>%
      with(sum(k > 1))
  )
})) %>%
  arrange(within_share)

save_tab(variation, "t16_within_beach_variation")
cat("\n== 4.6 Share of each predictor's variance that is within-beach ==\n")
print(variation, n = Inf)

plot_variation <- variation %>%
  mutate(var    = factor(var, levels = var),
         is_lever = ifelse(var %in% levers, "Safety lever", "Other driver")) %>%
  ggplot(aes(x = within_share, y = var, fill = is_lever)) +
  geom_col() +
  geom_text(aes(label = paste0(beaches_changing, " beaches")), hjust = -0.1, size = 2.6) +
  scale_x_continuous(labels = percent, limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = c("Safety lever" = "#2471A3", "Other driver" = "grey60")) +
  labs(title = "Figure 4.6: Share of each predictor's variation that is within beaches",
       subtitle = "Near 0% = effect only estimable by comparing beaches. Label = beaches where it changes",
       x = NULL, y = NULL, fill = NULL)

save_fig(plot_variation, "f4_6_within_beach_variation", height = 6.5)

######################################################################
# 4.7 time and beach effects the drivers do not capture
# base model: Poisson rate model with the reduced drivers + month

# (a) time: likelihood-ratio tests for a fiscal-year effect and the FY2020-21 dip flag
lrt_time <- bind_rows(lapply(c("Total", "Severe"), function(y) {
  base_vars <- c(keep_vars, "month_lab")
  base      <- glm(rate_formula(y, base_vars), family = poisson, data = model_df)
  add_fy    <- glm(rate_formula(y, c(base_vars, "factor(fy)")), family = poisson, data = model_df)
  add_dip   <- glm(rate_formula(y, c(base_vars, "dip_fy2020_21")), family = poisson, data = model_df)
  lr <- function(m1, label) {
    a <- anova(base, m1, test = "LRT")
    tibble(sev = y, added = label, df = a$Df[2], deviance_drop = a$Deviance[2], p_value = a$`Pr(>Chi)`[2])
  }
  bind_rows(lr(add_fy, "fiscal year"), lr(add_dip, "FY2020-21 dip flag"))
}))
save_tab(lrt_time, "t17_lrt_time_effects")
cat("\n== 4.7a LR tests: does a year effect remain after the drivers? ==\n")
print(lrt_time)

# a year effect is easier to judge from the actual vs expected ratio for each year
year_ae <- bind_rows(lapply(c("Total", "Severe"), function(y) {
  base <- glm(rate_formula(y, c(keep_vars, "month_lab")), family = poisson, data = model_df)
  model_df %>%
    mutate(expected = fitted(base)) %>%
    group_by(fy) %>%
    summarise(actual = sum(.data[[y]]), expected = sum(expected), .groups = "drop") %>%
    mutate(sev = y, a_over_e = actual / expected,
           lower = qchisq(0.025, 2 * actual) / 2 / expected,          # exact Poisson 95% CI
           upper = qchisq(0.975, 2 * (actual + 1)) / 2 / expected)
}))

lrt_label <- lrt_time %>%
  filter(added == "fiscal year") %>%
  transmute(sev, label = sprintf("LR test for year: p = %s", format.pval(p_value, digits = 2)))

plot_year <- year_ae %>%
  left_join(lrt_label, by = "sev") %>%
  mutate(sev = paste0(sev, " (", label, ")")) %>%
  ggplot(aes(x = fy, y = a_over_e)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  annotate("rect", xmin = 2019.5, xmax = 2021.5, ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_pointrange(aes(ymin = lower, ymax = upper), colour = "#2471A3") +
  facet_wrap(~ sev, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = 2017:2026) +
  labs(title = "Figure 4.7a: Actual / expected by fiscal year after the drivers",
       subtitle = "Shaded = FY2020-21 dip. Points away from 1 = a year effect the drivers miss",
       x = "Fiscal year", y = "Actual / expected")

save_fig(plot_year, "f4_7a_year_effects", height = 6)

# (b) beaches: actual vs expected by beach. A fixed beach factor cannot be fitted for
# Severe (many beaches have zero deaths, so their coefficients run off to -Inf), so test
# beach heterogeneity with a Pearson chi-square on beach-level A vs E instead.
# Large X^2 / df => unmeasured beach risk => beach random effects or credibility weighting
beach_ae <- bind_rows(lapply(c("Total", "Severe"), function(y) {
  base <- glm(rate_formula(y, c(keep_vars, "month_lab")), family = poisson, data = model_df)
  model_df %>%
    mutate(expected = fitted(base)) %>%
    group_by(beach_nm) %>%
    summarise(actual = sum(.data[[y]]), expected = sum(expected), .groups = "drop") %>%
    mutate(sev = y, a_over_e = actual / expected, .before = 1)
}))

beach_het <- beach_ae %>%
  group_by(sev) %>%
  summarise(
    beaches        = n(),
    pearson_x2     = sum((actual - expected)^2 / expected),
    x2_over_df     = pearson_x2 / (beaches - 1),
    p_value        = pchisq(pearson_x2, df = beaches - 1, lower.tail = FALSE),
    worst_beach    = beach_nm[which.max(abs(actual - expected) / sqrt(expected))],
    worst_a_over_e = a_over_e[which.max(abs(actual - expected) / sqrt(expected))],
    .groups        = "drop"
  )

save_tab(beach_ae,  "t18_beach_actual_vs_expected")
save_tab(beach_het, "t19_beach_heterogeneity")
cat("\n== 4.7b Beach heterogeneity after the drivers (beach-level A vs E) ==\n")
print(beach_het, width = Inf)

# funnel plot: with no beach effect, beaches fall inside the 95% band around 1
het_label <- beach_het %>%
  transmute(sev, label = sprintf("%s (X^2/df = %.1f)", sev, x2_over_df))

beach_ae_plot <- beach_ae %>%
  left_join(het_label, by = "sev") %>%
  mutate(z = (actual - expected) / sqrt(expected))

funnel <- beach_ae_plot %>%
  group_by(label) %>%
  reframe(expected = exp(seq(log(min(expected)), log(max(expected)), length.out = 100))) %>%
  mutate(lower = pmax(1 - 1.96 / sqrt(expected), 0),
         upper = 1 + 1.96 / sqrt(expected))

plot_beach_ae <- ggplot(beach_ae_plot, aes(x = expected, y = a_over_e)) +
  geom_ribbon(data = funnel, aes(x = expected, ymin = lower, ymax = upper),
              inherit.aes = FALSE, fill = "grey85") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(colour = abs(z) > 1.96), alpha = 0.8) +
  geom_text(data = beach_ae_plot %>% group_by(label) %>% slice_max(abs(z), n = 4),
            aes(label = beach_nm), size = 2.5, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~ label, scales = "free") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.25))) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#C0392B"), guide = "none") +
  labs(title = "Figure 4.7b: Actual / expected by beach after the drivers",
       subtitle = "Grey band = 95% range if there were no beach effect; red = outside it",
       x = "Expected people, 10 years (log scale)", y = "Actual / expected")

save_fig(plot_beach_ae, "f4_7b_beach_actual_vs_expected", width = 10)

######################################################################
# 4.8 sample size for the severe model
# rule of thumb: >= 10 events per estimated parameter
n_severe  <- sum(model_df$Severe)
n_params  <- length(coef(glm(rate_formula("Severe", c(keep_vars, "month_lab")), family = poisson, data = model_df)))
cat(sprintf("\n== 4.8 Severe model: %d deaths / %d parameters (reduced set + month) = %.1f deaths per parameter ==\n",
            n_severe, n_params, n_severe / n_params))
if (n_severe / n_params < 10) {
  cat("Below 10: the severe model needs fewer terms, or should borrow strength from the\n")
  cat("more frequent CloseCall / Total counts (e.g. severe = total rate x severe share).\n")
}

######################################################################
# 4.9 variables that must NOT be used as predictors
cat("\n== 4.9 Excluded from modelling ==\n")
cat("- response_min, n_lifeguards, drone_present and all person fields: recorded after\n")
cat("  the incident, not known in advance (leakage). Descriptive only (Figure 3.7).\n")
cat("- icecream_visits, n_bins: proxies for crowd size and heat, not drivers of risk\n")
cat("  (correlations in 4.3). Dropped: visitors is already the offset and temperature\n")
cat("  is measured directly, so they would only blur the real effects.\n")
cat("- visitors: the exposure. Enters as offset(log(visitors)), not as a predictor.\n")

cat("\nSaved EDA and pre-modelling outputs to", eda_dir, "\n")


######################################################################
# PRE-MODELLING NOTES (from the run on data to Aug 2026 - re-check if the data change)
######################################################################
# 1. Distribution (Fig 4.1): Mild / Total over-dispersed after the drivers (Pearson
#    ~2.1-2.3, NB AIC ~3,500-4,100 lower) -> negative binomial. Severe is close to
#    Poisson (Pearson ~1.2, NB improves AIC by only ~2). Observed zeros are within ~3%
#    of what the fitted models expect -> no zero-inflation needed.
# 2. Exposure (Fig 4.2): with all drivers the log(visitors) slope is ~0.86 for Total and
#    1 is inside the CI for Severe -> offset(log(visitors)) is reasonable for Severe; a
#    sensitivity with log(visitors) as a free covariate is worth reporting for Mild / Total.
# 3. Multicollinearity (Fig 4.3, 4.4): with patrolled and the crowd proxies removed,
#    the only VIF over 5 is the swell pair (rho 0.97) -> the automatic pruning drops
#    swell_avg_m. Every remaining VIF is <= 3.3, incl. log_lg_hrs at ~1.7.
# 4. Confounding (Fig 4.5): patrols run Sep-Apr and campaigns Nov-Feb only, so crude
#    effects include the season. Always adjust for month / season; within season the
#    patrolled and campaign gaps shrink but remain (t14, t15).
# 5. Identification (Fig 4.6): pa_system and flags change at only 6 beaches each, and
#    the static beach features never change -> their effects come from comparing beaches
#    and may absorb unmeasured beach risk. Treat as weaker evidence than the lifeguard /
#    drone / campaign effects, which vary within beaches.
# 6. Beach and time effects (Fig 4.7): strong unexplained beach heterogeneity for Total
#    (X^2/df ~22) but none for Severe (~1.0) -> beach random effects or credibility
#    weighting for 'where'. A year effect remains for Total but not Severe; the
#    FY2020-21 dip is explained by the drivers.
# 7. Severe has ~6.5 deaths per parameter (reduced set + month) -> keep the severe model lean.