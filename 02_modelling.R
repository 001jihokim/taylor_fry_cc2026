# =====================================================================
# 02_modelling.R
# UNSW x TAYLOR FRY T3 2026 Case Competition - Omega Consulting
# Members - Jiho Kim, Becky Thai, Christina Lou, Timothy Bae
#
# Part 2 of 2: model validation, final models, drivers, 2026/27 forecast and scenarios.
# Reads:  temp_output/clean/*.rds (from 01_setup_eda_cleaning.R)
# Writes: temp_output/models/
#
# Derived from main.R (kept unchanged as the original). Differences from main.R:
#   1. Data come from temp_output/clean/model_df.rds instead of being rebuilt here,
#      so EDA and models use one definition of the data.
#   2. Output csv is written to temp_output/models/ instead of the working directory.
# Everything from section 2 onward is main.R's logic, unchanged.
# NB theta is estimated by ML; the boosted model is xgboost.
# =====================================================================

# ---------------------------------------------------------------------
# 0. SETUP (same paths as 01_setup_eda_cleaning.R)
# ---------------------------------------------------------------------
out_root  <- "temp_output"
clean_dir <- file.path(out_root, "clean")
model_dir <- file.path(out_root, "models")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(clean_dir, "model_df.rds")))
  stop("Cleaned data not found. Run 01_setup_eda_cleaning.R first.", call. = FALSE)

# Packages are installed by 01_setup_eda_cleaning.R. MASS is called as MASS:: (it masks dplyr::select).
library(dplyr); library(tidyr); library(xgboost)
set.seed(2026)

# ---------------------------------------------------------------------
# 1. DATA: one row per beach-month with counts of PEOPLE by severity (built in 01_setup_eda_cleaning.R)
# ---------------------------------------------------------------------
dat <- readRDS(file.path(clean_dir, "model_df.rds")) %>%
  mutate(marine_health = as.character(marine_health),   # plain factor in formulas (as in main.R)
         postcode      = as.integer(as.character(postcode)),
         mh            = c(L = 0, M = 1, H = 2)[marine_health])
inc <- readRDS(file.path(clean_dir, "incidents_clean.rds"))
stopifnot(nrow(dat) == 9600, sum(dat$Severe) == sum(inc$inc_type == "Severe"))

# NOT used as predictors (leakage / not known at forecast time): response_min and all
# person-level fields in incidents.csv (age, sex, swim_level, alone, drone_present, n_lifeguards).

# ---------------------------------------------------------------------
# 2. QUICK CHECKS (full EDA is in 01_setup_eda_cleaning.R)
# ---------------------------------------------------------------------
cat("\nPeople by fiscal year\n")
print(dat %>% group_by(fy) %>% summarise(visitors_m = sum(visitors)/1e6, Mild = sum(Mild),
                                         CloseCall = sum(CloseCall), Severe = sum(Severe)))
cat("\nSevere deaths per million visits, patrolled vs not\n")
print(dat %>% group_by(patrolled) %>% summarise(sev_per_m = 1e6*sum(Severe)/sum(visitors),
                                                close_per_m = 1e6*sum(CloseCall)/sum(visitors)))
cat("\nBy calendar month\n")
print(dat %>% group_by(cal_m) %>% summarise(visitors_m = sum(visitors)/1e6, Mild = sum(Mild), Severe = sum(Severe)))

# ---------------------------------------------------------------------
# 3. MODEL FORMULAS (rate model: counts ~ drivers with log(visitors) offset)
# ---------------------------------------------------------------------
RED  <- "patrolled + log1p(lg_hrs) + log1p(drone_hrs) + campaign_idx + I(first_aid_cap/10) +
         pa_system + flags + I(equip_rating/10) + swell_avg_m + water_temp + tourist_ratio +
         surfer_prop + n_rips + fishing + I(n_bins/10) + rockpool"
FULL <- paste(RED, "+ swell_max_m + wind_days + air_temp_hi + I(rain_mm/100) + I(turbidity/10) +
         icecream_visits + factor(marine_health) + I(shore_len_m/1000) + rocky")
f_of <- function(y, rhs) as.formula(paste(y, "~", rhs, "+ offset(log(visitors))"))
feats <- c("visitors","shore_len_m","n_rips","air_temp_hi","swell_avg_m","swell_max_m","wind_days","water_temp",
           "rain_mm","turbidity","tourist_ratio","lg_hrs","drone_hrs","campaign_idx","first_aid_cap","pa_system",
           "flags","patrolled","equip_rating","rocky","fishing","n_bins","rockpool","icecream_visits",
           "surfer_prop","mh","cal_m")

pdev <- function(y, mu) { mu <- pmax(mu, 1e-9); 2 * sum(ifelse(y > 0, y * log(pmax(y, 1) / mu), 0) - (y - mu)) }

fit_xgb <- function(tr, y) {
  d <- xgb.DMatrix(as.matrix(tr[, feats]), label = tr[[y]], base_margin = log(tr$visitors))
  xgb.train(params = list(objective = "count:poisson", max_depth = 3, eta = 0.05,
                          subsample = 0.8, colsample_bytree = 0.8, seed = 1), data = d, nrounds = 200, verbose = 0)
}
pred_xgb <- function(m, te) predict(m, xgb.DMatrix(as.matrix(te[, feats]), base_margin = log(te$visitors)))

# ---------------------------------------------------------------------
# 4. FORWARD (time-based) VALIDATION: train on years < t, test on t = 2024, 2025, 2026
# ---------------------------------------------------------------------
res <- list()
for (y in c("Mild", "CloseCall", "Severe")) for (ty in 2024:2026) {
  tr <- dat %>% filter(fy < ty); te <- dat %>% filter(fy == ty)
  P <- list()
  rate <- tr %>% group_by(beach_nm) %>% summarise(r = sum(.data[[y]]) / sum(visitors))
  P[["Beach rate (naive)"]] <- rate$r[match(te$beach_nm, rate$beach_nm)] * te$visitors
  P[["GLM full"]]    <- predict(glm(f_of(y, FULL), poisson, tr), te, type = "response")
  P[["GLM reduced"]] <- predict(glm(f_of(y, RED),  poisson, tr), te, type = "response")
  P[["NB reduced"]]  <- predict(MASS::glm.nb(f_of(y, RED), data = tr), te, type = "response")
  P[["GBM (xgboost Poisson)"]] <- pred_xgb(fit_xgb(tr, y), te)
  for (k in names(P)) {
    mu <- P[[k]]; by <- data.frame(a = te[[y]], p = mu, b = te$beach_nm) %>% group_by(b) %>% summarise(a = sum(a), p = sum(p))
    res[[length(res) + 1]] <- data.frame(sev = y, model = k, test_year = ty, deviance = pdev(te[[y]], mu) / nrow(te),
                                         rmse_bm = sqrt(mean((te[[y]] - mu)^2)), rmse_beach_year = sqrt(mean((by$a - by$p)^2)),
                                         actual = sum(te[[y]]), pred = sum(mu))
  }
}
cvres <- bind_rows(res)
cat("\nOut-of-sample scores (mean of 3 test years; lower deviance/RMSE is better)\n")
print(cvres %>% group_by(sev, model) %>% summarise(deviance = mean(deviance), rmse_bm = mean(rmse_bm),
                                                   rmse_beach_year = mean(rmse_beach_year), A_over_E = sum(actual) / sum(pred), .groups = "drop") %>%
        arrange(sev, deviance), n = 20)

# ---------------------------------------------------------------------
# 5. FINAL MODELS on all data + diagnostics + AIC/BIC
# ---------------------------------------------------------------------
M <- list(); for (y in c("Mild", "CloseCall", "Severe")) {
  M[[y]] <- list(pois = glm(f_of(y, RED), poisson, dat), full = glm(f_of(y, FULL), poisson, dat),
                 nb = MASS::glm.nb(f_of(y, RED), data = dat))
  cat(sprintf("\n%s | Pearson dispersion %.2f | AIC full %.0f, reduced %.0f, NB %.0f | BIC reduced %.0f, NB %.0f\n", y,
              sum(residuals(M[[y]]$pois, "pearson")^2) / df.residual(M[[y]]$pois),
              AIC(M[[y]]$full), AIC(M[[y]]$pois), AIC(M[[y]]$nb), BIC(M[[y]]$pois), BIC(M[[y]]$nb)))
}
# Actual vs fitted by fiscal year (check for bias)
print(cbind(dat %>% group_by(fy) %>% summarise(Mild = sum(Mild), CC = sum(CloseCall), Sev = sum(Severe)),
            fit_Mild = tapply(fitted(M$Mild$pois), dat$fy, sum), fit_Sev = tapply(fitted(M$Severe$pois), dat$fy, sum)))

# ---------------------------------------------------------------------
# 6. DRIVERS (severe model): rate ratios with 95% CI
# ---------------------------------------------------------------------
S <- M$Severe$pois; b <- coef(S); V <- vcov(S)
rr <- function(contrast) { # contrast = named vector over coefficient names
  v <- setNames(rep(0, length(b)), names(b)); v[names(contrast)] <- contrast
  est <- sum(v * b); se <- sqrt(drop(t(v) %*% V %*% v)); exp(c(rr = est, lo = est - 1.96 * se, hi = est + 1.96 * se)) }
drivers <- rbind(
  "Non-local share (+1 SD)"      = rr(c(tourist_ratio = sd(dat$tourist_ratio))),
  "Average swell (+1 SD)"        = rr(c(swell_avg_m = sd(dat$swell_avg_m))),
  "Fishing spot (yes vs no)"     = rr(c(fishing = 1)),
  "Rip currents (+1 SD)"         = rr(c(n_rips = sd(dat$n_rips))),
  "Surfer share (+1 SD)"         = rr(c(surfer_prop = sd(dat$surfer_prop))),
  "Rip flags (yes vs no)"        = rr(c(flags = 1)),
  "Campaign 0.7 vs none"         = rr(c(campaign_idx = 0.7)),
  "Drones 60 h/month vs none"    = rr(setNames(log(61), "log1p(drone_hrs)")),
  "Lifeguards 600 h/month vs none" = rr(c(patrolled = 1, setNames(log(601), "log1p(lg_hrs)"))))
cat("\nSevere-incident rate ratios\n"); print(round(drivers, 2))

# ---------------------------------------------------------------------
# 7. FORECAST FRAME for 2026/27: last year's safety levels, 2024-26 average weather, visitors +2.5%
# ---------------------------------------------------------------------
wx <- c("air_temp_hi","swell_avg_m","swell_max_m","wind_days","water_temp","rain_mm","turbidity",
        "tourist_ratio","surfer_prop","icecream_visits")
avg <- dat %>% filter(fy >= 2024) %>% group_by(beach_nm, cal_m) %>% summarise(across(all_of(wx), mean), .groups = "drop")
fc <- dat %>% filter(fy == 2026) %>% select(-all_of(wx)) %>% left_join(avg, by = c("beach_nm", "cal_m")) %>%
  mutate(visitors = visitors * 1.025, fy = 2027)

# ---------------------------------------------------------------------
# 8. POINT FORECAST + INTERVALS (coefficient uncertainty + negative binomial counts)
# ---------------------------------------------------------------------
fc$gbm_Mild      <- pred_xgb(fit_xgb(dat, "Mild"), fc)
fc$gbm_CloseCall <- pred_xgb(fit_xgb(dat, "CloseCall"), fc)
sim_total <- function(m, newd, nsim = 1000, theta = NULL) {
  X <- model.matrix(delete.response(terms(m)), newd); off <- log(newd$visitors)
  B <- MASS::mvrnorm(nsim, coef(m), vcov(m))
  sapply(seq_len(nsim), function(j) { mu <- exp(drop(X %*% B[j, ]) + off)
  sum(if (is.null(theta)) rpois(length(mu), mu) else rnbinom(length(mu), size = theta, mu = mu)) })
}
out <- list()
for (y in c("Mild", "CloseCall", "Severe")) {
  m <- if (y == "Severe") M[[y]]$pois else M[[y]]$nb
  th <- if (y == "Severe") NULL else m$theta
  glm_pt <- sum(predict(m, fc, type = "response"))
  pt <- if (y == "Severe") glm_pt else sum(fc[[paste0("gbm_", y)]])      # GBM point for mild/close call
  sims <- sim_total(m, fc, theta = th) * pt / glm_pt                      # rescale range to the chosen point forecast
  out[[y]] <- c(point = pt, glm_point = glm_pt, p05 = quantile(sims, .05), p95 = quantile(sims, .95),
                p99 = quantile(sims, .99), TVaR95 = mean(sims[sims >= quantile(sims, .95)]))
}
cat("\n2026/27 forecast (people)\n"); print(round(do.call(rbind, out), 1))

# ---------------------------------------------------------------------
# 9. WHERE AND WHEN
# ---------------------------------------------------------------------
fc$mu_Severe <- predict(M$Severe$pois, fc, type = "response")
cat("\nTop 10 beaches by expected deaths\n")
print(fc %>% group_by(beach_nm) %>% summarise(exp_deaths = sum(mu_Severe), visits_m = sum(visitors)/1e6) %>%
        arrange(desc(exp_deaths)) %>% head(10))
cat("\nBy calendar month\n")
print(fc %>% group_by(cal_m) %>% summarise(mild = sum(gbm_Mild), close = sum(gbm_CloseCall), deaths = sum(mu_Severe)))

# ---------------------------------------------------------------------
# 10. SCENARIOS (GLM totals; relative to base)
# ---------------------------------------------------------------------
tot <- function(d) sapply(c("Mild", "CloseCall", "Severe"), function(y)
  sum(predict(if (y == "Severe") M[[y]]$pois else M[[y]]$nb, d, type = "response")))
b25 <- dat %>% filter(fy == 2025) %>% select(beach_nm, cal_m, camp25 = campaign_idx)
scen <- list(
  "Base"                         = fc,
  "Visitors +5%"                 = mutate(fc, visitors = visitors / 1.025 * 1.05),
  "Visitors flat"                = mutate(fc, visitors = visitors / 1.025),
  "Visitors -10%"                = mutate(fc, visitors = visitors / 1.025 * 0.9),
  "Non-local share +10%"         = mutate(fc, tourist_ratio = pmin(tourist_ratio * 1.1, 1)),
  "Swell +20%"                   = mutate(fc, swell_avg_m = swell_avg_m * 1.2),
  "Water temp +1C"               = mutate(fc, water_temp = water_temp + 1),
  "Campaigns at 2024/25 level"   = left_join(fc, b25, by = c("beach_nm", "cal_m")) %>% mutate(campaign_idx = camp25),
  "No campaigns"                 = mutate(fc, campaign_idx = 0),
  "Lifeguard hours -10%"         = mutate(fc, lg_hrs = lg_hrs * 0.9))
cat("\nScenarios (GLM expected people)\n"); print(round(t(sapply(scen, tot)), 1))

write.csv(fc[, c("beach_nm", "month", "visitors", "gbm_Mild", "gbm_CloseCall", "mu_Severe")],
          file.path(model_dir, "forecast_2026_27_by_beach_month.csv"), row.names = FALSE)
cat("\nSaved forecast to", model_dir, "\n")