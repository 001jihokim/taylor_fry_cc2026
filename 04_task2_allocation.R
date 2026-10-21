######################################################################
# 04_task2_allocation.R | Taylor Fry 2026 Case Competition - Task 2
# Uses Task 1 FY2027 future inputs + fitted severe model, and pricing.csv.
# Selects a greedy, budget-feasible portfolio, plus restricted alternatives.
# IMPORTANT: predictive associations are NOT causal treatment effects.
# Run: source('02_modelling_updated.R') once, then source('04_task2_allocation.R')
######################################################################

stopifnot(getRversion() >= '4.1.0')
needed <- c('dplyr', 'tidyr', 'readr', 'ggplot2')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install packages: ', paste(missing, collapse=', '))
library(dplyr)
library(tidyr)
library(readr)

BUDGET <- 2000000
model_dir <- 'temp_output/modelling'
out_dir <- 'temp_output/task2'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

files <- c('05_FY2027_assumed_inputs.csv', '06_fitted_models.rds',
           '07_FY2027_beach_month_forecast.csv')
missing_paths <- file.path(model_dir, files)[!file.exists(file.path(model_dir,files))]
if (length(missing_paths)) stop('Run 02_modelling_updated.R first; missing: ',
                               paste(missing_paths, collapse=', '))
price_candidates <- c('data/pricing.csv', 'pricing.csv', '/mnt/data/pricing.csv')
price_path <- price_candidates[file.exists(price_candidates)][1]
if (is.na(price_path)) stop('Pricing CSV missing: put pricing.csv in data/ or project root.')
pricing <- readr::read_csv(price_path, show_col_types=FALSE) %>%
  filter(!is.na(item), nzchar(trimws(item))) %>% mutate(item=trimws(item))
price_of <- function(item) {
  x <- pricing$cost_fd[match(item, pricing$item)]
  if (length(x) != 1 || is.na(x)) stop('Missing pricing item: ', item)
  as.numeric(x)
}
COST <- c(
  drone_hardware = price_of('Drone hardware'),
  drone_training = price_of('Drone operator training'),
  drone_hour = price_of('Drone operator (casual)'),
  casual_lg_hour = price_of('Casual SLSC lifeguard'),
  fulltime_lg_year = price_of('Full-time SLSC lifeguard'),
  fulltime_recruit = price_of('Full-time recruitment'),
  campaign = price_of('Beach safety education campaign'),
  community_training = price_of('Community first aid & rescue training'),
  pa_install = price_of('PA system'),
  signage = price_of('Rip-warning flags/signage')
)
# Fixed physical effects and implementation conditions from pricing.csv:
# full-time: +1750 hours/year, spread evenly; drone capacity 120 hours/month;
# campaign Nov-Feb sets index to MAX(current,0.7); first aid +10, capped at
# 100, effective next month; PA maintenance $2,000/year; drone hardware 3yr.
PA_MAINTENANCE <- 2000
MAX_DRONE_HOURS_PER_UNIT <- 120

future <- readr::read_csv(
  file.path(model_dir, files[1]), show_col_types=FALSE,
  col_types=readr::cols(.default=readr::col_guess(), postcode=readr::col_character())
) %>%
  mutate(date=as.Date(date), beach_nm=as.character(beach_nm),
         postcode=as.character(postcode),
         season=factor(as.character(season),levels=c('Spring','Summer','Autumn','Winter'))) %>%
  arrange(beach_nm,date)
if (nrow(future) != 960 || anyNA(future$lg_hrs) || anyNA(future$drone_hrs))
  stop('Expected 960 FY2027 beach-months with complete safety-resource inputs.')
final_models <- readRDS(file.path(model_dir,files[2]))
if (is.null(final_models$Severe)) stop('No final Severe model in 06_fitted_models.rds')
severe_obj <- final_models$Severe
if (!severe_obj$kind %in% c('Poisson','NegBin','NegBin_GLMM','TwoStage'))
  stop('Unsupported Severe model type: ',severe_obj$kind)

predict_severe <- function(df) {
  # The same prediction logic as Task 1, without sourcing / refitting Task 1.
  df$log_lg_hrs <- log1p(df$lg_hrs)
  df$log_drone_hrs <- log1p(df$drone_hrs)
  df$log_visitors <- log(df$visitors)
  if (severe_obj$kind == 'TwoStage') {
    freq <- predict(severe_obj$fit$frequency,newdata=df,type='response')
    share <- predict(severe_obj$fit$share,newdata=df,type='response')
    result <- as.numeric(freq) * as.numeric(share)
  } else {
    result <- as.numeric(predict(severe_obj$fit,newdata=df,type='response'))
  }
  if (any(!is.finite(result)) || any(result<0))
    stop('Nonfinite/negative severe predictions')
  result
}

baseline <- predict_severe(future)
reference <- readr::read_csv(file.path(model_dir,files[3]),show_col_types=FALSE) %>%
  arrange(beach_nm,date)
if (length(baseline) != nrow(reference) ||
    max(abs(baseline-reference$pred_Severe)) > 1e-5)
  stop('Task 2 baseline does not match Task 1 predictions. Check input/model alignment.')

# Every candidate is a purchasable BUNDLE with its exact months & cost.
# A new drone is paired with one trained operator before any flight takes place.
# No ownership is assumed for beaches without historical/baseline drone hours.
# For beaches with current drone hours we assume existing hardware + trained
# operator and enough free capacity up to 120 hours per month (a limitation).
# Summer = Nov-Feb. Monthly casual blocks = 40 additional hours, Sep-Apr.
# Campaign purchase applies to ALL beaches within its postcode.
# Community training purchased Sep 2026, benefits Oct 2026 onward.
# No more than one drone package (peak / year-round) per beach.

catalog <- list()
append_bundle <- function(category, name, idx, variable, mode, increment,
                          annual_cost=0, oneoff_cost=0, beach=NA_character_,
                          postcode=NA_character_, exclusive=NA_character_) {
  idx <- sort(unique(as.integer(idx)))
  if (!length(idx) || annual_cost+oneoff_cost <= 0) return(invisible(NULL))
  catalog[[length(catalog)+1L]] <<- list(
    id=sprintf('A%04d', length(catalog)+1L), category=category,
    name=name, beach=beach, postcode=postcode,
    idx=idx, variable=variable, mode=mode, increment=increment,
    annual_cost=as.numeric(annual_cost), oneoff_cost=as.numeric(oneoff_cost),
    cost=as.numeric(annual_cost+oneoff_cost), exclusive=exclusive
  )
  invisible(NULL)
}
for (b in unique(future$beach_nm)) {
  all_i <- which(future$beach_nm == b)
  pc <- future$postcode[all_i[1]]
  peak_i <- all_i[future$cal_m[all_i] %in% c(11,12,1,2)]
  patrol_i <- all_i[future$cal_m[all_i] %in% c(9,10,11,12,1,2,3,4)]
  for (i in patrol_i) {
    append_bundle('Casual lifeguards',paste0(b, ' | ',format(future$date[i],'%b %Y'), ' +40h'),
                  i,'lg_hrs','add',40, annual_cost=40*COST['casual_lg_hour'],
                  beach=b,postcode=pc)
  }
  append_bundle('Casual lifeguards',paste0(b,' | Nov-Feb +40h/month'),
                peak_i,'lg_hrs','add',40,
                annual_cost=length(peak_i)*40*COST['casual_lg_hour'],beach=b,postcode=pc)
  append_bundle('Full-time lifeguards',paste0(b,' | +1 FTE year-round'),
                all_i,'lg_hrs','add',1750/12,
                annual_cost=COST['fulltime_lg_year'],oneoff_cost=COST['fulltime_recruit'],
                beach=b,postcode=pc)
  if (future$pa_system[all_i[1]] == 0 && all(future$pa_system[all_i] == 0)) {
    append_bundle('PA systems',paste0(b,' | install PA'),all_i,
                  'pa_system','max',1,annual_cost=PA_MAINTENANCE,
                  oneoff_cost=COST['pa_install'],beach=b,postcode=pc)
  }
  if (all(future$flags[all_i]==0)) {
    append_bundle('Rip-warning signage',paste0(b,' | install warning flags'),all_i,
                  'flags','max',1,oneoff_cost=COST['signage'],beach=b,postcode=pc)
  }
  existing <- any(future$drone_hrs[all_i]>0)
  fixed_cost <- if (existing) 0 else COST['drone_hardware']+COST['drone_training']
  # Conservative one-drone-capacity simplification: do not assume spare
  # capacity when observed projected hours exceed the documented cap.
  for (pattern in c('Peak','Year-round')) {
    i <- if (pattern=='Peak') peak_i else all_i
    hours <- if (pattern=='Peak') 20 else 10
    if (any(future$drone_hrs[i]+hours > MAX_DRONE_HOURS_PER_UNIT)) next
    append_bundle('Drone surveillance',paste0(b,' | ',pattern,' +',hours,'h/month'),
                  i,'drone_hrs','add',hours,
                  annual_cost=length(i)*hours*COST['drone_hour'],
                  oneoff_cost=fixed_cost,beach=b,postcode=pc,
                  exclusive=paste0('Drone:',b))
  }
}
for (pc in unique(future$postcode)) {
  all_i <- which(future$postcode==pc)
  campaign_i <- all_i[future$cal_m[all_i] %in% c(11,12,1,2)]
  if (any(future$campaign_idx[campaign_i] < 0.7-1e-8)) {
    append_bundle('Education campaigns',paste0('Postcode ',pc,' | Nov-Feb education'),
                  campaign_i,'campaign_idx','max',0.7,
                  annual_cost=COST['campaign'],postcode=pc)
  }
  firstaid_i <- all_i[future$date[all_i] >= as.Date('2026-10-01')]
  if (any(future$first_aid_cap[firstaid_i] < 100-1e-8)) {
    append_bundle('Community first aid',paste0('Postcode ',pc,' | training from Oct'),
                  firstaid_i,'first_aid_cap','add_capped',10,
                  annual_cost=COST['community_training'],postcode=pc)
  }
}
stopifnot(length(catalog)>0)

apply_bundle <- function(df,a) {
  ii <- a$idx
  old <- df[[a$variable]][ii]
  if (a$mode=='add') new <- old+a$increment
  else if (a$mode=='max') new <- pmax(old,a$increment)
  else if (a$mode=='add_capped') new <- pmin(100,old+a$increment)
  else stop('Unrecognised change type')
  df[[a$variable]][ii] <- new
  df
}

# A sequential, budget-constrained greedy optimizer.
# Cache candidate benefits, and re-evaluate only actions whose beach-months
# overlap an earlier purchase. This respects nonlinear effects and interactions
# of the fitted models. It is NOT a proof of a globally optimal solution.
score_bundle <- function(df, predictions, a) {
  changed <- apply_bundle(df[a$idx,,drop=FALSE],modifyList(a,list(idx=seq_along(a$idx))))
  gained <- sum(predictions[a$idx]-predict_severe(changed))
  max(0, as.numeric(gained))
}

run_greedy <- function(allowed_categories, cap=BUDGET, keep_trace=FALSE) {
  allowed <- which(vapply(catalog,function(a) a$category %in% allowed_categories, logical(1)))
  df <- future
  p <- baseline
  scores <- rep(NA_real_,length(catalog))
  active <- rep(FALSE,length(catalog)); active[allowed] <- TRUE
  used_groups <- character()
  trace <- list()
  spent <- 0
  for (j in allowed) scores[j] <- score_bundle(df,p,catalog[[j]])
  repeat {
    available <- which(active & vapply(catalog,function(a) a$cost <= cap-spent+1e-7,logical(1)))
    if (!length(available)) break
    ratios <- vapply(available,function(j) scores[j]/catalog[[j]]$cost,numeric(1))
    best <- available[which.max(ratios)]
    if (!length(best) || !is.finite(scores[best]) || scores[best] <= 1e-10) break
    a <- catalog[[best]]
    previous <- sum(p[a$idx])
    df <- apply_bundle(df,a)
    p[a$idx] <- predict_severe(df[a$idx,,drop=FALSE])
    benefit <- previous-sum(p[a$idx])
    spent <- spent+a$cost
    trace[[length(trace)+1L]] <- tibble(
      action_id=a$id, category=a$category, action=a$name,
      beach_nm=a$beach, postcode=a$postcode,
      annual_cost=a$annual_cost, oneoff_cost=a$oneoff_cost,
      total_cost=a$cost, incremental_deaths_avoided=benefit,
      cumulative_cost=spent,
      cumulative_deaths_avoided=sum(baseline)-sum(p)
    )
    active[best] <- FALSE
    if (!is.na(a$exclusive)) {
      for (j in which(active)) {
        if (!is.na(catalog[[j]]$exclusive) && catalog[[j]]$exclusive==a$exclusive)
          active[j] <- FALSE
      }
    }
    # Only rescore bundles touching the changed future beach-months.
    for (j in which(active)) {
      if (any(catalog[[j]]$idx %in% a$idx))
        scores[j] <- score_bundle(df,p,catalog[[j]])
    }
  }
  list(actions=bind_rows(trace),future=df,pred=p,cost=spent,
       baseline=sum(baseline),after=sum(p),avoided=sum(baseline)-sum(p))
}

categories <- sort(unique(vapply(catalog,`[[`,character(1),'category')))
cat('\nTask 2: evaluating ',length(catalog),' purchasable bundles...\n',sep='')
# Evaluate the unrestricted greedy run and restricted approaches. A greedy
# method is not globally optimal; on occasion a restricted run can actually
# beat an unrestricted run because greedy decisions interact. Recommend the
# BEST of the tested feasible portfolios rather than presuming mix is best.
strategy_candidates <- list(
  'Unrestricted mix' = run_greedy(categories),
  'Lifeguards only' = run_greedy(c('Casual lifeguards','Full-time lifeguards')),
  'Drones only' = run_greedy('Drone surveillance'),
  'Infrastructure & education only' = run_greedy(c('PA systems','Rip-warning signage',
                                                    'Education campaigns','Community first aid'))
)
benefits <- vapply(strategy_candidates, function(x) x$avoided, numeric(1))
winning_idx <- which.max(benefits)
winning_strategy <- names(strategy_candidates)[winning_idx]
portfolio <- strategy_candidates[[winning_idx]]
names(strategy_candidates)[winning_idx] <- paste0('Recommended: ',winning_strategy)
cat('Recommended: ',winning_strategy, ' | ',nrow(portfolio$actions),
    ' actions | FD$',format(portfolio$cost,big.mark=','),
    ' allocated | model-implied avoided severe outcomes ',
    round(portfolio$avoided,3),'\n',sep='')

strategy_summary <- bind_rows(lapply(names(strategy_candidates), function(n) {
  x <- strategy_candidates[[n]]
  tibble(strategy=n, spending=x$cost, unspent=BUDGET-x$cost,
         baseline_severe=x$baseline, predicted_severe=x$after,
         model_implied_deaths_avoided=x$avoided,
         cost_per_model_implied_death_avoided=ifelse(x$avoided>1e-8,x$cost/x$avoided,NA_real_),
         actions=nrow(x$actions))
}))

# These REALISATION assumptions are illustrative scenario haircuts, NOT
# statistically identified causal treatment-effect estimates.
realisation <- c(0,0.25,0.5,1.0)
sensitivity <- bind_rows(lapply(realisation,function(frac) {
  tibble(realisation_fraction=frac,baseline_severe=sum(baseline),
         estimated_deaths_avoided=portfolio$avoided*frac,
         remaining_severe=sum(baseline)-portfolio$avoided*frac,
         spending=portfolio$cost,
         fd_per_avoided_death=ifelse(frac>0 && portfolio$avoided>0,
                                    portfolio$cost/(portfolio$avoided*frac),NA_real_))
}))

# Annual first-year costing: full purchase prices + 2027 maintenance.
# Drone hardware has a 3-year life; full acquisition cost is booked now.
actions <- portfolio$actions
if (!nrow(actions)) {
  stop('No candidate reduced predicted Severe outcomes under these fitted models. ',
       'Do not claim an evidence-based reduction. Review coefficient directions ',
       'and consider non-modelled policy criteria before allocating funding.')
}
category_spend <- actions %>% group_by(category) %>%
  summarise(actions=n(),annual_cost=sum(annual_cost),oneoff_cost=sum(oneoff_cost),
            total_cost=sum(total_cost),marginal_modelled_deaths_avoided=
              sum(incremental_deaths_avoided),.groups='drop') %>%
  arrange(desc(total_cost))
beach_alloc <- actions %>% filter(!is.na(beach_nm)) %>%
  group_by(beach_nm) %>% summarise(total_cost=sum(total_cost),
         actions=n(),marginal_modelled_deaths_avoided=sum(incremental_deaths_avoided),
         .groups='drop')
postcode_alloc <- actions %>% filter(!is.na(postcode)) %>%
  group_by(postcode) %>% summarise(total_cost=sum(total_cost),
         actions=n(),marginal_modelled_deaths_avoided=sum(incremental_deaths_avoided),
         .groups='drop')
scenario <- future %>%
  transmute(beach_nm,postcode,date,cal_m,visitors,
            baseline_severe=baseline,post_investment_severe=portfolio$pred,
            reduction=baseline_severe-post_investment_severe,
            reduction_at_50pct_realisation=0.5*reduction)

# Output tables used by the figures.
write_csv(actions,file.path(out_dir,'01_recommended_actions.csv'))
write_csv(category_spend,file.path(out_dir,'02_spending_by_initiative.csv'))
write_csv(beach_alloc,file.path(out_dir,'03_spending_by_beach.csv'))
write_csv(postcode_alloc,file.path(out_dir,'04_spending_by_postcode.csv'))
write_csv(scenario,file.path(out_dir,'05_monthly_beach_counterfactual.csv'))
write_csv(strategy_summary,file.path(out_dir,'06_alternative_strategies.csv'))
write_csv(sensitivity,file.path(out_dir,'07_causal_realisation_sensitivity.csv'))
write_csv(actions %>% select(action_id,category,action,total_cost,cumulative_cost,
                              cumulative_deaths_avoided),
          file.path(out_dir,'08_greedy_spend_path.csv'))
write_csv(tibble(budget=BUDGET,spent=portfolio$cost,remaining=BUDGET-portfolio$cost,
                 baseline_deaths=sum(baseline),modelled_after=portfolio$after,
                 model_implied_avoided=portfolio$avoided,
                 half_realisation_avoided=portfolio$avoided*0.5,
                 selected_severe_model=severe_obj$kind,
                 selected_strategy=winning_strategy,
                 decision_method='Best of unrestricted and 3 restricted greedy portfolios'),
          file.path(out_dir,'09_task2_executive_summary.csv'))

# Write a complete, data-driven recommendation note with actual numbers, not
# numbers invented in advance of running the Task 1 models.
fmt <- function(x,dp=2) format(round(x,dp),big.mark=',',nsmall=dp,trim=TRUE)
top <- head(category_spend,4)
lines <- c(
  '# Task 2 — additional FD$2 million allocation (FY2027)',
  '',
  paste0('**Recommended allocation:** FD$',fmt(portfolio$cost,0),' across ',
         nrow(actions),' actions; unallocated FD$',fmt(BUDGET-portfolio$cost,0),'.'),
  paste0('**Baseline forecast:** ',fmt(sum(baseline)), ' severe outcomes; **post-investment model scenario:** ',
         fmt(portfolio$after),'.'),
  paste0('**Model-implied reduction:** ',fmt(portfolio$avoided),
         ' severe outcomes. This is NOT a validated causal estimate.'),
  paste0('**Illustrative 50% effect realisation:** ',fmt(0.5*portfolio$avoided),
         ' avoided severe outcomes.'),
  '',
  '## Funding mix',
  if(nrow(top)) paste0('- ',top$category,': FD$',fmt(top$total_cost,0),
                     ' (',top$actions,' actions)') else '- No positive modeled-benefit actions selected.',
  '',
  '## Comparison with alternatives',
  paste0('- ',strategy_summary$strategy,': FD$',fmt(strategy_summary$spending,0),
         ' spent; ',fmt(strategy_summary$model_implied_deaths_avoided),
         ' model-implied avoided severe outcomes.'),
  '',
  '## Why this allocation?',
  paste0('The selected strategy is **',winning_strategy,'**. Four portfolios were compared. Each greedily purchases eligible interventions in descending order of the *current* modelled marginal reduction in severe incidents per FD$, recalculating overlapping benefits after each purchase. The highest-benefit feasible result is recommended. Constraints include costs, postcode application, effective dates, drone prerequisites and the first-year FD$2m cap. This is a transparent heuristic, not proof of a global optimum.'),
  '',
  '## Assumptions and limitations',
  '- Existing FY2026 coverage is continued in the Task 1 FY2027 baseline; future visitation and weather are assumed, not observed.',
  '- Extra casual patrols are offered in 40-hour monthly blocks Sep–Apr, with a Nov–Feb block option; only one FTE per beach and one additional drone package per beach are considered. This candidate menu is not exhaustive.',
  '- Full-time patrol hours are spread uniformly over 12 months, as specified in pricing; real rosters may differ.',
  '- Drone hardware + operator training are purchased for any beach with no baseline drone hours. Existing drone flights are taken as evidence of existing hardware and trained staffing; 120-hour/month per-drone capacity is assumed.',
  '- First-aid training is purchased in September and only affects October onward; campaigns affect November–February at postcode level.',
  '- PA installs include FD$2,000 maintenance in FY2027; drone hardware is charged in full in its purchase year, despite a three-year life.',
  '- Static PA and signage coefficients may be biased by differences across beaches. Likewise, patrol/drone hours may be targeted toward riskier beaches. Modelled counterfactuals must not be described as proven lives saved.',
  '- The zero, 25%, 50% and 100% effect-realisation cases are illustrative stress tests, not confidence intervals or formally estimated causal treatment effects.',
  '- Greedy cost-effectiveness ranking may differ from a true global mixed-integer optimum. Intervention impacts may not generalise outside the historical range.',
  '- The 2 million is treated as first-year cash funding; ongoing cost in later years needs a separate sustainability check.',
  '- Individual incident severity is person-level: a multi-person event can contribute multiple outcomes.',
  '',
  '## Implementation and measurement',
  'Confirm operational rosters, local permissions and trained drone personnel; implement contracted and fixed-cost initiatives before their active months; track intervention delivery, visitors and incident severity monthly; compare outcomes with matched control sites where feasible, and update allocations based on measured results.',
  '',
  'See 01_recommended_actions.csv, 06_alternative_strategies.csv and figures for detailed evidence.'
)
writeLines(lines,file.path(out_dir,'Task2_recommendation.md'))
cat('\nSaved Task 2 outputs to: ',normalizePath(out_dir),'\n',sep='')
cat('NOTE: model-implied reduction != causal number of lives saved.\n')
