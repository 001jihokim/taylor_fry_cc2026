######################################################################
# 05_task2_presentation_figures.R | 9 slide-ready ggplot figures
# Run AFTER: 04_task2_allocation.R
# Outputs: temp_output/task2/presentation_figures/*.png (16:9, 300 dpi)
# Blue-only palette, consistent with Task 1's 03_presentation_figures.R.
######################################################################

packages <- c('dplyr','tidyr','readr','ggplot2','scales')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install: ',paste(missing,collapse=', '))
library(dplyr); library(tidyr); library(readr); library(ggplot2); library(scales)

data_dir <- 'temp_output/task2'
fig_dir <- file.path(data_dir,'presentation_figures')
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
read_t2 <- function(name) {
  path <- file.path(data_dir,name)
  if(!file.exists(path)) stop('Missing ',path,'; run 04_task2_allocation.R first.')
  readr::read_csv(path,show_col_types=FALSE)
}
portfolio <- read_t2('01_recommended_actions.csv')
spend <- read_t2('02_spending_by_initiative.csv')
beach_spend <- read_t2('03_spending_by_beach.csv')
pc_spend <- read_t2('04_spending_by_postcode.csv')
months <- read_t2('05_monthly_beach_counterfactual.csv') %>% mutate(date=as.Date(date))
alts <- read_t2('06_alternative_strategies.csv')
sensitivity <- read_t2('07_causal_realisation_sensitivity.csv')
path <- read_t2('08_greedy_spend_path.csv')
summary <- read_t2('09_task2_executive_summary.csv')

BLUE <- list(navy='#103C68',dark='#155E9A',mid='#318AC5',sky='#71BCE6',
             light='#A8DBF4',ice='#DDF2FC',pale='#F2F9FD',ink='#213E59',
             muted='#64819A',grid='#E1EEF5',white='#FFFFFF')
colours <- c(BLUE$navy,BLUE$dark,BLUE$mid,BLUE$sky,BLUE$light,BLUE$ice)
slide_theme <- function() {
  theme_minimal(base_size=15,base_family='sans') + theme(
    plot.background=element_rect(fill='white',colour=NA),
    panel.background=element_rect(fill='white',colour=NA),
    plot.title=element_text(size=22,face='bold',colour=BLUE$navy,margin=margin(b=9)),
    plot.subtitle=element_text(size=13,colour=BLUE$muted,margin=margin(b=16)),
    plot.caption=element_text(size=10,colour=BLUE$muted,hjust=0,margin=margin(t=13)),
    axis.title=element_text(size=13,colour=BLUE$ink),
    axis.text=element_text(size=12,colour=BLUE$ink),
    axis.ticks=element_blank(),panel.grid.minor=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.major.x=element_line(colour=BLUE$grid,linewidth=0.4),
    legend.position='bottom',legend.title=element_blank(),
    legend.text=element_text(size=12,colour=BLUE$ink),
    plot.margin=margin(20,27,20,20))
}
save_slide <- function(p,filename) {
  ggplot2::ggsave(file.path(fig_dir,filename),p,width=13.333,height=7.5,
                 dpi=300,units='in',bg='white')
  cat('Saved ',filename,'\n',sep='')
}
fd <- scales::label_dollar(prefix='FD$',accuracy=1,scale_cut=cut_short_scale())
model_note <- 'Conditional model scenario, not proven causal reductions in deaths.'

# Figure 01 - spending mix
p1 <- ggplot(spend %>% mutate(category=reorder(category,total_cost)),
             aes(x=total_cost,y=category,fill=category)) +
  geom_col(width=.67,show.legend=FALSE) +
  geom_text(aes(label=paste0('FD$',scales::comma(round(total_cost)))),
            hjust=-.12,fontface='bold',size=4,colour=BLUE$ink) +
  scale_fill_manual(values=rep_len(colours,nrow(spend))) +
  scale_x_continuous(labels=fd,expand=expansion(mult=c(0,.25))) +
  labs(title='Where the additional budget goes',
       subtitle='FY2027 selected interventions | Acquisition and annual charges included',
       x='Allocated first-year funding',y=NULL,
       caption=paste0('Budget FD$2,000,000 | Allocated FD$',comma(round(summary$spent[1])),
                      ' | Unspent FD$',comma(round(summary$remaining[1])))) + slide_theme()
save_slide(p1,'01_budget_by_intervention.png')

# Figure 02 - expected deaths scenario with expressly hypothetical effectiveness
bars <- tibble(
  case=factor(c('Baseline','Illustrative 50%\nrealisation','Full model-implied\neffect'),
              levels=c('Baseline','Illustrative 50%\nrealisation','Full model-implied\neffect')),
  deaths=c(summary$baseline_deaths[1],
           summary$baseline_deaths[1]-summary$half_realisation_avoided[1],
           summary$modelled_after[1]))
p2 <- ggplot(bars,aes(x=case,y=deaths,fill=case)) +
  geom_col(width=.53,show.legend=FALSE) +
  geom_text(aes(label=sprintf('%.2f',deaths)),vjust=-.6,size=6,
            fontface='bold',colour=BLUE$navy) +
  scale_fill_manual(values=c(BLUE$ice,BLUE$sky,BLUE$navy)) +
  scale_y_continuous(limits=c(0,max(bars$deaths)*1.18),expand=c(0,0)) +
  labs(title='How much could the portfolio change severe outcomes?',
       subtitle='FY2027 severe incident expectations | Effect-realisation sensitivity',
       x=NULL,y='Expected severe outcomes (people)',
       caption='50% realisation is an assumed scenario, NOT a model confidence interval or causal estimate.') +
  slide_theme() + theme(axis.text.x=element_text(face='bold'))
save_slide(p2,'02_expected_deaths_scenarios.png')

# Figure 03 - marginal expected reductions by action vs cost (evidence for choices)
p3 <- ggplot(portfolio,aes(x=total_cost,y=incremental_deaths_avoided,colour=category)) +
  geom_point(size=3,alpha=.83) +
  scale_colour_manual(values=setNames(rep_len(colours,length(unique(portfolio$category))),
                                      sort(unique(portfolio$category)))) +
  scale_x_continuous(labels=fd) +
  scale_y_continuous(labels=number_format(accuracy=.001)) +
  labs(title='Which selected purchases carry the most modelled benefit?',
       subtitle='Each point is one purchased bundle; benefit is recalculated after earlier choices',
       x='First-year bundle cost',y='Marginal model-implied severe outcomes avoided',
       caption=model_note) + slide_theme()
save_slide(p3,'03_action_cost_vs_benefit.png')

# Figure 04 - 12 beaches with highest baseline expected severe incidents,
# contrasted with their post-intervention model scenario.
beach_rank <- months %>% group_by(beach_nm) %>%
  summarise(baseline=sum(baseline_severe),after=sum(post_investment_severe),
            reduction=sum(reduction),.groups='drop') %>%
  slice_max(baseline,n=12,with_ties=FALSE) %>%
  mutate(beach_nm=factor(beach_nm,levels=rev(beach_nm)))
p4 <- ggplot(beach_rank,aes(y=beach_nm)) +
  geom_segment(aes(x=after,xend=baseline,yend=beach_nm),linewidth=3,
               colour=BLUE$light) +
  geom_point(aes(x=baseline,colour='Baseline'),size=4.2) +
  geom_point(aes(x=after,colour='With investment'),size=4.2) +
  scale_colour_manual(values=c('Baseline'=BLUE$navy,'With investment'=BLUE$mid)) +
  scale_x_continuous(labels=number_format(accuracy=.01),
                     expand=expansion(mult=c(.03,.12))) +
  labs(title='How does funding change risk at priority beaches?',
       subtitle='12 beaches with the highest baseline expected severe-incident counts',
       x='FY2027 expected severe outcomes',y=NULL,caption=model_note) + slide_theme()
save_slide(p4,'04_priority_beaches_before_after.png')

# Figure 05 - timing of model-implied improvements
by_month <- months %>% mutate(calendar_month=factor(month.abb[cal_m],
                      levels=month.abb[c(9:12,1:8)])) %>%
  group_by(calendar_month) %>%
  summarise(baseline=sum(baseline_severe),after=sum(post_investment_severe),
            reduction=sum(reduction),.groups='drop') %>%
  pivot_longer(c(baseline,after),names_to='phase',values_to='expected') %>%
  mutate(phase=factor(phase,levels=c('baseline','after'),
                      labels=c('Baseline','With investment')))
p5 <- ggplot(by_month,aes(calendar_month,expected,fill=phase)) +
  geom_col(position=position_dodge(width=.75),width=.68) +
  scale_fill_manual(values=c('Baseline'=BLUE$light,'With investment'=BLUE$navy)) +
  scale_y_continuous(labels=number_format(accuracy=.1),expand=expansion(mult=c(0,.08))) +
  labs(title='When should benefits of the investment appear?',
       subtitle='FY2027 monthly forecast | Peak demand and campaign months are visible',
       x='Month (September 2026 – August 2027)',y='Expected severe outcomes',
       caption=model_note) + slide_theme()
save_slide(p5,'05_monthly_baseline_vs_investment.png')

# Figure 06 - sequential spend versus reduction (greedy selection path)
frontier <- bind_rows(tibble(cumulative_cost=0,cumulative_deaths_avoided=0),
                      path %>% select(cumulative_cost,cumulative_deaths_avoided))
p6 <- ggplot(frontier,aes(cumulative_cost,cumulative_deaths_avoided)) +
  geom_area(fill=BLUE$ice,alpha=.9) +
  geom_line(colour=BLUE$navy,linewidth=1.5) +
  geom_point(data=tail(frontier,1),colour=BLUE$dark,size=4) +
  geom_vline(xintercept=2000000,linetype='dashed',colour=BLUE$mid,linewidth=.8) +
  scale_x_continuous(labels=fd,limits=c(0,2200000)) +
  scale_y_continuous(labels=number_format(accuracy=.01)) +
  labs(title='How the proposed portfolio uses additional funding',
       subtitle='Cumulative spending versus model-implied reduction in severe outcomes',
       x='Cumulative FY2027 investment',y='Cumulative model-implied severe outcomes avoided',
       caption='Sequential greedy purchase path; not a globally optimal efficient frontier. ' %>%
         paste0(model_note)) + slide_theme()
save_slide(p6,'06_budget_vs_benefit_path.png')

# Figure 07 - alternatives, directly required by Task 2 question 2
alt_plot <- alts %>%
  mutate(strategy=factor(strategy,levels=rev(strategy)))
p7 <- ggplot(alt_plot,aes(x=model_implied_deaths_avoided,y=strategy,
                         fill=grepl('^Recommended:',as.character(strategy)))) +
  geom_col(width=.64,show.legend=FALSE) +
  geom_text(aes(label=paste0(sprintf('%.2f',model_implied_deaths_avoided),
                             '  |  FD$',comma(round(spending)))),
            hjust=-.1,size=4.1,fontface='bold',colour=BLUE$ink) +
  scale_fill_manual(values=c('TRUE'=BLUE$navy,'FALSE'=BLUE$light)) +
  scale_x_continuous(expand=expansion(mult=c(0,.38)),labels=number_format(accuracy=.1)) +
  labs(title='Why choose the recommended investment strategy?',
       subtitle='Four feasible strategies compared using the same severe model and FD$2m budget',
       x='Model-implied severe outcomes avoided',y=NULL,caption=model_note) + slide_theme()
save_slide(p7,'07_alternative_strategy_comparison.png')

# Figure 08 - investment geographic concentration by postcode
p8data <- pc_spend %>% arrange(desc(total_cost)) %>% slice_head(n=15) %>%
  mutate(postcode=factor(postcode,levels=rev(postcode)))
p8 <- ggplot(p8data,aes(total_cost,postcode)) +
  geom_col(fill=BLUE$mid,width=.64) +
  geom_text(aes(label=paste0('FD$',comma(round(total_cost)))),
            hjust=-.12,size=4,colour=BLUE$navy,fontface='bold') +
  scale_x_continuous(labels=fd,expand=expansion(mult=c(0,.25))) +
  labs(title='Which postcodes receive the most investment?',
       subtitle='15 highest-funded postcodes | Includes direct beach costs and postcode programmes',
       x='FY2027 investment',y='Postcode',
       caption='A postcode can contain more than one beach; campaign and community-training charges are not duplicated.') +
  slide_theme()
save_slide(p8,'08_investment_by_postcode.png')

# Figure 09 - sensitivity: no proven causality, show realistic uncertainty
p9 <- ggplot(sensitivity,aes(x=realisation_fraction,y=estimated_deaths_avoided)) +
  geom_area(fill=BLUE$ice) + geom_line(colour=BLUE$navy,linewidth=1.5) +
  geom_point(colour=BLUE$mid,size=3.5) +
  geom_text(aes(label=sprintf('%.2f',estimated_deaths_avoided)),
            vjust=-1.0,colour=BLUE$dark,fontface='bold',size=4) +
  scale_x_continuous(labels=percent_format(accuracy=1),breaks=c(0,.25,.5,.75,1)) +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  labs(title='How sensitive are claimed benefits to causal assumptions?',
       subtitle='From no realised intervention effect to the full model-implied association',
       x='Illustrative fraction of modelled effect realised',
       y='Severe outcomes avoided (scenario)',
       caption='Illustrative policy scenarios only; these percentages are not probabilities or confidence intervals.') +
  slide_theme()
save_slide(p9,'09_effect_realisation_sensitivity.png')

cat('\nAll nine slide-ready figures saved to: ',normalizePath(fig_dir),'\n',sep='')
