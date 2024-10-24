## ---------------------------
##
## Program: 8a. Incidence rates and hazard ratios
##
## Purpose: 
## 1. Calculate incidence rates, rate ratios, and rate differences with and without weights for the outcome by exposure group.
## 2. Calculate cox proportional hazard ratios with and without rates for the outcome by exposure group. 
##    - For clinical outcomes, competing risk of death is handled using a Fine-Gray Cox regression model. 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-10-15
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive, antidiabetic
exposure <- 'antihypertensive'

# outcome: all-cause mortality, suicidal ideation
outcome <- 'all-cause mortality'

# analysis: main, flexible_grace_period, 90_day_grace_period, male, female
# young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed
analysis <- 'main'

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(ggplot2)
library(writexl)
library(lubridate)

#### DEFINE PATHS ####

path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_final_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'final', sep = '/')

path_comorb_cprd <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/data', exposure, 'comorbidities', 'Aurum code comorbidities', sep = '/')
comorb_cprd_files <- list.files(path_comorb_cprd, pattern = '.xlsx', all.files = TRUE, full.names = TRUE)

setwd(path_intermediate_res)

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))
cohort_long <- readRDS(file = paste(path_intermediate_res, 'cohort_ipcw.rds', sep = '/'))
covariates <- readRDS(file = paste(path_intermediate_res, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_intermediate_res, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_intermediate_res, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_intermediate_res, 'dec_comorb.rds', sep = '/'))

if (analysis == 'depressed' | analysis == 'not_depressed') {
  comorbidities <- comorbidities[!comorbidities %in% c('depression')]
  base_comorb <- base_comorb[!base_comorb %in% c('depression_base')]
  dec_comorb <- dec_comorb[!dec_comorb %in% c('depression_d1', 'depression_d2', 'depression_d3',
                                                  'depression_d4', 'depression_d5', 'depression_d6', 
                                                  'depression_d7', 'depression_d8', 'depression_d9')]
}

times_dec <- readRDS(paste(path_intermediate_res, 'times_dec.rds', sep = '/'))
times_dec <- times_dec[1:9] # remove last decile (100%) for splitting later on
times_dec

table(cohort$trt)

#### ITT ANALYSIS: INCIDENCE RATES ####

## Compute weighted person-time and events to calculate weighted IR
## for unstabilized and stabilized IPTW using ITT_cohort_analytic df
## note: IR calculated per YEAR

## Summary of incidence over entire study period

# quick overview table of events
table(cohort$itt_event, cohort$trt_dummy)
cohort$itt_event <- as.numeric(cohort$itt_event)

cohort_analytic_itt <- cohort
cohort_analytic_itt$itt_event <- as.numeric(cohort_analytic_itt$itt_event)

cohort_analytic_itt %<>%
  mutate (
    person_time = as.numeric(itt_exit_date - entry_date),
    iptw_person_time = person_time * iptw,
    siptw_person_time = person_time * siptw,
    iptw_event = itt_event * iptw,
    siptw_event = itt_event * siptw,
  )

# no weights
itt_ir <- cohort_analytic_itt %>%
  summarise (
    n = length(unique(id)),
    n_ref = n_distinct(id[trt_dummy == 0]), # number of patients in reference trt group
    n_comp= n_distinct(id[trt_dummy == 1]), # number of patients in comparator trt group
    n_events = sum(itt_event), # total events
    n_events_ref = n_distinct(id[itt_event == 1 & trt_dummy == 0]), # number of events in ref trt group
    n_events_comp = n_distinct(id[itt_event == 1 & trt_dummy == 1]), # number of events in comparator trt group
    total_follow_up_yrs = sum(person_time)/365, # total person-time (years)
    IR = n_events / total_follow_up_yrs, # overall incidence rate (per year)
    IR_ref = n_events_ref/sum(person_time[trt_dummy == 0]) * 365, # incidence rate in reference group (per year)
    IR_comp = n_events_comp/sum(person_time[trt_dummy == 1]) * 365, # incidence rate in comparator group (per year)
    IRR = IR_comp / IR_ref, # incidence rate ratio
    IRD = IR_comp - IR_ref # incidence rate difference
  )

itt_ir$model <- 'itt_ir'

# iptw and siptw weights
itt_ir_iptw <- cohort_analytic_itt %>%
  summarise (
    n_events = sum(iptw_event),
    n_events_ref = sum(iptw_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_event),
    total_follow_up_stab_yrs = sum(siptw_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

itt_ir_iptw$model <- 'itt_ir_iptw'

#### ITT ANALYSIS: HAZARD RATIOS ####

## no weights
cox_itt <- coxph(
  Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ trt_dummy,
  data = cohort_analytic_itt,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_itt, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) # check proportional hazard assumption, all p-values should be large (>0.05)
ph

png(filename = 'cox_itt_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_itt)
saveRDS(cox_itt, file = paste(path_final_res, 'cox_itt.rds', sep = '/'))

## IPTW weights

# unstab
cox_itt_iptw <- coxph(
  Surv(as.numeric(itt_exit_date, entry_date), itt_event) ~ trt_dummy,
  data = cohort_analytic_itt,
  weights = iptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_itt_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_itt_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_itt_iptw)
saveRDS(cox_itt_iptw, file = paste(path_final_res, 'cox_itt_iptw.rds', sep = '/'))

# stab
cox_itt_siptw <- coxph(
  Surv(as.numeric(itt_exit_date, entry_date), itt_event) ~ trt_dummy,
  data = cohort_analytic_itt,
  weights = siptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_itt_siptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_itt_siptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_itt_siptw)
saveRDS(cox_itt_siptw, file = paste(path_final_res, 'cox_itt_siptw.rds', sep = '/'))
rm(cox_itt, cox_itt_iptw, cox_itt_siptw)

#### AT ANALYSIS: ANALYTIC COHORT DATAFRAME ####

cohort_analytic_at <- cohort

## Create analytic dataframe with long cohort split by event or at exit time

# time: time to event or cohort at cohort exit
# counting_time: time in days between cohort entry and time at exit
# event_at_time: whether patient experienced the event at counting time 

cohort_analytic_at$time <- pmin(cohort_analytic_at$at_exit_date, cohort_analytic_at$at_event_date, na.rm = TRUE)
cohort_analytic_at$event_at_time <- ifelse(is.na(cohort_analytic_at$at_event_date), 0, ifelse(cohort_analytic_at$time == cohort_analytic_at$at_event_date, 1, 0))
cohort_analytic_at$event_counting_time <- as.numeric(difftime(cohort_analytic_at$at_event_date, cohort_analytic_at$entry_date, units = 'days'))
cohort_analytic_at$counting_time <- as.numeric(difftime(cohort_analytic_at$time, cohort_analytic_at$entry_date, units = "days"))
cohort_analytic_at <- cohort_analytic_at[order(cohort_analytic_at$counting_time),]

cohort_analytic_at$Tstart <- 0
cohort_analytic_at <- survSplit(cohort_analytic_at, cut=times_dec, end="counting_time",
                                start="Tstart", event="event_at_time")

names(cohort_analytic_at)[names(cohort_analytic_at) == 'counting_time'] <- 'Tstop'
cohort_analytic_at$event_at_tstop <- if_else(is.na(cohort_analytic_at$event_counting_time), 0, # indicator for having event for given time interval
                                             if_else(cohort_analytic_at$Tstop == cohort_analytic_at$event_counting_time, 1, 0)) 

# retrieve covariate values at deciles
for (i in 1:length(comorbidities)) {
  comorb_name <- comorbidities[i]
  cohort_analytic_at$comorb <- if_else(
    cohort_analytic_at$Tstart == 0,
    cohort_analytic_at[[paste(comorb_name, 'base', sep = '_')]],
    if_else(
      cohort_analytic_at$Tstart == times_dec[[1]],
      cohort_analytic_at[[paste(comorb_name, 'd1', sep = '_')]],
      if_else(
        cohort_analytic_at$Tstart == times_dec[[2]],
        cohort_analytic_at[[paste(comorb_name, 'd2', sep = '_')]],
        if_else(
          cohort_analytic_at$Tstart == times_dec[[3]],
          cohort_analytic_at[[paste(comorb_name, 'd3', sep = '_')]],
          if_else(
            cohort_analytic_at$Tstart == times_dec[[4]],
            cohort_analytic_at[[paste(comorb_name, 'd4', sep = '_')]],
            if_else(
              cohort_analytic_at$Tstart == times_dec[[5]],
              cohort_analytic_at[[paste(comorb_name, 'd5', sep = '_')]],
              if_else(
                cohort_analytic_at$Tstart == times_dec[[6]],
                cohort_analytic_at[[paste(comorb_name, 'd6', sep = '_')]],
                if_else(
                  cohort_analytic_at$Tstart == times_dec[[7]],
                  cohort_analytic_at[[paste(comorb_name, 'd7', sep = '_')]],
                  if_else(
                    cohort_analytic_at$Tstart == times_dec[[8]],
                    cohort_analytic_at[[paste(comorb_name, 'd8', sep = '_')]],
                    cohort_analytic_at[[paste(comorb_name, 'd9', sep = '_')]])
                )
              )
            )
          )
        )
      )
    )
  )
  names(cohort_analytic_at)[names(cohort_analytic_at) == 'comorb'] <- comorb_name
}

table(cohort_analytic_at$event_at_tstop, cohort_analytic_at$trt)

# add IPCW corresponding to time interval
cohort_analytic_at$decTstart <- if_else(
  cohort_analytic_at$Tstart >= 0 & cohort_analytic_at$Tstart < times_dec[1],
  0,
  if_else(
    cohort_analytic_at$Tstart >= times_dec[1] & cohort_analytic_at$Tstart < times_dec[2],
    times_dec[1],
    if_else(
      cohort_analytic_at$Tstart >= times_dec[2] & cohort_analytic_at$Tstart < times_dec[3],
      times_dec[2],
      if_else(
        cohort_analytic_at$Tstart >= times_dec[3] & cohort_analytic_at$Tstart < times_dec[4],
        times_dec[3],
        if_else(
          cohort_analytic_at$Tstart >= times_dec[4] & cohort_analytic_at$Tstart < times_dec[5],
          times_dec[4],
          if_else(
            cohort_analytic_at$Tstart >= times_dec[5] & cohort_analytic_at$Tstart < times_dec[6],
            times_dec[5],
            if_else(
              cohort_analytic_at$Tstart >= times_dec[6] & cohort_analytic_at$Tstart < times_dec[7],
              times_dec[6],
              if_else(
                cohort_analytic_at$Tstart >= times_dec[7] & cohort_analytic_at$Tstart < times_dec[8],
                times_dec[7],
                if_else(
                  cohort_analytic_at$Tstart >= times_dec[8] & cohort_analytic_at$Tstart < times_dec[9],
                  times_dec[8], 
                  times_dec[9]
                )
              )
            )
          )
        )
      )
    )
  )
)

ipcw_weights <- cohort_long %>% dplyr::select(id, Tstart, ipcw_str_lag, sipcw_str_lag, ipcw_str_nonlag, sipcw_str_nonlag, ipcw_pl_nonlag, sipcw_pl_nonlag, ipcw_pl_lag, sipcw_pl_lag)

cohort_analytic_at <- cohort_analytic_at %>%
  dplyr::left_join(ipcw_weights, by = c('id', 'decTstart' = 'Tstart'), relationship = 'many-to-one')

cohort_analytic_at$at_event <- as.numeric(cohort_analytic_at$at_event)

# for plots later
saveRDS(cohort_analytic_at, file = paste(path_intermediate_res, 'cohort_analytic_at.rds', sep = '/'))

#### AT ANALYSIS: INCIDENCE RATES ####

## Compute weighted person-time and events to calculate weighted IR
## for unstabilized and stabilized IPTW, IPCW, and combination weights
## for the AT analytic dataframe.

# convert event to a numeric for calculations

cohort_analytic_at %<>%
  mutate (
    person_time = (Tstop - Tstart),
    
    # IPTW and sIPTW
    iptw_person_time = person_time * iptw,
    iptw_event = event_at_tstop * iptw,
    
    siptw_person_time = person_time * siptw,
    siptw_event = event_at_tstop * siptw,
    
    
    # IPCW and sIPCW - stratified lagged
    ipcw_str_lag_person_time = person_time * ipcw_str_lag,
    ipcw_str_lag_event = event_at_tstop * ipcw_str_lag,
    
    sipcw_str_lag_person_time = person_time * sipcw_str_lag,
    sipcw_str_lag_event = event_at_tstop * sipcw_str_lag,
    
    iptw_ipcw_str_lag_person_time = person_time * iptw * ipcw_str_lag,
    iptw_ipcw_str_lag_event = event_at_tstop * iptw * ipcw_str_lag,
    
    siptw_sipcw_str_lag_person_time = person_time * siptw * sipcw_str_lag,
    siptw_sipcw_str_lag_event = event_at_tstop * siptw * sipcw_str_lag,
    
    
    # ipcw_str_nonlag and sipcw_str_nonlag - stratified non-lagged
    ipcw_str_nonlag_person_time = person_time * ipcw_str_nonlag,
    ipcw_str_nonlag_event = event_at_tstop * ipcw_str_nonlag,
    
    sipcw_str_nonlag_person_time = person_time * sipcw_str_nonlag,
    sipcw_str_nonlag_event = event_at_tstop * sipcw_str_nonlag,
    
    iptw_ipcw_str_nonlag_person_time = person_time * iptw * ipcw_str_nonlag,
    iptw_ipcw_str_nonlag_event = event_at_tstop * iptw * ipcw_str_nonlag,
    
    siptw_sipcw_str_nonlag_person_time = person_time * siptw * sipcw_str_nonlag,
    siptw_sipcw_str_nonlag_event = event_at_tstop * siptw * sipcw_str_nonlag,
    
    
    # ipcw_pl_nonlag and sipcw_pl_nonlag - pooled non-lagged
    ipcw_pl_nonlag_person_time = person_time * ipcw_pl_nonlag,
    ipcw_pl_nonlag_event = event_at_tstop * ipcw_pl_nonlag,
    
    sipcw_pl_nonlag_person_time = person_time * sipcw_pl_nonlag,
    sipcw_pl_nonlag_event = event_at_tstop * sipcw_pl_nonlag,
    
    iptw_ipcw_pl_nonlag_person_time = person_time * iptw * ipcw_pl_nonlag,
    iptw_ipcw_pl_nonlag_event = event_at_tstop * iptw * ipcw_pl_nonlag,
    
    siptw_sipcw_pl_nonlag_person_time = person_time * siptw * sipcw_pl_nonlag,
    siptw_sipcw_pl_nonlag_event = event_at_tstop * siptw * sipcw_pl_nonlag,
    
    # ipcw_pl_lag and sipcw_pl_lag - pooled lagged
    ipcw_pl_lag_person_time = person_time * ipcw_pl_lag,
    ipcw_pl_lag_event = event_at_tstop * ipcw_pl_lag,
    
    sipcw_pl_lag_person_time = person_time * sipcw_pl_lag,
    sipcw_pl_lag_event = event_at_tstop * sipcw_pl_lag,
    
    iptw_ipcw_pl_lag_person_time = person_time * iptw * ipcw_pl_lag,
    iptw_ipcw_pl_lag_event = event_at_tstop * iptw * ipcw_pl_lag,
    
    siptw_sipcw_pl_lag_person_time = person_time * siptw * sipcw_pl_lag,
    siptw_sipcw_pl_lag_event = event_at_tstop * siptw * sipcw_pl_lag
  )

## Summary of incidence over entire study period

# no weights
at_ir <- cohort_analytic_at %>%
  summarise (
    n = length(unique(id)),
    n_ref = n_distinct(id[trt_dummy == 0]), # number of patients in reference trt group
    n_comp= n_distinct(id[trt_dummy == 1]), # number of patients in comparator trt group
    n_events = sum(event_at_tstop), # total events
    n_events_ref = n_distinct(id[event_at_tstop == 1 & trt_dummy == 0]), # number of events in ref trt group
    n_events_comp = n_distinct(id[event_at_tstop == 1 & trt_dummy == 1]), # number of events in comparator trt group
    total_follow_up_yrs = sum(person_time)/365, # total person-time (years)
    IR = n_events / total_follow_up_yrs, # overall incidence rate (per year)
    IR_ref = n_events_ref/sum(person_time[trt_dummy == 0]) * 365, # incidence rate in reference group (per year)
    IR_comp = n_events_comp/sum(person_time[trt_dummy == 1]) * 365, # incidence rate in comparator group (per year)
    IRR = IR_comp / IR_ref, # incidence rate ratio
    IRD = IR_comp - IR_ref # incidence rate difference
  )

at_ir$model <- 'at_ir'

# iptw and siptw weights
at_ir_iptw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_event),
    n_events_ref = sum(iptw_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_event),
    total_follow_up_stab_yrs = sum(siptw_person_time) / 365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw$model <- 'at_ir_iptw'

# stratified lagged ipcw weights
at_ir_ipcw_str_lag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_str_lag_event),
    n_events_ref = sum(ipcw_str_lag_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_str_lag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_str_lag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_str_lag_event),
    total_follow_up_stab_yrs = sum(sipcw_str_lag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_str_lag_event[trt_dummy == 0])/sum(sipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_str_lag_event[trt_dummy == 1])/sum(sipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_str_lag$model <- 'at_ir_ipcw_str_lag'

at_ir_iptw_ipcw_str_lag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_str_lag_event),
    n_events_ref = sum(iptw_ipcw_str_lag_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_str_lag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_str_lag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_str_lag_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_str_lag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_str_lag_event[trt_dummy == 0])/sum(siptw_sipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_str_lag_event[trt_dummy == 1])/sum(siptw_sipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_str_lag$model <- 'at_ir_iptw_ipcw_str_lag'

# stratified non-lagged ipcw weights
at_ir_ipcw_str_nonlag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_str_nonlag_event),
    n_events_ref = sum(ipcw_str_nonlag_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_str_nonlag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_str_nonlag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_str_nonlag_event),
    total_follow_up_stab_yrs = sum(sipcw_str_nonlag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_str_nonlag_event[trt_dummy == 0])/sum(sipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_str_nonlag_event[trt_dummy == 1])/sum(sipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_str_nonlag$model <- 'at_ir_ipcw_str_nonlag'

at_ir_iptw_ipcw_str_nonlag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_str_nonlag_event),
    n_events_ref = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_str_nonlag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_str_nonlag_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_str_nonlag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_str_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_str_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_str_nonlag$model <- 'at_ir_iptw_ipcw_str_nonlag'


# pooled non-lagged ipcw weights
at_ir_ipcw_pl_nonlag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_pl_nonlag_event),
    n_events_ref = sum(ipcw_pl_nonlag_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_pl_nonlag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_pl_nonlag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_pl_nonlag_event),
    total_follow_up_stab_yrs = sum(sipcw_pl_nonlag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_pl_nonlag_event[trt_dummy == 0])/sum(sipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pl_nonlag_event[trt_dummy == 1])/sum(sipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_pl_nonlag$model <- 'at_ir_ipcw_pl_nonlag'

at_ir_iptw_ipcw_pl_nonlag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_pl_nonlag_event),
    n_events_ref = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_pl_nonlag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_pl_nonlag_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_pl_nonlag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_pl_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pl_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_pl_nonlag$model <- 'at_ir_iptw_ipcw_pl_nonlag'


# pooled lagged ipcw weights
at_ir_ipcw_pl_lag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_pl_lag_event),
    n_events_ref = sum(ipcw_pl_lag_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_pl_lag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_pl_lag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_pl_lag_event),
    total_follow_up_stab_yrs = sum(sipcw_pl_lag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_pl_lag_event[trt_dummy == 0])/sum(sipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pl_lag_event[trt_dummy == 1])/sum(sipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_pl_lag$model <- 'at_ir_ipcw_pl_lag'

at_ir_iptw_ipcw_pl_lag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_pl_lag_event),
    n_events_ref = sum(iptw_ipcw_pl_lag_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_pl_lag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_pl_lag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_pl_lag_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_pl_lag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_pl_lag_event[trt_dummy == 0])/sum(siptw_sipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pl_lag_event[trt_dummy == 1])/sum(siptw_sipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_pl_lag$model <- 'at_ir_iptw_ipcw_pl_lag'

#### AT ANALYSIS: HAZARD RATIOS ####

## NO WEIGHTS
cox_at <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at)
saveRDS(cox_at, file = paste(path_final_res, 'cox_at.rds', sep = '/'))

## IPTW weights
# unstab
cox_at_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw)
saveRDS(cox_at_iptw, file = paste(path_final_res, 'cox_at_iptw.rds', sep = '/'))

# stab
cox_at_siptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw)
saveRDS(cox_at_siptw, file = paste(path_final_res, 'cox_at_siptw.rds', sep = '/'))

## STRATIFIED LAGGED WEIGHTS
# unstab
cox_at_ipcw_str_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_str_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_ipcw_str_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ipcw_str_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_ipcw_str_lag)
saveRDS(cox_at_ipcw_str_lag, file = paste(path_final_res, 'cox_at_ipcw_str_lag.rds', sep = '/'))

# stab
cox_at_sipcw_str_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_str_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_str_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_str_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_str_lag)
saveRDS(cox_at_sipcw_str_lag, file = paste(path_final_res, 'cox_at_sipcw_str_lag.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_str_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_str_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_str_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_str_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_str_lag)
saveRDS(cox_at_iptw_ipcw_str_lag, file = paste(path_final_res, 'cox_at_iptw_ipcw_str_lag.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_str_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_str_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_str_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_str_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_str_lag)
saveRDS(cox_at_siptw_sipcw_str_lag, file = paste(path_final_res, 'cox_at_siptw_sipcw_str_lag.rds', sep = '/'))

## STRATIFIED NON-LAGGED WEIGHTS
# unstab
cox_at_icpw_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_str_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_icpw_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_icpw_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_icpw_nonlag)
saveRDS(cox_at_icpw_nonlag, file = paste(path_final_res, 'cox_at_ipcw_str_nonlag.rds', sep = '/'))

# stab
cox_at_sipcw_str_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_str_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_str_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_str_nonlag.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_str_nonlag)
saveRDS(cox_at_sipcw_str_nonlag, file = paste(path_final_res, 'cox_at_sipcw_str_nonlag.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_str_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_str_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_str_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_str_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_str_nonlag)
saveRDS(cox_at_iptw_ipcw_str_nonlag, file = paste(path_final_res, 'cox_at_iptw_ipcw_str_nonlag.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_str_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_str_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_str_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_str_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_str_nonlag)
saveRDS(cox_at_siptw_sipcw_str_nonlag, file = paste(path_final_res, 'cox_at_siptw_sipcw_str_nonlag.rds', sep = '/'))

## POOLED NON-LAGGED WEIGHTS 
# unstab
cox_at_ipcw_pl_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_pl_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_ipcw_pl_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ipcw_pl_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_ipcw_pl_nonlag)
saveRDS(cox_at_ipcw_pl_nonlag, file = paste(path_final_res, 'cox_at_ipcw_pl_nonlag.rds', sep = '/'))

# stab
cox_at_sipcw_pl_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_pl_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_pl_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_pl_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_pl_nonlag)
saveRDS(cox_at_sipcw_pl_nonlag, file = paste(path_final_res, 'cox_at_sipcw_pl_nonlag.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_pl_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_pl_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_pl_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_pl_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_pl_nonlag)
saveRDS(cox_at_iptw_ipcw_pl_nonlag, file = paste(path_final_res, 'cox_at_iptw_ipcw_pl_nonlag.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_pl_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_pl_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_pl_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_pl_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_pl_nonlag)
saveRDS(cox_at_siptw_sipcw_pl_nonlag, file = paste(path_final_res, 'cox_at_siptw_sipcw_pl_nonlag.rds', sep = '/'))


## POOLED LAGGED WEIGHTS
# unstab
cox_at_ipcw_pl_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_pl_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_ipcw_pl_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ipcw_pl_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_ipcw_pl_lag)
saveRDS(cox_at_ipcw_pl_lag, file = paste(path_final_res, 'cox_at_ipcw_pl_lag.rds', sep = '/'))

# stab
cox_at_sipcw_pl_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_pl_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_pl_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_pl_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_pl_lag)
saveRDS(cox_at_sipcw_pl_lag, file = paste(path_final_res, 'cox_at_sipcw_pl_lag.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_pl_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_pl_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_pl_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_pl_lag.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_pl_lag)
saveRDS(cox_at_iptw_ipcw_pl_lag, file = paste(path_final_res, 'cox_at_iptw_ipcw_pl_lag.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_pl_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_pl_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_pl_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_pl_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_pl_lag)
saveRDS(cox_at_siptw_sipcw_pl_lag, file = paste(path_final_res, 'cox_at_siptw_sipcw_pl_lag.rds', sep = '/'))

#### DATAFRAME WITH INCIDENCE RATES FOR EACH MODEL ####

incidence_rates <- bind_rows(
  itt_ir,
  itt_ir_iptw,
  at_ir,
  at_ir_iptw,
  at_ir_ipcw_str_lag,
  at_ir_iptw_ipcw_str_lag,
  at_ir_ipcw_str_nonlag,
  at_ir_iptw_ipcw_str_nonlag,
  at_ir_ipcw_pl_nonlag,
  at_ir_iptw_ipcw_pl_nonlag,
  at_ir_ipcw_pl_lag,
  at_ir_iptw_ipcw_pl_lag
)

incidence_rates %<>% relocate(model)

write_xlsx(incidence_rates, paste(path_final_res, 'incidence_rates.xlsx', sep ='/'))

