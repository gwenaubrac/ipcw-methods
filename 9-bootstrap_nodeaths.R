## ---------------------------
##
## Program: 9. Bootstrapping 
##
## Purpose: Bootstrap cohort to obtain percentile confidence intervals for IR, IRR, and IRD.
## The same steps from program 4b (defining outcome) to 8 (running analyses) are replicated.
## For efficiency, the code was rewritten using data.table. 
##
## Follow-up intervals in which patients died are excluded from the IPCW model fitting. 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-11-01
##
## ---------------------------
##
## Notes: Less comments/documentation because the same steps from previous programs
## are simply repeated inside the bootstrap function. 
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive
exposure <- 'antihypertensive'
outcome <- 'all-cause mortality'

# analysis: main, flexible_grace_period, 90_day_grace_period, male, female, young, old
analysis <- '90_day_grace_period'

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(lubridate)
library(purrr)
library(data.table)
library(fastglm)
library(boot)
library(parallel)
library(writexl)
library(tidyr)

#### DEFINE PATHS ####

path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, 'main', 'intermediate', sep = '/')
path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_final_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, analysis, 'final', sep = '/')

## Variables and dataframes used within bootstrap

cohort <- readRDS(file = paste(path_intermediate_res_main, 'cohort_outcome.rds', sep = '/'))

cohort <- switch(
  analysis,
  "male" = cohort %<>% filter(sex == 'Male'),
  "female" = cohort %<>% filter(sex == 'Female'),
  "young" = cohort %<>% filter(age_at_entry < 65),
  "old" = cohort %<>% filter(age_at_entry >= 65),
  cohort
)

# for sensitivity analyses, take cohort from censored analyses
if (analysis == 'flexible_grace_period' | analysis == '90_day_grace_period') {
  cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_outcome.rds', sep = '/'))
}

covariates <- readRDS(file = paste(path_intermediate_res_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_intermediate_res_main, 'comorbidities.rds', sep = '/'))
breaks <- readRDS(file = paste(path_intermediate_res_main, 'interval_breaks.rds', sep = '/'))
model <- readRDS(file = paste(path_final_res, 'iptw_model.rds', sep = '/'))

first_comorb <- readRDS(paste(path_intermediate_res, 'first_comorb.rds', sep = '/'))
first_comorb <- first_comorb %>%
  pivot_wider(id_cols = id, names_from = comorb, values_from = comorb_date)

base_comorb <- readRDS(paste(path_intermediate_res_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(paste(path_intermediate_res_main, 'dec_comorb.rds', sep = '/'))
cat_variables <- c('sex', 'year', 'ethnicity', 'deprivation', base_comorb, dec_comorb, 'disc', 'switch')

# remove variables for subgroup analyses
variables <- c(covariates, comorbidities)

if (analysis == 'male' | analysis == 'female') {
  variables <- variables[!variables %in% c('sex')]
} else if (analysis == 'young' | analysis == 'old') {
  cohort$age_group <- droplevels(cohort$age_group)
}
variables

# remove variables that violate positivity in IPCW

if (exposure == 'antidepressant' & analysis == 'female') {variables <- variables[!variables %in% c('pacemaker')]}
if (exposure == 'antidepressant' & analysis == 'young') {variables <- variables[!variables %in% c('pacemaker', 'valvular_heart_disease', 'lvh')]}
# if (exposure == 'antidepressant' &
#     analysis == 'old') {
#   variables <- variables[!variables %in% c(
#     'suicidal_ideation_self_harm',
#     'deprivation',
#     'epilepsy',
#     'pacemaker'
#   )]
# }

## Functions used within bootstrap ##

# IPTW model
predictors <- as.formula(paste('~', as.character(model)[-1][2], sep = ''))
predictors

# IPCW model - no interactions
p_uncens_predictors <- as.formula(paste('uncensored_at_tstop', "~", paste(variables, collapse = ' + ')))
p_uncens_predictors_pooled <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'dec'), collapse = " + ")))

p_uncens_predictors
p_uncens_predictors_pooled


fit_and_predict <- function(data_subset) {
  x <- model.matrix(p_uncens_predictors, data = data_subset)
  y <- data_subset$uncensored_at_tstop
  model <- fastglm(x=x, y=y, family = binomial(link = "logit"))
  data_subset$p_uncens <- predict(model, newdata = x, type = "response")
  return(data_subset)
}

fit_and_predict_pooled <- function(data_subset) {
  x <- model.matrix(p_uncens_predictors_pooled, data = data_subset)
  y <- data_subset$uncensored_at_tstop
  model <- fastglm(x=x, y=y, family = binomial(link = "logit"))
  data_subset$p_uncens <- predict(model, newdata = x, type = "response")
  return(data_subset)
}

rm(variables, covariates, model)

setDT(cohort)
setDT(first_comorb)

#### BOOTSTRAP FUNCTION ####

# note: run one time inside fx to get result names

bs <- function(data, indices) {
  
  # for initial run inside fx:
  # indices <- sample(x = 1:length(cohort$id), size = length(cohort$id), replace = TRUE)
  # data <- cohort
  # or d <- cohort for original sample
  
  d <- data[indices,]
  d[, boot_id := .I]
  
  #### 0. GET COVARIATE VALUES AT INTERVALS ####
  
  d[, `:=`(
    censor = as.integer(switch == 1 | disc == 1)
  )]
  
  d[, `:=`(
    censor_date = if_else(censor == 1, pmin(switch_date, disc_date, na.rm = TRUE), NA)
  )]
  
  d[, censor_date := as.Date(censor_date)]
  
  censored_only <- d[!is.na(censor_date)]
  
  censoring_times <- as.numeric(censored_only$censor_date - censored_only$entry_date)
  times_dec <- quantile(censoring_times, probs = breaks)
  
  if (exposure == 'antidepressant') {
    # as in original sample, push start of 2nd interval by 1 day
    if (times_dec[1] != times_dec[2]) {
      times_dec[[1]] <- times_dec[[1]] + 1
      # but if start of interval 2 = start of interval 3, push start of 3rd interval by 1 day
    } else if (times_dec[1] == times_dec[2]) {
      times_dec[[2]] <- times_dec[[2]] + 1
    }
  }
  
  times_dec <- times_dec[1:9]
  
  # get relative date of censoring interval for each patient
  d[, paste0("cens_d", 1:9) := lapply(1:9, function(i) entry_date + times_dec[[i]])]
  
  for (i in 2:ncol(first_comorb)) {
    comorb_name <- names(first_comorb)[i]
    first_comorb_i <- first_comorb[, .(id, comorb_date = get(comorb_name))]
    
    d <- merge(d, first_comorb_i, by = "id", all.x = TRUE)
    
    d[, `:=`(
      comorb_base = as.integer(!is.na(comorb_date) & comorb_date < entry_date),
      comorb_d1 = as.integer(!is.na(comorb_date) & comorb_date < cens_d1),
      comorb_d2 = as.integer(!is.na(comorb_date) & comorb_date < cens_d2),
      comorb_d3 = as.integer(!is.na(comorb_date) & comorb_date < cens_d3),
      comorb_d4 = as.integer(!is.na(comorb_date) & comorb_date < cens_d4),
      comorb_d5 = as.integer(!is.na(comorb_date) & comorb_date < cens_d5),
      comorb_d6 = as.integer(!is.na(comorb_date) & comorb_date < cens_d6),
      comorb_d7 = as.integer(!is.na(comorb_date) & comorb_date < cens_d7),
      comorb_d8 = as.integer(!is.na(comorb_date) & comorb_date < cens_d8),
      comorb_d9 = as.integer(!is.na(comorb_date) & comorb_date < cens_d9)
    )]
    
    setnames(d, "comorb_date", paste(comorb_name, "date", sep = "_"))
    setnames(d, "comorb_base", paste(comorb_name, "base", sep = "_"))
    setnames(d, "comorb_d1", paste(comorb_name, "d1", sep = "_"))
    setnames(d, "comorb_d2", paste(comorb_name, "d2", sep = "_"))
    setnames(d, "comorb_d3", paste(comorb_name, "d3", sep = "_"))
    setnames(d, "comorb_d4", paste(comorb_name, "d4", sep = "_"))
    setnames(d, "comorb_d5", paste(comorb_name, "d5", sep = "_"))
    setnames(d, "comorb_d6", paste(comorb_name, "d6", sep = "_"))
    setnames(d, "comorb_d7", paste(comorb_name, "d7", sep = "_"))
    setnames(d, "comorb_d8", paste(comorb_name, "d8", sep = "_"))
    setnames(d, "comorb_d9", paste(comorb_name, "d9", sep = "_"))
  }
  
  d[, (cat_variables) := lapply(.SD, as.factor), .SDcols = cat_variables]
  
  if (exposure == 'antidepressant') {
    d[, depression_base := depression_d1]
  }
  
  if (exposure == 'antihypertensive') {
    d[, hypertension_base := hypertension_d1]
  }
  
  #### 1. CALCUALTE IPTW ####
  
  y <- as.numeric(as.character(d$trt_dummy))
  x_denom <- model.matrix(predictors, data = d)

  denom_fit <- fastglm(x=x_denom, y=y, family = binomial(link = 'logit'))
  pd_trt <- predict(denom_fit, newdata = x_denom, type = 'response')
  d$iptw <- 1 / if_else(d$trt_dummy == 0, 1 - pd_trt, pd_trt)
  
  rm(y, x_denom, denom_fit, pd_trt)
  
  #### 2. CALCULATE IPCW ####
  
  ## Prepare data ##
  d_long <- d
  
  d_long$counting_time <- as.numeric(difftime(d_long$at_exit_date, d_long$entry_date, units = 'days'))
  d_long$censor_counting_time <- as.numeric(difftime(d_long$censor_date, d_long$entry_date, units = 'days'))
  d_long$uncensored <- if_else(d_long$censor == 0, 1, 0)
  d_long <- d_long[order(d_long$counting_time),]
  
  # split up data into intervals
  d_long$Tstart <- 0
  d_long <- survSplit(d_long, cut = times_dec, end = 'counting_time', start = 'Tstart', event = 'uncensored')
  names(d_long)[names(d_long) == 'counting_time'] <- 'Tstop' 
  d_long$uncensored_at_tstop <- if_else(is.na(d_long$censor_counting_time), 1,
                                        if_else(d_long$Tstop == d_long$censor_counting_time, 0, 1)) 
  
  rm(censored_only, censoring_times)
  
  # retrieve covariate values at deciles
  setDT(d_long)
  
  for (comorb_name in comorbidities) {
    base_col <- paste(comorb_name, 'base', sep = '_')
    dec_cols <- paste(comorb_name, paste0('d', 1:9), sep = '_')
    
    d_long[, (comorb_name) := fcase(
      Tstart == 0, get(base_col),
      Tstart == times_dec[[1]], get(dec_cols[1]),
      Tstart == times_dec[[2]], get(dec_cols[2]),
      Tstart == times_dec[[3]], get(dec_cols[3]),
      Tstart == times_dec[[4]], get(dec_cols[4]),
      Tstart == times_dec[[5]], get(dec_cols[5]),
      Tstart == times_dec[[6]], get(dec_cols[6]),
      Tstart == times_dec[[7]], get(dec_cols[7]),
      Tstart == times_dec[[8]], get(dec_cols[8]),
      Tstart == times_dec[[9]], get(dec_cols[9])
    )]
  }
  
  d_long %<>%
    group_by(boot_id) %>%
    mutate(dec = as.factor(row_number())) %>% 
    ungroup()
  
  ## Remove patients who died during interval
  d_long_nodeaths <- d_long %>% 
    mutate(death_fu = as.numeric(dod - entry_date)) %>% 
    mutate(death_in_interval = if_else (!is.na(dod) & death_fu >= Tstart & death_fu <= Tstop, 1, 0)) %>% 
    filter(death_in_interval == 0) 
  
  ## Stratified lagged and non-lagged weights ##
  
  setDT(d_long_nodeaths)
  d_long_nodeaths <- d_long_nodeaths[, fit_and_predict(.SD), by = .(dec, trt_dummy)]

  d_long_nodeaths[, `:=`(
    weight = 1 / p_uncens
  ), by = boot_id]
  
  str_weights <- d_long_nodeaths[, .(boot_id, dec, weight)]
  
  d_long <- merge(d_long, str_weights, all.x = TRUE, by = c('boot_id', 'dec'))
  
  setDT(d_long)
  d_long[is.na(weight), weight := 1]
  
  d_long[, `:=`(
    ipcw_str_nonlag = cumprod(weight)
  ), by = boot_id]
  
  d_long[, `:=`(
    ipcw_str_lag = shift(ipcw_str_nonlag, n = 1, type = "lag", fill = 1)
  ), by = boot_id]
  
  d_long[, `:=`(weight = NULL)]
  
  gc()
  rm(str_weights)
  
  ## Pooled lagged and non-lagged weights ##
  
  d_long_nodeaths <- d_long_nodeaths[, fit_and_predict_pooled(.SD), by = .(trt_dummy)]
  
  d_long_nodeaths[, `:=`(
    weight = 1 / p_uncens
  ), by = boot_id]
  
  pl_weights <- d_long_nodeaths[, .(boot_id, dec, weight)]
  
  d_long <- merge(d_long, pl_weights, all.x = TRUE, by = c('boot_id', 'dec'))
  
  d_long[is.na(weight), weight := 1]
  
  d_long[, `:=`(
    ipcw_pl_nonlag = cumprod(weight)
  ), by = boot_id]
  
  d_long[, `:=`(
    ipcw_pl_lag = shift(ipcw_pl_nonlag, n = 1, type = "lag", fill = 1)
  ), by = boot_id]
  
  d_long[, `:=`(weight = NULL)]
  
  gc()
  rm(pl_weights)
  
  #### 3. CALCULATE INCIDENCE RATES ####
  
  ## A. ITT analyses ##
  
  setDT(d)
  
  itt_d <- d[, .(boot_id, trt, trt_dummy, itt_event, itt_exit_date, entry_date, iptw, month_year)]
  itt_d[, itt_event := as.numeric(itt_event)]
  
  itt_d[, `:=`(
    person_time = as.numeric(itt_exit_date - entry_date),
    iptw_person_time = NA_real_,
    iptw_event = NA_real_
  )]
  
  itt_d[, `:=`(
    iptw_person_time = person_time * iptw,
    iptw_event = itt_event * iptw
  )]
  
  setDT(itt_d)
  
  # no weights
  ir_itt <- itt_d[, .(
    IR = sum(itt_event) / (sum(person_time) / 365),
    IR_ref = sum(itt_event == 1 & trt_dummy == 0) / sum(person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(itt_event == 1 & trt_dummy == 1) / sum(person_time[trt_dummy == 1]) * 365
  )]
  
  ir_itt[, `:=`(
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_itt <- data.frame(t(ir_itt))
  colnames(ir_itt) <- 'itt'
  ir_itt$row_names <- rownames(ir_itt)
  ir_itt <- ir_itt %>% relocate(row_names)
  
  # iptw 
  ir_itt_iptw <- itt_d[, .(
    IR = sum(iptw_event) / sum(iptw_person_time) * 365,
    IR_ref = sum(iptw_event[trt_dummy == 0]) / sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_event[trt_dummy == 1]) / sum(iptw_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_itt_iptw <- ir_itt_iptw[, .(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_itt_iptw <- data.frame(t(ir_itt_iptw))
  colnames(ir_itt_iptw) <- 'itt.iptw'
  ir_itt_iptw$row_names <- rownames(ir_itt_iptw)
  ir_itt_iptw <- ir_itt_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_itt, ir_itt_iptw, by = 'row_names')
  rm(ir_itt, ir_itt_iptw)
  
  ## B. AT ANALYSES ##
  
  at_d <- copy(d)
  at_d[, time := pmin(at_exit_date, at_event_date, na.rm = TRUE)]
  at_d[, `:=`(
    event_at_time = fifelse(is.na(at_event_date), 0, fifelse(time == at_event_date, 1, 0)),
    event_counting_time = as.numeric(difftime(at_event_date, entry_date, units = 'days')),
    counting_time = as.numeric(difftime(time, entry_date, units = "days")),
    Tstart = 0
  )]
  at_d <- at_d[order(counting_time)]
  
  at_d <- survSplit(at_d, cut=times_dec, end="counting_time",
                    start="Tstart", event="event_at_time")
  
  setnames(at_d, old = "counting_time", new = "Tstop")
  
  setDT(at_d)
  at_d[, event_at_tstop := fifelse(
    is.na(event_counting_time), 
    0, 
    fifelse(Tstop == event_counting_time, 1, 0)
  )]
  
  # retrieve covariate values at deciles
  for (comorb_name in comorbidities) {
    base_col <- paste(comorb_name, 'base', sep = '_')
    dec_cols <- paste(comorb_name, paste0('d', 1:9), sep = '_')
    
    at_d[, (comorb_name) := fcase(
      Tstart == 0, get(base_col),
      Tstart == times_dec[[1]], get(dec_cols[1]),
      Tstart == times_dec[[2]], get(dec_cols[2]),
      Tstart == times_dec[[3]], get(dec_cols[3]),
      Tstart == times_dec[[4]], get(dec_cols[4]),
      Tstart == times_dec[[5]], get(dec_cols[5]),
      Tstart == times_dec[[6]], get(dec_cols[6]),
      Tstart == times_dec[[7]], get(dec_cols[7]),
      Tstart == times_dec[[8]], get(dec_cols[8]),
      Tstart == times_dec[[9]], get(dec_cols[9])
    )]
  }
  
  # add IPCW corresponding to time interval
  times_dec_bins <- c(0, times_dec)
  at_d[, decTstart := times_dec_bins[findInterval(Tstart, times_dec_bins)]]
  ipcw_weights <- d_long[, .(boot_id, Tstart, dec, ipcw_str_lag, ipcw_str_nonlag, ipcw_pl_nonlag, ipcw_pl_lag)]
  setnames(ipcw_weights, old = "Tstart", new = "decTstart")
  at_d <- merge(at_d, ipcw_weights, by = c('boot_id', 'decTstart'), all.x = TRUE)
  at_d <- at_d[, .(boot_id, Tstart, Tstop, event_at_tstop, iptw, ipcw_str_lag, ipcw_str_nonlag,
                   ipcw_pl_nonlag, ipcw_pl_lag, dec, trt, trt_dummy, month_year)]
  
  at_d[, person_time := Tstop - Tstart]
  
  at_d[, `:=`(
    # IPTW
    iptw_person_time = person_time * iptw,
    iptw_event = event_at_tstop * iptw,
    
    # stratified lagged IPCW
    ipcw_str_lag_person_time = person_time * ipcw_str_lag,
    ipcw_str_lag_event = event_at_tstop * ipcw_str_lag,
    iptw_ipcw_str_lag_person_time = person_time * iptw * ipcw_str_lag,
    iptw_ipcw_str_lag_event = event_at_tstop * iptw * ipcw_str_lag,
    
    # stratified non-lagged IPCW
    ipcw_str_nonlag_person_time = person_time * ipcw_str_nonlag,
    ipcw_str_nonlag_event = event_at_tstop * ipcw_str_nonlag,
    iptw_ipcw_str_nonlag_person_time = person_time * iptw * ipcw_str_nonlag,
    iptw_ipcw_str_nonlag_event = event_at_tstop * iptw * ipcw_str_nonlag,

    # pooled non-lagged IPCW
    ipcw_pl_nonlag_person_time = person_time * ipcw_pl_nonlag,
    ipcw_pl_nonlag_event = event_at_tstop * ipcw_pl_nonlag,
    iptw_ipcw_pl_nonlag_person_time = person_time * iptw * ipcw_pl_nonlag,
    iptw_ipcw_pl_nonlag_event = event_at_tstop * iptw * ipcw_pl_nonlag,

    # pooled lagged IPCW
    ipcw_pl_lag_person_time = person_time * ipcw_pl_lag,
    ipcw_pl_lag_event = event_at_tstop * ipcw_pl_lag,
    iptw_ipcw_pl_lag_person_time = person_time * iptw * ipcw_pl_lag,
    iptw_ipcw_pl_lag_event = event_at_tstop * iptw * ipcw_pl_lag
  )]
  
  rm(d_long, ipcw_weights)
  
  # no weights
  ir_at <- at_d[, .(
    IR = sum(event_at_tstop) / (sum(person_time)/365), 
    IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
    IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at[, `:=` (
    IRR = IR_comp / IR_ref, 
    IRD = IR_comp - IR_ref
  )]
  
  ir_at <- data.frame(t(ir_at))
  colnames(ir_at) <- 'at'
  ir_at$row_names <- rownames(ir_at)
  ir_at <- ir_at %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at, by = 'row_names')
  rm(ir_at)
  
  # iptw 
  ir_at_iptw <- at_d[, .(
    IR = (sum(iptw_event) / sum(iptw_person_time)) * 365, 
    IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_iptw <- data.frame(t(ir_at_iptw))
  colnames(ir_at_iptw) <- 'at.iptw'
  ir_at_iptw$row_names <- rownames(ir_at_iptw)
  ir_at_iptw <- ir_at_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw, by = 'row_names')
  rm(ir_at_iptw)
  
  # stratified lagged weights
  ir_at_ipcw_str_lag <- at_d[, .(
    IR = (sum(ipcw_str_lag_event) / sum(ipcw_str_lag_person_time)) * 365, 
    IR_ref = sum(ipcw_str_lag_event[trt_dummy == 0])/sum(ipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_str_lag_event[trt_dummy == 1])/sum(ipcw_str_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_str_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_ipcw_str_lag <- data.frame(t(ir_at_ipcw_str_lag))
  colnames(ir_at_ipcw_str_lag) <- 'at.ipcw_str_lag'
  ir_at_ipcw_str_lag$row_names <- rownames(ir_at_ipcw_str_lag)
  ir_at_ipcw_str_lag <- ir_at_ipcw_str_lag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_str_lag, by = 'row_names')
  rm(ir_at_ipcw_str_lag)
  
  ir_at_iptw_ipcw_str_lag <- at_d[, .(
    IR = (sum(iptw_ipcw_str_lag_event) / sum(iptw_ipcw_str_lag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_str_lag_event[trt_dummy == 0])/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_str_lag_event[trt_dummy == 1])/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_str_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_iptw_ipcw_str_lag <- data.frame(t(ir_at_iptw_ipcw_str_lag))
  colnames(ir_at_iptw_ipcw_str_lag) <- 'at.iptw.ipcw_str_lag'
  ir_at_iptw_ipcw_str_lag$row_names <- rownames(ir_at_iptw_ipcw_str_lag)
  ir_at_iptw_ipcw_str_lag <- ir_at_iptw_ipcw_str_lag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_str_lag, by = 'row_names')
  rm(ir_at_iptw_ipcw_str_lag)
  
  # stratified non-lagged weights
  ir_at_ipcw_str_nonlag <- at_d[, .(
    IR = (sum(ipcw_str_nonlag_event) / sum(ipcw_str_nonlag_person_time)) * 365, 
    IR_ref = sum(ipcw_str_nonlag_event[trt_dummy == 0])/sum(ipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_str_nonlag_event[trt_dummy == 1])/sum(ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_str_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_ipcw_str_nonlag <- data.frame(t(ir_at_ipcw_str_nonlag))
  colnames(ir_at_ipcw_str_nonlag) <- 'at.ipcw_str_nonlag'
  ir_at_ipcw_str_nonlag$row_names <- rownames(ir_at_ipcw_str_nonlag)
  ir_at_ipcw_str_nonlag <- ir_at_ipcw_str_nonlag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_str_nonlag, by = 'row_names')
  rm(ir_at_ipcw_str_nonlag)
  
  ir_at_iptw_ipcw_str_nonlag <- at_d[, .(
    IR = (sum(iptw_ipcw_str_nonlag_event) / sum(iptw_ipcw_str_nonlag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 0])/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 1])/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_str_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_iptw_ipcw_str_nonlag <- data.frame(t(ir_at_iptw_ipcw_str_nonlag))
  colnames(ir_at_iptw_ipcw_str_nonlag) <- 'at.iptw.ipcw_str_nonlag'
  ir_at_iptw_ipcw_str_nonlag$row_names <- rownames(ir_at_iptw_ipcw_str_nonlag)
  ir_at_iptw_ipcw_str_nonlag <- ir_at_iptw_ipcw_str_nonlag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_str_nonlag, by = 'row_names')
  rm(ir_at_iptw_ipcw_str_nonlag)
  
  # pooled non-lagged weights
  ir_at_ipcw_pl_nonlag <- at_d[, .(
    IR = (sum(ipcw_pl_nonlag_event) / sum(ipcw_pl_nonlag_person_time)) * 365, 
    IR_ref = sum(ipcw_pl_nonlag_event[trt_dummy == 0])/sum(ipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_pl_nonlag_event[trt_dummy == 1])/sum(ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_pl_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_ipcw_pl_nonlag <- data.frame(t(ir_at_ipcw_pl_nonlag))
  colnames(ir_at_ipcw_pl_nonlag) <- 'at.ipcw_pl_nonlag'
  ir_at_ipcw_pl_nonlag$row_names <- rownames(ir_at_ipcw_pl_nonlag)
  ir_at_ipcw_pl_nonlag <- ir_at_ipcw_pl_nonlag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_pl_nonlag, by = 'row_names')
  rm(ir_at_ipcw_pl_nonlag)
  
  ir_at_iptw_ipcw_pl_nonlag <- at_d[, .(
    IR = (sum(iptw_ipcw_pl_nonlag_event) / sum(iptw_ipcw_pl_nonlag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 0])/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 1])/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_pl_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_iptw_ipcw_pl_nonlag <- data.frame(t(ir_at_iptw_ipcw_pl_nonlag))
  colnames(ir_at_iptw_ipcw_pl_nonlag) <- 'at.iptw.ipcw_pl_nonlag'
  ir_at_iptw_ipcw_pl_nonlag$row_names <- rownames(ir_at_iptw_ipcw_pl_nonlag)
  ir_at_iptw_ipcw_pl_nonlag <- ir_at_iptw_ipcw_pl_nonlag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_pl_nonlag, by = 'row_names')
  rm(ir_at_iptw_ipcw_pl_nonlag)
  
  # pooled lagged weights
  ir_at_ipcw_pl_lag <- at_d[, .(
    IR = (sum(ipcw_pl_lag_event) / sum(ipcw_pl_lag_person_time)) * 365, 
    IR_ref = sum(ipcw_pl_lag_event[trt_dummy == 0])/sum(ipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_pl_lag_event[trt_dummy == 1])/sum(ipcw_pl_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_pl_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_ipcw_pl_lag <- data.frame(t(ir_at_ipcw_pl_lag))
  colnames(ir_at_ipcw_pl_lag) <- 'at.ipcw_pl_lag'
  ir_at_ipcw_pl_lag$row_names <- rownames(ir_at_ipcw_pl_lag)
  ir_at_ipcw_pl_lag <- ir_at_ipcw_pl_lag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_pl_lag, by = 'row_names')
  rm(ir_at_ipcw_pl_lag)
  
  ir_at_iptw_ipcw_pl_lag <- at_d[, .(
    IR = (sum(iptw_ipcw_pl_lag_event) / sum(iptw_ipcw_pl_lag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_pl_lag_event[trt_dummy == 0])/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_pl_lag_event[trt_dummy == 1])/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_pl_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_iptw_ipcw_pl_lag <- data.frame(t(ir_at_iptw_ipcw_pl_lag))
  colnames(ir_at_iptw_ipcw_pl_lag) <- 'at.iptw.ipcw_pl_lag'
  ir_at_iptw_ipcw_pl_lag$row_names <- rownames(ir_at_iptw_ipcw_pl_lag)
  ir_at_iptw_ipcw_pl_lag <- ir_at_iptw_ipcw_pl_lag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_pl_lag, by = 'row_names')
  rm(ir_at_iptw_ipcw_pl_lag)
  gc()
  
  ### 4. RESULTS ###
  
  ir_results <- ir_results[,-1]
  
  result <- c(
    IR = (ir_results[1,]),
    IR_ref = (ir_results[2,]),
    IR_comp = (ir_results[3,]),
    IRR = (ir_results[4,]),
    IRD = (ir_results[5,])
  )
  
  result <- unlist(result)
  
  # if use boot::boot, save manually along the way (in case of crash)
  bs_results_manual <<- bind_rows(bs_results_manual, result)
  saveRDS(bs_results_manual, file = paste(path_final_res, 'bs_results_manual_nodeaths.rds', sep='/'))
  
  # if extracting CIs using manual save, can get results for main sample
  #bs_results_manual_t0 <<- rbind(bs_results_manual_t0, result)
  #saveRDS(bs_results_manual_t0, file = paste(path_final_res, 'bs_results_manual_t0_nodeaths.rds', sep='/'))
  
  return(result)
  
}

# get column names to store results
# incidence rates(rows): 5 (IR, IR_ref, IR_comp, ...)
# anlaysis types: 12 (itt, itt.iptw, at, ...)
# total columns per measure: 12*5 = 60

# col_names <- names(ir_results)
# tags <- c(
#   rep('IR.', each = 12),
#   rep('IR_ref.', each = 12),
#   rep('IR_comp.', each = 12),
#   rep('IRR.', each = 12),
#   rep('IRD.', each = 12)
# )
# 
# bs_col_names <- rep(col_names, times = 5)
# bs_col_names <- paste(tags, bs_col_names, sep = '')
# 
# saveRDS(bs_col_names, file = paste(path_intermediate_res_main, 'bs_col_names_nodeaths.rds', sep = '/'))

#### BOOTSTRAP EXECUTION USING BOOT::BOOT ####

set.seed(1)
R = 1000

bs_col_names <- readRDS(file = paste(path_intermediate_res_main, 'bs_col_names_nodeaths.rds', sep='/'))
bs_results_manual <- data.frame(matrix(nrow = 1, ncol = length(bs_col_names)))
names(bs_results_manual) <- bs_col_names

# time_to_start <- force_tz(ymd_hms("2024-11-28 00:30:00"), tzone = "America/Montreal")
# while (Sys.time() < time_to_start) {}

# run either (1) or (2)
# (1) without clusters
system.time({
  bs_results_boot <- boot(
    data = cohort,
    statistic = bs,
    R = R,
    parallel = "snow"
  )
})

bs_results_manual <- bs_results_manual[-1,]
saveRDS(bs_results_manual, file = paste(path_final_res, 'bs_results_manual_nodeaths.rds', sep='/'))
saveRDS(bs_results_boot, file = paste(path_final_res, 'bs_results_boot_nodeaths.rds', sep='/'))

# (2) with clusters
num_cores <- 2
cluster <- makeCluster(num_cores)

clusterExport(
  cluster,
  list(
    "bs",
    "cohort",
    "comorbidities",
    "fit_and_predict",
    "fit_and_predict_pooled",
    "breaks",
    "p_uncens_predictors",
    "p_uncens_predictors_pooled",
    "predictors",
    "bs_results_manual",
    "path_final_res",
    "exposure",
    "first_comorb",
    "cat_variables"
  )
)

clusterEvalQ(cluster, {
 library(dplyr)
 library(magrittr)
 library(survival)
 library(lubridate)
 library(purrr)
 library(data.table)
 library(fastglm)
})

system.time({
 bs_results_boot <- boot(
   data = cohort,
   statistic = bs,
   R = R,
   parallel = "snow",
   ncpus = num_cores,
   cl = cluster
 )

})

stopCluster(cluster)

# bs_results_manual <- bs_results_manual[-1,]
# saveRDS(bs_results_manual, file = paste(path_final_res, 'bs_results_manual_nodeaths.rds', sep='/'))
saveRDS(bs_results_boot, file = paste(path_final_res, 'bs_results_boot_nodeaths.rds', sep='/'))

# just to get t0 (run only inside bs fx once using cohort instead of sample)
# bs_results_manual_t0 <- data.frame(matrix(nrow = 1, ncol = length(bs_col_names)))
# names(bs_results_manual_t0) <- (bs_col_names)

#### BUILD BOOTSTRAPPED CIs - USING BOOT::BOOT.CI ####

bs_results_boot <- readRDS (file = paste(path_final_res, 'bs_results_boot_nodeaths.rds', sep = '/'))
bs_ci_boot <- data.frame(matrix(nrow = 0, ncol = 3))
bs_col_names <- readRDS(file = paste(path_intermediate_res_main, 'bs_col_names_nodeaths.rds', sep='/'))
colnames(bs_ci_boot) <- c('estimate', 'lower_ci', 'upper_ci')

for (i in 1:length(bs_results_boot$t0)) {
  stat <- names(bs_results_boot[[1]][i])
  stat
  
  stat_ci <- data.frame(estimate = NA,
                        lower_ci = NA,
                        upper_ci = NA)
  
  conf.int <- boot.ci(bs_results_boot, type = 'perc', index = i)
  CI <- conf.int$percent
  
  stat_ci$estimate <- bs_results_boot$t0[i]
  stat_ci$lower_ci <- t(CI[, 4])[[1]]
  stat_ci$upper_ci <- t(CI[, 5])[[1]]
  bs_ci_boot <- rbind(bs_ci_boot, stat_ci)
  rownames(bs_ci_boot)[i] <- stat
}
row.names(bs_ci_boot) <- bs_col_names
saveRDS(bs_ci_boot, file = paste(path_final_res, 'bs_ci_boot_nodeaths.rds', sep='/'))

## Extract CIs for overall IRR 

rows_to_keep <- c(
  
  'IRR.itt',
  'IRD.itt',
  'IRR.itt.iptw',
  'IRD.itt.iptw',
  
  'IRR.at',
  'IRD.at',
  'IRR.at.iptw',
  'IRD.at.iptw',
  
  # stratified lagged weights
  'IRR.at.ipcw_str_lag',
  'IRD.at.ipcw_str_lag',
  'IRR_stab.at.ipcw_str_lag',
  'IRD_stab.at.ipcw_str_lag',
  'IRR.at.iptw.ipcw_str_lag',
  'IRD.at.iptw.ipcw_str_lag',
  
  # stratified non-lagged weights
  'IRR.at.ipcw_str_nonlag',
  'IRD.at.ipcw_str_nonlag',
  'IRR_stab.at.ipcw_str_nonlag',
  'IRD_stab.at.ipcw_str_nonlag',
  'IRR.at.iptw.ipcw_str_nonlag',
  'IRD.at.iptw.ipcw_str_nonlag',
  
  # pooled non-lagged weights
  'IRR.at.ipcw_pl_nonlag',
  'IRD.at.ipcw_pl_nonlag',
  'IRR_stab.at.ipcw_pl_nonlag',
  'IRD_stab.at.ipcw_pl_nonlag',
  'IRR.at.iptw.ipcw_pl_nonlag',
  'IRD.at.iptw.ipcw_pl_nonlag',
  
  # pooled lagged weights
  'IRR.at.ipcw_pl_lag',
  'IRD.at.ipcw_pl_lag',
  'IRR_stab.at.ipcw_pl_lag',
  'IRD_stab.at.ipcw_pl_lag',
  'IRR.at.iptw.ipcw_pl_lag',
  'IRD.at.iptw.ipcw_pl_lag'
)

ir_ci_boot <- bs_ci_boot[row.names(bs_ci_boot) %in% rows_to_keep,]
ir_ci_boot$variable <- row.names(ir_ci_boot)
ir_ci_boot <- ir_ci_boot %>%
  select(variable, everything())

saveRDS(ir_ci_boot, paste(path_final_res, 'ir_ci_boot_nodeaths.rds', sep ='/'))
write_xlsx(ir_ci_boot, paste(path_final_res, 'ir_ci_boot_nodeaths.xlsx', sep ='/'))

#### BUILD BOOTSTRAPPED CIs - USING MANUAL SAVE ####

# if results saved in different chunks can combine them
save1 <- readRDS(paste(path_final_res, 'bs_results_manual_nodeaths1.rds', sep='/'))
save2 <- readRDS(paste(path_final_res, 'bs_results_manual_nodeaths2.rds', sep='/'))
results <- rbind(save1, save2)

results %<>% filter(!is.na(IR.itt))
saveRDS(results, file = paste(path_final_res, 'bs_results_manual_nodeaths_final.rds', sep = '/'))
# results <- results[-1,]

# results <- readRDS(paste(path_final_res, 'bs_results_manual_nodeaths.rds', sep='/'))
results_t0 <- readRDS(paste(path_final_res, 'bs_results_manual_t0_nodeaths.rds', sep='/'))

bs_ci_manual <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(bs_ci_manual) <- c('estimate', 'lower_ci', 'upper_ci')

get_stat_ci <- function(stat_name) {
  
  stat_ci <- data.frame(
    estimate = NA,
    lower_ci = NA,
    upper_ci = NA
  )
  
  stat_values <- results[stat_name]
  
  stat_ci$estimate <- results_t0[2, stat_name] # to get estimate from main sample
  
  #stat_ci$estimate <- quantile(stat_values, prob = c(0.5), na.rm = TRUE) # to get estimate from 50%
  stat_ci$lower_ci <- quantile(stat_values, prob = c(0.025), na.rm = TRUE)
  stat_ci$upper_ci <- quantile(stat_values, prob = c(0.975), na.rm = TRUE)
  
  rownames(stat_ci)[1] <- stat_name
  bs_ci_manual <<- rbind(bs_ci_manual, stat_ci)
  
}

stat_names <- names(results)
invisible(lapply(stat_names, get_stat_ci))

bs_ci_manual$variable <- row.names(bs_ci_manual)

bs_ci_manual <- bs_ci_manual %>%
  select(variable, everything())

saveRDS(bs_ci_manual, file = paste(path_final_res, 'bs_ci_manual_nodeaths.rds', sep='/'))

## extract CIs for overall IRR ##

rows_to_keep <- c(
  
  'IRR.itt',
  'IRD.itt',
  'IRR.itt.iptw',
  'IRD.itt.iptw',
  'IRR_stab.itt.iptw',
  'IRD_stab.itt.iptw',
  
  'IRR.at',
  'IRD.at',
  'IRR.at.iptw',
  'IRD.at.iptw',
  'IRR_stab.at.iptw',
  'IRD_stab.at.iptw',
  
  # stratified lagged weights
  'IRR.at.ipcw_str_lag',
  'IRD.at.ipcw_str_lag',
  'IRR_stab.at.ipcw_str_lag',
  'IRD_stab.at.ipcw_str_lag',
  'IRR.at.iptw.ipcw_str_lag',
  'IRD.at.iptw.ipcw_str_lag',
  'IRR_stab.at.iptw.ipcw_str_lag',
  'IRD_stab.at.iptw.ipcw_str_lag',
  
  # stratified non-lagged weights
  'IRR.at.ipcw_str_nonlag',
  'IRD.at.ipcw_str_nonlag',
  'IRR_stab.at.ipcw_str_nonlag',
  'IRD_stab.at.ipcw_str_nonlag',
  'IRR.at.iptw.ipcw_str_nonlag',
  'IRD.at.iptw.ipcw_str_nonlag',
  'IRR_stab.at.iptw.ipcw_str_nonlag',
  'IRD_stab.at.iptw.ipcw_str_nonlag',
  
  # pooled non-lagged weights
  'IRR.at.ipcw_pl_nonlag',
  'IRD.at.ipcw_pl_nonlag',
  'IRR_stab.at.ipcw_pl_nonlag',
  'IRD_stab.at.ipcw_pl_nonlag',
  'IRR.at.iptw.ipcw_pl_nonlag',
  'IRD.at.iptw.ipcw_pl_nonlag',
  'IRR_stab.at.iptw.ipcw_pl_nonlag',
  'IRD_stab.at.iptw.ipcw_pl_nonlag',
  
  # pooled lagged weights
  'IRR.at.ipcw_pl_lag',
  'IRD.at.ipcw_pl_lag',
  'IRR_stab.at.ipcw_pl_lag',
  'IRD_stab.at.ipcw_pl_lag',
  'IRR.at.iptw.ipcw_pl_lag',
  'IRD.at.iptw.ipcw_pl_lag',
  'IRR_stab.at.iptw.ipcw_pl_lag',
  'IRD_stab.at.iptw.ipcw_pl_lag'
)

ir_ci_manual <- bs_ci_manual[row.names(bs_ci_manual) %in% rows_to_keep,]
ir_ci_manual$variable <- row.names(ir_ci_manual)
ir_ci_manual <- ir_ci_manual %>%
  select(variable, everything())

write_xlsx(ir_ci_manual, paste(path_final_res, 'ir_ci_manual_nodeaths.xlsx', sep ='/'))

#### TEST FOR POSITIVITY VIOLATIONS ####
# see if within each decile there is a reasonable count
# of patients with certain covariates

cohort_long <- readRDS(paste(path_intermediate_res, 'cohort_ipcw.rds', sep = '/'))

# ref group
trt0 <- cohort_long %>% filter(trt_dummy == 0)
trt0_1 <- trt0 %>% filter(dec == 1)
summary(trt0_1[variables])

trt0_2 <- trt0 %>% filter(dec == 2)
summary(trt0_2[variables])

trt0_3 <- trt0 %>% filter(dec == 3)
summary(trt0_3[variables])

trt0_4 <- trt0 %>% filter(dec == 4)
summary(trt0_4[variables])

trt0_5 <- trt0 %>% filter(dec == 5)
summary(trt0_5[variables])

trt0_6 <- trt0 %>% filter(dec == 6)
summary(trt0_6[variables])

trt0_7 <- trt0 %>% filter(dec == 7)
summary(trt0_7[variables])

trt0_8 <- trt0 %>% filter(dec == 8)
summary(trt0_8[variables])

trt0_9 <- trt0 %>% filter(dec == 9)
summary(trt0_9[variables])

trt0_10 <- trt0 %>% filter(dec == 10)
summary(trt0_10[variables])

# comparator group
trt1 <- cohort_long %>% filter(trt_dummy == 1)
trt1_1 <- trt1 %>% filter(dec == 1)
summary(trt1_1[variables])

trt1_2 <- trt1 %>% filter(dec == 2)
summary(trt1_2[variables])

trt1_3 <- trt1 %>% filter(dec == 3)
summary(trt1_3[variables])

trt1_4 <- trt1 %>% filter(dec == 4)
summary(trt1_4[variables])

trt1_5 <- trt1 %>% filter(dec == 5)
summary(trt1_5[variables])

trt1_6 <- trt1 %>% filter(dec == 6)
summary(trt1_6[variables])

trt1_7 <- trt1 %>% filter(dec == 7)
summary(trt1_7[variables])

trt1_8 <- trt1 %>% filter(dec == 8)
summary(trt1_8[variables])

trt1_9 <- trt1 %>% filter(dec == 9)
summary(trt1_9[variables])

trt1_10 <- trt1 %>% filter(dec == 10)
summary(trt1_10[variables])
