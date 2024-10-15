## ---------------------------
##
## Program: 9. Bootstrapping 
##
## Purpose: Bootstrap cohort to obtain percentile confidence intervals for IR, IRR, and IRD.
##
## Author: Gwen Aubrac
##
## Date Created: 2024-08-30
##
## ---------------------------
##
## Notes: Confidence intervals constructed using boot (bootstrap_results_boot, bootstrap_results_manual to save along iterations).
## Computation speed was compared with future (bootstrap_results_future).
##
## ---------------------------

# analysis: flex_grace_period, 90_day_grace_period
# male, female, young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed

analysis <- ''

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

#### DEFINE PATHS ####

path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 

if (analysis == 'main' | analysis == '') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/main" 
} else if (analysis == 'male') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/male" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/male" 
} else if (analysis == 'female') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/female" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/female" 
} else if (analysis == 'young') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/young" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/young" 
} else if (analysis == 'old') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/old" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/old" 
} else if (analysis == '2019') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2019" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2019" 
} else if (analysis == '2020') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2020" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2020" 
} else if (analysis == '2021') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2021" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2021" 
} else if (analysis == '2022') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2022" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2022" 
} else if (analysis == 'flex_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/flex_grace_period" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/sensitivity/flex_grace_period" 
} else if (analysis == '90_day_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/90_day_grace_period" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/sensitivity/90_day_grace_period" 
} else if (analysis == 'depressed') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/depressed" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/depressed" 
} else if (analysis == 'not_depressed') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/not_depressed" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/not_depressed" 
} 

## Variables and dataframes used within boostrap

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))
covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))

# times_dec <- readRDS(paste(path_cohort, 'times_dec.rds', sep = '/'))
# times_dec <- times_dec[1:9] 
# times_dec

# breaks for the time intervals
breaks <- readRDS(file = paste(path_cohort, 'interval_breaks.rds', sep = '/'))

# remove variables that violate positivity assumption
variables <- c(covariates, comorbidities)
# variables <- variables[!variables %in% c('hypocalcemia', 'hypomagnesemia', 'acute_renal_disease')] 
variables <- variables[!variables %in% c('hypocalcemia', 'hypomagnesemia', 'acute_renal_disease', 'hypokalemia', 'cardiomyopathy')] 
# variables <- variables[!variables %in% c('hypocalcemia', 'hypomagnesemia', 'acute_renal_disease', 'ethnicity')] 
variables

model <- readRDS(file = paste(path_results, 'iptw_model.rds', sep = '/'))

# remove variables for subgroup analyses
if (analysis == 'male' | analysis == 'female') {
  variables <- variables[!variables %in% c('sex')]
} else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
  variables <- variables[!variables %in% c('year', 'month_year')]
} else if (analysis == 'depressed' | analysis == 'not_depressed') {
  variables <- variables[!variables %in% c('depression')]
} else if (analysis == 'young' | analysis == 'old') {
  cohort$age_group <- droplevels(cohort$age_group)
}
variables

## Functions used within bootstrap ##

predictors <- as.formula(paste('~', as.character(model)[-1][2], sep = ''))
predictors

p_uncens_predictors <- as.formula(paste("~", paste(c(variables, 'age_group*depression'), collapse = " + ")))
if (analysis == 'depressed' | analysis == 'not_depressed') {
  p_uncens_predictors <- as.formula(paste("~", paste(c(variables), collapse = " + ")))
}
p_uncens_predictors

p_uncens_predictors_pooled <- as.formula(paste("~", paste(c(variables, 'age_group*depression', 'dec'), collapse = " + ")))
if (analysis == 'depressed' | analysis == 'not_depressed') {
  p_uncens_predictors_pooled <- as.formula(paste("~", paste(c(variables, 'dec'), collapse = " + ")))
}
p_uncens_predictors_pooled

fit_and_predict <- function(data_subset) {
  x <- model.matrix(p_uncens_predictors, data = data_subset)
  y <- data_subset$uncensored_at_tstop
  model <- fastglm(x=x, y=y, family = binomial(link = "logit"))
  data_subset$p_uncens <- predict(model, newdata = x, type = "response")
  return(data_subset)
}

fit_and_predict_base <- function(data_subset) {
  x <- model.matrix(~one, data = data_subset)
  y <- data_subset$uncensored_at_tstop
  model <- fastglm(x=x, y=y, family = binomial(link = "logit"))
  data_subset$p_uncens_base <- predict(model, newdata = x, type = "response")
  return(data_subset)
}

fit_and_predict_pooled <- function(data_subset) {
  x <- model.matrix(p_uncens_predictors_pooled, data = data_subset)
  y <- data_subset$uncensored_at_tstop
  model <- fastglm(x=x, y=y, family = binomial(link = "logit"))
  data_subset$p_uncens <- predict(model, newdata = x, type = "response")
  return(data_subset)
}

rm(variables, covariates, model, analysis)

setDT(cohort)

#### BOOTSTRAP FUNCTION ####

# note: run one time inside fx to get result names

bs <- function(data, indices) {
  
  d <- data[indices,]
  d[, boot_id := .I]
  
  #### 1. CALCUALTE IPTW ####
  
  d$one <- 1
  y <- as.numeric(as.character(d$trt_dummy))
  x_denom <- model.matrix(predictors, data = d)
  x_num <- model.matrix(~one, data = d)
  
  denom_fit <- fastglm(x=x_denom, y=y, family = binomial(link = 'logit'))
  numer_fit <- fastglm(x=x_num, y=y, family = binomial(link = 'logit'))
  
  pn_trt <- predict(numer_fit, newdata = x_num, type = 'response')
  pd_trt <- predict(denom_fit, newdata = x_denom, type = 'response')
  
  d$iptw <- 1 / if_else(d$trt_dummy == 0, 1 - pd_trt, pd_trt)
  
  d$siptw <- if_else(d$trt_dummy == 0, 
                     ((1-pn_trt) / (1-pd_trt)),
                     pn_trt/pd_trt)
  
  rm(y, x_denom, x_num, denom_fit, numer_fit, pn_trt, pd_trt)
  
  #### 2. CALCULATE IPCW ####
  
  ## Prepare data ##
  d_long <- d
  
  d_long$counting_time <- as.numeric(difftime(d_long$at_exit_date, d_long$entry_date, units = 'days'))
  d_long$censor_counting_time <- as.numeric(difftime(d_long$censor_date, d_long$entry_date, units = 'days'))
  d_long$uncensored <- if_else(d_long$censor == 0, 1, 0)
  d_long <- d_long[order(d_long$counting_time),]
  
  # define interval breaks
  censored_only <- d %>%
    filter(!is.na(censor_date))
  
  censoring_times <- as.numeric(censored_only$censor_date - censored_only$entry_date)
  times_dec <- quantile(censoring_times, probs = breaks)
  
  if (analysis == 'main') {
    times_dec[[1]] <- 59
  } else if (analysis == '90_day_grace_period') {
    times_dec[[1]] <- 119
  } else if (analysis == 'flex_grace_period') {
    times_dec[[1]] <- 57
  }
  
  times_dec <- times_dec[1:9]
  
  d_long$Tstart <- 0
  d_long <- survSplit(d_long, cut = times_dec, end = 'counting_time', start = 'Tstart', event = 'uncensored')
  names(d_long)[names(d_long) == 'counting_time'] <- 'Tstop' 
  d_long$uncensored_at_tstop <- if_else(is.na(d_long$censor_counting_time), 1,
                                        if_else(d_long$Tstop == d_long$censor_counting_time, 0, 1)) 
  
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
  
  ## Stratified lagged and non-lagged weights ##
  
  d_long$one <- 1
  setDT(d_long)
  d_long <- d_long[, fit_and_predict(.SD), by = .(dec, trt_dummy)]
  d_long <- d_long[, fit_and_predict_base(.SD), by = .(dec, trt_dummy)]
  
  d_long[, `:=`(
    weight = 1 / p_uncens
  ), by = boot_id]
  
  d_long[, `:=`(
    ipcw_str_nonlag = cumprod(weight),
    sipcw_str_nonlag = ipcw_str_nonlag * p_uncens_base
  ), by = boot_id]
  
  d_long[, `:=`(
    ipcw_str_lag = shift(ipcw_str_nonlag, n = 1, type = "lag", fill = 1),
    sipcw_str_lag = ipcw_str_lag * p_uncens_base
  ), by = boot_id]
  
  d_long[, `:=`(weight = NULL, stab_weight = NULL, p_uncens = NULL)]
  
  gc()
  
  ## Pooled lagged and non-lagged weights ##
  
  d_long <- d_long[, fit_and_predict_pooled(.SD), by = .(trt_dummy)]
  
  d_long[, `:=`(
    weight = 1 / p_uncens,
  ), by = boot_id]
  
  d_long[, `:=`(
    ipcw_pl_nonlag = cumprod(weight),
    sipcw_pl_nonlag = ipcw_pl_nonlag * p_uncens_base
  ), by = boot_id]
  
  d_long[, `:=`(
    ipcw_pl_lag = shift(ipcw_pl_nonlag, n = 1, type = "lag", fill = 1),
    sipcw_pl_lag = ipcw_pl_lag * p_uncens_base
  ), by = boot_id]
  
  d_long %<>% select(-weight, -stab_weight, -p_uncens, -p_uncens_base)
  
  gc()
  
  #### 3. CALCULATE INCIDENCE RATES ####
  
  ## A. ITT analyses ##
  
  setDT(d)
  
  itt_d <- d[, .(boot_id, trt, trt_dummy, itt_event, itt_exit_date, entry_date, iptw, siptw, month_year)]
  itt_d[, itt_event := as.numeric(itt_event)]
  
  itt_d[, `:=`(
    person_time = as.numeric(itt_exit_date - entry_date),
    iptw_person_time = NA_real_,
    siptw_person_time = NA_real_, 
    iptw_event = NA_real_,
    siptw_event = NA_real_
  )]
  
  itt_d[, `:=`(
    iptw_person_time = person_time * iptw,
    siptw_person_time = person_time * siptw,
    iptw_event = itt_event * iptw,
    siptw_event = itt_event * siptw
  )]
  
  setDT(itt_d)
  
  ## over entire study
  
  # no weights
  ir_itt <- itt_d[, .(
    IR = sum(itt_event) / (sum(person_time) / 365),
    IR_ref = sum(itt_event == 1 & trt_dummy == 0) / sum(person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(itt_event == 1 & trt_dummy == 1) / sum(person_time[trt_dummy == 1]) * 365
  )]
  
  ir_itt[, `:=`(
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = 0,
    IR_ref_stab = 0,
    IR_comp_stab = 0,
    IRR_stab = 0,
    IRD_stab = 0
  )]
  
  ir_itt <- data.frame(t(ir_itt))
  colnames(ir_itt) <- 'itt'
  ir_itt$row_names <- rownames(ir_itt)
  ir_itt <- ir_itt %>% relocate(row_names)
  
  # iptw and siptw weights
  ir_itt_iptw <- itt_d[, .(
    IR = sum(iptw_event) / sum(iptw_person_time) * 365,
    IR_ref = sum(iptw_event[trt_dummy == 0]) / sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_event[trt_dummy == 1]) / sum(iptw_person_time[trt_dummy == 1]) * 365,
    IR_stab = sum(siptw_event) / sum(siptw_person_time) * 365,
    IR_ref_stab = sum(siptw_event[trt_dummy == 0]) / sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1]) / sum(siptw_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_itt_iptw <- ir_itt_iptw[, .(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_itt_iptw <- data.frame(t(ir_itt_iptw))
  colnames(ir_itt_iptw) <- 'itt.iptw'
  ir_itt_iptw$row_names <- rownames(ir_itt_iptw)
  ir_itt_iptw <- ir_itt_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_itt, ir_itt_iptw, by = 'row_names')
  rm(ir_itt, ir_itt_iptw)
  
  ## over time
  
  # no weights
  ir_itt_time <- itt_d[, .(
    IR = sum(itt_event) / (sum(person_time) / 365),
    IR_ref = sum(itt_event[trt_dummy == 0]) / sum(person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(itt_event[trt_dummy == 1]) / sum(person_time[trt_dummy == 1]) * 365,
    IR_stab = 0,
    IR_ref_stab = 0,
    IR_comp_stab = 0,
    IRR = NA_real_,
    IRD = NA_real_,
    IRR_stab = 0,
    IRD_stab = 0
  ), by = month_year]
  
  ir_itt_time[, `:=`(
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_itt_time <- data.frame(t(ir_itt_time))
  colnames(ir_itt_time) <- ir_itt_time[1,]
  ir_itt_time <- ir_itt_time[-1,]
  colnames(ir_itt_time) <- paste('itt', colnames(ir_itt_time), sep = '.')
  ir_itt_time$row_names <- rownames(ir_itt_time)
  ir_itt_time <- ir_itt_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_itt_time, by = 'row_names')
  rm(ir_itt_time)
  
  # iptw and siptw weights
  ir_itt_iptw_time <- itt_d [, .(
    IR = (sum(iptw_event) / sum(iptw_person_time)) * 365, 
    IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_event) / sum(siptw_person_time)) * 365, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_itt_iptw_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_itt_iptw_time <- data.frame(t(ir_itt_iptw_time))
  colnames(ir_itt_iptw_time) <- ir_itt_iptw_time[1,]
  ir_itt_iptw_time <- ir_itt_iptw_time[-1,]
  colnames(ir_itt_iptw_time) <- paste('itt', colnames(ir_itt_iptw_time), sep = '.')
  ir_itt_iptw_time$row_names <- rownames(ir_itt_iptw_time)
  ir_itt_iptw_time <- ir_itt_iptw_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_itt_iptw_time, by = 'row_names', suffix = c('', '.iptw'))
  rm(ir_itt_iptw_time)
  
  rm(itt_d)
  gc()
  
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
  ipcw_weights <- d_long[, .(boot_id, Tstart, dec, ipcw_str_lag, sipcw_str_lag, ipcw_str_nonlag, sipcw_str_nonlag, ipcw_pl_nonlag, sipcw_pl_nonlag, ipcw_pl_lag, sipcw_pl_lag)]
  setnames(ipcw_weights, old = "Tstart", new = "decTstart")
  at_d <- at_d[ipcw_weights, on = .(boot_id, decTstart), nomatch = 0]
  at_d <- at_d[, .(boot_id, Tstart, Tstop, event_at_tstop, iptw, siptw, ipcw_str_lag, sipcw_str_lag, ipcw_str_nonlag, sipcw_str_nonlag,
                   ipcw_pl_nonlag, sipcw_pl_nonlag, ipcw_pl_lag, sipcw_pl_lag, dec, trt, trt_dummy, month_year)]
  
  at_d[, person_time := Tstop - Tstart]
  
  at_d[, `:=`(
    # IPTW
    iptw_person_time = person_time * iptw,
    iptw_event = event_at_tstop * iptw,
    siptw_person_time = person_time * siptw,
    siptw_event = event_at_tstop * siptw,
    
    # stratified lagged IPCW
    ipcw_str_lag_person_time = person_time * ipcw_str_lag,
    ipcw_str_lag_event = event_at_tstop * ipcw_str_lag,
    sipcw_str_lag_person_time = person_time * sipcw_str_lag,
    sipcw_str_lag_event = event_at_tstop * sipcw_str_lag,
    iptw_ipcw_str_lag_person_time = person_time * iptw * ipcw_str_lag,
    iptw_ipcw_str_lag_event = event_at_tstop * iptw * ipcw_str_lag,
    siptw_sipcw_str_lag_person_time = person_time * siptw * sipcw_str_lag,
    siptw_sipcw_str_lag_event = event_at_tstop * siptw * sipcw_str_lag,
    
    # stratified non-lagged IPCW
    ipcw_str_nonlag_person_time = person_time * ipcw_str_nonlag,
    ipcw_str_nonlag_event = event_at_tstop * ipcw_str_nonlag,
    sipcw_str_nonlag_person_time = person_time * sipcw_str_nonlag,
    sipcw_str_nonlag_event = event_at_tstop * sipcw_str_nonlag,
    iptw_ipcw_str_nonlag_person_time = person_time * iptw * ipcw_str_nonlag,
    iptw_ipcw_str_nonlag_event = event_at_tstop * iptw * ipcw_str_nonlag,
    siptw_sipcw_str_nonlag_person_time = person_time * siptw * sipcw_str_nonlag,
    siptw_sipcw_str_nonlag_event = event_at_tstop * siptw * sipcw_str_nonlag,
    
    # pooled non-lagged IPCW
    ipcw_pl_nonlag_person_time = person_time * ipcw_pl_nonlag,
    ipcw_pl_nonlag_event = event_at_tstop * ipcw_pl_nonlag,
    sipcw_pl_nonlag_person_time = person_time * sipcw_pl_nonlag,
    sipcw_pl_nonlag_event = event_at_tstop * sipcw_pl_nonlag,
    iptw_ipcw_pl_nonlag_person_time = person_time * iptw * ipcw_pl_nonlag,
    iptw_ipcw_pl_nonlag_event = event_at_tstop * iptw * ipcw_pl_nonlag,
    siptw_sipcw_pl_nonlag_person_time = person_time * siptw * sipcw_pl_nonlag,
    siptw_sipcw_pl_nonlag_event = event_at_tstop * siptw * sipcw_pl_nonlag,
    
    # pooled lagged IPCW
    ipcw_pl_lag_person_time = person_time * ipcw_pl_lag,
    ipcw_pl_lag_event = event_at_tstop * ipcw_pl_lag,
    sipcw_pl_lag_person_time = person_time * sipcw_pl_lag,
    sipcw_pl_lag_event = event_at_tstop * sipcw_pl_lag,
    iptw_ipcw_pl_lag_person_time = person_time * iptw * ipcw_pl_lag,
    iptw_ipcw_pl_lag_event = event_at_tstop * iptw * ipcw_pl_lag,
    siptw_sipcw_pl_lag_person_time = person_time * siptw * sipcw_pl_lag,
    siptw_sipcw_pl_lag_event = event_at_tstop * siptw * sipcw_pl_lag
  )]
  
  rm(d_long, ipcw_weights)
  
  ## over entire study
  
  # no weights
  ir_at <- at_d[, .(
    IR = sum(event_at_tstop) / (sum(person_time)/365), 
    IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
    IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at[, `:=` (
    IRR = IR_comp / IR_ref, 
    IRD = IR_comp - IR_ref,
    IR_stab = 0, 
    IR_ref_stab = 0,
    IR_comp_stab = 0,
    IRR_stab = 0,
    IRD_stab = 0
  )]
  
  ir_at <- data.frame(t(ir_at))
  colnames(ir_at) <- 'at'
  ir_at$row_names <- rownames(ir_at)
  ir_at <- ir_at %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at, by = 'row_names')
  rm(ir_at)
  
  # iptw and siptw weights
  ir_at_iptw <- at_d[, .(
    IR = (sum(iptw_event) / sum(iptw_person_time)) * 365, 
    IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_event) / sum(siptw_person_time)) * 365, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(ipcw_str_lag_event[trt_dummy == 1])/sum(ipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_str_lag_event) / sum(sipcw_str_lag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_str_lag_event[trt_dummy == 0])/sum(sipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_str_lag_event[trt_dummy == 1])/sum(sipcw_str_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_str_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(iptw_ipcw_str_lag_event[trt_dummy == 1])/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_str_lag_event) / sum(siptw_sipcw_str_lag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_str_lag_event[trt_dummy == 0])/sum(siptw_sipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_str_lag_event[trt_dummy == 1])/sum(siptw_sipcw_str_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_str_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(ipcw_str_nonlag_event[trt_dummy == 1])/sum(ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_str_nonlag_event) / sum(sipcw_str_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_str_nonlag_event[trt_dummy == 0])/sum(sipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_str_nonlag_event[trt_dummy == 1])/sum(sipcw_str_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_str_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 1])/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_str_nonlag_event) / sum(siptw_sipcw_str_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_str_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_str_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_str_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_str_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(ipcw_pl_nonlag_event[trt_dummy == 1])/sum(ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_pl_nonlag_event) / sum(sipcw_pl_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_pl_nonlag_event[trt_dummy == 0])/sum(sipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pl_nonlag_event[trt_dummy == 1])/sum(sipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_pl_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 1])/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_pl_nonlag_event) / sum(siptw_sipcw_pl_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_pl_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pl_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_pl_nonlag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(ipcw_pl_lag_event[trt_dummy == 1])/sum(ipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_pl_lag_event) / sum(sipcw_pl_lag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_pl_lag_event[trt_dummy == 0])/sum(sipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pl_lag_event[trt_dummy == 1])/sum(sipcw_pl_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_ipcw_pl_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
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
    IR_comp = sum(iptw_ipcw_pl_lag_event[trt_dummy == 1])/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_pl_lag_event) / sum(siptw_sipcw_pl_lag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_pl_lag_event[trt_dummy == 0])/sum(siptw_sipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pl_lag_event[trt_dummy == 1])/sum(siptw_sipcw_pl_lag_person_time[trt_dummy == 1]) * 365
  )]
  
  ir_at_iptw_ipcw_pl_lag[, `:=` (
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_iptw_ipcw_pl_lag <- data.frame(t(ir_at_iptw_ipcw_pl_lag))
  colnames(ir_at_iptw_ipcw_pl_lag) <- 'at.iptw.ipcw_pl_lag'
  ir_at_iptw_ipcw_pl_lag$row_names <- rownames(ir_at_iptw_ipcw_pl_lag)
  ir_at_iptw_ipcw_pl_lag <- ir_at_iptw_ipcw_pl_lag %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_pl_lag, by = 'row_names')
  rm(ir_at_iptw_ipcw_pl_lag)
  
  ## over time
  
  # no weights
  ir_at_time <- at_d[, .(
    IR = sum(event_at_tstop) / (sum(person_time) / 365),
    IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
    IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365, 
    IR_stab = 0,
    IR_ref_stab = 0,
    IR_comp_stab = 0,
    IRR = NA_real_,
    IRD = NA_real_,
    IRR_stab = 0,
    IRD_stab = 0
  ), by = month_year]
  
  ir_at_time[, `:=`(
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref
  )]
  
  ir_at_time <- data.frame(t(ir_at_time))
  colnames(ir_at_time) <- ir_at_time[1,]
  ir_at_time <- ir_at_time[-1,]
  colnames(ir_at_time) <- paste('at', colnames(ir_at_time), sep = '.')
  ir_at_time$row_names <- rownames(ir_at_time)
  ir_at_time <- ir_at_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_time, by = 'row_names')
  rm(ir_at_time)
  
  # iptw and siptw weights
  ir_at_iptw_time <- at_d [, .(
    IR = (sum(iptw_event) / sum(iptw_person_time)) * 365, 
    IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_event) / sum(siptw_person_time)) * 365, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_iptw_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_iptw_time <- data.frame(t(ir_at_iptw_time))
  colnames(ir_at_iptw_time) <- ir_at_iptw_time[1,]
  ir_at_iptw_time <- ir_at_iptw_time[-1,]
  colnames(ir_at_iptw_time) <- paste('at', colnames(ir_at_iptw_time), sep = '.')
  ir_at_iptw_time$row_names <- rownames(ir_at_iptw_time)
  ir_at_iptw_time <- ir_at_iptw_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_time, by = 'row_names', suffix = c('', '.iptw'))
  rm(ir_at_iptw_time)
  
  # stratified lagged weights
  ir_at_ipcw_str_lag_time <- at_d [, .(
    IR = (sum(ipcw_str_lag_event) / sum(ipcw_str_lag_person_time)) * 365, 
    IR_ref = sum(ipcw_str_lag_event[trt_dummy == 0])/sum(ipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_str_lag_event[trt_dummy == 1])/sum(ipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_str_lag_event) / sum(sipcw_str_lag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_str_lag_event[trt_dummy == 0])/sum(sipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_str_lag_event[trt_dummy == 1])/sum(sipcw_str_lag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_ipcw_str_lag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_ipcw_str_lag_time <- data.frame(t(ir_at_ipcw_str_lag_time))
  colnames(ir_at_ipcw_str_lag_time) <- ir_at_ipcw_str_lag_time[1,]
  ir_at_ipcw_str_lag_time <- ir_at_ipcw_str_lag_time[-1,]
  colnames(ir_at_ipcw_str_lag_time) <- paste('at', colnames(ir_at_ipcw_str_lag_time), sep = '.')
  ir_at_ipcw_str_lag_time$row_names <- rownames(ir_at_ipcw_str_lag_time)
  ir_at_ipcw_str_lag_time <- ir_at_ipcw_str_lag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_str_lag_time, by = 'row_names', suffix = c('', '.ipcw_str_lag'))
  rm(ir_at_ipcw_str_lag_time)
  
  ir_at_iptw_ipcw_str_lag_time <- at_d [, .(
    IR = (sum(iptw_ipcw_str_lag_event) / sum(iptw_ipcw_str_lag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_str_lag_event[trt_dummy == 0])/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_str_lag_event[trt_dummy == 1])/sum(iptw_ipcw_str_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_str_lag_event) / sum(siptw_sipcw_str_lag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_str_lag_event[trt_dummy == 0])/sum(siptw_sipcw_str_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_str_lag_event[trt_dummy == 1])/sum(siptw_sipcw_str_lag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_iptw_ipcw_str_lag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_iptw_ipcw_str_lag_time <- data.frame(t(ir_at_iptw_ipcw_str_lag_time))
  colnames(ir_at_iptw_ipcw_str_lag_time) <- ir_at_iptw_ipcw_str_lag_time[1,]
  ir_at_iptw_ipcw_str_lag_time <- ir_at_iptw_ipcw_str_lag_time[-1,]
  colnames(ir_at_iptw_ipcw_str_lag_time) <- paste('at', colnames(ir_at_iptw_ipcw_str_lag_time), sep = '.')
  ir_at_iptw_ipcw_str_lag_time$row_names <- rownames(ir_at_iptw_ipcw_str_lag_time)
  ir_at_iptw_ipcw_str_lag_time <- ir_at_iptw_ipcw_str_lag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_str_lag_time, by = 'row_names', suffix = c('', '.iptw.ipcw_str_lag'))
  rm(ir_at_iptw_ipcw_str_lag_time)
  
  # stratified non-lagged weights
  ir_at_ipcw_str_nonlag_time <- at_d [, .(
    IR = (sum(ipcw_str_nonlag_event) / sum(ipcw_str_nonlag_person_time)) * 365, 
    IR_ref = sum(ipcw_str_nonlag_event[trt_dummy == 0])/sum(ipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_str_nonlag_event[trt_dummy == 1])/sum(ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_str_nonlag_event) / sum(sipcw_str_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_str_nonlag_event[trt_dummy == 0])/sum(sipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_str_nonlag_event[trt_dummy == 1])/sum(sipcw_str_nonlag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_ipcw_str_nonlag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_ipcw_str_nonlag_time <- data.frame(t(ir_at_ipcw_str_nonlag_time))
  colnames(ir_at_ipcw_str_nonlag_time) <- ir_at_ipcw_str_nonlag_time[1,]
  ir_at_ipcw_str_nonlag_time <- ir_at_ipcw_str_nonlag_time[-1,]
  colnames(ir_at_ipcw_str_nonlag_time) <- paste('at', colnames(ir_at_ipcw_str_nonlag_time), sep = '.')
  ir_at_ipcw_str_nonlag_time$row_names <- rownames(ir_at_ipcw_str_nonlag_time)
  ir_at_ipcw_str_nonlag_time <- ir_at_ipcw_str_nonlag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_str_nonlag_time, by = 'row_names', suffix = c('', '.ipcw_str_nonlag'))
  rm(ir_at_ipcw_str_nonlag_time)
  
  ir_at_iptw_ipcw_str_nonlag_time <- at_d [, .(
    IR = (sum(iptw_ipcw_str_nonlag_event) / sum(iptw_ipcw_str_nonlag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 0])/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_str_nonlag_event[trt_dummy == 1])/sum(iptw_ipcw_str_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_str_nonlag_event) / sum(siptw_sipcw_str_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_str_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_str_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_str_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_str_nonlag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_iptw_ipcw_str_nonlag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_iptw_ipcw_str_nonlag_time <- data.frame(t(ir_at_iptw_ipcw_str_nonlag_time))
  colnames(ir_at_iptw_ipcw_str_nonlag_time) <- ir_at_iptw_ipcw_str_nonlag_time[1,]
  ir_at_iptw_ipcw_str_nonlag_time <- ir_at_iptw_ipcw_str_nonlag_time[-1,]
  colnames(ir_at_iptw_ipcw_str_nonlag_time) <- paste('at', colnames(ir_at_iptw_ipcw_str_nonlag_time), sep = '.')
  ir_at_iptw_ipcw_str_nonlag_time$row_names <- rownames(ir_at_iptw_ipcw_str_nonlag_time)
  ir_at_iptw_ipcw_str_nonlag_time <- ir_at_iptw_ipcw_str_nonlag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_str_nonlag_time, by = 'row_names', suffix = c('', '.iptw.ipcw_str_nonlag'))
  rm(ir_at_iptw_ipcw_str_nonlag_time)
  
  # pooled non-lagged weights
  ir_at_ipcw_pl_nonlag_time <- at_d [, .(
    IR = (sum(ipcw_pl_nonlag_event) / sum(ipcw_pl_nonlag_person_time)) * 365, 
    IR_ref = sum(ipcw_pl_nonlag_event[trt_dummy == 0])/sum(ipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_pl_nonlag_event[trt_dummy == 1])/sum(ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_pl_nonlag_event) / sum(sipcw_pl_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_pl_nonlag_event[trt_dummy == 0])/sum(sipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pl_nonlag_event[trt_dummy == 1])/sum(sipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_ipcw_pl_nonlag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_ipcw_pl_nonlag_time <- data.frame(t(ir_at_ipcw_pl_nonlag_time))
  colnames(ir_at_ipcw_pl_nonlag_time) <- ir_at_ipcw_pl_nonlag_time[1,]
  ir_at_ipcw_pl_nonlag_time <- ir_at_ipcw_pl_nonlag_time[-1,]
  colnames(ir_at_ipcw_pl_nonlag_time) <- paste('at', colnames(ir_at_ipcw_pl_nonlag_time), sep = '.')
  ir_at_ipcw_pl_nonlag_time$row_names <- rownames(ir_at_ipcw_pl_nonlag_time)
  ir_at_ipcw_pl_nonlag_time <- ir_at_ipcw_pl_nonlag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_pl_nonlag_time, by = 'row_names', suffix = c('', '.ipcw_pl_nonlag'))
  rm(ir_at_ipcw_pl_nonlag_time)
  
  ir_at_iptw_ipcw_pl_nonlag_time <- at_d [, .(
    IR = (sum(iptw_ipcw_pl_nonlag_event) / sum(iptw_ipcw_pl_nonlag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 0])/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_pl_nonlag_event[trt_dummy == 1])/sum(iptw_ipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_pl_nonlag_event) / sum(siptw_sipcw_pl_nonlag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_pl_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_pl_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pl_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_pl_nonlag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_iptw_ipcw_pl_nonlag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_iptw_ipcw_pl_nonlag_time <- data.frame(t(ir_at_iptw_ipcw_pl_nonlag_time))
  colnames(ir_at_iptw_ipcw_pl_nonlag_time) <- ir_at_iptw_ipcw_pl_nonlag_time[1,]
  ir_at_iptw_ipcw_pl_nonlag_time <- ir_at_iptw_ipcw_pl_nonlag_time[-1,]
  colnames(ir_at_iptw_ipcw_pl_nonlag_time) <- paste('at', colnames(ir_at_iptw_ipcw_pl_nonlag_time), sep = '.')
  ir_at_iptw_ipcw_pl_nonlag_time$row_names <- rownames(ir_at_iptw_ipcw_pl_nonlag_time)
  ir_at_iptw_ipcw_pl_nonlag_time <- ir_at_iptw_ipcw_pl_nonlag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_pl_nonlag_time, by = 'row_names', suffix = c('', '.iptw.ipcw_pl_nonlag'))
  rm(ir_at_iptw_ipcw_pl_nonlag_time)
  
  # pooled lagged weights
  ir_at_ipcw_pl_lag_time <- at_d [, .(
    IR = (sum(ipcw_pl_lag_event) / sum(ipcw_pl_lag_person_time)) * 365, 
    IR_ref = sum(ipcw_pl_lag_event[trt_dummy == 0])/sum(ipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(ipcw_pl_lag_event[trt_dummy == 1])/sum(ipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(sipcw_pl_lag_event) / sum(sipcw_pl_lag_person_time)) * 365, 
    IR_ref_stab = sum(sipcw_pl_lag_event[trt_dummy == 0])/sum(sipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pl_lag_event[trt_dummy == 1])/sum(sipcw_pl_lag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_ipcw_pl_lag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_ipcw_pl_lag_time <- data.frame(t(ir_at_ipcw_pl_lag_time))
  colnames(ir_at_ipcw_pl_lag_time) <- ir_at_ipcw_pl_lag_time[1,]
  ir_at_ipcw_pl_lag_time <- ir_at_ipcw_pl_lag_time[-1,]
  colnames(ir_at_ipcw_pl_lag_time) <- paste('at', colnames(ir_at_ipcw_pl_lag_time), sep = '.')
  ir_at_ipcw_pl_lag_time$row_names <- rownames(ir_at_ipcw_pl_lag_time)
  ir_at_ipcw_pl_lag_time <- ir_at_ipcw_pl_lag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_ipcw_pl_lag_time, by = 'row_names', suffix = c('', '.ipcw_pl_lag'))
  rm(ir_at_ipcw_pl_lag_time)
  
  ir_at_iptw_ipcw_pl_lag_time <- at_d [, .(
    IR = (sum(iptw_ipcw_pl_lag_event) / sum(iptw_ipcw_pl_lag_person_time)) * 365, 
    IR_ref = sum(iptw_ipcw_pl_lag_event[trt_dummy == 0])/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = sum(iptw_ipcw_pl_lag_event[trt_dummy == 1])/sum(iptw_ipcw_pl_lag_person_time[trt_dummy == 1]) * 365,
    IR_stab = (sum(siptw_sipcw_pl_lag_event) / sum(siptw_sipcw_pl_lag_person_time)) * 365, 
    IR_ref_stab = sum(siptw_sipcw_pl_lag_event[trt_dummy == 0])/sum(siptw_sipcw_pl_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pl_lag_event[trt_dummy == 1])/sum(siptw_sipcw_pl_lag_person_time[trt_dummy == 1]) * 365
  ), by = month_year]
  
  ir_at_iptw_ipcw_pl_lag_time[, `:=`(
    IR = IR,
    IR_ref = IR_ref,
    IR_comp = IR_comp,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    IR_stab = IR_stab,
    IR_ref_stab = IR_ref_stab,
    IR_comp_stab = IR_comp_stab,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )]
  
  ir_at_iptw_ipcw_pl_lag_time <- data.frame(t(ir_at_iptw_ipcw_pl_lag_time))
  colnames(ir_at_iptw_ipcw_pl_lag_time) <- ir_at_iptw_ipcw_pl_lag_time[1,]
  ir_at_iptw_ipcw_pl_lag_time <- ir_at_iptw_ipcw_pl_lag_time[-1,]
  colnames(ir_at_iptw_ipcw_pl_lag_time) <- paste('at', colnames(ir_at_iptw_ipcw_pl_lag_time), sep = '.')
  ir_at_iptw_ipcw_pl_lag_time$row_names <- rownames(ir_at_iptw_ipcw_pl_lag_time)
  ir_at_iptw_ipcw_pl_lag_time <- ir_at_iptw_ipcw_pl_lag_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, ir_at_iptw_ipcw_pl_lag_time, by = 'row_names', suffix = c('', '.iptw.ipcw_pl_lag'))
  rm(ir_at_iptw_ipcw_pl_lag_time)
  
  rm(at_d)
  gc()
  
  ### 4. RESULTS ###
  
  ir_results <- ir_results[,-1]
  
  result <- c(
    IR = as.numeric(ir_results[1,]),
    IR_ref = as.numeric(ir_results[2,]),
    IR_comp = as.numeric(ir_results[3,]),
    IRR = as.numeric(ir_results[4,]),
    IRD = as.numeric(ir_results[5,]),
    IR_stab = as.numeric(ir_results[6,]),
    IR_ref_stab = as.numeric(ir_results[7,]),
    IR_comp_stab = as.numeric(ir_results[8,]),
    IRR_stab = as.numeric(ir_results[9,]),
    IRD_stab = as.numeric(ir_results[10,])
  )
  
  # if use boot:
  bootstrap_results_manual <<- rbind(bootstrap_results_manual, result)
  saveRDS(bootstrap_results_manual, file = paste(path_results, 'bootstrap_results_manual_save.rds', sep='/'))
  
  #bootstrap_results_manual_t0 <<- rbind(bootstrap_results_manual_t0, result)
  #saveRDS(bootstrap_results_manual_t0, file = paste(path_results, 'bootstrap_results_manual_t0.rds', sep='/'))
  
  return(result)
  
}

# get result names
# anlaysis types: 12 (itt, itt.iptw, at, etc...)
# time results: 12 * 48 (itt_time, itt.iptw_time, etc...)
# total columns per measure: 12 + 12*48 = 588
# number of measures: 10 (IR, IR_ref, IR_comp, etc...)

# names <- names(ir_results)
# endings <- c(
#   rep('.IR', each = 588),
#   rep('.IR_ref', each = 588),
#   rep('.IR_comp', each = 588),
#   rep('.IRR', each = 588),
#   rep('.IRD', each = 588),
#   rep('.IR_stab', each = 588),
#   rep('.IR_ref_stab', each = 588),
#   rep('.IR_comp_stab', each = 588),
#   rep('.IRR_stab', each = 588),
#   rep('.IRD_stab', each = 588)
#   )
# 
# new_names <- rep(names, times = 10)
# new_names <- paste(new_names, endings, sep = '')
# 
# saveRDS(new_names, file = paste(path_main, 'result_names.rds', sep = '/'))

#### BOOTSTRAP EXECUTION USING BOOT ####

set.seed(1)
R = 10

result_names <- readRDS(file = paste(path_main, 'result_names.rds', sep='/'))
bootstrap_results_manual <- data.frame(matrix(nrow = 1, ncol = length(result_names)))
names(bootstrap_results_manual) <- result_names

# run either of (1) or (2)
# (1) without clusters - best
system.time({
  bootstrap_results_boot <- boot(
    data = cohort,
    statistic = bs,
    R = R,
    parallel = "snow"
  )
})

test <- bootstrap_results_manual %>% select(at.iptw.ipcw_str_nonlag.IRR_stab)
summary(test)

# (2) with clusters - does NOT work on large cohort (>50,000 obs)
# num_cores <- 2
# cluster <- makeCluster(num_cores)
# 
# clusterExport(cluster, list("bs", "cohort", "comorbidities", "fit_and_predict", "fit_and_predict_base", "fit_and_predict_pooled",
#                       "times_dec", "p_uncens_predictors", "p_uncens_predictors_pooled", "predictors", "bootstrap_results_manual", "path_results"))
# 
# clusterEvalQ(cluster, {
#  library(dplyr)
#  library(magrittr)
#  library(survival)
#  library(lubridate)
#  library(purrr)
#  library(data.table)
#  library(fastglm)
# })
# 
# system.time({
#  bootstrap_results_boot <- boot(
#    data = cohort,
#    statistic = bs,
#    R = R,
#    parallel = "snow",
#    ncpus = num_cores,
#    cl = cluster
#  )
#  
# })
# 
# stopCluster(cluster)

# test <- bootstrap_results_manual %>% select(at.iptw.ipcw_str_nonlag.IRR_stab)
# summary(test)

bootstrap_results_manual <- bootstrap_results_manual[-1,]

test <- bootstrap_results_manual %>%
 select(itt.IRR, itt.iptw.IRR, at.IRR, at.iptw.IRR, at.iptw.ipcw_str_lag.IRR, at.iptw.ipcw_str_nonlag.IRR, at.iptw.ipcw_pl_nonlag.IRR, at.iptw.ipcw_pl_lag.IRR)

saveRDS(bootstrap_results_manual, file = paste(path_results, 'bootstrap_results_manual.rds', sep='/'))
saveRDS(bootstrap_results_boot, file = paste(path_results, 'bootstrap_results_boot.rds', sep='/'))

# just to get t0 (run only inside bs fx once using cohort instead of sample)
# bootstrap_results_manual_t0 <- data.frame(matrix(nrow = 1, ncol = length(result_names)))
# names(bootstrap_results_manual_t0) <- result_names

#### BUILD BOOTSTRAPPED CIs - USING BOOT.CI ####

bootstrap_results_boot <- readRDS (file = paste(path_results, 'bootstrap_results_boot.rds', sep = '/'))
bootstrap_ci <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(bootstrap_ci) <- c('estimate', 'lower_ci', 'upper_ci')

for (i in 1:length(bootstrap_results_boot$t0)) {
   stat <- names(bootstrap_results_boot[[1]][i])
   stat

   stat_ci <- data.frame(estimate = NA,
                         lower_ci = NA,
                         upper_ci = NA)

 conf.int <- boot.ci(bootstrap_results_boot, type = 'perc', index = i)
 CI <- conf.int$percent

 stat_ci$estimate <- bootstrap_results_boot$t0[i]
 stat_ci$lower_ci <- t(CI[, 4])[[1]]
 stat_ci$upper_ci <- t(CI[, 5])[[1]]
 bootstrap_ci <- rbind(bootstrap_ci, stat_ci)
 rownames(bootstrap_ci)[i] <- stat
}
row.names(bootstrap_ci) <- result_names
saveRDS(bootstrap_ci, file = paste(path_results, 'bootstrap_ci_boot.rds', sep='/'))

## Extract CIs for overall IRR 

# start: itt. or at.
# middle: nothing (overall) or .y-m (ex: 2019-01)
# middle: .iptw, .iptw_ipcw
# end: .IRR, .IRD, .IRR_stab, .IRD_stab

rows_to_keep <- c(
 ## ITT
  
 'itt.IRR',
 'itt.IRD',
 'itt.iptw.IRR',
 'itt.iptw.IRD',
 'itt.iptw.IRR_stab',
 'itt.iptw.IRD_stab',
 'at.IRR',
 'at.IRD',
 'at.iptw.IRR',
 'at.iptw.IRD',
 'at.iptw.IRR_stab',
 'at.iptw.IRD_stab',
 
 ## AT
 # lagged weights
 'at.ipcw_str_lag.IRR',
 'at.ipcw_str_lag.IRD',
 'at.ipcw_str_lag.IRR_stab',
 'at.ipcw_str_lag.IRD_stab',
 'at.iptw.ipcw_str_lag.IRR',
 'at.iptw.ipcw_str_lag.IRD',
'at.iptw.ipcw_str_lag.IRR_stab',
 'at.iptw.ipcw_str_lag.IRD_stab',

 # non-lagged weights
 'at.ipcw_str_nonlag.IRR',
 'at.ipcw_str_nonlag.IRD',
 'at.ipcw_str_nonlag.IRR_stab',
 'at.ipcw_str_nonlag.IRD_stab',
 'at.iptw.ipcw_str_nonlag.IRR',
 'at.iptw.ipcw_str_nonlag.IRD',
 'at.iptw.ipcw_str_nonlag.IRR_stab',
 'at.iptw.ipcw_str_nonlag.IRD_stab',

 # pooled weights
 'at.ipcw_pl_nonlag.IRR',
 'at.ipcw_pl_nonlag.IRD',
 'at.ipcw_pl_nonlag.IRR_stab',
 'at.ipcw_pl_nonlag.IRD_stab',
 'at.iptw.ipcw_pl_nonlag.IRR',
 'at.iptw.ipcw_pl_nonlag.IRD',
 'at.iptw.ipcw_pl_nonlag.IRR_stab',
 'at.iptw.ipcw_pl_nonlag.IRD_stab',

# modified non-lagged weights
 'at.ipcw_pl_lag.IRR',
 'at.ipcw_pl_lag.IRD',
 'at.ipcw_pl_lag.IRR_stab',
 'at.ipcw_pl_lag.IRD_stab',
 'at.iptw.ipcw_pl_lag.IRR',
  'at.iptw.ipcw_pl_lag.IRD',
  'at.iptw.ipcw_pl_lag.IRR_stab',
  'at.iptw.ipcw_pl_lag.IRD_stab'
)

summary_results <- bootstrap_ci[row.names(bootstrap_ci) %in% rows_to_keep,]
summary_results$variable <- row.names(summary_results)
summary_results <- summary_results %>%
  select(variable, everything())

write_xlsx(summary_results, paste(path_results, 'incidence_rates_ci_fx.xlsx', sep ='/'))
saveRDS(summary_results, paste(path_results, 'incidence_rates_ci_fx.rds', sep ='/'))

#### BOOTSTRAP EXECUTION USING FUTURE ####
# 
# library(doFuture)
# library(foreach)
# 
# set.seed(1)
# R <- 2
# 
# doFuture::registerDoFuture()
# future::plan(multisession, workers = 2)
# options(future.globals.maxSize= 891289600)
# 
# system.time({
#   bootstrap_results_future <- foreach(
#     i = 1:R,
#     .combine = rbind,
#     .options.future = list(
#       packages = c(
#         "dplyr",
#         "magrittr",
#         "survival",
#         "lubridate",
#         "purrr",
#         "data.table",
#         "fastglm"
#       ),
#       seed = TRUE
#     )
#   ) %dofuture% {
#     bs(cohort)
#   }
# })
# 
# result_names <- readRDS(file = paste(path_results, 'result_names.rds', sep='/'))
# bootstrap_results_future <- as.data.frame(bootstrap_results_future)
# names(bootstrap_results_future) <- result_names
# saveRDS(bootstrap_results_future, paste(path_results, 'bootstrap_results_future.rds', sep ='/'))

#### BUILD BOOTSTRAPPED CIs - USING FUTURE OR MANUAL ####

results <- readRDS(paste(path_results, 'bootstrap_results_manual.rds', sep='/'))
results_t0 <- readRDS(paste(path_results, 'bootstrap_results_manual_t0.rds', sep='/'))
#results <- readRDS(paste(path_results, 'bootstrap_results_future.rds', sep='/'))

results %<>% filter(at.ipcw_str_nonlag.IRR_stab <= 3)

bootstrap_ci <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(bootstrap_ci) <- c('estimate', 'lower_ci', 'upper_ci')

get_stat_ci <- function(stat_name) {

  stat_ci <- data.frame(
    estimate = NA,
    lower_ci = NA,
    upper_ci = NA
  )

  stat_values <- results[stat_name]
  
  stat_ci$estimate <- results_t0[2, stat_name]

  #stat_ci$estimate <- quantile(stat_values, prob = c(0.5), na.rm = TRUE)
  stat_ci$lower_ci <- quantile(stat_values, prob = c(0.025), na.rm = TRUE)
  stat_ci$upper_ci <- quantile(stat_values, prob = c(0.975), na.rm = TRUE)

  rownames(stat_ci)[1] <- stat_name
  bootstrap_ci <<- rbind(bootstrap_ci, stat_ci)

}

stat_names <- names(results)
invisible(lapply(stat_names, get_stat_ci))

bootstrap_ci$variable <- row.names(bootstrap_ci)

bootstrap_ci <- bootstrap_ci %>%
  select(variable, everything())

saveRDS(bootstrap_ci, file = paste(path_results, 'bootstrap_ci.rds', sep='/'))

## Extract CIs for overall IRR 

# start: itt. or at. 
# middle: nothing (overall) or .y-m (ex: 2019-01)
# middle: .iptw, .iptw_ipcw
# end: .IRR, .IRD, .IRR_stab, .IRD_stab

rows_to_keep <- c(
  ## ITT
  'itt.IRR',
  'itt.IRD',
  'itt.iptw.IRR',
  'itt.iptw.IRD',
  'itt.iptw.IRR_stab',
  'itt.iptw.IRD_stab',
  'at.IRR',
  'at.IRD',
  'at.iptw.IRR',
  'at.iptw.IRD',
  'at.iptw.IRR_stab',
  'at.iptw.IRD_stab',

  ## AT
  # sratified lagged weights
  'at.ipcw_str_lag.IRR',
  'at.ipcw_str_lag.IRD',
  'at.ipcw_str_lag.IRR_stab',
  'at.ipcw_str_lag.IRD_stab',
  'at.iptw.ipcw_str_lag.IRR',
  'at.iptw.ipcw_str_lag.IRD',
  'at.iptw.ipcw_str_lag.IRR_stab',
  'at.iptw.ipcw_str_lag.IRD_stab',

  # stratified non-lagged weights
  'at.ipcw_str_nonlag.IRR',
  'at.ipcw_str_nonlag.IRD',
  'at.ipcw_str_nonlag.IRR_stab',
  'at.ipcw_str_nonlag.IRD_stab',
  'at.iptw.ipcw_str_nonlag.IRR',
  'at.iptw.ipcw_str_nonlag.IRD',
  'at.iptw.ipcw_str_nonlag.IRR_stab',
  'at.iptw.ipcw_str_nonlag.IRD_stab',


  # pooled non-lagged weights
  'at.ipcw_pl_nonlag.IRR',
  'at.ipcw_pl_nonlag.IRD',
  'at.ipcw_pl_nonlag.IRR_stab',
  'at.ipcw_pl_nonlag.IRD_stab',
  'at.iptw.ipcw_pl_nonlag.IRR',
  'at.iptw.ipcw_pl_nonlag.IRD',
  'at.iptw.ipcw_pl_nonlag.IRR_stab',
  'at.iptw.ipcw_pl_nonlag.IRD_stab',

  # pooled lagged weights
  'at.ipcw_pl_lag.IRR',
  'at.ipcw_pl_lag.IRD',
  'at.ipcw_pl_lag.IRR_stab',
  'at.ipcw_pl_lag.IRD_stab',
  'at.iptw.ipcw_pl_lag.IRR',
  'at.iptw.ipcw_pl_lag.IRD',
  'at.iptw.ipcw_pl_lag.IRR_stab',
  'at.iptw.ipcw_pl_lag.IRD_stab'
)

summary_results <- bootstrap_ci[row.names(bootstrap_ci) %in% rows_to_keep,]
summary_results$variable <- row.names(summary_results)
summary_results <- summary_results %>%
  select(variable, everything())

write_xlsx(summary_results, paste(path_results, 'incidence_rates_ci.xlsx', sep ='/'))
saveRDS(summary_results, paste(path_results, 'incidence_rates_ci.rds', sep ='/'))


