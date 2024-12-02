## ---------------------------
##
## Program: 7. IPCW
##
## Purpose: Calculate inverse probability of censoring weights at each decile of the censoring
## distribution, with updated values for the covariates. Four different methods are compared:
## 1) Stratified nonlagged weights
## 2) Stratified lagged weights
## 3) Pooled nonlagged weights
## 4) Pooled lagged weights
##
##
## Author: Gwen Aubrac
##
## Date Created: 2024-10-15
##
## ---------------------------
##
## Notes:
##
## 1. Censoring models are fitted separately for each exposure group. 
## 2. A long format of the cohort dataframe ('cohort_long') is created and split by the deciles of censoring,
## with an indicator for whether or not patient was censored at the Tstop for that decile
## and an indicator for whether or not the patient experienced the event at the Tstop for that decile. 
## 3. There were too few counts for some of the covariates (hypocalcemia, hypomagnesemia, acute renal disease)
## when splitting by quartile, which violates the positivity assumption for inverse weighting.
## These covariates were removed from the IPTW and IPCW models. 
## 4. Patients who died during the follow-up interval are remove from the dataset used to fit and estimate
## the IPCW weights and receive a weight of 1 for that interval.
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive
exposure <- 'antihypertensive'
outcome <- 'all-cause mortality'

# analysis: main, flexible_grace_period, 90_day_grace_period, male, female, young, old
analysis <- 'main'

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(ggplot2)
library(writexl)
library(broom)
library(purrr)

#### DEFINE PATHS ####

path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_final_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, analysis, 'final', sep = '/')
path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, 'main', 'intermediate', sep = '/')

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))

# setwd(path_intermediate_res) 
# ipcw_desc <- "ipcw_desc.txt"
# writeLines("IPCW description:", ipcw_desc)
# cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = ipcw_desc, append= TRUE)

covariates <- readRDS(file = paste(path_intermediate_res_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_intermediate_res_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_intermediate_res_main, 'base_comorb.rds', sep = '/'))

#### PREPARE DATA ####

# remove cardiomyopathy from comorbidities because too few counts in intervals
comorbidities <- comorbidities[!comorbidities %in% c('cardiomyopathy')]
base_comorb <- base_comorb[!base_comorb %in% c('cardiomyopathy_base')]

cohort_long <- cohort
variables <- c(covariates, comorbidities)

if (analysis == 'male' | analysis == 'female') {
  variables <- variables[!variables %in% c('sex')]
}  
variables

base_variables <- c(covariates, base_comorb)

if (analysis == 'male' | analysis == 'female') {
  base_variables <- base_variables[!base_variables %in% c('sex')]
} 
base_variables

#### EXAMINE PREDICTORS OF CENSORING ####

# set reference group to largest group
predictors_switch <- glm(reformulate(base_variables, 'switch'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_switch), paste(path_final_res, 'predictors_switch.xlsx', sep ='/'))

predictors_disc <- glm(reformulate(base_variables, 'disc'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_disc), paste(path_final_res, 'predictors_disc.xlsx', sep ='/'))

predictors_cens <- glm(reformulate(base_variables, 'censor'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_cens), paste(path_final_res, 'predictors_cens.xlsx', sep ='/'))

rm(predictors_switch, predictors_disc, predictors_cens)

# treatment-specific probability of censoring
cohort_ref <- cohort %>% filter(trt_dummy == 0)
predictors_cens_ref <- glm(reformulate(base_variables, 'censor'), data = cohort_ref, family = binomial())
write_xlsx(tidy(predictors_cens_ref), paste(path_final_res, 'predictors_cens_ref.xlsx', sep ='/'))

cohort_comp <- cohort %>% filter(trt_dummy == 1)
predictors_cens_comp <- glm(reformulate(base_variables, 'censor'), data = cohort_comp, family = binomial())
write_xlsx(tidy(predictors_cens_comp), paste(path_final_res, 'predictors_cens_comp.xlsx', sep ='/'))

rm(cohort_ref, cohort_comp, predictors_cens_ref, predictors_cens_comp)

#### CONVERT TO COUNTING TIME ####

# counting_time: time in days between cohort entry and censoring or cohort exit, whichever occurred first
# times_dec: intervals of the censoring distribution for any reason (based on distribution of censoring times due to disc or switch)

cohort_long$counting_time <- as.numeric(difftime(cohort_long$at_exit_date, cohort_long$entry_date, units = 'days'))
cohort_long$censor_counting_time <- as.numeric(difftime(cohort_long$censor_date, cohort_long$entry_date, units = 'days'))
cohort_long$uncensored <- if_else(cohort_long$censor == 0, 1, 0)
cohort_long <- cohort_long[order(cohort_long$counting_time),]

times_dec <- readRDS(paste(path_intermediate_res_main, 'times_dec.rds', sep = '/'))
times_dec <- times_dec[1:9] # remove last decile (100%)
times_dec

#### SPLIT DATA INTO INTERVALS OF CENSORING TIMES ####

cohort_long$Tstart <- 0
cohort_long <- survSplit(cohort_long, cut = times_dec, end = 'counting_time', start = 'Tstart', event = 'uncensored')
names(cohort_long)[names(cohort_long) == 'counting_time'] <- 'Tstop' 
cohort_long$uncensored_at_tstop <- if_else(is.na(cohort_long$censor_counting_time), 1, # indicator for being uncensored by start of next interval
                                           if_else(cohort_long$Tstop == cohort_long$censor_counting_time, 0, 1)) 

# retrieve covariate values at intervals
for (i in 1:length(comorbidities)) {
  comorb_name <- comorbidities[i]
  cohort_long$comorb <- if_else(
    cohort_long$Tstart == 0,
    cohort_long[[paste(comorb_name, 'base', sep = '_')]],
    if_else(
      cohort_long$Tstart == times_dec[[1]],
      cohort_long[[paste(comorb_name, 'd1', sep = '_')]],
      if_else(
        cohort_long$Tstart == times_dec[[2]],
        cohort_long[[paste(comorb_name, 'd2', sep = '_')]],
        if_else(
          cohort_long$Tstart == times_dec[[3]],
          cohort_long[[paste(comorb_name, 'd3', sep = '_')]],
          if_else(
            cohort_long$Tstart == times_dec[[4]],
            cohort_long[[paste(comorb_name, 'd4', sep = '_')]],
            if_else(
              cohort_long$Tstart == times_dec[[5]],
              cohort_long[[paste(comorb_name, 'd5', sep = '_')]],
              if_else(
                cohort_long$Tstart == times_dec[[6]],
                cohort_long[[paste(comorb_name, 'd6', sep = '_')]],
                if_else(
                  cohort_long$Tstart == times_dec[[7]],
                  cohort_long[[paste(comorb_name, 'd7', sep = '_')]],
                  if_else(
                    cohort_long$Tstart == times_dec[[8]],
                    cohort_long[[paste(comorb_name, 'd8', sep = '_')]], cohort_long[[paste(comorb_name, 'd9', sep = '_')]])
                )
              )
            )
          )
        )
      )
    )
  )
  names(cohort_long)[names(cohort_long) == 'comorb'] <- comorb_name
}

cohort_long %<>% 
  group_by(id) %>% 
  mutate(dec = as.factor(row_number())) %>% 
  arrange(id, dec) %>% 
  ungroup()

# p_uncens_formula <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables), collapse = " + ")))
# p_uncens_formula_pooled <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'dec'), collapse = " + ")))

# add interaction between anxiety and age for antidepressant group
if (exposure == 'antidepressant') {
  p_uncens_formula <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'age_group*anxiety'), collapse = " + ")))
  p_uncens_formula_pooled <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'age_group*anxiety', 'dec'), collapse = " + ")))
}

# add interactions for heart failure/ischemic heart disease/MI for antihypertensive group
if (exposure == 'antihypertensive') {
  p_uncens_formula <- as.formula(paste('uncensored_at_tstop', "~", 
                                       paste(c(variables, 
                                               'age_group*heart_failure',
                                               'age_group*ischemic_heart_disease',
                                               'age_group*myocardial_infarction',
                                               'age_group*lvh',
                                               'age_group*valvular_heart_disease',
                                               'age_group*arrhythmia',
                                               'heart_failure*heart_failure',
                                               'heart_failure*sex',
                                               'heart_failure*year',
                                               'heart_failure*deprivation',
                                               'heart_failure*ethnicity',
                                               'heart_failure*myocardial_infarction'), 
                                             collapse = " + ")))
  
  p_uncens_formula_pooled <- as.formula(paste('uncensored_at_tstop', "~", 
                                              paste(c(variables, 
                                                      'age_group*heart_failure',
                                                      'age_group*ischemic_heart_disease',
                                                      'age_group*myocardial_infarction',
                                                      'age_group*lvh',
                                                      'age_group*valvular_heart_disease',
                                                      'age_group*arrhythmia',
                                                      'heart_failure*heart_failure',
                                                      'heart_failure*sex',
                                                      'heart_failure*year',
                                                      'heart_failure*deprivation',
                                                      'heart_failure*ethnicity',
                                                      'heart_failure*myocardial_infarction',
                                                      'dec'), 
                                                    collapse = " + ")))
}

if (exposure == 'antihypertensive' && (analysis == 'male' | analysis == 'female')) {
  p_uncens_formula <- as.formula(paste('uncensored_at_tstop', "~", 
                                       paste(c(variables, 
                                               'age_group*heart_failure',
                                               'age_group*ischemic_heart_disease',
                                               'age_group*myocardial_infarction',
                                               'age_group*lvh',
                                               'age_group*valvular_heart_disease',
                                               'age_group*arrhythmia',
                                               'heart_failure*heart_failure',
                                               'heart_failure*year',
                                               'heart_failure*deprivation',
                                               'heart_failure*ethnicity',
                                               'heart_failure*myocardial_infarction'), 
                                             collapse = " + ")))
  
  p_uncens_formula_pooled <- as.formula(paste('uncensored_at_tstop', "~", 
                                              paste(c(variables, 
                                                      'age_group*heart_failure',
                                                      'age_group*ischemic_heart_disease',
                                                      'age_group*myocardial_infarction',
                                                      'age_group*lvh',
                                                      'age_group*valvular_heart_disease',
                                                      'age_group*arrhythmia',
                                                      'heart_failure*heart_failure',
                                                      'heart_failure*year',
                                                      'heart_failure*deprivation',
                                                      'heart_failure*ethnicity',
                                                      'heart_failure*myocardial_infarction',
                                                      'dec'), 
                                                    collapse = " + ")))
}

p_uncens_formula
p_uncens_formula_pooled

#### STRATIFIED: NON-LAGGED AND LAGGED WEIGHTS ####

# flag and remove patients who died during follow-up interval
cohort_long_nodeaths <- cohort_long %>% 
  mutate(death_fu = as.numeric(dod - entry_date)) %>% 
  mutate(death_in_interval = if_else (!is.na(dod) & death_fu >= Tstart & death_fu <= Tstop, 1, 0)) 

cohort_long_nodeaths <- cohort_long_nodeaths %>% 
  filter(death_in_interval == 0) 

# function to get stratified parametric probability of uncensored by next interval
fit_and_predict <- function(data_subset) {
  model <- glm(p_uncens_formula,
               family = binomial(link = "logit"),
               data = data_subset)
  
  data_subset$p_uncens <- predict(model, type = "response")
  return(data_subset)
}

# function to get marginal probability of uncensored by next interval
# (for stabilization)
fit_and_predict_base <- function(data_subset) {
  model <- glm(uncensored_at_tstop ~ 1,
               family = binomial(link = "logit"),
               data = data_subset)
  
  data_subset$p_uncens_base <- predict(model, type = "response")
  return(data_subset)
}

# apply functions to cohort dataframe
getweights <- cohort_long_nodeaths %>%
  group_by(dec, trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict) %>% 
  list_rbind()

getweights %<>%
  group_by(dec, trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict_base) %>% 
  list_rbind()

# calculate stratified weights
getweights %<>%
  group_by(id) %>% 
  arrange(id, dec) %>% 
  mutate(weight = 1/p_uncens)

getweights %<>% select(id, dec, weight, p_uncens_base)
cohort_long <- merge(cohort_long, getweights, by = c('id', 'dec'), all.x = TRUE)

# replace NAs (patiens who died in interval) with 1
cohort_long %<>% mutate(
  weight = if_else(is.na(weight), 1, weight)
)

cohort_long %<>%
  group_by(id) %>% 
  arrange(id, dec) %>% 
  mutate(
    ipcw_str_nonlag = cumprod(weight),
    ipcw_str_lag = if_else(row_number() == 1, 1, lag(ipcw_str_nonlag))
  )

cat(paste('Unstabilized stratified non-lagged IPCW weights: \n min:', min(cohort_long$ipcw_str_nonlag), '\n max:', max(cohort_long$ipcw_str_nonlag), '\n mean:', mean(cohort_long$ipcw_str_nonlag), '\n sd:', sd(cohort_long$ipcw_str_nonlag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized stratified non-lagged IPCW weights: \n min:', min(cohort_long$sipcw_str_nonlag), '\n max:', max(cohort_long$sipcw_str_nonlag), '\n mean:', mean(cohort_long$sipcw_str_nonlag), '\n sd:', sd(cohort_long$sipcw_str_nonlag), '\n'), file = ipcw_desc, append = TRUE)

cat(paste('Unstabilized stratified lagged IPCW weights: \n min:', min(cohort_long$ipcw_str_lag), '\n max:', max(cohort_long$ipcw_str_lag), '\n mean:', mean(cohort_long$ipcw_str_lag), '\n sd:', sd(cohort_long$ipcw_str_lag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized stratified lagged IPCW weights: \n min:', min(cohort_long$sipcw_str_lag), '\n max:', max(cohort_long$sipcw_str_lag), '\n mean:', mean(cohort_long$sipcw_str_lag), '\n sd:', sd(cohort_long$sipcw_str_lag), '\n'), file = ipcw_desc, append = TRUE)

#### POOLED: NON-LAGGED AND LAGGED WEIGHTS ####

# function to get pooled parametric probability of uncensored by next interval
fit_and_predict_pooled <- function(data_subset) {
  model <- glm(p_uncens_formula_pooled,
               family = binomial(link = "logit"),
               data = data_subset)
  
  data_subset$p_uncens <- predict(model, type = "response")
  return(data_subset)
}

# apply functions to cohort dataframe
getpooledweights <- cohort_long_nodeaths %>%
  group_by(trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict_pooled) %>% 
  list_rbind()

getpooledweights %<>%
  group_by(trt_dummy, dec) %>%
  group_split() %>%
  map(fit_and_predict_base) %>% 
  list_rbind()

# calculate pooled weights
getpooledweights %<>%
  arrange(id, dec) %>% 
  group_by(id) %>% 
  mutate(weight_pl = 1/p_uncens
  )

getpooledweights %<>% select(id, dec, weight_pl)
cohort_long <- merge(cohort_long, getpooledweights, by = c('id', 'dec'), all.x = TRUE)

# replace NAs with 1
cohort_long %<>% mutate(
  weight_pl = if_else(is.na(weight_pl), 1, weight_pl))

cohort_long %<>%
  arrange(id, dec) %>% 
  group_by(id) %>% 
  mutate(
    ipcw_pl_nonlag = cumprod(weight_pl),
    ipcw_pl_lag = if_else(row_number() == 1, 1, lag(ipcw_pl_nonlag))
  )

cat(paste('Unstabilized pooled nonlagged IPCW weights: \n min:', min(cohort_long$ipcw_pl_nonlag), '\n max:', max(cohort_long$ipcw_pl_nonlag), '\n mean:', mean(cohort_long$ipcw_pl_nonlag), '\n sd:', sd(cohort_long$ipcw_pl_nonlag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized pooled nonlagged IPCW weights: \n min:', min(cohort_long$sipcw_pl_nonlag), '\n max:', max(cohort_long$sipcw_pl_nonlag), '\n mean:', mean(cohort_long$sipcw_pl_nonlag), '\n sd:', sd(cohort_long$sipcw_pl_nonlag), '\n'), file = ipcw_desc, append = TRUE)

cat(paste('Unstabilized pooled lagged IPCW weights: \n min:', min(cohort_long$ipcw_pl_lag), '\n max:', max(cohort_long$ipcw_pl_lag), '\n mean:', mean(cohort_long$ipcw_pl_lag), '\n sd:', sd(cohort_long$ipcw_pl_lag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized pooled lagged IPCW weights: \n min:', min(cohort_long$sipcw_pl_lag), '\n max:', max(cohort_long$sipcw_pl_lag), '\n mean:', mean(cohort_long$sipcw_pl_lag), '\n sd:', sd(cohort_long$sipcw_pl_lag), '\n'), file = ipcw_desc, append = TRUE)

saveRDS(cohort_long, file = paste(path_intermediate_res, 'cohort_ipcw.rds', sep='/'))

#### CHECKING WEIGHTS AND MODEL ####

# pooled IPCW model coefficients
ipcw_pl_nonlagged_fit_ref <- glm(p_uncens_formula_pooled, 
                family = binomial(), 
                data = cohort_long, 
                subset = trt_dummy == 0)

write_xlsx(tidy(ipcw_pl_nonlagged_fit_ref), paste(path_final_res, 'ipcw_pl_nonlagged_fit_ref.xlsx', sep ='/'))

ipcw_pl_nonlagged_fit_comp <- glm(p_uncens_formula_pooled, 
                            family = binomial(), 
                            data = cohort_long, 
                            subset = trt_dummy == 1)

write_xlsx(tidy(ipcw_pl_nonlagged_fit_comp), paste(path_final_res, 'ipcw_pl_nonlagged_fit_comp.xlsx', sep ='/'))

# check if censored patients have higher IPCW weights in one treatment group
censored_ref <- cohort_long %>% 
  filter(censor == 1, trt_dummy == 0)

summary(censored_ref$iptw)
sd(censored_ref$iptw)
summary(censored_ref$ipcw_str_lag)
sd(censored_ref$ipcw_str_lag)
summary(censored_ref$ipcw_str_nonlag)
sd(censored_ref$ipcw_str_nonlag)
summary(censored_ref$ipcw_pl_nonlag)
sd(censored_ref$ipcw_pl_nonlag)
summary(censored_ref$ipcw_pl_lag)
sd(censored_ref$ipcw_pl_lag)

censored_comp <- cohort_long %>% 
  filter(censor == 1, trt_dummy == 1)

summary(censored_comp$iptw)
sd(censored_comp$iptw)
summary(censored_comp$ipcw_str_lag)
sd(censored_comp$ipcw_str_lag)
summary(censored_comp$ipcw_str_nonlag)
sd(censored_comp$ipcw_str_nonlag)
summary(censored_comp$ipcw_pl_nonlag)
sd(censored_comp$ipcw_pl_nonlag)
summary(censored_comp$ipcw_pl_lag)
sd(censored_comp$ipcw_pl_lag)
