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
## Date Created: 2024-07-22
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
## ---------------------------

# analysis: flex_grace_period, 90_day_grace_period
# male, female, young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed

analysis <- ''

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(ggplot2)
library(writexl)
library(broom)
library(purrr)

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


path_comorb_cprd <- "Z:/EPI/Protocol 24_004042/Gwen/data/comorbidities/Aurum codes comorbidities" 

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))

setwd(path_results) 
ipcw_desc <- "ipcw_desc.txt"
writeLines("IPCW description:", ipcw_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = ipcw_desc, append= TRUE)

covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_main, 'dec_comorb.rds', sep = '/'))

#### PREPARE DATA ####

cohort_long <- cohort
variables <- c(covariates, comorbidities)
variables <- variables[!variables %in% c('hypocalcemia', 'hypomagnesemia', 'acute_renal_disease')] # too few counts/violate positivity assumption

if (analysis == 'male' | analysis == 'female') {
  variables <- variables[!variables %in% c('sex')]
} else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
  variables <- variables[!variables %in% c('year', 'month_year')]
} else if (analysis == 'depressed' | analysis == 'not_depressed') {
  variables <- variables[!variables %in% c('depression')]
} 

variables

base_variables <- c(covariates, base_comorb)
base_variables <- base_variables[!base_variables %in% c('hypocalcemia_d1', 'hypomagnesemia_d1', 'acute_renal_disease_d1',
                                            'hypocalcemia_d2', 'hypomagnesemia_d2', 'acute_renal_disease_d1',
                                            'hypocalcemia_d3', 'hypomagnesemia_d3', 'acute_renal_disease_d3',
                                            'hypocalcemia_d4', 'hypomagnesemia_d4', 'acute_renal_disease_d4',
                                            'hypocalcemia_d5', 'hypomagnesemia_d5', 'acute_renal_disease_d5',
                                            'hypocalcemia_d6', 'hypomagnesemia_d6', 'acute_renal_disease_d6',
                                            'hypocalcemia_d7', 'hypomagnesemia_d7', 'acute_renal_disease_d7',
                                            'hypocalcemia_d8', 'hypomagnesemia_d8', 'acute_renal_disease_d8',
                                            'hypocalcemia_d9', 'hypomagnesemia_d9', 'acute_renal_disease_d9')]

if (analysis == 'male' | analysis == 'female') {
  base_variables <- base_variables[!base_variables %in% c('sex')]
} else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
  base_variables <- base_variables[!base_variables %in% c('year')]
} else if (analysis == 'depressed' | analysis == 'not_depressed') {
  base_variables <- base_variables[!base_variables %in% c('depression_base')]
  }


base_variables

#### EXAMINE PREDICTORS OF CENSORING ####

# set reference group to largest group
cohort$age_group <- relevel(cohort$age_group, ref = "25-34")
cohort$deprivation <- relevel(cohort$deprivation, ref = "5")
cohort$ethnicity <- relevel(cohort$ethnicity, ref = "White")

predictors_switch <- glm(reformulate(base_variables, 'switch'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_switch), paste(path_results, 'predictors_switch.xlsx', sep ='/'))

predictors_disc <- glm(reformulate(base_variables, 'disc'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_disc), paste(path_results, 'predictors_disc.xlsx', sep ='/'))

predictors_cens <- glm(reformulate(base_variables, 'censor'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_cens), paste(path_results, 'predictors_cens.xlsx', sep ='/'))

# treatment-specific probability of censoring
cohort_snri <- cohort %>% filter(trt == 'snri')
predictors_cens_snri <- glm(reformulate(base_variables, 'censor'), data = cohort_snri, family = binomial())
write_xlsx(tidy(predictors_cens_snri), paste(path_results, 'predictors_cens_snri.xlsx', sep ='/'))

cohort_ssri <- cohort %>% filter(trt == 'ssri')
predictors_cens_ssri <- glm(reformulate(base_variables, 'censor'), data = cohort_ssri, family = binomial())
write_xlsx(tidy(predictors_cens_ssri), paste(path_results, 'predictors_cens_ssri.xlsx', sep ='/'))

#### CONVERT TO COUNTING TIME ####

# counting_time: time in days between cohort entry and censoring or cohort exit, whichever occurred first
# times_dec: deciles of the censoring distribution for any reason (based on distribution of at_exit_date)

cohort_long$counting_time <- as.numeric(difftime(cohort_long$at_exit_date, cohort_long$entry_date, units = 'days'))
cohort_long$censor_counting_time <- as.numeric(difftime(cohort_long$censor_date, cohort_long$entry_date, units = 'days'))
cohort_long$uncensored <- if_else(cohort_long$censor == 0, 1, 0)
cohort_long <- cohort_long[order(cohort_long$counting_time),]

times_dec <- readRDS(paste(path_cohort, 'times_dec.rds', sep = '/'))
times_dec <- times_dec[1:9] # remove last decile (100%)
times_dec

#### SPLIT DATA INTO INTERVALS OF CENSORING TIMES ####

cohort_long$Tstart <- 0
cohort_long <- survSplit(cohort_long, cut = times_dec, end = 'counting_time', start = 'Tstart', event = 'uncensored')
names(cohort_long)[names(cohort_long) == 'counting_time'] <- 'Tstop' 
cohort_long$uncensored_at_tstop <- if_else(is.na(cohort_long$censor_counting_time), 1, # indicator for being uncensored by start of next interval
                                           if_else(cohort_long$Tstop == cohort_long$censor_counting_time, 0, 1)) 

# retrieve covariate values at deciles
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

p_uncens_formula <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'age_group*depression'), collapse = " + ")))
p_uncens_formula_pooled <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'age_group*depression', 'dec'), collapse = " + ")))

if (analysis == 'depressed' | analysis == 'not_depressed') {
  p_uncens_formula <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables), collapse = " + ")))
  p_uncens_formula_pooled <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'dec'), collapse = " + ")))
}

p_uncens_formula
p_uncens_formula_pooled

#### STRATIFIED: NON-LAGGED AND LAGGED WEIGHTS ####

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
getweights <- cohort_long %>%
  group_by(dec, trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict) %>% 
  list_rbind()

getweights %<>%
  group_by(trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict_base) %>% 
  list_rbind()

# calculate stratified weights
getweights %<>%
  group_by(id) %>% 
  arrange(id, dec) %>% 
  mutate(weight = 1/p_uncens,
         ipcw_str_nonlag = cumprod(weight),
         ipcw_str_lag = if_else(row_number() == 1, 1, lag(ipcw_str_nonlag)),
         sipcw_str_nonlag = ipcw_str_nonlag * p_uncens_base,
         sipcw_str_lag = ipcw_str_lag * p_uncens_base
  )

getweights %<>% select(id, dec, ipcw_str_lag, ipcw_str_nonlag, sipcw_str_lag, sipcw_str_nonlag)
summary(getweights)

cohort_long <- merge(cohort_long, getweights, by = c('id', 'dec'), all.x = TRUE)
cohort_long %<>% arrange(id, dec)

cat(paste('Unstabilized lagged IPCW weights: \n min:', min(cohort_long$ipcw_str_lag), '\n max:', max(cohort_long$ipcw_str_lag), '\n mean:', mean(cohort_long$ipcw_str_lag), '\n sd:', sd(cohort_long$ipcw_str_lag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized lagged IPCW weights: \n min:', min(cohort_long$sipcw_str_lag), '\n max:', max(cohort_long$sipcw_str_lag), '\n mean:', mean(cohort_long$sipcw_str_lag), '\n sd:', sd(cohort_long$sipcw_str_lag), '\n'), file = ipcw_desc, append = TRUE)

cat(paste('Unstabilized non-lagged IPCW weights: \n min:', min(cohort_long$ipcw_str_nonlag), '\n max:', max(cohort_long$ipcw_str_nonlag), '\n mean:', mean(cohort_long$ipcw_str_nonlag), '\n sd:', sd(cohort_long$ipcw_str_nonlag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized non-lagged IPCW weights: \n min:', min(cohort_long$sipcw_str_nonlag), '\n max:', max(cohort_long$sipcw_str_nonlag), '\n mean:', mean(cohort_long$sipcw_str_nonlag), '\n sd:', sd(cohort_long$sipcw_str_nonlag), '\n'), file = ipcw_desc, append = TRUE)

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
getpooledweights <- cohort_long %>%
  group_by(trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict_pooled) %>% 
  list_rbind()

getpooledweights %<>%
  group_by(trt_dummy) %>%
  group_split() %>%
  map(fit_and_predict_base) %>% 
  list_rbind()

# calculate pooled weights
getpooledweights %<>%
  arrange(id, dec) %>% 
  group_by(id) %>% 
  mutate(weight = 1/p_uncens,
         ipcw_pl_nonlag = cumprod(weight),
         ipcw_pl_lag = if_else(row_number() == 1, 1, lag(ipcw_pl_nonlag)),
         sipcw_pl_nonlag = ipcw_pl_nonlag * p_uncens_base,
         sipcw_pl_lag = ipcw_pl_lag * p_uncens_base
  )

getpooledweights %<>% select(id, dec, ipcw_pl_nonlag, sipcw_pl_nonlag, weight)
summary(getpooledweights)

cohort_long <- merge(cohort_long, getpooledweights, by = c('id', 'dec'), all.x = TRUE)

cat(paste('Unstabilized pooled IPCW weights: \n min:', min(cohort_long$ipcw_pl_nonlag), '\n max:', max(cohort_long$ipcw_pl_nonlag), '\n mean:', mean(cohort_long$ipcw_pl_nonlag), '\n sd:', sd(cohort_long$ipcw_pl_nonlag), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized pooled IPCW weights: \n min:', min(cohort_long$sipcw_pl_nonlag), '\n max:', max(cohort_long$sipcw_pl_nonlag), '\n mean:', mean(cohort_long$sipcw_pl_nonlag), '\n sd:', sd(cohort_long$sipcw_pl_nonlag), '\n'), file = ipcw_desc, append = TRUE)

saveRDS(cohort_long, file = paste(path_cohort, 'antidepressant_cohort_ipcw.rds', sep='/'))

#### CHECKING WEIGHTS AND MODEL ####

# pooled IPCW model coefficients
ipcw_pl_nonlaged_fit_snri <- glm(p_uncens_formula_pooled, 
                family = binomial(), 
                data = cohort_long, 
                subset = trt_dummy == 0)

write_xlsx(tidy(ipcw_pl_nonlaged_fit_snri), paste(path_results, 'ipcw_pl_nonlaged_fit_snri.xlsx', sep ='/'))

ipcw_pl_nonlaged_fit_ssri <- glm(p_uncens_formula_pooled, 
                            family = binomial(), 
                            data = cohort_long, 
                            subset = trt_dummy == 1)

write_xlsx(tidy(ipcw_pl_nonlaged_fit_ssri), paste(path_results, 'ipcw_pl_nonlaged_fit_ssri.xlsx', sep ='/'))

# check if censored patients have higher IPCW weights in one treatment group
censored_snri <- cohort_long %>% 
  filter(censor == 1, trt_dummy == 0)

summary(censored_snri$siptw)
sd(censored_snri$siptw)
summary(censored_snri$sipcw_str_lag)
sd(censored_snri$sipcw_str_lag)
summary(censored_snri$sipcw_str_nonlag)
sd(censored_snri$sipcw_str_nonlag)
summary(censored_snri$sipcw_pl_nonlag)
sd(censored_snri$sipcw_pl_nonlag)
summary(censored_snri$sipcw_pl_lag)
sd(censored_snri$sipcw_pl_lag)

censored_ssri <- cohort_long %>% 
  filter(censor == 1, trt_dummy == 1)

summary(censored_ssri$siptw)
sd(censored_ssri$siptw)
summary(censored_ssri$sipcw_str_lag)
sd(censored_ssri$sipcw_str_lag)
summary(censored_ssri$sipcw_str_nonlag)
sd(censored_ssri$sipcw_str_nonlag)
summary(censored_ssri$sipcw_pl_nonlag)
sd(censored_ssri$sipcw_pl_nonlag)
summary(censored_ssri$sipcw_pl_lag)
sd(censored_ssri$sipcw_pl_lag)
