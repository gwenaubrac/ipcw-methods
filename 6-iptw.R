## ---------------------------
##
## Program: 6. IPTW weights
##
## Purpose: Calculate inverse probability of treatment weights.
##
## Author: Gwen Aubrac
##
## Date Created: 2024-15-15
##
## ---------------------------
##
## Notes: Different models can be specified and tested for IPTW to achieve balance. Here, an interaction term was added
## between baseline depression and age group.
## Hypomagnesemia, hypocalcemia, and acute renal disease were removed from the IPTW due to positivity assumption violation
## (too few counts resulting in no contrast in some bootstrap samples)
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
# hypertensive, not_hypertensive
analysis <- 'main'

#### LOAD PACKAGES ####

library(lubridate)
library(dplyr)
library(magrittr)
library(fastDummies)
library(ggplot2)
library(tidyr)
library(cobalt)
library(survival)
library(survminer)
library(splines)
library(writexl)
library(broom)

#### DEFINE PATHS ####

path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_final_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'final', sep = '/')

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_update_cov.rds', sep = '/'))
setwd(path_intermediate_res)
iptw_desc <- "iptw_desc.txt"
writeLines("IPTW description:", iptw_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = iptw_desc, append= TRUE)

#### SET UP ####

covariates <- readRDS(file = paste(path_intermediate_res, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_intermediate_res, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_intermediate_res, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_intermediate_res, 'dec_comorb.rds', sep = '/'))

base_variables <- c(covariates, base_comorb)

if (analysis == 'male' | analysis == 'female') {
  base_variables <- base_variables[!base_variables %in% c('sex')]
} else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
  base_variables <- base_variables[!base_variables %in% c('year', 'month_year')]
} else if (analysis == 'depressed' | analysis == 'not_depressed') {
  base_variables <- base_variables[!base_variables %in% c('depression_base')]
} else if (analysis == 'young' | analysis == 'old') {
  cohort$age_group <- droplevels(cohort$age_group)
}

summary(cohort[base_variables])

cat(paste('Covariates included in IPTW are:', '\n'), file = iptw_desc, append = TRUE)
cat(paste(base_variables, '\n'), file = iptw_desc, append = TRUE)

covs <- subset(cohort, select = base_variables)
saveRDS(covs, file = paste(path_intermediate_res, 'iptw_covs.rds', sep ='/')) # for plots later

#### CALCULATE IPTW WEIGHTS ####

# set reference group to largest group for multi-category variables
most_common_age_cat <- cohort %>%
  count(age_group) %>%
  slice_max(n) %>%
  pull(age_group)

cohort$age_group <- relevel(cohort$age_group, ref = as.character(most_common_age_cat))
levels(cohort$age_group)

most_common_depr_cat <- cohort %>%
  count(deprivation) %>%
  slice_max(n) %>%
  pull(deprivation)

cohort$deprivation <- relevel(cohort$deprivation, ref = as.character(most_common_depr_cat))
levels(cohort$deprivation)

most_common_ethn_cat <- cohort %>%
  count(ethnicity) %>%
  slice_max(n) %>%
  pull(ethnicity)

cohort$ethnicity <- relevel(cohort$ethnicity, ref = as.character(most_common_ethn_cat))
levels(cohort$ethnicity)

## Define IPTW model

base_model <- reformulate(base_variables, 'trt_dummy')

# add  interaction between age & anxiety for antidepressant group
if (exposure == 'antidepressant') {
  inter_model <- as.formula(paste('trt_dummy', "~", paste(c(base_variables, 'age_group*anxiety_base'), collapse = " + ")))
  model <- base_model
  
  if (analysis == 'depressed' | analysis == 'not_depressed') {
    model <- base_model
  }
}

# add interaction between age & heart failure/ischemic heart disease/MI for antihypertensive group
if (exposure == 'antihypertensive') {
  inter_model <- as.formula(paste('trt_dummy', "~", 
                                  paste(c(base_variables, 
                                          'age_group*heart_failure_base', 
                                          'age_group*ischemic_heart_disease_base', 
                                          'age_group*myocardial_infarction_base',
                                          'age_group*lvh_base',
                                          'age_group*valvular_heart_disease_base',
                                          'age_group*arrhythmia_base',
                                          'heart_failure_base*heart_failure_base',
                                          'heart_failure_base*sex',
                                          'heart_failure_base*deprivation',
                                          'heart_failure_base*ethnicity',
                                          'heart_failure_base*year',
                                          'heart_failure_base*myocardial_infarction_base'), 
                                        collapse = " + ")))
  model <- inter_model
}

model
saveRDS(model, file = paste(path_final_res, 'iptw_model.rds', sep = '/')) # for later

## Unstabilized weights

iptw_fit <- glm(model, 
                family = binomial(), 
                data = cohort)
summary(iptw_fit)

write_xlsx(tidy(iptw_fit), paste(path_final_res, 'iptw_fit.xlsx', sep ='/'))

prop_score <- if_else(cohort$trt_dummy == 0, 
                      1 - predict(iptw_fit, type = 'response'),
                      predict(iptw_fit, type = 'response'))

cohort$ps <- prop_score
cohort$iptw <- 1/prop_score
summary(cohort$iptw)
sd(cohort$iptw)

cat('Unstabilized IPTW weights: \n min:', min(cohort$iptw), '\n max:', max(cohort$iptw), '\n mean:', mean(cohort$iptw), '\n sd:', sd(cohort$iptw))
cat(paste('Unstabilized IPTW weights: \n min:', min(cohort$iptw), '\n max:', max(cohort$iptw), '\n mean:', mean(cohort$iptw), '\n sd:', sd(cohort$iptw), '\n'), file = iptw_desc, append = TRUE)

## Stabilized weights

# numerator: marginal probability of treatment
numer_fit <- glm(trt_dummy~1, family = binomial(), data = cohort)
summary(numer_fit)
pn_trt <- predict(numer_fit, type = 'response')

# denominator: conditional probability of treatment
denom_fit <- glm(model,
                 family = binomial(),
                 data = cohort)

summary(denom_fit)
pd_trt <- predict(denom_fit, type = 'response')

# calculate weights
cohort$siptw <- if_else(cohort$trt_dummy == 0, 
                        ((1-pn_trt) / (1-pd_trt)),
                        pn_trt/pd_trt)
summary(cohort$siptw)
sd(cohort$siptw)

cat('Stabilized IPTW weights: \n min:', min(cohort$siptw), '\n max:', max(cohort$siptw), '\n mean:', mean(cohort$siptw), '\n sd:', sd(cohort$siptw))
cat(paste('Stabilized IPTW weights: \n min:', min(cohort$siptw), '\n max:', max(cohort$siptw), '\n mean:', mean(cohort$siptw), '\n sd:', sd(cohort$siptw), '\n'), file = iptw_desc, append = TRUE)

rm(iptw_fit, numer_fit, denom_fit)
rm(base_model, inter_model)

saveRDS(cohort, file = paste(path_intermediate_res, 'cohort_iptw.rds', sep='/'))
