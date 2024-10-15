## ---------------------------
##
## Program: 6. IPTW weights
##
## Purpose: Calculate inverse probability of treatment weights.
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes: Different models can be specified and tested for IPTW to achieve balance. Here, an interaction term was added
## between baseline depression and age group.
## Hypomagnesemia, hypocalcemia, and acute renal disease were removed from the IPTW due to positivity assumption violation
## (too few counts resulting in no contrast in some bootstrap samples)
##
## ---------------------------

# analysis: flex_grace_period, 90_day_grace_period
# male, female, young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed

analysis <- ''

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

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))
setwd(path_results)
iptw_desc <- "iptw_desc.txt"
writeLines("IPTW description:", iptw_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = iptw_desc, append= TRUE)

#### SET UP ####

covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_main, 'dec_comorb.rds', sep = '/'))

comorbidities <- comorbidities[!comorbidities %in% c('hypocalcemia', 'hypomagnesemia', 'acute_renal_disease')]
base_comorb <- base_comorb[!base_comorb %in% c('hypocalcemia_base', 'hypomagnesemia_base', 'acute_renal_disease_base')]
dec_comorb <- dec_comorb[!dec_comorb %in% c('hypocalcemia_d1', 'hypomagnesemia_d1', 'acute_renal_disease_d1',
                  'hypocalcemia_d2', 'hypomagnesemia_d2', 'acute_renal_disease_d1',
                  'hypocalcemia_d3', 'hypomagnesemia_d3', 'acute_renal_disease_d3',
                  'hypocalcemia_d4', 'hypomagnesemia_d4', 'acute_renal_disease_d4',
                  'hypocalcemia_d5', 'hypomagnesemia_d5', 'acute_renal_disease_d5',
                  'hypocalcemia_d6', 'hypomagnesemia_d6', 'acute_renal_disease_d6',
                  'hypocalcemia_d7', 'hypomagnesemia_d7', 'acute_renal_disease_d7',
                  'hypocalcemia_d8', 'hypomagnesemia_d8', 'acute_renal_disease_d8',
                  'hypocalcemia_d9', 'hypomagnesemia_d9', 'acute_renal_disease_d9')]

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

cat(paste('Covariates included are:', base_variables, '\n'), file = iptw_desc, append = TRUE)

covs <- subset(cohort, select = base_variables)
saveRDS(covs, file = paste(path_cohort, 'iptw_covs.rds', sep ='/')) # for plots later

#### CALCULATE IPTW WEIGHTS ####

# set reference group to largest group
cohort$age_group <- relevel(cohort$age_group, ref = "25-34")
cohort$deprivation <- relevel(cohort$deprivation, ref = "5")
cohort$ethnicity <- relevel(cohort$ethnicity, ref = "White")

## Define IPTW model

base_model <- reformulate(base_variables, 'trt_dummy')
inter_model <- as.formula(paste('trt_dummy', "~", paste(c(base_variables, 'age_group*depression_base'), collapse = " + ")))

model <- inter_model

if (analysis == 'depressed' | analysis == 'not_depressed') {
  model <- base_model
}

saveRDS(model, file = paste(path_results, 'iptw_model.rds', sep = '/')) # for later
model

## Unstabilized weights

iptw_fit <- glm(model, 
                family = binomial(), 
                data = cohort)
summary(iptw_fit)

write_xlsx(tidy(iptw_fit), paste(path_results, 'iptw_fit.xlsx', sep ='/'))

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

saveRDS(cohort, file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep='/'))
