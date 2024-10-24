## ---------------------------
##
## Program: 5. Subgroup Analyses
##
## Purpose: Create the following subgroups of patients for subgroup analyses:
## 1. Male vs Female
## 2. Year of cohort entry (2019, 2020, 2021, 2022)
## 3. Age group (<65, >=65)
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
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive, antidiabetic
exposure <- 'antihypertensive'

# outcome: all-cause mortality, suicidal ideation
outcome <- 'all-cause mortality'

#### LOAD PACKAGES ####

library(dplyr)

#### DEFINE PATHS ####

path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, 'main', 'intermediate', sep = '/')
cohort <- readRDS(file = paste(path_intermediate_res_main, 'cohort_update_cov.rds', sep = '/'))

path_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, sep = '/')

subgroup_male <- cohort %>% 
  filter (sex == 'Male')
saveRDS(subgroup_male, file = paste(path_res, 'male', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_female <- cohort %>% 
  filter (sex == 'Female')
saveRDS(subgroup_female, file = paste(path_res, 'female', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_young <- cohort %>% 
  filter (age_at_entry < 65)
saveRDS(subgroup_young, file = paste(path_res, 'young', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_old <- cohort %>% 
  filter (age_at_entry >= 65)
saveRDS(subgroup_old, file = paste(path_res, 'old', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_2019 <- cohort %>% 
  filter (year == 2019)
saveRDS(subgroup_2019, file = paste(path_res, '2019', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_2020 <- cohort %>% 
  filter (year == 2020)
saveRDS(subgroup_2020, file = paste(path_res, '2020', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_2021 <- cohort %>% 
  filter (year == 2021)
saveRDS(subgroup_2021, file = paste(path_res, '2021', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

subgroup_2022 <- cohort %>% 
  filter (year == 2022)
saveRDS(subgroup_2022, file = paste(path_res, '2022', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

if (exposure == 'antidepressant') {
  subgroup_depressed <- cohort %>% 
    filter (depression_base == 1)
  saveRDS(subgroup_depressed, file = paste(path_res, 'depressed', 'intermediate', 'cohort_update_cov.rds', sep= '/'))
  
  subgroup_not_depressed <- cohort %>% 
    filter (depression_base == 0)
  saveRDS(subgroup_not_depressed, file = paste(path_res, 'not_depressed', 'intermediate', 'cohort_update_cov.rds', sep= '/'))

  rm(subgroup_depressed, subgroup_not_depressed)
  }

if (exposure == 'antihypertensive') {
  subgroup_hypertensive <- cohort %>% 
    filter (hypertension_base == 1)
  saveRDS(subgroup_hypertensive, file = paste(path_res, 'hypertensive', 'cohort_update_cov.rds', sep= '/'))
  
  subgroup_not_hypertensive <- cohort %>% 
    filter (hypertension_base == 0)
  saveRDS(subgroup_not_hypertensive, file = paste(path_res, 'not_hypertensive', 'cohort_update_cov.rds', sep= '/'))

  rm(subgroup_hypertensive, subgroup_not_hypertensive)
  }

rm(subgroup_male, subgroup_female, subgroup_young, subgroup_old, subgroup_2019, subgroup_2020, subgroup_2021, subgroup_2022)

