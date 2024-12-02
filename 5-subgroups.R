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

# cohort: antidepressant, antihypertensive
exposure <- 'antihypertensive'
outcome <- 'all-cause mortality'

#### LOAD PACKAGES ####

library(dplyr)

#### DEFINE PATHS ####

path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, 'main', 'intermediate', sep = '/')
cohort <- readRDS(file = paste(path_intermediate_res_main, 'cohort_update_cov.rds', sep = '/'))

path_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, sep = '/')

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

rm(subgroup_male, subgroup_female, subgroup_young, subgroup_old)

