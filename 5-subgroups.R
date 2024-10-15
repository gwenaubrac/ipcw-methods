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
## Date Created: 2024-07-09
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

#### LOAD PACKAGES ####

library(dplyr)

#### DEFINE PATHS ####

path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort"
cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_outcome.rds', sep = '/'))

subgroup_male <- cohort %>% 
  filter (sex == 'Male')
saveRDS(subgroup_male, file = paste(path_cohort, 'subgroup', 'male', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_female <- cohort %>% 
  filter (sex == 'Female')
saveRDS(subgroup_female, file = paste(path_cohort, 'subgroup', 'female', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_young <- cohort %>% 
  filter (age_at_entry < 65)
saveRDS(subgroup_young, file = paste(path_cohort, 'subgroup', 'young', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_old <- cohort %>% 
  filter (age_at_entry >= 65)
saveRDS(subgroup_old, file = paste(path_cohort, 'subgroup', 'old', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_2019 <- cohort %>% 
  filter (year == 2019)
saveRDS(subgroup_2019, file = paste(path_cohort, 'subgroup', '2019', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_2020 <- cohort %>% 
  filter (year == 2020)
saveRDS(subgroup_2020, file = paste(path_cohort, 'subgroup', '2020', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_2021 <- cohort %>% 
  filter (year == 2021)
saveRDS(subgroup_2021, file = paste(path_cohort, 'subgroup', '2021', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_2022 <- cohort %>% 
  filter (year == 2022)
saveRDS(subgroup_2022, file = paste(path_cohort, 'subgroup', '2022', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_depressed <- cohort %>% 
  filter (depression_base == 1)
saveRDS(subgroup_depressed, file = paste(path_cohort, 'subgroup', 'depressed', 'antidepressant_cohort_covariates.rds', sep= '/'))

subgroup_not_depressed <- cohort %>% 
  filter (depression_base == 0)
saveRDS(subgroup_not_depressed, file = paste(path_cohort, 'subgroup', 'not_depressed', 'antidepressant_cohort_covariates.rds', sep= '/'))

rm(subgroup_male, subgroup_female, subgroup_young, subgroup_old, subgroup_2019, subgroup_2020, subgroup_2021, subgroup_2022)
rm(subgroup_depressed, subgroup_not_depressed)
