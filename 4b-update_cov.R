## ---------------------------
##
## Program: 4c. Define censoring and update covariates
##
## Purpose: Define censoring and follow-up intervals and 
## update occurrence of baseline covariates at start of each interval. 
## 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-10-22
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive
exposure <- 'antihypertensive'
outcome <- 'all-cause mortality'

# analysis: main, flexible_grace_period, 90_day_grace_period
analysis <- 'flexible_grace_period'

#### LOAD PACKAGES ####

library(lubridate)
library(readxl)
library(dplyr)
library(magrittr)
library(haven)
library(parallel)
library(data.table)
library(tidyr)

#### DEFINE PATHS ####

path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_linkage_1 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt1'
path_linkage_2 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt2'

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_outcome.rds', sep = '/'))
switched_to <- readRDS(file = paste('Z:/EPI/Protocol 24_004042/Gwen - ipcw_methods/results', exposure, 'all-cause mortality', 'main', 'intermediate', 'switched_to.rds', sep = '/'))

setwd(path_intermediate_res)

update_cov_desc <- "update_cov_desc.txt"
writeLines("Outcome description:", update_cov_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = update_cov_desc, append= TRUE)

if (outcome == 'suicidal ideation') {
  outcome <- 'suicidal_ideation_self_harm'
}

#### DEFINE OVERALL CENSOR DATE FROM TRT DISCONTINUATION OR SWITCH ####

cohort %<>%
  mutate (censor = if_else (switch == 1 | disc == 1, 1, 0),
          censor_date = if_else(censor == 1, pmin(switch_date, disc_date, na.rm = TRUE), NA))

table(cohort$censor)
length(which(is.na(cohort$censor_date)))

cat('Number of patients censored due to trt switch or discontinuation:', sum(cohort$censor))
cat(paste('Number of patients censored due to trt switch or discontinuation:', sum(cohort$censor), '\n'), file = update_cov_desc, append = TRUE)

cat('Number of patients censored due to trt switch:', sum(cohort$switch))
cat(paste('Number of patients censored due to trt switch:', sum(cohort$switch), '\n'), file = update_cov_desc, append = TRUE)

cat('Number of patients censored due to trt discontinuation:', sum(cohort$disc))
cat(paste('Number of patients censored due to trt discontinuation:', sum(cohort$disc), '\n'), file = update_cov_desc, append = TRUE)

cat('Number of patients not censored:', sum(cohort$censor == 0))
cat(paste('Number of patients not censored:', sum(cohort$censor == 0), '\n'), file = update_cov_desc, append = TRUE)

# describe trt switches
switched_to %<>% select(id, trt_seq)
cohort <- merge(cohort, switched_to, by='id', all.x = TRUE)

switch_tbl <- table(cohort[cohort$switch ==1, 'trt_seq'])
cat(paste('Treatment switches:', '\n'), file = update_cov_desc, append = TRUE)
write.table(switch_tbl, file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

cat(paste('Number of ITT events:', '\n'), file = update_cov_desc, append = TRUE)
write.table(table(cohort$trt, cohort$itt_event), file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

cat(paste('Number of AT events:', '\n'), file = update_cov_desc, append = TRUE)
write.table(table(cohort$trt, cohort$at_event), file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)


#### UPDATE PATIENT COVARIATE VALUES AT INTERVALS ####

## Set up

# define deciles of censoring distribution in follow-up counting time
# based on distribution of censoring times among censored patients
censored_only <- cohort %>%
  filter(!is.na(censor_date))

censoring_times <- as.numeric(censored_only$censor_date - censored_only$entry_date)
quantile(censoring_times, probs = seq(0, 1, by = 0.1))
cat(paste('The original distribution of censoring times is:', '\n'), file = update_cov_desc, append = TRUE)
cat(paste(quantile(censoring_times, probs = seq(0, 1, by = 0.1))), file = update_cov_desc, append = TRUE)

## for ANTIDEPRESSANT cohort:
# the start of the 1st, 2nd, and 3rd interval are the same at day 58
# so let us split data into deciles starting from the 3rd interval (~35%)
# to be more flexible

if (exposure == 'antidepressant') {
  breaks <- round(seq(
    from = 0.35,
    to = 1,
    length.out = 10
  ), 2)
  breaks
  
  times_dec <- quantile(censoring_times, probs = breaks)
  times_dec
  
  saveRDS(breaks, file = paste(path_intermediate_res, 'interval_breaks.rds', sep = '/'))
  
  length(which(censoring_times < times_dec[[1]]))
  length(which(censoring_times >= times_dec[[1]] &
                 censoring_times < times_dec[[2]]))
  length(which(censoring_times == times_dec[[1]]))
  
  # let us shift the concentration of censoring to the first decile
  # by setting the second decile at +1 days
  times_dec[[1]] <- times_dec[[1]] + 1
  times_dec
  
  saveRDS(times_dec, file = paste(path_intermediate_res, 'times_dec.RDS', sep = '/'))
  
  cat(paste('The new distribution of censoring times is:', '\n'), file = update_cov_desc, append = TRUE)
  cat(paste(times_dec), file = update_cov_desc, append = TRUE)
  
}

### for ANTIHYPERTENSIVE cohort
# the start of the 2nd and 3rd interval are the same at day 58
# so let us split the data in 11 even intervals
# and remove the duplicate interval starts
# (end up with 10 intervals)
# to be more flexible

if (exposure == 'antihypertensive') {
  
  breaks <- round(seq(
    from = 0.1,
    to = 1,
    length.out = 11
  ), 2)
  breaks
  
  times_dec <- quantile(censoring_times, probs = breaks)
  times_dec
  
  breaks <- breaks[-3]
  breaks
  times_dec <- quantile(censoring_times, probs = breaks)
  times_dec
  
  saveRDS(breaks, file = paste(path_intermediate_res, 'interval_breaks.rds', sep = '/'))
  
  length(which(censoring_times < times_dec[[1]])) # patients censored in 1st interval
  length(which(censoring_times >= times_dec[[1]] & # patients censored in 2nd interval
                 censoring_times < times_dec[[2]]))
  length(which(censoring_times == times_dec[[1]])) # patients censored on start of 2nd interval
  
  saveRDS(times_dec, file = paste(path_intermediate_res, 'times_dec.RDS', sep = '/'))
  
  cat(paste('The new distribution of censoring times is:', '\n'), file = update_cov_desc, append = TRUE)
  cat(paste(quantile(censoring_times, probs = breaks)), file = update_cov_desc, append = TRUE)
  
}

# get relative date of censoring interval for each patient
cohort <- cohort %>% 
  mutate (cens_d1 = entry_date + times_dec[[1]],
          cens_d2 = entry_date + times_dec[[2]],
          cens_d3 = entry_date + times_dec[[3]], 
          cens_d4 = entry_date + times_dec[[4]], 
          cens_d5 = entry_date + times_dec[[5]], 
          cens_d6 = entry_date + times_dec[[6]], 
          cens_d7 = entry_date + times_dec[[7]], 
          cens_d8 = entry_date + times_dec[[8]], 
          cens_d9 = entry_date + times_dec[[9]]
  )

# create a column for each comorbidity with the date of first occurrence in each cell;
# for each comorb, define whether it occurred at baseline and quartiles of censoring
# and create corresponding columns in cohort dataframe

first_comorb <- readRDS(paste(path_intermediate_res, 'first_comorb.rds', sep = '/'))

first_comorb <- first_comorb %>%
  pivot_wider(id_cols = id, names_from = comorb, values_from = comorb_date)

for (i in 2:ncol(first_comorb)) {
  comorb_name <- names(first_comorb[i])
  first_comorb_i <- first_comorb[, c('id', comorb_name)]
  cohort <- merge(cohort, first_comorb_i, by = 'id', all.x = TRUE)
  names(cohort)[names(cohort) == comorb_name] <- 'comorb_date'
  cohort %<>%
    mutate (comorb_base = as.factor(if_else (!is.na(comorb_date) & comorb_date<entry_date, 1, 0)),
            comorb_d1 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d1, 1, 0)),
            comorb_d2 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d2, 1, 0)),
            comorb_d3 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d3, 1, 0)),
            comorb_d4 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d4, 1, 0)),
            comorb_d5 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d5, 1, 0)),
            comorb_d6 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d6, 1, 0)),
            comorb_d7 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d7, 1, 0)),
            comorb_d8 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d8, 1, 0)),
            comorb_d9 = as.factor(if_else (!is.na(comorb_date) & comorb_date<cens_d9, 1, 0))
    )
  names(cohort)[names(cohort) == 'comorb_date'] <- paste(comorb_name, 'date', sep = '_')
  names(cohort)[names(cohort) == 'comorb_base'] <- paste(comorb_name, 'base', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d1'] <- paste(comorb_name, 'd1', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d2'] <- paste(comorb_name, 'd2', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d3'] <- paste(comorb_name, 'd3', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d4'] <- paste(comorb_name, 'd4', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d5'] <- paste(comorb_name, 'd5', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d6'] <- paste(comorb_name, 'd6', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d7'] <- paste(comorb_name, 'd7', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d8'] <- paste(comorb_name, 'd8', sep = '_')
  names(cohort)[names(cohort) == 'comorb_d9'] <- paste(comorb_name, 'd9', sep = '_')
  
}

# quickly convert all binary/categorical variables to factors
base_comorb <- readRDS(paste(path_intermediate_res, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(paste(path_intermediate_res, 'dec_comorb.rds', sep = '/'))
cat_variables <- c('sex', 'year', 'ethnicity', 'deprivation', base_comorb, dec_comorb, 'disc', 'switch')

cohort <- cohort %>% 
  mutate(across(cat_variables, as.factor)) 

rm(cat_variables)

# check depression at baseline: 
table(cohort$depression_base, cohort$trt)
table(cohort$depression_d1, cohort$trt)

cat(paste('Depression at baseline:', '\n'), file = update_cov_desc, append = TRUE)
write.table(table(cohort$depression_base, cohort$trt), file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

cat(paste('Depression at d1:', '\n'), file = update_cov_desc, append = TRUE)
write.table(table(cohort$depression_d1, cohort$trt), file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

# check anxiety at baseline
table(cohort$anxiety_base, cohort$trt)
table(cohort$anxiety_d1, cohort$trt)

# check hypertension at baseline
table(cohort$hypertension_base, cohort$trt)
table(cohort$hypertension_d1, cohort$trt)

cat(paste('Hypertension at baseline:', '\n'), file = update_cov_desc, append = TRUE)
write.table(table(cohort$hypertension_base, cohort$trt), file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

cat(paste('Hypertension at d1:', '\n'), file = update_cov_desc, append = TRUE)
write.table(table(cohort$hypertension_d1, cohort$trt), file = "update_cov_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

# because conditions tend to be undercoded in EHR data
# let us use the condition at the first interval
# as baseline, for each respective exposure of interest

if (exposure == 'antidepressant' & outcome == 'all-cause mortality') {
  cohort <- cohort %>% 
    mutate(depression_base = depression_d1)
}

if (exposure == 'antihypertensive' & outcome == 'all-cause mortality') {
  cohort <- cohort %>% 
    mutate(hypertension_base = hypertension_d1)
}

# add code for diabetes...

saveRDS(cohort, file = paste(path_intermediate_res, 'cohort_update_cov.rds', sep='/'))

