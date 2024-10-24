## ---------------------------
##
## Program: 2c. Censoring - Flexible Grace Period
##
## Purpose: Change the grace period to be flexible for sensitivity analyses.
##
## Author: Gwen Aubrac
##
## Date Created: 2024-10-15
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive, antidiabetic
exposure <- 'antidepressant'

#### LOAD PACKAGES ####

library(lubridate)
library(dplyr)
library(magrittr)
library(haven) 
library(parallel)
library(survival) 
library(ggplot2)

#### DEFINE PATHS ####

path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, 'all-cause mortality', 'main', 'intermediate', sep = '/')
path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, 'all-cause mortality', 'flexible_grace_period', 'intermediate', sep = '/')
cohort <- readRDS(file = paste(path_intermediate_res_main, 'cohort_creation.rds', sep = '/'))

study_start = ymd(20190101)
study_end = ymd(20221231) 
study_follow_up_end = ymd(20240331) 

setwd(path_intermediate_res)

censor_desc <- "censor_desc_flex.txt"
writeLines("Censoring description:", censor_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = censor_desc, append= TRUE)

#### SENSITIVITY ANALYSIS: FLEXIBLE GRACE PERIOD ####
# The grace period for each patient is set to be the duration of the supply for each rx

trt_supply_flex <- readRDS(file = paste(path_intermediate_res_main, 'trt_supply_raw.rds', sep = '/'))
summary(trt_supply_flex$date)

# define three scenarios for discontinuation: single refill for trt, gap in refill, or last refill
trt_supply_flex$disc_date_single <- NA
trt_supply_flex$disc_date_gap <- NA

# for patients with only 1 refill, disc_date is the first supply_end date + rx duration
gc()
trt_supply_flex %<>%
  mutate (disc_date_single = if_else (max_refills == 1, supply_end+duration, disc_date_single))

groups(trt_supply_flex) # check that grouped by id

# for patients with 1+ refills, disc_date is the supply_end date prior to gap + rx duration
trt_supply_flex %<>%
  mutate (prior_supply_end = lag(supply_end),
          prior_rx_duration = lag(duration),
          gap = if_else (!is.na(gap_in_supply) & gap_in_supply>prior_rx_duration, 1, 0),
          disc_date_gap = if_else (gap == 1, prior_supply_end+prior_rx_duration, disc_date_gap))

trt_supply_flex %<>% # set earliest disc date, produces warning but still runs (returns inf for patients with no gap)
  mutate (disc_date_gap = min(disc_date_gap, na.rm = TRUE))

# for patients with no gap, disc_date is the latest supply_end + rx duration
# so we identify patients with NA (or infinite) for single and gap disc dates
# and compute last supply disc date as last supply end + duration for those patients
# then merge dates back with df

last_disc_dates_flex <- trt_supply_flex %>%
  filter(all(is.na(disc_date_single)) & all(is.na(disc_date_gap) | is.infinite(disc_date_gap))) %>%
  arrange(desc(supply_end)) %>% 
  slice(1) %>% 
  mutate(disc_date_last = supply_end + duration)

last_disc_dates_flex %<>%
  select(id, disc_date_last)

trt_supply_flex <- merge(trt_supply_flex, last_disc_dates_flex, by = 'id', all.x = TRUE)

# extract earliest disc date from all three scenarios (should be only one date per patient, but to be sure)
# and if disc_date is after study_follow_up_end, patient is not considered to have discontinued

groups(trt_supply_flex)
gc()

disc_dates_flex <- trt_supply_flex %>% 
  group_by(id) %>% 
  mutate (disc_date = min(disc_date_single, disc_date_gap, disc_date_last, na.rm = TRUE),
          disc_date = if_else(disc_date>study_follow_up_end, NA, disc_date),
          disc_date_gap = if_else(is.infinite(disc_date_gap), NA, disc_date_gap)) %>% 
  slice(1) 

length(unique(disc_dates_flex$id))
length(disc_dates_flex$id)

# describe number of patients who discontinued based on 3 scenarios
# note: not representative of final cohort, since we are only considering treatment supply here
# (and not other reasons for censoring that might precede trt discontinuation)

cat('Number of patients who discontinued their trt during study:', length(which(!is.na(disc_dates_flex$disc_date))))
cat(paste('Number of patients who discontinued their trt during study:', length(which(!is.na(disc_dates_flex$disc_date))), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients who did not discontinue their trt during study:', length(which(is.na(disc_dates_flex$disc_date))))
cat(paste('Number of patients who did not discontinue their trt during study:', length(which(is.na(disc_dates_flex$disc_date))), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients who had a single refill:', length(which(disc_dates_flex$disc_date == disc_dates_flex$disc_date_single)))
cat(paste('Number of patients who had a single refill:', length(which(disc_dates_flex$disc_date == disc_dates_flex$disc_date_single)), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients whose last refill supply ended prior to study end:', length(which(disc_dates_flex$disc_date == disc_dates_flex$disc_date_last)))
cat(paste('Number of patients whose last refill supply ended prior to study end:', length(which(disc_dates_flex$disc_date == disc_dates_flex$disc_date_last)), '\n'), file = censor_desc, append = TRUE)

disc_dates_flex %<>% select (id, disc_date) 
cohort <- merge(cohort, disc_dates_flex, by = 'id', all.x = TRUE)

rm(trt_supply_flex, disc_dates_flex, last_disc_dates_flex)
saveRDS(cohort, file = paste(path_intermediate_res, 'cohort_censored.rds', sep='/'))
