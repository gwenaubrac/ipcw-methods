## ---------------------------
##
## Program: 2b. Censoring - 90-Day Grace Period
##
## Purpose: Change the grace period to be 90-days for sensitivity analyses.
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes:
##
## ---------------------------

#### LOAD PACKAGES ####

library(lubridate)
library(dplyr)
library(magrittr)
library(haven) 
library(parallel)
library(survival) 
library(ggplot2)

#### DEFINE PATHS ####

path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main"
path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/90_day_grace_period"
cohort <- readRDS(file = paste(path_main, 'antidepressant_cohort.rds', sep = '/'))
study_follow_up_end = ymd(20240331) # end of follow-up

setwd(path_cohort)

censor_desc <- "censor_desc.txt"
writeLines("Censoring description:", censor_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = censor_desc, append= TRUE)

#### DEFINE GRACE PERIOD ####

# grace period: number of days for gap to be considered trt discontinuation
grace_period <- 90

#### CENSORING DUE TO TREATMENT DISCONTINUATION ####

# dataframe with number of days between end of rx supply and date of next rx
trt_supply <- readRDS(file = paste(path_main, 'trt_supply_raw.rds', sep = '/'))

## Define scenarios for treatment discontinuation
# 1. A patient had only 1 refill (single refill)
# 2. A patient had a gap in their refills (gap in refill)
# 3. A patient's last refill did not last beyond study end (last refill)

# 1. for patients with only 1 refill, disc_date is the first supply_end date + grace_period
trt_supply$disc_date_single <- NA
gc()
trt_supply %<>%
  mutate (disc_date_single = if_else (max_refills == 1, supply_end+grace_period, disc_date_single))

groups(trt_supply) # check that grouped by id

# 2. for patients with 1+ refills, disc_date is the supply_end date prior to gap + grace_period
trt_supply$disc_date_gap <- NA
gc()
trt_supply %<>%
  mutate (prior_supply_end = lag(supply_end),
          gap = if_else (!is.na(gap_in_supply) & gap_in_supply>grace_period, 1, 0),
          disc_date_gap = if_else (gap == 1, prior_supply_end+grace_period, disc_date_gap))

trt_supply %<>% # set earliest disc date from gap, produces warning but still runs (returns inf for patients with no gap)
  mutate (disc_date_gap = min(disc_date_gap, na.rm = TRUE))

# 3. for patients with no gap, disc_date is the latest supply_end + grace_period
# so we identify patients with NA (or infinite) for single and gap disc dates
# and compute their last supply disc date as last supply end + grace period
# then merge dates back with df

gc()
last_disc_dates <- trt_supply %>%
  filter(all(is.na(disc_date_single)) & all(is.na(disc_date_gap) | is.infinite(disc_date_gap))) %>%
  arrange(desc(supply_end)) %>% # already arranged by id
  slice(1) %>% 
  mutate(disc_date_last = supply_end + grace_period)

last_disc_dates %<>%
  select(id, disc_date_last)

trt_supply <- merge(trt_supply, last_disc_dates, by = 'id', all.x = TRUE)

# extract earliest disc date from all three scenarios (should be only one date per patient, but to be sure)
# if disc_date is after study_follow_up_end, patient is not considered to have discontinued
# and if the gap disc date is infinite (meaning no gap), we replace it with NA

groups(trt_supply)

gc()
disc_dates <- trt_supply %>% 
  group_by(id) %>% 
  mutate (disc_date = min(disc_date_single, disc_date_gap, disc_date_last, na.rm = TRUE),
          disc_date = if_else(disc_date>study_follow_up_end, NA, disc_date),
          disc_date_gap = if_else(is.infinite(disc_date_gap), NA, disc_date_gap)) %>% 
  slice(1) 

length(unique(disc_dates$id))
length(disc_dates$id)

# describe number of patients who discontinued based on 3 scenarios
# note: not representative of final cohort, since we are only considering treatment supply here
# (and not other reasons for censoring that might precede trt discontinuation)

cat('Number of patients who discontinued their trt during study:', length(which(!is.na(disc_dates$disc_date))))
cat(paste('Number of patients who discontinued their trt during study:', length(which(!is.na(disc_dates$disc_date))), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients who did not discontinue their trt during study:', length(which(is.na(disc_dates$disc_date))))
cat(paste('Number of patients who did not discontinue their trt during study:', length(which(is.na(disc_dates$disc_date))), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients who had a single refill:', length(which(disc_dates$disc_date == disc_dates$disc_date_single)))
cat(paste('Number of patients who had a single refill:', length(which(disc_dates$disc_date == disc_dates$disc_date_single)), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients who had a', grace_period, 'day gap in refills:', length(which(disc_dates$disc_date == disc_dates$disc_date_gap)))
cat(paste('Number of patients who had a', grace_period, 'day gap in refills:', length(which(disc_dates$disc_date == disc_dates$disc_date_gap)), '\n'), file = censor_desc, append = TRUE)

cat('Number of patients whose last refill supply ended prior to study end:', length(which(disc_dates$disc_date == disc_dates$disc_date_last)))
cat(paste('Number of patients whose last refill supply ended prior to study end:', length(which(disc_dates$disc_date == disc_dates$disc_date_last)), '\n'), file = censor_desc, append = TRUE)

disc_dates %<>% select (id, disc_date) 
cohort <- merge(cohort, disc_dates, by = 'id', all.x = TRUE)

rm(trt_supply, disc_dates, last_disc_dates)
saveRDS(cohort, file = paste(path_cohort, 'antidepressant_cohort_censored.rds', sep='/'))
