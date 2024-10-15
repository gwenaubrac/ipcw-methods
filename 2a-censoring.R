## ---------------------------
##
## Program: 2a. Censoring (main)
##
## Purpose: Identify the date for the following censoring events in the as-treated analysis:
## 1. Patients are censored when there is a 30+ day gap (or otherwise defined grace period) in the treatment class that led to cohort entry.
## 2. Patients are censored at the date they switch to another treatment for the same indication (identified in previous program).
##
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes: The grace period defines the number of days for what is considered a gap in treatment supply
## and can be adapted to be more or less conservative.
##
## Algorithm used to define treatment duration was developed by Pauline Reynier at the Lady Davis Institute. 
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
library(data.table)
library(readxl)
library(writexl)

#### DEFINE PATHS ####

path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort.rds', sep = '/'))
rx_for_exposure <- readRDS(file = paste(path_cohort, 'rx_for_exposure.rds', sep ='/'))

study_start = ymd(20190101)
study_end = ymd(20221231)
study_follow_up_end = ymd(20240331)

setwd(path_cohort)
censor_desc <- "censor_desc.txt"
writeLines("Censoring description:", censor_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = censor_desc, append= TRUE)

#### DEFINE PRESCRIPTION DURATION ####

# This algorithm to define treatment duration in the CPRD was developed by Pauline Reynier at the Lady Davis Institute

common_dosages <- fread('Z:/EPI/Protocol 24_004042/202406_lookup/202406_Lookups_CPRDAurum/202406_Lookups_CPRDAurum/common_dosages.txt')
common_dosages %<>% select(dosageid, dosage_text, daily_dose)
summary(rx_for_exposure$duration == 0)

## Get daily dose

# retrieve dose information
rx_for_exposure <- merge(rx_for_exposure, common_dosages, by = 'dosageid', all.x = TRUE)
length(which(rx_for_exposure$daily_dose == 0))

zero <- rx_for_exposure %>% 
  filter(daily_dose == 0 & dosage_text != '')

dosage_text <- unique(zero$dosage_text) 
dosage_text_sum <- zero %>% 
  group_by(dosage_text) %>% 
  summarize(count = n()) %>% 
  arrange(desc(count))

# replace 0 or empty doses based on dosage text
rx_for_exposure <- rx_for_exposure %>%
  mutate(daily_dose = if_else(
    (daily_dose == 0 | daily_dose == '') &
      dosage_text %in% c(
        '1',
        '1 AS DIRECTED',
        '1 AS DIRECTED WHEN REQUIRED',
        '1 AS REQUIRED',
        '1 CAPSULE A DAY',
        '1 D',
        '1 WHEN NEEDED',
        '1D',
        'APPLY AS DIRECTED',
        'AS ADV',
        'AS DIR’,’AS DIRECTED',
        'AS DIRECTED BY CONSULTANT',
        'AS DIRECTED BY HOSPITAL',
        'AS DIRECTED BY PAEDIATRICIAN',
        'AS DIRECTED BY PAEDS',
        'AS DIRECTED BY SPECIALIST',
        'AS DIRECTED BY THE HOSPITAL',
        'AS DIRECTED FROM HOSPITAL',
        'AS NEEDED',
        'AS PER INSTRUCTION',
        'ASD',
        'ASD BY HOSPITAL',
        'ASD BY SPECIALIST',
        'ASDIR',
        'D',
        'ID',
        'MDU',
        'ONE AS DIRECTED',
        'ONE TO BE TAKEN AS REQUIRED',
        'ONE WHEN NEEDED',
        'TAKE AS DIRECTED',
        'TAKE ONE AS DIRECTED',
        'TAKE ONE AS NEEDED',
        'TAKE ONE CAPSULE A DAY',
        'TAKE ONE WHEN NEEDED',
        'TAKE ONE WHEN REQUIRED',
        'TO BE TAKEN AS DIRECTED',
        'TO BE USED AS DIRECTED',
        'USE AS DIECTED',
        'USE AS DIRECTED',
        'USE AS DIRECTED BY HOSPITAL',
        'WHEN REQUIRED',
        'AS DIRECTED',
        'AS DIR',
        "ONE TABLET AS DIRECTED",
        "TAKE ONE",
        "1 TABLET AS DIRECTED",
        "TAKE ONE TABLET",
        "ONE TO BE TAKEN",
        "USE ONE A DAY",
        "ONCE NIGHTLY",
        "ONE A TNIGHT",
        "AS DIRECTED.",
        "AS DIRETED",
        "TAKEN AS DIRECTED",
        "AS DISCUSSED"
      ),
    1,
    daily_dose
  ))

rx_for_exposure <- rx_for_exposure %>%
  mutate(daily_dose = if_else(
    (daily_dose == 0 | daily_dose == '') &
      dosage_text %in% c('HALF A TABLET AS REQUIRED'),
    0.5,
    daily_dose
  ))


rx_for_exposure <- rx_for_exposure %>%
  mutate(daily_dose = if_else(
    (daily_dose == 0 | daily_dose == '') &
      dosage_text %in% c(
        '1 -2 AS REQUIRED',
        'ONE OR TWO TO BE TAKEN AS DIRECTED',
        '1-2 AS REQUIRED',
        'ONE OR TWO TO BE TAKEN AS REQUIRED',
        'TAKE 1 OR 2 AS DIRECTED',
        'TAKE 1-2 WHEN REQUIRED',
        'TAKE ONE OR TWO AS DIRECTED',
        "1-2 AS DIRECTED",
        "1-2 AS NEEDED",
        "1-2"
      ),
    1.5,
    daily_dose
  ))


rx_for_exposure <- rx_for_exposure %>%
  mutate(daily_dose = if_else(
    (daily_dose == 0 | daily_dose == '') &
      dosage_text %in% c('2 D','2D', 'TAKE TWO CAPSULES A DAY', "2 CAPSULES A DAY"),
    2,
    daily_dose
  ))

# still have some patients with missing daily dose
# with following dosage text
length(which(rx_for_exposure$daily_dose == 0))
zero <- rx_for_exposure %>% 
  filter(daily_dose == 0 & dosage_text != '')
dosage_text <- unique(zero$dosage_text) 
dosage_text

length(which(rx_for_exposure$qty == 0))
length(which(rx_for_exposure$daily_dose == 0))
length(which(cohort$id %in% rx_for_exposure$id))

## Choose between duration1 (from CPRD) and duration2 (calculated from qty/daily_dose)

rx_for_exposure$duration1 <- NA
rx_for_exposure %<>% mutate(duration1 = if_else(1<=duration & duration<=183, duration, NA))
rx_for_exposure %<>% mutate(duration2 = if_else(qty>0 & daily_dose>0, round(qty/daily_dose, 1), NA))
rx_for_exposure %<>% mutate(duration2 = if_else(duration2<1 | duration2>183, NA, duration2))

rx_for_exposure %<>% mutate(same_duration = if_else(!is.na(duration1) & !is.na(duration2) & duration1 == duration2, 1, 0))
length(which(is.na(rx_for_exposure$duration1)))
length(which(is.na(rx_for_exposure$duration2)))
table(rx_for_exposure$same_duration)
cat(paste("Same duration:", table(rx_for_exposure$same_duration)), file = censor_desc, append= TRUE)

# same duration: use CPRD duration
rx_for_exposure$duration <- NA
rx_for_exposure %<>% mutate(duration = if_else(duration1>0 & duration2>0 & same_duration == 1, duration1, duration))

# different duration: use standard 28 days for a duration for a prescription
rx_for_exposure %<>% mutate(duration = if_else(
  is.na(duration) &
    same_duration == 0 &
    duration1 == 28 &
    (duration2 > 28 & duration2 <= 30),
  duration1,
  duration
))

rx_for_exposure %<>% mutate(duration = if_else(
  is.na(duration) &
    same_duration == 0 &
    duration2 == 28 &
    (duration1 > 28 & duration1 <= 30),
  duration2,
  duration
))

# different duration: durations are within 2 days then take max
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      !is.na(duration1) & !is.na(duration2) &
      duration1 - duration2 >= -2 | duration1 - duration2 <= 2,
    max(duration1, duration2),
    duration
  )
)

# different duration: duration1 = 1 and there is a value for duration 2 then take duration2
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      duration1 == 1 & !is.na(duration2),
    duration2,
    duration
  )
)

# different duration: duration2 = 1 and there is a value for duration 1 then take duration1
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      duration2 == 1 & !is.na(duration1),
    duration1,
    duration
  )
)

# daily dose != 1 and duration1 = qty then take duration2
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      daily_dose != 1 & qty == duration1 & daily_dose %in% c(0.5, 1.5, 2),
    duration2, 
    duration
  )
)

# daily dose != 1 and duration1 = 28 then take duration2
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      !is.na(duration2) & duration1 == 28 & daily_dose %in% c(0.5, 1.5, 2),
    duration2, 
    duration
  )
)


# one of the durations is equal to 28 days then use 28 days
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      !is.na(duration1) &
      !is.na(duration2) &
      (duration1 == 28 | duration2 == 28),
    28,
    duration
  )
)

# duration1 and duration2 are different with no specific criteria then take duration2
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      !is.na(duration1) &
      !is.na(duration2),
    duration2,
    duration
  )
)

# only duration1 available then take duration1
rx_for_exposure %<>% mutate(duration = if_else(
  is.na(duration) &
    same_duration == 0 &
    !is.na(duration1) &
    is.na(duration2),
  duration1,
  duration
))

# only duration2 available then take duration2
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      is.na(duration1) &
      !is.na(duration2),
    duration2,
    duration
  )
)

# duration1 and duration2 are not available but qty = 28 or 56 use qty
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      is.na(duration1) & 
      is.na(duration2) &
      daily_dose %in% c(0, '') &
      qty %in% c(28, 56),
    qty,
    duration
  )
)

# still no duration but valid qty then use qty
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      is.na(duration1) & 
      is.na(duration2) &
      qty>=1 & qty<=183,
    qty,
    duration
  )
)


# qty > 183
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration) &
      same_duration == 0 &
      daily_dose %in% c(0, '') &
      qty > 183, 
    28,
    duration
  )
)

# all other scenarios use 28
rx_for_exposure %<>% mutate(
  duration = if_else(
    is.na(duration),
    28,
    duration
  )
)

## Describe duration for the prescription

summary_rx_duration <- rx_for_exposure %>% 
  group_by(duration) %>% 
  summarize(count = n()) %>% 
  arrange(desc(count))

write_xlsx(summary_rx_duration, paste(path_cohort, 'summary_rx_duration.xlsx', sep ='/'))
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_for_exposure_with_duration.rds', sep='/'))
rx_for_exposure <- readRDS(file = paste(path_cohort, 'rx_for_exposure_with_duration.rds', sep='/'))

#### DEFINE GRACE PERIOD ####

# grace period: number of days for gap to be considered trt discontinuation
#readRDS(paste(path_cohort, 'rx_for_exposure_with_duration.rds', sep='/'))

grace_period <- 30

#### CENSORING DUE TO TREATMENT DISCONTINUATION ####

# retrieve all prescriptions for the trt that led to cohort entry
trt_supply <- rx_for_exposure %>% 
  filter(id %in% cohort$id) %>% 
  select(-dosageid, -dosage_text, -daily_dose, -duration1, -duration2, -same_duration, -ProductName, -product_code, -qty)

id_exposures <- cohort %>% select (id, trt, entry_date)

trt_supply <- merge(trt_supply, id_exposures, by.x = c('id', 'trt', 'date'), by.y = c('id', 'trt', 'entry_date'), all.x = TRUE)

length(unique(trt_supply$id)) # should be same as number of patients in cohort

# calculate number of days between end of rx supply and date of next rx
gc()
trt_supply <- trt_supply %>% 
  group_by (id) %>% 
  arrange (id, date) %>% 
  dplyr::select(id, trt, date, duration, first_rx) %>% 
  mutate (refill_number = row_number(),
          supply_end = date+duration,
          gap_in_supply = date - lag(supply_end), # gap between current rx and end of supply from previous rx
          max_refills = max(refill_number, na.rm = TRUE)) # total number of refills for each patient

saveRDS(trt_supply, file = paste(path_cohort, 'trt_supply_raw.rds', sep='/'))
# trt_supply <- readRDS(file = paste(path_cohort, 'trt_supply_raw.rds', sep = '/'))

# examine distribution of rx dates
summary(trt_supply$date)

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
