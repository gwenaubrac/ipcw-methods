## ---------------------------
##
## Program: 4b. Outcomes and follow-up (clinical outcome)
##
## Purpose: Flag the occurrence and date of first occurrence of the outcome of interest since cohort entry.
## Here, the outcome is the first occurrence of a clinical outcome. 
## 
## Define follow-up for the two following analyses types: 
## - intention-to-treat analysis (ITT): patients are analyzed according to their original 
## exposure group. Follow-up ends at the earliest among: death ('dod'), occurrence of the event 
## of interest ('event_date'), departure from the CPRD ('regend'), end of linkage to HES or ONS, 
## last available data ('lcd'), or study follow-up end ('study_follow_up_end')
## - as-treated: patients are analyzed according to their actual exposure group, and are thus
## censored when they switch or discontinue their treatment. Follow-up ends at the earliest among
## the above plus treatment switch ('switch_date') or treatment discontinuation ('disc_date'). 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-10-15
##
## ---------------------------
##
## Notes: For this project, both CPRD and HES are used to identify outcomes. The date of first occurrence
## of the outcome is first identified in CPRD and HES separately, and then the earliest between both
## is used to define first outcome occurrence date. 
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive, antidiabetic
exposure <- 'antihypertensive'

# outcome: suicidal ideation
outcome <- 'suicidal ideation'

# analysis: main, flexible_grace_period, 90_day_grace_period
analysis <- 'main'

#### LOAD PACKAGES ####

library(lubridate)
library(readxl)
library(dplyr)
library(magrittr)
library(haven)
library(parallel)
library(data.table)
library(tidyr)

#### SPECIFY STUDY DESIGN ####

study_start = ymd(20190101)
study_end = ymd(20221231)
study_follow_up_end = ymd(20240331)

#### DEFINE PATHS ####

path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_linkage_1 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt1'
path_linkage_2 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt2'

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_covariates.rds', sep = '/'))

setwd(path_intermediate_res)

outcome_desc <- "outcome_desc.txt"
writeLines("Outcome description:", outcome_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = outcome_desc, append= TRUE)

if (outcome == 'suicidal ideation') {
  outcome <- 'suicidal_ideation_self_harm'
}

#### DEFINE ALL-CAUSE MORTALITY/DATE OF DEATH ####

# merge with ONS dod
col_classes <- c("character", "character", "character", "character", "character", "character", "character", "character")
death_ons1 <- fread(paste(path_linkage_1, 'death_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
death_ons2 <- fread(paste(path_linkage_2, 'death_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
death_ons <- bind_rows(death_ons1, death_ons2)

rm(death_ons1, death_ons2, col_classes)

length(which(cohort$id %in% death_ons$patid))
cat(paste('Patients with ONS information:', length(which(cohort$id %in% death_ons$patid)), '\n'), file = outcome_desc, append= TRUE)

death_ons %<>% 
  select(patid, dod) %>% 
  rename(ons_dod = dod) %>% 
  mutate(ons_dod = as.Date(ons_dod, format = '%d/%m/%Y'))

cohort <- merge(cohort, death_ons, by.x = 'id', by.y = 'patid', all.x = TRUE)

dod_match <- cohort %>% 
  select(id, emis_dod, cprd_dod, ons_dod) %>% 
  mutate(same_dod = if_else((!is.na(emis_dod) & !is.na(ons_dod) & !is.na(cprd_dod)) & # patient has dod for all sources
                              emis_dod == ons_dod & ons_dod == cprd_dod, # and all dod are the same
                            1, 0),
         no_ons_dod = if_else( (!is.na(emis_dod) | !(is.na(cprd_dod))) & is.na(ons_dod), 1, 0)
  )

length(which(dod_match$same_dod == 1))  
cat(paste('Matching CPRD/EMIS/ONS dod:', length(which(dod_match$same_dod == 1)) , '\n'), file = outcome_desc, append= TRUE)

length(which(dod_match$no_ons_dod == 1))  
cat(paste('Missing ONS dod with EMIS or CPRD dod:', length(which(dod_match$no_ons_dod == 1)) , '\n'), file = outcome_desc, append= TRUE)


# set dod as earliest among CPRD, and ONS
cohort$dod <- pmin(
  cohort$emis_dod,
  cohort$cprd_dod,
  cohort$ons_dod,
  na.rm = TRUE
)

# remove patients with issues in dod encoding
length(which(cohort$dod<cohort$birthdate))

cohort <- cohort %>% 
  mutate(dod_issue = if_else(!is.na(dod), # if patient died
                             if_else(dod<birthdate, 1, 0), # ensure date of death is after date of birth
                             0))

cat('Number of patients with issues in date of death encoding (removed):', length(which(cohort$dod_issue==1)))
cat(paste('Number of patients with issues in date of death encoding (removed):', length(which(cohort$dod_issue==1)), '\n'), file = outcome_desc, append = TRUE)

length(which(cohort$dod_issue==0))

cohort %<>% filter (dod_issue == 0) %>% select (-dod_issue)
length(unique(cohort$id))

summary(cohort$dod)

cat('Number of patients who died on day of cohort entry (removed):', length(which(cohort$dod == cohort$entry_date)))
cat(paste('Number of patients who died on day of cohort entry (removed):', length(which(cohort$dod == cohort$entry_date)), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who died prior to cohort entry (removed):', length(which(cohort$dod < cohort$entry_date)))
cat(paste('Number of patients who died prior to cohort entry (removed):', length(which(cohort$dod < cohort$entry_date)), '\n'), file = outcome_desc, append = TRUE)

cohort <- cohort %>% 
  filter (entry_date<dod | is.na(dod))

cat('Number of patients who died during study period:', length(which(!is.na(cohort$dod))))
cat(paste('Number of patients who died during study period:', length(which(!is.na(cohort$dod))), '\n'), file = outcome_desc, append = TRUE)

trt_death_tbl <- table(!is.na(cohort$dod), cohort$trt)

cat('Number of patients who died during study period on each treatment:', table(!is.na(cohort$dod), cohort$trt))
write.table(trt_death_tbl, file = "outcome_desc.txt", append = TRUE, sep = "\t", col.names = TRUE)

rm(death_ons, dod_match)

#### EXCLUDE PATIENTS WITH EVENT PRIOR TO ENTRY ####

first_comorb <- readRDS(paste(path_intermediate_res_mort_main, 'first_comorb.rds', sep = '/'))

first_comorb <- first_comorb %>%
  pivot_wider(id_cols = id, names_from = comorb, values_from = comorb_date)

first_comorb_outcome <- first_comorb %>% 
  mutate(event_date = first_comorb[[outcome]]) %>% 
  select(id, event_date)

cohort <- merge(cohort, first_comorb_outcome, by = 'id')

event_before_entry <- cohort %>% 
  filter(event_date <= entry_date)

cat('Number of patients with even prior to or on entry:', length(unique(event_before_entry$id)))
cat(paste('Number of patients with even prior to or on entry:', length(unique(event_before_entry$id)), '\n'), file = outcome_desc, append = TRUE)

cohort %<>%
  filter(is.na(event_date) | event_date > entry_date)

#### DEFINING FOLLOW-UP ####

## ITT
# define ITT exit
cohort$itt_exit_date <- pmin(
  cohort$lcd, 
  cohort$regend,
  cohort$dod,
  cohort$event_date,
  study_follow_up_end,
  na.rm = TRUE
)

# cap ITT follow-up at 730 days
cohort %<>%
  mutate(
    itt_exit_date = if_else(time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$itt_exit_date)), 'days') > 730,
                            entry_date + 730,
                            itt_exit_date)
  )

# define outcome of interest
# event if patient had outcome prior to or on ITT exit
cohort %<>%
  mutate (itt_event = if_else(!is.na(event_date) & event_date <= itt_exit_date, 1, 0),
          itt_event_date = if_else(itt_event == 1, event_date, NA))

table(cohort$itt_event)
summary(cohort$itt_event_date)

cohort$itt_follow_up <- time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$itt_exit_date)), 'days')
summary(cohort$itt_follow_up)

## AT
# define AT exit
# same as ITT with additional censoring when switch or disc
cohort$at_exit_date <- pmin(
  cohort$itt_exit_date,
  cohort$disc_date,
  cohort$switch_date,
  na.rm = TRUE
)

# cap AT follow-up at 730 days
cohort %<>%
  mutate(
    at_exit_date = if_else(time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$at_exit_date)), 'days') > 730, 
                           entry_date + 730,
                           at_exit_date)
  )

# define AT events
# AT event if patient experience outcome prior to or on AT exit
cohort <- cohort %>% 
  mutate (at_event = if_else((!is.na(event_date) & event_date <= at_exit_date), 1, 0),
          at_event_date = if_else(at_event == 1, event_date, NA))

table(cohort$at_event)
summary(cohort$at_event_date)

cohort$at_follow_up <- time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$at_exit_date)), 'days')
summary(cohort$at_follow_up)

# define discontinuation and switch flags
# disc: patients discontinued if they exited the cohort because their trt disc occurred first
# patients cannot discontinue after they die
# so ensure disc date occurs before death, if death occurred
# similarly, if patients experience event and discontinue on the same day
# they will be considered to have exited due to the event and not the disc
cohort <- cohort %>% 
  mutate(
    disc = if_else(
      !is.na(disc_date) & at_exit_date == disc_date &
        ((!is.na(dod) & disc_date < dod) | is.na(dod)) &
        ((!is.na(event_date) & disc_date != event_date) | is.na(event_date)), 
      1, 
      0
    )
  )
table(cohort$disc)

# switch: patients switched if they exited the cohort because their switch occurred first
# patients cannot switch after they die
# so ensure switch date occurs before death, if death occurred
# similarly, if patients experience event and switch on the same day
# they will be considered to have exited due to the event and not the switch
cohort <- cohort %>% 
  mutate(
    switch = if_else(
      !is.na(switch_date) & at_exit_date == switch_date & 
        ((!is.na(dod) & switch_date < dod) | is.na(dod)) &
        ((!is.na(event_date) & switch_date != event_date) | is.na(event_date)), 
      
      1, 
      0
    )
  )
table(cohort$switch)

# examine follow up times
summary(cohort$itt_follow_up)
summary(cohort$at_follow_up)

test <- cohort %>% 
  filter(itt_follow_up <= 0) %>% 
  mutate(regend_before_exit = if_else(regend <= itt_exit_date, 1, 0)) %>% 
  select(id, entry_date, regend, itt_exit_date, at_exit_date, dod, itt_follow_up, at_follow_up, regend_before_exit)

table(test$regend_before_exit)

# some patients have a negative ITT follow-up
# because their regimen ended prior to their cohort exit
length(which(cohort$itt_follow_up<0))
length(which(cohort$itt_follow_up==0)) 

cat('Number of patients with 0 or negative ITT follow-up (removed):', length(test$id))
cat(paste('Number of patients with 0 or negative ITT follow-up (removed):', length(test$id), '\n'), file = outcome_desc, append = TRUE)

cohort %<>% filter (itt_follow_up > 0)

summary(cohort$itt_follow_up)
summary(cohort$at_follow_up)

length(which(cohort$at_follow_up<=0)) 

test <- cohort %>% 
  filter(at_follow_up <= 0) %>% 
  mutate(switch_on_entry = if_else (switch_date == entry_date, 1, 0)) %>% 
  select(id, entry_date, regend, itt_exit_date, at_exit_date, dod, itt_follow_up, at_follow_up, switch_on_entry)

table(test$switch_on_entry)

# some patients have 0 AT follow-up because switched on the same day they entered
cat('Number of patients with 0 or negative AT follow-up (removed):', length(which(cohort$at_follow_up<=0)) )
cat(paste('Number of patients with 0 or negative AT follow-up (removed):', length(which(cohort$at_follow_up<=0)) , '\n'), file = outcome_desc, append = TRUE)

cohort %<>% filter (at_follow_up > 0)

summary(cohort$itt_follow_up)
summary(cohort$at_follow_up)

# ensure no patients who discontinued have a follow-up less than the grace period
# or have the event (since censored at disc and event, so cannot have both)
test <- cohort %>% 
  filter(disc == 1)

summary(test$at_follow_up)
summary(test$at_event)

test <- cohort %>% # should be empty
  filter(disc == 1 & at_event == 1) %>% 
  select(id, entry_date, at_event_date, disc_date)

saveRDS(cohort, file = paste(path_intermediate_res, 'cohort_outcome.rds', sep='/'))
