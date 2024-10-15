## ---------------------------
##
## Program: 4a. Outcomes and follow-up
##
## Purpose: Flag the occurrence and date of first occurrence of the outcome of interest since cohort entry.
## Here, the outcome is all-cause mortality, where date of death (dod) is the earliest among:
## Emis dod, CPRD dod, HES dod, and ONS dod. 
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
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes: For this project, both CPRD and HES are used to identify outcomes. The date of first occurrence
## of the outcome is first identified in CPRD and HES separately, and then the earliest between both
## is used to define first outcome occurrence date. 
##
## ---------------------------

# analysis: main, flex_grace_period, or 90_day_grace_period

analysis <- ''

#### LOAD PACKAGES ####

library(lubridate)
library(readxl)
library(dplyr)
library(magrittr)
library(haven)
library(parallel)
library(data.table)

#### SPECIFY STUDY DESIGN ####

study_start = ymd(20190101)
study_end = ymd(20221231)
study_follow_up_end = ymd(20240331)

#### DEFINE PATHS ####

if (analysis == 'main' | analysis == '') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
} else if (analysis == 'flex_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/flex_grace_period" 
} else if (analysis == '90_day_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/90_day_grace_period" 
} 

path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
path_linkage_1 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt1'
path_linkage_2 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt2'

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))
switched_to <- readRDS(file = paste(path_main, 'switched_to.rds', sep = '/'))

setwd(path_cohort)

outcome_desc <- "outcome_desc.txt"
writeLines("Outcome description:", outcome_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = outcome_desc, append= TRUE)

#### A. DEFINE ALL-CAUSE MORTALITY/DATE OF DEATH ####

# merge with ONS dod
col_classes <- c("character", "character", "character", "character", "character", "character", "character", "character")
death_ons1 <- fread(paste(path_linkage_1, 'death_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
death_ons2 <- fread(paste(path_linkage_2, 'death_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
death_ons <- bind_rows(death_ons1, death_ons2)

rm(death_ons1, death_ons2)

length(which(cohort$id %in% death_ons$patid))
cat(paste('Patients with ONS information:', length(which(cohort$id %in% death_ons$patid)), '\n'), file = outcome_desc, append= TRUE)

death_ons %<>% 
  select(patid, dod) %>% 
  rename(ons_dod = dod)

cohort <- merge(cohort, death_ons, by.x = 'id', by.y = 'patid', all.x = TRUE)

cohort$ons_dod <- dmy(cohort$ons_dod)

dod_match <- cohort %>% 
  select(id, emis_dod, cprd_dod, ons_dod) %>% 
  mutate(same_dod = if_else(!is.na(emis_dod) & emis_dod == cprd_dod &
                            !is.na(cprd_dod) & cprd_dod == ons_dod &
                            !is.na(ons_dod), 1, 0),
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

cat('Number of patients with issues in date of death encoding:', length(which(cohort$dod_issue==1)))
cat(paste('Number of patients with issues in date of death encoding:', length(which(cohort$dod_issue==1)), '\n'), file = outcome_desc, append = TRUE)

length(which(cohort$dod_issue==0))

cohort %<>% filter (dod_issue == 0) %>% select (-dod_issue)
length(unique(cohort$id))

summary(cohort$dod)

cat('Number of patients who experienced event on day of cohort entry (removed):', length(which(cohort$dod == cohort$entry_date)))
cat(paste('Number of patients who experienced event on day of cohort entry (removed):', length(which(cohort$dod == cohort$entry_date)), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who experienced event prior to cohort entry (removed):', length(which(cohort$dod < cohort$entry_date)))
cat(paste('Number of patients who experienced event prior to cohort entry (removed):', length(which(cohort$dod < cohort$entry_date)), '\n'), file = outcome_desc, append = TRUE)

cohort <- cohort %>% 
  filter (entry_date<dod | is.na(dod))

cat('Number of patients who died during study period:', length(which(!is.na(cohort$dod))))
cat(paste('Number of patients who died during study period:', length(which(!is.na(cohort$dod))), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who died during study period on SNRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'snri')))
cat(paste('Number of patients who died during study period on SNRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'snri')), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who died during study period on SSRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'ssri')))
cat(paste('Number of patients who died during study period on SSRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'ssri')), '\n'), file = outcome_desc, append = TRUE)


#### DEFINING FOLLOW-UP ####

## ITT
# define ITT exit excluding event timing
cohort$itt_exit_date <- pmin(
  cohort$lcd, 
  cohort$regend,
  cohort$dod,
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

# define ITT event
# no event if no date of death or if died after ITT exit
cohort %<>%
  mutate (itt_event = if_else(is.na(dod) | (!is.na(dod) & dod > itt_exit_date), 0, 1),
          itt_event_date = if_else(itt_event == 1, dod, NA))

cohort$itt_follow_up <- time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$itt_exit_date)), 'days')

## AT
# define AT exit excluding event timing
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
cohort <- cohort %>% 
  mutate (at_event = if_else(itt_event == 1 & itt_event_date <= at_exit_date, 1, 0),
          at_event_date = if_else(itt_event_date <= at_exit_date, itt_event_date, NA))

cohort$at_follow_up <- time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$at_exit_date)), 'days')

# define discontinuation and switch flags
# disc: patients discontinued if they exited the cohort because their trt disc occurred first
cohort <- cohort %>% 
  mutate(
    disc = if_else(
      !is.na(disc_date) & at_exit_date == disc_date & # if there is a disc date and it led to cohort exit
        ((!is.na(itt_event_date) & at_exit_date != itt_event_date) | is.na(itt_event_date)), # and if there was an event then disc date was not on the day of the event
      1, 
      0
    )
  )

# switch: patients switched if they exited the cohort because their switch occurred first
cohort <- cohort %>% 
  mutate(
    switch = if_else(
      !is.na(switch_date) & at_exit_date == switch_date & 
        ((!is.na(itt_event_date) & at_exit_date != itt_event_date) | is.na(itt_event_date)), 
      1, 
      0
    )
  )


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

# ensure no patients who discontinued have a follow-up less than the grace period or have an event
test <- cohort %>% 
  filter(disc == 1)

summary(test$at_follow_up)
summary(test$at_event)

test <- cohort %>% # should be empty
  filter(disc == 1 & at_event == 1) %>% 
  select(id, entry_date, at_event_date, disc_date)

#### DEFINE OVERALL CENSOR DATE FROM TRT DISCONTINUATION OR SWITCH ####

cohort %<>%
  mutate (censor = if_else (switch == 1 | disc == 1, 1, 0),
          censor_date = if_else(censor == 1, pmin(switch_date, disc_date, na.rm = TRUE), NA))

table(cohort$censor)
length(which(is.na(cohort$censor_date)))

cat('Number of patients censored due to trt switch or discontinuation:', sum(cohort$censor))
cat(paste('Number of patients censored due to trt switch or discontinuation:', sum(cohort$censor), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients censored due to trt switch:', sum(cohort$switch))
cat(paste('Number of patients censored due to trt switch:', sum(cohort$switch), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients censored due to trt discontinuation:', sum(cohort$disc))
cat(paste('Number of patients censored due to trt discontinuation:', sum(cohort$disc), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients not censored:', sum(cohort$censor == 0))
cat(paste('Number of patients not censored:', sum(cohort$censor == 0), '\n'), file = outcome_desc, append = TRUE)

# describe trt switches
switched_to %<>% select(id, trt_seq)
cohort <- merge(cohort, switched_to, by='id', all.x = TRUE)
table(cohort[cohort$switch ==1, 'trt_seq'])
cat(paste('Switch SNRI to SSRI:', table(cohort[cohort$switch ==1, 'trt_seq']), '\n')[[1]], file = outcome_desc, append = TRUE)
cat(paste('Switch SSRI to SNRI:', table(cohort[cohort$switch ==1, 'trt_seq']), '\n')[[2]], file = outcome_desc, append = TRUE)

#### UPDATE PATIENT COVARIATE VALUES AT INTERVALS ####

## Set up

# define deciles of censoring distribution in follow-up counting time
# based on distribution of censoring times among censored patients
censored_only <- cohort %>%
  filter(!is.na(censor_date))

censoring_times <- as.numeric(censored_only$censor_date - censored_only$entry_date)
quantile(censoring_times, probs = seq(0, 1, by = 0.05))
cat(paste('The distribution of censoring times is:', quantile(censoring_times, probs = seq(0, 1, by = 0.05)), '\n'), file = cov_desc, append = TRUE)

# there is a concentration of censoring at 58 days (d1 = d2 = d3)
# so let us split data into deciles starting from the 3rd decile (~35%)
# to be more flexible
breaks <- round(seq(from = 0.35, to = 1, length.out = 10), 2)
breaks
saveRDS(breaks, file = paste(path_cohort, 'interval_breaks.rds', sep = '/'))

times_dec <- quantile(censoring_times, probs = breaks)
times_dec

length(which(censoring_times < times_dec[[1]]))
length(which(censoring_times >= times_dec[[1]] & censoring_times < times_dec[[2]]))
length(which(censoring_times == times_dec[[1]]))

# let us shift the concentration of censoring to the first decile
# by setting the second decile at 59 days (for main)
times_dec[[1]] <- 59
times_dec

if (analysis == '90_day_grace_period') {
  times_dec[[1]] <- 119
} else if (analysis == 'flex_grace_period') {
  times_dec[[1]] <- 57
}

saveRDS(times_dec, file = paste(path_cohort, 'times_dec.RDS', sep = '/'))

cat('The deciles of censoring distribution are:', paste(times_dec, sep = ','))
cat(paste('The deciles of censoring distribution are:', paste(times_dec, sep = ','), '\n'), file = cov_desc, append = TRUE)

# get relative date of censoring quartile for each patient
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

first_comorb <- readRDS(paste(path_cohort, 'first_comorb.rds', sep = '/'))

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
base_comorb <- readRDS(paste(path_cohort, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(paste(path_cohort, 'dec_comorb.rds', sep = '/'))
names <- c('sex', 'year', 'ethnicity', 'deprivation', base_comorb, dec_comorb, 'disc', 'switch')

cohort <- cohort %>% 
  mutate(across(names, as.factor)) 

# depression tends to be undercoded in EHR data
# so to have more representative baseline depression
# let us change baseline depression to d1 depression
table(cohort$depression_base, cohort$trt)
table(cohort$depression_d1, cohort$trt)

cohort <- cohort %>% 
  mutate(depression_base = depression_d1) 

# check how many patients have anxiety at baseline
table(cohort$anxiety_base, cohort$trt)
table(cohort$anxiety_d1, cohort$trt)

saveRDS(cohort, file = paste(path_cohort, 'antidepressant_cohort_outcome.rds', sep='/'))
