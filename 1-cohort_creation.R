## ---------------------------
##
## Program: 1. Cohort creation
##
## Purposes:
## - Specify study design, including study period (start, recruitment end, follow-up end)
## and codes for exposures of interest in CPRD Aurum. 
## - Create a 'cohort' dataframe of first-time users for specified exposure meeting inclusion and exclusion criteria.
## - Create a 'rx_for_exposure' dataframe containing the prescriptions of all CPRD Aurum patients for the exposures of interest 
## (to be used later to assess treatment discontinuation).
##
## Inclusion criteria:
## - Prescription for exposure of interest during specified study period. 
##
## Exclusion criteria:
## - Less than 18 years old at cohort entry
## - Less than 365 days of look-back available in the database from cohort entry (based on regstart - patient regiment start)
## - Less than 365 days of follow-up available in the database since cohort entry (based on lcd - last practice collection date)

## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes: 
## 1. The folder for the exposure codes should contain separate excel files for each 
## exposure of interest (ex: 'ssri.xlsx' contains SSRI codes and 'snri.xlsx' contains SNRI codes).
## Please adapt 'ProdCodeId' to the column name corresponding to the exposure codes in excel file. 
##
## 2. The data used for this project was converted from SAS to R. 
##
## ---------------------------

#### LOAD PACKAGES ####

library(lubridate)
library(readxl)
library(dplyr)
library(magrittr)
library(haven)
library(parallel)
library(data.table)

options(scipen = 999)

#### DEFINE PATHS ####

path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main"
path_cprdA <- "Z:/EPI/Protocol 24_004042/dataA"
path_cprdB <- "Z:/EPI/Protocol 24_004042/dataB" 
path_cprdC <- "Z:/EPI/Protocol 24_004042/dataC (no followup)" 
path_exposure <- "Z:/EPI/Protocol 24_004042/Gwen/data/exposures/antidepressants" # folder containing separate excel for exposures
path_linkage_1 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt1'
path_linkage_2 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt2'
setwd(path_cohort)

# create note to keep track of cohort creation numbers
cohort_creation_desc <- "cohort_creation_desc.txt"
writeLines("Cohort creation:", cohort_creation_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = cohort_creation_desc, append= TRUE)

# get linkable patients
col_classes <- c("character", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")
linkage_eligibility1 <- fread(paste(path_linkage_1, '24_004042_linkage_eligibility_aurum.txt', sep = '/'), colClasses = col_classes)
linkage_eligibility2 <- fread(paste(path_linkage_2, '24_004042_linkage_eligibility_aurum.txt', sep = '/'), colClasses = col_classes)
linkage_eligibility <- bind_rows(linkage_eligibility1, linkage_eligibility2)

#### SPECIFY STUDY DESIGN ####

study_start = ymd(20190101)
study_end = ymd(20221231) # end of recruitment
study_follow_up_end = ymd(20240331) # end of follow-up

cat(paste("Study start:", format(study_start, "%Y-%m-%d"), '\n'), file = cohort_creation_desc, append= TRUE)
cat(paste("Study end:", format(study_end, "%Y-%m-%d"), '\n'), file = cohort_creation_desc, append= TRUE)
cat(paste("Follow-up end:", format(study_follow_up_end, "%Y-%m-%d"), '\n'), file = cohort_creation_desc, append= TRUE)

# dataframe containing all exposure codes and their corresponding treatment
exposure <- data.frame()
exposure_files <- list.files(path_exposure, pattern = '.xlsx', all.files = TRUE, full.names = TRUE)

for (file in exposure_files) {
  data <- read_excel(file, trim_ws = TRUE, col_type = 'text')
  data$trt <- tools::file_path_sans_ext(basename(file))
  exposure <- bind_rows(exposure, data)
}

rm(data, file)

#### INCLUSION CRITERIA ####

## Identify first prescription for exposure of interest occurring during the study period.
## Cohort entry date is defined as the date of the first prescription for exposure of interest.
## If there are 1+ first prescriptions (meaning patient received rx for different treatment classes), 
## the second one indicates a treatment switch date. 
## We end up with 'first_rx_study' that contains information on the 1st prescription during study period + switch date.

## Iterate through data tables and retain those containing product codes for exposure

filter_prescriptions <- function (file_path) {
  file <- read_sas(file_path) 
  file_filtered <- file %>%
    dplyr::select(id, product_code, date, dosageid, qty, duration) %>% 
    filter (product_code %in% exposure$ProdCodeId) %>%
    group_by (id) %>% 
    mutate (first_rx = min(date)) %>%
    filter (first_rx >= study_start & first_rx <= study_follow_up_end) %>% 
    arrange (id, date)
  rx_for_exposure <<- bind_rows(rx_for_exposure, file_filtered)
}


therapy_filesA <- list.files(
  path_cprdA,
  pattern = 'therapy',
  all.files = TRUE,
  full.names = TRUE
)

therapy_filesA1 <- therapy_filesA[1:60]
therapy_filesA2 <- therapy_filesA[61:113]

gc()
rx_for_exposure <- data.frame()
mclapply(therapy_filesA1, filter_prescriptions)
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_A1.rds', sep='/'))
rx_A1 <- readRDS(file = paste(path_cohort, 'rx_A1.rds', sep = '/'))
rm(rx_for_exposure)

gc()
rx_for_exposure <- data.frame()
mclapply(therapy_filesA2, filter_prescriptions)
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_A2.rds', sep='/'))
rm(rx_for_exposure)

therapy_filesB <- list.files(
  path_cprdB,
  pattern = 'therapy',
  all.files = TRUE,
  full.names = TRUE
)

therapy_filesB1 <- therapy_filesB[1:60]
therapy_filesB2 <- therapy_filesB[61:120]

gc()
rx_for_exposure <- data.frame()
mclapply(therapy_filesB1, filter_prescriptions)
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_B1.rds', sep='/'))
rm(rx_for_exposure)

gc()
rx_for_exposure <- data.frame()
mclapply(therapy_filesB2, filter_prescriptions)
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_B2.rds', sep='/'))
rm(rx_for_exposure)

therapy_filesC <- list.files(
  path_cprdC,
  pattern = 'therapy',
  all.files = TRUE,
  full.names = TRUE
)

gc()
rx_for_exposure <- data.frame()
mclapply(therapy_filesC, filter_prescriptions)
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_C.rds', sep='/'))
rm(rx_for_exposure)

rx_A1 <- readRDS(file = paste(path_cohort, 'rx_A1.rds', sep = '/'))
rx_A2 <- readRDS(file = paste(path_cohort, 'rx_A2.rds', sep = '/'))
rx_B1 <- readRDS(file = paste(path_cohort, 'rx_B1.rds', sep = '/'))
rx_B2 <- readRDS(file = paste(path_cohort, 'rx_B2.rds', sep = '/'))
rx_C <- readRDS(file = paste(path_cohort, 'rx_C.rds', sep = '/'))

rx_for_exposure <- rbind(rx_A1, rx_A2, rx_B1, rx_B2, rx_C)
rx_for_exposure <- merge(rx_for_exposure, exposure, by.x = 'product_code', by.y = 'ProdCodeId', all.x = TRUE)
saveRDS(rx_for_exposure, file = paste(path_cohort, 'rx_for_exposure.rds', sep ='/'))

length(unique(rx_for_exposure$id))
length(rx_for_exposure$id)

## Select earliest prescription for each exposure group

# split prescriptions by treatment
rx_groups <- split(x = rx_for_exposure, f = rx_for_exposure$trt, drop = FALSE)
first_rx <- data.frame()

for (i in 1:length(rx_groups)) {
  rx_name <- names(rx_groups[i])
  rx_df <- rx_groups[[i]]
  rx_df %<>%
    group_by (id) %>% 
    arrange (date) %>% 
    slice (1)
  first_rx <- bind_rows(first_rx, rx_df)
}

length(unique(first_rx$id))
length(first_rx$id)

switched_to <- first_rx %>% # for later use, to determine what trt patients switched to
  group_by(id) %>% 
  arrange(date) %>% 
  summarize(trt_seq_list = list(trt))

switched_to <- switched_to %>% 
  group_by(id) %>% 
  mutate(trt_seq = paste(trt_seq_list[[1]][1], trt_seq_list[[1]][2], sep =' to '))

# identify patients with a first rx for 1+ treatments
multiple_first_rx <- first_rx %>% 
  group_by(id) %>% 
  arrange(date) %>% 
  mutate (trt_count = n()) %>% 
  filter(trt_count>1) %>% 
  arrange(id) %>% 
  slice(1)

cat('Number of first rx for any trt:', length(first_rx$id))
cat(paste('Number of first rx for any trt:', length(first_rx$id), '\n'), file = cohort_creation_desc, append = TRUE)

for (trt in unique(exposure$trt)) {
  cat('Number of first rx for', trt, ':', length(which(first_rx$trt == trt)), '\n')
  cat(paste('Number of first rx for', trt, ':', length(which(first_rx$trt == trt)), '\n'), file = cohort_creation_desc, append = TRUE)
}

cat('Number of patients with multiple first rx:', length(multiple_first_rx$id))
cat(paste('Number of patients with multiple first rx:', length(multiple_first_rx$id), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of unique patients with first rx for either trt:', length(unique(first_rx$id)))
cat(paste('Number of unique patients with first rx for either trt:', length(unique(first_rx$id)), '\n'), file = cohort_creation_desc, append = TRUE)

rm(rx_name, rx_df, rx_groups, i, trt, multiple_first_rx)

# if a patient has 1+ first rx, the first one designates cohort entry
# and the second one designates the date of trt switch
first_rx %<>% 
  group_by (id) %>% 
  arrange (date) %>% 
  mutate (trt_count = n()) %>% 
  mutate (switch_date = if_else (trt_count>1, date[2], NA)) %>% 
  arrange (id)

groups(first_rx)

first_rx %<>% 
  slice (1) %>% 
  rename (entry_date = date,
          trt_code = product_code)

length(unique(first_rx$id))
length(first_rx$id)

# keep only first rx occurring during study period
length(which(first_rx$entry_date < study_start))
length(which(first_rx$entry_date > study_end))
cat('Number of patients who initiated trt outside of study period:', length(which(first_rx$entry_date < study_start | first_rx$entry_date > study_end)))
cat(paste('Number of patients who initiated trt outside of study period:', length(which(first_rx$entry_date < study_start | first_rx$entry_date > study_end)), '\n'), file = cohort_creation_desc, append = TRUE)

first_rx_study <- first_rx %>%
  filter (study_start <= entry_date & entry_date <= study_end) %>% 
  dplyr::select (id, trt, trt_code, entry_date, switch_date)

cat('Number of patients who initiated trt during study period:', length(first_rx_study$id))
cat(paste('Number of patients who initiated trt during study period:', length(first_rx_study$id), '\n'), file = cohort_creation_desc, append = TRUE)

summary(first_rx_study$entry_date)

rm(rx_A1, rx_A2, rx_B1, rx_B2, rx_C)

#### EXCLUSION CRITERIA ####

## Exclude patients who are:
## 1. >18 years old
## 2. Less than 365 days look-back available since cohort entry
## 3. Less than 365 days of follow-up available from cohort entry. 
## 4. No linkage to HES or ONS

## Extract patient information

# read patient data
patients1 <- read_sas(paste(path_cprdA, 'patient1.sas7bdat', sep = '/'))
patients2 <- read_sas(paste(path_cprdB, 'patient1.sas7bdat', sep = '/'))
patients3 <- read_sas(paste(path_cprdC, 'patient1.sas7bdat', sep = '/'))
patients <- bind_rows(patients1, patients2, patients3)

length(unique(patients$id))
length(patients$id)
cat(paste('Number of patients in CPRD files:', length(unique(patients$id)), '\n'), file = cohort_creation_desc, append = TRUE)

# join with practice data
practice1 <- read_sas(paste(path_cprdA, 'practice1.sas7bdat', sep = '/'))
practice2 <- read_sas(paste(path_cprdB, 'practice1.sas7bdat', sep = '/'))
practice3 <- read_sas(paste(path_cprdB, 'practice1.sas7bdat', sep = '/'))
practice <- bind_rows(practice1, practice2, practice3)
length(unique(practice$practice))
length(practice$practice)

# may have duplicate practices since merged data from males and females
# so remove duplicates
practice %<>%
  group_by(practice) %>% 
  slice(1)

patients_first_rx <- merge(first_rx_study, patients, by = 'id', all.x = TRUE)
length(unique(patients_first_rx$id))
length(patients_first_rx$id)

patients_first_rx <- patients_first_rx %>%
  left_join (practice, by = 'practice', relationship = 'many-to-one')
length(unique(patients_first_rx$id))
length(patients_first_rx$id)

table(patients_first_rx$acceptable)
summary(patients_first_rx$lcd)
length(which(!is.na(patients_first_rx$lcd)))

## 1. Apply exclusion criteria (1): >18 at cohort entry

# define age at cohort entry as latest possible based on yob/mob
length(which(is.na(patients_first_rx$yob)))

cohort <- patients_first_rx %>%
  mutate (birthdate = if_else(is.na(mob), 
                              ymd(paste(yob, '12', '31' , sep = '')), 
                              ceiling_date(ym(paste(yob, mob, sep = '')), "month") - 1),
          age_at_entry = time_length(difftime(entry_date, birthdate), 'years'))

length(unique(cohort$id))
length(cohort$id)

# exclude patients <18
cat('Number of patients <18 at entry:', length(which(cohort$age_at_entry<18)))
cat(paste('Number of patients <18 at entry:', length(which(cohort$age_at_entry<18)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients >=18 at cohort entry:', length(which(cohort$age_at_entry>=18)))
cat(paste('Number of patients >=18 at cohort entry:', length(which(cohort$age_at_entry>=18)), '\n'), file = cohort_creation_desc, append = TRUE)

cohort %<>% filter (age_at_entry>=18)
length(cohort$id)

## 2. Apply exclusion criteria (2): >365 days of look-back available from cohort entry

# exclude patients with less than 1 year medical history
cat('Number of patients with <365 days lookback (regstart):', length(which(cohort$entry_date - cohort$regstart < 365)))
cat(paste('Number of patients with <365 days lookback (regstart):', length(which(cohort$entry_date - cohort$regstart < 365)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients with >= 365 days lookback (regstart):', length(which(cohort$entry_date - cohort$regstart >= 365)))
cat(paste('Number of patients with >= 365 days lookback (regstart):', length(which(cohort$entry_date - cohort$regstart >= 365)), '\n'), file = cohort_creation_desc, append = TRUE)

cohort %<>% filter (cohort$entry_date - cohort$regstart >= 365)
length(unique(cohort$id))

## 3. Apply exclusion criteria (3): >365 days of follow-up available since cohort entry
## based on last collection date of practice

cat('Number of patients with <365 days follow-up (lcd):', length(which(cohort$lcd - cohort$entry_date < 365)))
cat(paste('Number of patients with <365 days follow-up (lcd):', length(which(cohort$lcd - cohort$entry_date < 365)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients with >=365 days follow-up (lcd):', length(which(is.na(cohort$lcd) | cohort$lcd - cohort$entry_date >= 365)))
cat(paste('Number of patients with >=365 days follow-up (lcd):', length(which(is.na(cohort$lcd) | cohort$lcd - cohort$entry_date >= 365)), '\n'), file = cohort_creation_desc, append = TRUE)

cohort %<>% filter (is.na(lcd) | lcd - entry_date >= 365)

## Investigate death within one year of entry by regimen duration

# some patients' regimen ends when they die
# so defining follow-up based on regimen end could exclude such patients
# and introduce survival bias

cohort$dod <- pmin(cohort$emis_dod, cohort$cprd_dod)

cohort %<>%
  mutate (death_within_365 = if_else(!is.na(dod) & entry_date - dod < 365, 1, 0))

table(cohort$death_within_365)

cat('Number of patients with <365 days follow-up who die prior to 365 days (regend):', length(which(cohort$death_within_365 == 1 & cohort$regend - cohort$entry_date < 365)))
cat(paste('Number of patients with <365 days follow-up who die prior to 365 days (regend):', length(which(cohort$death_within_365 == 1 & cohort$regend - cohort$entry_date < 365)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients with <365 days follow-up who did not die prior to 365 days (regend):', length(which(cohort$death_within_365 == 0 & cohort$regend - cohort$entry_date < 365)))
cat(paste('Number of patients with <365 days follow-up who did not die prior to 365 days (regend):', length(which(cohort$death_within_365 == 0 & cohort$regend - cohort$entry_date < 365)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients with >=365 days follow-up who die prior to 365 days (regend):', length(which(cohort$death_within_365 == 1 & (is.na(cohort$regend) | cohort$regend - cohort$entry_date >= 365))))
cat(paste('Number of patients with >=365 days follow-up who die prior to 365 days (regend):', length(which(cohort$death_within_365 == 1 & (is.na(cohort$regend) | cohort$regend - cohort$entry_date >= 365))), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients with >=365 days follow-up who did not die prior to 365 days (regend):', length(which(cohort$death_within_365 == 0 & (is.na(cohort$regend) | cohort$regend - cohort$entry_date >= 365))))
cat(paste('Number of patients with >=365 days follow-up who did not die prior to 365 days (regend):', length(which(cohort$death_within_365 == 0 & (is.na(cohort$regend) | cohort$regend - cohort$entry_date >= 365))), '\n'), file = cohort_creation_desc, append = TRUE)

cohort %<>% select(-dod)

## 4. Apply exclusion criteria (4): linkage to HES and ONS

linkage_eligible_ids <- linkage_eligibility %>% 
  filter(ons_death_e == 1 &
           hes_apc_e == 1 &
           lsoa_e ==1)

ons_e_ids <- linkage_eligibility %>% filter(ons_death_e == 1)
hes_e_ids <- linkage_eligibility %>% filter(hes_apc_e == 1)
lsoa_e_ids <- linkage_eligibility %>% filter(lsoa_e == 1)

cat('Number of patients not linkable:', length(which(!cohort$id %in% linkage_eligible_ids$patid)))
cat(paste('Number of patients not linkable:', length(which(!cohort$id %in% linkage_eligible_ids$patid)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients not linkable to ONS:', length(which(!cohort$id %in% ons_e_ids$patid)))
cat(paste('Number of patients not linkable to ONS:', length(which(!cohort$id %in% ons_e_ids$patid)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients not linkable to HES:', length(which(!cohort$id %in% hes_e_ids$patid)))
cat(paste('Number of patients not linkable to HES:', length(which(!cohort$id %in% hes_e_ids$patid)), '\n'), file = cohort_creation_desc, append = TRUE)

cat('Number of patients not linkable to small area data:', length(which(!cohort$id %in% lsoa_e_ids$patid)))
cat(paste('Number of patients not linkable to small area data:', length(which(!cohort$id %in% lsoa_e_ids$patid)), '\n'), file = cohort_creation_desc, append = TRUE)


cohort %<>% filter (id %in% linkage_eligible_ids$patid)

length(cohort$id)
length(unique(cohort$id))

cat('Number of patients meeting inclusion and exclusion criteria:', length(unique(cohort$id)))
cat(paste('Number of patients meeting inclusion and exclusion criteria:', length(unique(cohort$id)), '\n'), file = cohort_creation_desc, append = TRUE)

# for clarity rename 'male' as 'sex'
cohort %<>%
  mutate (sex = as.factor(if_else(male == 0, 'Female', 'Male'))) %>% 
  select(-male)

rm(first_rx, first_rx_study, patients_first_rx)
rm(patients1, patients2, patients3, patients, practice1, practice2, practice3, practice)
rm(exposure_files, filter_prescriptions, exposure)
rm(therapy_filesA, therapy_filesB, therapy_filesA1, therapy_filesA2, therapy_filesB1, therapy_filesB2, therapy_filesC)

table(cohort$trt)

# save cohort
saveRDS(cohort, file = paste(path_cohort, 'antidepressant_cohort.rds', sep='/'))
saveRDS(switched_to, file = paste(path_cohort, 'switched_to.rds', sep = '/'))
