## ---------------------------
##
## Program: 3. Define covariates
##
## Purpose: Update the cohort with values of covariates and comorbidities at study entry 
## and at the quartiles of the censoring distribution ('cohort_with_covariates'). 
## Censoring for the distribution of censoring times can be for any reason, 
# i.e. event, death, switch, discontinuation, linkage end, departure from CPRD...
## which is captured in at_exit_date. 
##
## Comorbidities are assessed from both CPRD and HES. 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes: 
## 1. Baseline and quartile dates are relative to each patient (and not to the study timeline).
## Quartiles are based on the number of days of follow-up from the overall follow-up distribution.
## Relative dates are calculated based on each patient's entry date.
##
## 2. The folder for comorbidity codes should contain separate excel files for each comorbidity. 
## The comorbidity will be named according to the excel file name, and its flag value at each quartile
## will be named 'Q1_comorb', 'Q2_comorb' ... with the respective 'comorb' name. 
##
## 3.The occurrence of comorbidities is assessed based on CPRD and HES. The date of first
## occurrence of the comorbidity is used to define the presence of the comorbidity in the cohort. 
##
## ---------------------------

# analysis: main, flex_grace_period, or 90_day_grace_period

analysis <- ''

#### LOAD PACKAGES ####

library(tidyr)
library(lubridate)
library(readxl)
library(dplyr)
library(magrittr)
library(haven) 
library(parallel)
library(data.table)

#### DEFINE PATHS ####

if (analysis == 'main' | analysis == '') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
} else if (analysis == 'flex_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/flex_grace_period" 
} else if (analysis == '90_day_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/90_day_grace_period" 
} 

path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
path_cprdA <- "Z:/EPI/Protocol 24_004042/dataA" 
path_cprdB <- "Z:/EPI/Protocol 24_004042/dataB" 
path_cprdC <- "Z:/EPI/Protocol 24_004042/dataC (no followup)"
path_comorb_cprd <- "Z:/EPI/Protocol 24_004042/Gwen/data/comorbidities/Aurum codes comorbidities"
path_comorb_hes <- "Z:/EPI/Protocol 24_004042/Gwen/data/comorbidities/ICD codes comorbidities"
path_linkage_1 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt1'
path_linkage_2 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt2'

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_censored.rds', sep = '/'))

setwd(path_cohort)

cov_desc <- "cov_desc.txt"
writeLines("Covariate description:", cov_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = cov_desc, append= TRUE)

#### YEAR OF ENTRY, AGE GROUP, ETHNICITY, AND DEPRIVATION ####

# calendar month and year of cohort entry
cohort <- transform(cohort,
                    month = factor(format(entry_date, format = '%B'), levels = month.name),
                    year = format(entry_date, format = '%Y'))

cohort$month_year <- format(as.Date(cohort$entry_date), '%Y-%m')
str(cohort$year)
str(cohort$month_year)

cohort %<>%
  mutate(year = as.factor(year),
         month_year = as.factor(month_year))

# age group
cohort <- cohort %>% 
  mutate (
    age_group = dplyr::case_when(
      age_at_entry >= 18 & age_at_entry < 25 ~ '18-24',
      age_at_entry >= 25 & age_at_entry < 35 ~ '25-34',
      age_at_entry >= 35 & age_at_entry < 45 ~ '35-44',
      age_at_entry >= 45 & age_at_entry < 55 ~ '45-54',
      age_at_entry >= 55 & age_at_entry < 65 ~ '55-64',
      age_at_entry >= 65 & age_at_entry < 75 ~ '65-74',
      age_at_entry >= 75 & age_at_entry < 85 ~ '75-84',
      age_at_entry > 85 ~ '>85'
    ),
    age_group = factor(age_group, levels = c('18-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '>85'))
  )

summary(cohort$age_group)
length(which(is.na(cohort$age_group)))

# ethnicity
col_classes <- c("character", "character", "character", "character")
ethnicity1 <- fread(paste(path_linkage_1, 'hes_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
ethnicity2 <- fread(paste(path_linkage_2, 'hes_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
ethnicity <- bind_rows(ethnicity1, ethnicity2) %>% 
  select(patid, gen_ethnicity) %>% 
  rename (ethnicity = gen_ethnicity)

rm(ethnicity1, ethnicity2)

length(which(cohort$id %in% ethnicity$patid))
cat(paste('Patients with ethnicity information:', length(which(cohort$id %in% ethnicity$patid)), '\n'), file = cov_desc, append= TRUE)

cohort <- merge(cohort, ethnicity, by.x = 'id', by.y = 'patid', all.x = TRUE)

cohort %<>% mutate(ethnicity = as.factor(ethnicity))
summary(cohort$ethnicity)

cohort %<>% mutate(ethnicity = if_else(ethnicity == '', NA, ethnicity)) # replace blanks by NA
summary(cohort$ethnicity)

# deprivation index
col_classes <- c("character", "character", "character")
deprivation1 <- fread(paste(path_linkage_1, 'patient_2019_imd_24_004042.txt', sep = '/'), colClasses = col_classes)
deprivation2 <- fread(paste(path_linkage_2, 'patient_2019_imd_24_004042.txt', sep = '/'), colClasses = col_classes)
deprivation <- bind_rows(deprivation1, deprivation2) %>% 
  select(patid, e2019_imd_5) %>% 
  rename (deprivation = e2019_imd_5)

rm(deprivation1, deprivation2)

length(which(cohort$id %in% deprivation$patid))
cat(paste('Patients with deprivation information:', length(which(cohort$id %in% deprivation$patid)), '\n'), file = cov_desc, append= TRUE)

cohort <- merge(cohort, deprivation, by.x = 'id', by.y = 'patid', all.x = TRUE)

cohort %<>% mutate(deprivation = as.factor(deprivation))
summary(cohort$deprivation)

cohort %<>% mutate(deprivation = if_else(deprivation == '', NA, deprivation)) # replace blanks by NA
summary(cohort$deprivation)

## CPRD COMORBIDITIES

# read comorbidity codes from CPRD
comorb_cprd <- data.frame()
comorb_cprd_files <- list.files(path_comorb_cprd, pattern = '.xlsx', all.files = TRUE, full.names = TRUE)

for (file in comorb_cprd_files) {
  data <- read_excel(file, trim_ws = TRUE, col_type = 'text')
  data$comorb <- tools::file_path_sans_ext(basename(file))
  comorb_cprd <- bind_rows(comorb_cprd, data)
}

comorb_cprd %<>%
  select(MedCodeId, Term, comorb)

rm(data, file)

## CPRD: Identify comorbidities

cohort_ids <- as.list(cohort$id)

filter_comorb <- function (file_path) {
  file <- read_sas(file_path)
  file_filtered <- file %>%
    dplyr::select(id, medical_code, date) %>% 
    filter(id %in% cohort_ids & medical_code %in% comorb_cprd$MedCodeId) %>% 
    arrange (id, date)
  cohort_comorb_cprd <<- bind_rows(cohort_comorb_cprd, file_filtered)
}

observation_filesA <- list.files(
  path_cprdA,
  pattern = 'observation',
  all.files = TRUE,
  full.names = TRUE
)

observation_filesA1 <- observation_filesA[1:60]
observation_filesA2 <- observation_filesA[61:113]
observation_filesA3 <- observation_filesA[114:189]

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesA1, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort,  'cprd_comorb_A1.rds', sep='/'))
rm(cohort_comorb_cprd)

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesA2, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort,  'cprd_comorb_A2.rds', sep='/'))
rm(cohort_comorb_cprd)

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesA3, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort,  'cprd_comorb_A3.rds', sep='/'))
rm(cohort_comorb_cprd)

observation_filesB <- list.files(
  path_cprdB,
  pattern = 'observation',
  all.files = TRUE,
  full.names = TRUE
)

observation_filesB1 <- observation_filesB[1:60]
observation_filesB2 <- observation_filesB[61:113]
observation_filesB3 <- observation_filesB[114:197]

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesB1, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort, 'cprd_comorb_B1.rds', sep='/'))
rm(cohort_comorb_cprd)

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesB2, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort, 'cprd_comorb_B2.rds', sep='/'))
rm(cohort_comorb_cprd)

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesB3, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort, 'cprd_comorb_B3.rds', sep='/'))
rm(cohort_comorb_cprd)

observation_filesC <- list.files(
  path_cprdC,
  pattern = 'observation',
  all.files = TRUE,
  full.names = TRUE
)

gc()
cohort_comorb_cprd <- data.frame()
mclapply(observation_filesC, filter_comorb)
saveRDS(cohort_comorb_cprd, file = paste(path_cohort, 'cprd_comorb_C.rds', sep='/'))

rm(cohort_comorb_cprd)

cprd_comorb_A1 <- readRDS(file = paste(path_main, 'cprd_comorb_A1.rds', sep = '/'))
cprd_comorb_A2 <- readRDS(file = paste(path_main, 'cprd_comorb_A2.rds', sep = '/'))
cprd_comorb_A3 <- readRDS(file = paste(path_main, 'cprd_comorb_A3.rds', sep = '/'))
cprd_comorb_B1 <- readRDS(file = paste(path_main, 'cprd_comorb_B1.rds', sep = '/'))
cprd_comorb_B2 <- readRDS(file = paste(path_main, 'cprd_comorb_B2.rds', sep = '/'))
cprd_comorb_B3 <- readRDS(file = paste(path_main, 'cprd_comorb_B3.rds', sep = '/'))
cprd_comorb_C <- readRDS(file = paste(path_main, 'cprd_comorb_C.rds', sep = '/'))

cohort_comorb_cprd <- rbind(cprd_comorb_A1, cprd_comorb_A2, cprd_comorb_A3, cprd_comorb_B1, cprd_comorb_B2, cprd_comorb_B3, cprd_comorb_C)
cohort_comorb_cprd <- merge(cohort_comorb_cprd, comorb_cprd, by.x = 'medical_code', by.y = 'MedCodeId')

rm(observation_filesA, observation_filesB, observation_filesA1, observation_filesA2, observation_filesA3, 
   observation_filesB1, observation_filesB2, observation_filesB3, observation_filesC)

rm(cprd_comorb_A1, cprd_comorb_A2, cprd_comorb_A3, cprd_comorb_B1, cprd_comorb_B2, cprd_comorb_B3, cprd_comorb_C)

## HES comorbidities

# read comorbidity codes from hes
comorb_hes <- data.frame()

comorb_hes_files <- list.files(path_comorb_hes, pattern = '.xlsx', all.files = TRUE, full.names = TRUE)

for (file in comorb_hes_files) {
  data <- read_excel(file, trim_ws = TRUE, col_type = 'text')
  data$comorb <- tools::file_path_sans_ext(basename(file))
  comorb_hes <- bind_rows(comorb_hes, data)
}

rm(data, file)

# formatting
comorb_hes %<>% filter(!is.na(ICD10Code))
comorb_hes$comorb <- gsub('^icd-', '', comorb_hes$comorb)
comorb_hes$ICD10Code <- gsub('\\.x$', '', comorb_hes$ICD10Code)
comorb_hes$ICD10Code <- gsub('\\s+', '', comorb_hes$ICD10Code)

# split by ICD code length
summary(nchar(comorb_hes$ICD10Code))
comorb_hes_three <- comorb_hes %>% filter(nchar(ICD10Code) == 3)
comorb_hes_five <- comorb_hes %>% filter(nchar(ICD10Code) == 5)
comorb_hes_six <- comorb_hes %>% filter(nchar(ICD10Code) == 6)
comorb_hes_seven <- comorb_hes %>% filter(nchar(ICD10Code) == 7)

# 1. comorb from episodes
col_classes <- c("character", "character", "character", "character", "character", "character", "character", "character")
hes_dx_epi1 <- fread(paste(path_linkage_1, 'hes_diagnosis_epi_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
hes_dx_epi2 <- fread(paste(path_linkage_2, 'hes_diagnosis_epi_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
hes_dx_epi <- bind_rows(hes_dx_epi1, hes_dx_epi2)

rm(hes_dx_epi1, hes_dx_epi2)
cohort_ids <- as.list(cohort$id)

length(which(cohort$id %in% hes_dx_epi$patid))
cat(paste('Patients with HES episode information:', length(which(cohort$id %in% hes_dx_epi$patid)), '\n'), file = cov_desc, append= TRUE)

# HES ICD contains up to 5 characters (including .)
# and ICDx contains 4th or 5th ICD character (excluding .), which is 5th or 6th character (including .)
# so we concatenate 5-character long ICD with ICDx if there is an ICDx

cohort_comorb_hes_epi <- hes_dx_epi %>%
  select(patid, ICD, ICDx, epistart) %>%
  rename(id = patid) %>% 
  mutate(
    ICD = if_else(
      !is.na(ICDx) & ICDx != '-' & nchar(ICD) == 5, 
      paste(ICD, ICDx, sep=''), ICD),
    icd_code_length = nchar(ICD)
    ) 

cohort_comorb_hes_epi <- cohort_comorb_hes_epi %>% 
  filter(id %in% cohort_ids &
           (
           (icd_code_length == 3 & ICD %in% comorb_hes_three$ICD10Code) |
           (icd_code_length == 5 & ICD %in% comorb_hes_five$ICD10Code) |
           (icd_code_length == 6 & ICD %in% comorb_hes_six$ICD10Code) |
           (icd_code_length == 7 & ICD %in% comorb_hes_seven$ICD10Code)
           )
  ) %>%
  mutate(epistart = dmy(epistart)) %>% 
  arrange(id, epistart) %>% 
  rename (date = epistart)

length(unique(cohort_comorb_hes_epi$id))

cohort_comorb_hes_epi <- merge(cohort_comorb_hes_epi, comorb_hes, by.x = 'ICD', by.y = 'ICD10Code')

# 2. comorb from hospitalisations
col_classes <- c("character", "character", "character", "character", "character", "character")
hes_dx_hosp1 <- fread(paste(path_linkage_1, 'hes_diagnosis_hosp_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
hes_dx_hosp2 <- fread(paste(path_linkage_2, 'hes_diagnosis_hosp_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
hes_dx_hosp <- bind_rows(hes_dx_hosp1, hes_dx_hosp2)

rm(hes_dx_hosp1, hes_dx_hosp2)

length(which(cohort$id %in% hes_dx_hosp$patid))
cat(paste('Patients with HES hospitalization information:', length(which(cohort$id %in% hes_dx_hosp$patid)), '\n'), file = cov_desc, append= TRUE)

cohort_comorb_hes_hosp <- hes_dx_hosp %>%
  select(patid, ICD, ICDx, admidate) %>%
  rename(id = patid) %>% 
  mutate(
    ICD = if_else(
      !is.na(ICDx) & ICDx != '-' & nchar(ICD) == 5, 
      paste(ICD, ICDx, sep=''), ICD),
    icd_code_length = nchar(ICD)
  )

cohort_comorb_hes_hosp <- cohort_comorb_hes_hosp %>% 
  filter(id %in% cohort_ids &
           (
             (icd_code_length == 3 & ICD %in% comorb_hes_three$ICD10Code) |
               (icd_code_length == 5 & ICD %in% comorb_hes_five$ICD10Code) |
               (icd_code_length == 6 & ICD %in% comorb_hes_six$ICD10Code) |
               (icd_code_length == 7 & ICD %in% comorb_hes_seven$ICD10Code)
           )
  ) %>%
  mutate(admidate = dmy(admidate)) %>% 
  arrange(id, admidate) %>% 
  rename (date = admidate)

length(unique(cohort_comorb_hes_hosp$id))

cohort_comorb_hes_hosp <- merge(cohort_comorb_hes_hosp, comorb_hes, by.x = 'ICD', by.y = 'ICD10Code')

## Combine CPRD & HES comorb and select first one

first_comorb <- bind_rows(cohort_comorb_cprd, cohort_comorb_hes_epi, cohort_comorb_hes_hosp)

first_comorb <- first_comorb %>% 
  group_by(id, comorb) %>%
  summarise(comorb_date = min(date)) %>% 
  ungroup()

saveRDS(first_comorb, file = paste(path_cohort, 'first_comorb.rds', sep = '/'))

rm(first_comorb, comorb_cprd, comorb_hes, cohort_comorb_cprd, cohort_comorb_hes_epi, cohort_comorb_hes_hosp, cohort_ids)

#### FORMAT COVARIATE DATA ####

#trt_dummy = 0 (reference exposure group, here SNRI) and trt_dummy = 1 (comparator exposure group, here SSRI)
cohort$trt_dummy <- if_else (cohort$trt == 'ssri', 1, 0)
cohort %<>% mutate(trt_dummy = as.factor(trt_dummy), trt = as.factor(trt))

# create vectors with names for all user-specified covariates and comorbidities
covariates <- c('age_group', 'sex', 'year', 'ethnicity', 'deprivation') 
cat(paste('The baseline covariates are:', covariates, '\n'), file = cov_desc, append = TRUE)

comorbidities <- vector (mode = 'character', length = length(comorb_cprd_files))
base_comorb <- vector (mode = 'character', length = length(comorb_cprd_files))
dec_comorb <- c()

for (i in 1:length(comorb_cprd_files)) {
  comorb_name <- tools::file_path_sans_ext(basename(comorb_cprd_files)[i])
  base_comorb[i] <- paste(comorb_name, 'base', sep = '_')
  comorbidities[i] <- comorb_name
  decile_comorb <- vector (mode = 'character', length = 9)
  for (j in 1:9) {
    decile_comorb[j] <- paste(comorb_name, '_d', j, sep='')
  }
  dec_comorb <- c(dec_comorb, decile_comorb)
}

base_comorb
dec_comorb

saveRDS(covariates, file = paste(path_cohort, 'covariates.rds', sep = '/'))
saveRDS(comorbidities, file = paste(path_cohort, 'comorbidities.rds', sep = '/'))
saveRDS(base_comorb, file = paste(path_cohort, 'base_comorb.rds', sep = '/'))
saveRDS(dec_comorb, file = paste(path_cohort, 'dec_comorb.rds', sep = '/'))

# remove empty factor levels
cohort$deprivation <- droplevels(cohort$deprivation)
cohort$ethnicity <- droplevels(cohort$ethnicity)

# check for patients with missing covariate info
length(which(is.na(cohort$sex)))
cat(paste('Patients missing sex (removed from analyses):', length(which(is.na(cohort$sex))), '\n'), file = cov_desc, append = TRUE)

cohort <- cohort %>%
  filter (!is.na(sex))

length(which(is.na(cohort$ethnicity)))
cat(paste('Patients missing ethnicity (recategorized as unknown):', length(which(is.na(cohort$ethnicity))), '\n'), file = cov_desc, append = TRUE)

cohort <- cohort %>%
  mutate (ethnicity = as.factor(if_else(is.na(ethnicity), 'Unknown', ethnicity)))
summary(cohort$ethnicity)

length(which(is.na(cohort$deprivation)))
cat(paste('Patients missing deprivation (recategorized as unknown):', length(which(is.na(cohort$deprivation))), '\n'), file = cov_desc, append = TRUE)

cohort <- cohort %>%
  mutate (deprivation = as.factor(if_else(is.na(deprivation), 'Unknown', deprivation)))
summary(cohort$deprivation)

# recategorize ethnicity because too few counts of certain groups
cohort %<>%
  mutate(ethnicity = if_else(ethnicity == 'Unknown', 'Unknown', if_else(ethnicity == 'White', 'White', 'Non-white')),
         ethnicity = as.factor(ethnicity))

saveRDS(cohort, file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep='/'))

