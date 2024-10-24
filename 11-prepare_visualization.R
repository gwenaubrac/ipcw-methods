## ---------------------------
##
## Program: 11. Prepare visualization
##
## Purpose: Prepare the data for visualization. 
##
## Author: Gwen Aubrac
##
## Date Created: 2024/10/15
##
## ---------------------------
##
## Notes: Only have data for CPRD region. Will add results from other databases in the future. 
##
## ---------------------------

#### SPECIFY ANALYSIS ####

# cohort: antidepressant, antihypertensive, antidiabetic
exposure <- 'antihypertensive'

# outcome: all-cause mortality, suicidal ideation
outcome <- 'all-cause mortality'

# analysis: main, flexible_grace_period, 90_day_grace_period, male, female
# young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed
analysis <- 'main'

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(ggplot2)
library(writexl)
library(lubridate)
library(tidyr)
library(broom)

options(scipen = 999)

#### DEFINE PATHS ####

path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, 'main', 'intermediate', sep = '/')
path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_final_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'final', sep = '/')
path_vis <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'visualization', sep = '/')

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))
covariates <- readRDS(file = paste(path_intermediate_res_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_intermediate_res_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_intermediate_res_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_intermediate_res_main, 'dec_comorb.rds', sep = '/'))

# prepare data

cohort_vis <- cohort
cohort_vis$region <- 'CPRD'

#### All X BY YEAR: NUMBER OF PATIENTS ####

# prop_overall: no grouping
# prop_region: grouped by region (here, only CPRD)
# prop_region_year: grouped by region and year
# prop_region_trt: grouped by region and treatment

allxbyyear <- cohort_vis %>% 
  group_by(region, year, trt) %>% 
  summarize(
    num_patients = n()
  )

total_patients <- nrow(cohort_vis)

allxbyyear %<>%
  mutate(prop_overall = num_patients / total_patients) %>% 
  group_by(region) %>% 
  mutate(prop_region = num_patients / total_patients) %>%
  ungroup() %>% 
  group_by(region, year) %>%
  mutate(prop_region_year = num_patients / sum(num_patients)) %>%
  ungroup() %>%
  group_by(region, trt) %>%
  mutate(prop_region_trt = num_patients / sum(num_patients))

saveRDS(allxbyyear, file = paste(path_vis, 'allxbyyear.xlsx', sep = '/'))

#### ALL Y BY YEAR ####

cohort_vis <- cohort
cohort_vis$region <- 'CPRD'

# Get number of deaths per region and year

deaths_per_yr <- cohort_vis %>% 
  filter(itt_event == 1) %>% 
  mutate(event_yr = year(itt_event_date)) %>% 
  group_by(region, event_yr) %>% 
  summarize(num_deaths = n())

# Get days of follow-up contributed by patient in each year

get_yearly_followup <- function(year, year_start, year_end) {
  
  year_col <- paste('fu_days_', year, sep = '')
  cohort_vis[[year_col]] <- NA
  
  cohort_vis[[year_col]] <- case_when(
    cohort_vis$entry_date > year_end | cohort_vis$itt_exit_date <= year_start ~ 0,
    cohort_vis$entry_date <= year_start & cohort_vis$itt_exit_date > year_end ~ 365,
    cohort_vis$entry_date <= year_start & cohort_vis$itt_exit_date <= year_end ~ as.numeric(cohort_vis$itt_exit_date - year_start),
    cohort_vis$entry_date >= year_start & cohort_vis$itt_exit_date > year_end ~ as.numeric(year_end - cohort_vis$entry_date),
    cohort_vis$entry_date >= year_start & cohort_vis$itt_exit_date <= year_end ~ as.numeric(cohort_vis$itt_exit_date - cohort_vis$entry_date),
    TRUE ~ cohort_vis[[year_col]]
  )
  
  cohort_vis <<- cohort_vis
  
}

get_yearly_followup('2019', year_start = ymd(20190101), year_end = ymd(20191231))
get_yearly_followup('2020', year_start = ymd(20200101), year_end = ymd(20201231))
get_yearly_followup('2021', year_start = ymd(20210101), year_end = ymd(20211231))
get_yearly_followup('2022', year_start = ymd(20220101), year_end = ymd(20221231))
get_yearly_followup('2023', year_start = ymd(20230101), year_end = ymd(20231231))
get_yearly_followup('2024', year_start = ymd(20240101), year_end = ymd(20240331)) # study ended in March 2024

# Get days and years of follow up for each study year

followup_per_yr <- lapply(2019:2024, function(year) {
  cohort_vis %>%
    group_by(region) %>%
    summarize(
      pdays = sum(.data[[paste0("fu_days_", year)]]),  # Use .data for dynamic column access
      .groups = 'drop'
    ) %>%
    mutate(event_yr = as.character(year),
           pyears = pdays / 365)
})

followup_per_yr <- bind_rows(followup_per_yr)
  
# Create allybyyear dataframe

allybyyear <- merge(deaths_per_yr, followup_per_yr, by = c('region', 'event_yr')) %>% 
  mutate(IRper100 = num_deaths*100/pyears)
saveRDS(allybyyear, file = paste(path_vis, 'allybyyear.xlsx', sep = '/')) 

#### COVARIATES BY REGION ####

base_comorb <- base_comorb[!base_comorb %in% c('anxiety_base')]

# overall
covsbyreg <- cohort_vis %>%
  summarise(across(all_of(base_comorb), ~ sum(. == 1, na.rm = TRUE)))

covsbyreg <- covsbyreg %>% 
  pivot_longer(cols = -c(), names_to = "comorbidity", values_to = "count")

# for exposed
covsbyreg_exp1 <- cohort_vis %>%
  filter(trt_dummy == 1) %>% 
  summarise(across(all_of(base_comorb), ~ sum(. == 1, na.rm = TRUE)))

covsbyreg_exp1 <- covsbyreg_exp1 %>% 
  pivot_longer(cols = -c(), names_to = "comorbidity", values_to = "count_exp1")

# for unexposed
covsbyreg_exp0 <- cohort_vis %>%
  filter(trt_dummy == 0) %>% 
  summarise(across(all_of(base_comorb), ~ sum(. == 1, na.rm = TRUE)))

covsbyreg_exp0 <- covsbyreg_exp0 %>% 
  pivot_longer(cols = -c(), names_to = "comorbidity", values_to = "count_exp0")

# final covsbyreg
covsbyreg <- merge(covsbyreg, covsbyreg_exp1, by = 'comorbidity')
covsbyreg <- merge(covsbyreg, covsbyreg_exp0, by = 'comorbidity')
covsbyreg <- cbind(region = 'CPRD', covsbyreg, row.names = NULL)

covsbyreg %<>%
  group_by(region) %>% 
  mutate(comorb_freq = count / sum(count),
         comorb_freq_exp1 = count_exp1 / sum(count_exp1),
         comorb_freq_exp0 = count_exp0 / sum(count_exp0))

saveRDS(covscoef, file = paste(path_vis, 'covsbyreg.xlsx', sep = '/'))

#### PS COVARIATE COEFFICIENTS ####

ps_formula <- readRDS(file = paste(path_final_res, 'iptw_model.rds', sep = '/'))

ps_model <- glm(ps_formula, family = binomial(link = 'logit'), data = cohort)

ps_parameter_estimates <- tidy(ps_model)

ps_parameter_estimates$region <- 'CPRD'

ps_parameter_estimates %<>%
  filter(term != '(Intercept)') %>% 
  mutate(OR = exp(estimate)) %>% 
  rename(parameter = term) %>% 
  select(parameter, region, estimate, OR)

saveRDS(ps_parameter_estimates, file = paste(path_vis, 'pscofs.xlsx', sep = '/'))

#### MARGINAL BIAS TERMS ####

# can replace base_comorb by any list of variable names
# for binary variables for which want to calculate bias

marg_bias_results <- data.frame()

base_comorb <- base_comorb[!base_comorb %in% c('anxiety_base')]

for (cov in base_comorb) {
  
  # get proportioned of patients with covariate in exposed and unexposed groups
  prop_exposed <- cohort_vis %>%
    group_by(region, !!sym(cov), trt_dummy) %>%
    summarize(num_patients = n())
  
  prop_exposed %<>%
    group_by(trt_dummy) %>% 
    mutate(num_patients_trt = sum(num_patients)) %>% 
    ungroup () %>% 
    mutate(pct_col = num_patients / num_patients_trt)
  
  prop_exposed %<>%
    filter(!!sym(cov) == 1) %>%
    select(region, trt_dummy, pct_col) %>%
    pivot_wider(names_from = trt_dummy, values_from = pct_col, names_prefix = "Exp_")
  
  # get probability of death within 1 year associated with covariate
  model <- glm(reformulate(cov, 'death_within_365'),  data = cohort_vis, family = binomial(link = 'logit'))
  
  risk_ratio <- tidy(model)
  
  risk_ratio$region <- 'CPRD'
  
  risk_ratio %<>%
    filter(term != '(Intercept)') %>% 
    mutate(RR = exp(estimate)) %>% 
    rename(parameter = term) %>% 
    select(parameter, region, estimate, RR)
  
  # get final bias term dataframe by region
  merged <- merge(prop_exposed, risk_ratio, by = 'region')
  
  merged %<>%
    mutate(bias = if_else(
      estimate > 0,
      (Exp_1 * (RR - 1) + 1) / (Exp_0 * (RR - 1) + 1),
      (Exp_1 * (1 / RR - 1) + 1) / (Exp_0 * (1 / RR - 1) + 1)
    ),
    AbsBias = abs(log(bias)))
  
  marg_bias_results <- bind_rows(marg_bias_results, merged)
  
}

saveRDS(marg_bias_results, file = paste(path_vis, 'margbiasterms.xlsx', sep = '/'))

#### COX MODEL RESULTS ####

cox_itt <- readRDS(paste(path_final_res, 'cox_itt.rds', sep = '/'))
cox_itt_iptw <- readRDS(paste(path_final_res, 'cox_itt_iptw.rds', sep = '/'))
cox_at <- readRDS(paste(path_final_res, 'cox_at.rds', sep = '/'))
cox_at_iptw <- readRDS(paste(path_final_res, 'cox_at_iptw.rds', sep = '/'))
cox_at_ipcw_str_lag <- readRDS(paste(path_final_res, 'cox_at_ipcw_str_lag.rds', sep = '/'))
cox_at_iptw_ipcw_str_lag <- readRDS(paste(path_final_res, 'cox_at_iptw_ipcw_str_lag.rds', sep = '/'))
cox_at_ipcw_str_nonlag <- readRDS(paste(path_final_res, 'cox_at_ipcw_str_nonlag.rds', sep = '/'))
cox_at_iptw_ipcw_str_nonlag <- readRDS(paste(path_final_res, 'cox_at_iptw_ipcw_str_nonlag.rds', sep = '/'))
cox_at_ipcw_pl_nonlag <- readRDS(paste(path_final_res, 'cox_at_ipcw_pl_nonlag.rds', sep = '/'))
cox_at_iptw_ipcw_pl_nonlag <- readRDS(paste(path_final_res, 'cox_at_iptw_ipcw_pl_nonlag.rds', sep = '/'))
cox_at_ipcw_pl_lag <- readRDS(paste(path_final_res, 'cox_at_ipcw_pl_lag.rds', sep = '/'))
cox_at_iptw_ipcw_pl_lag <- readRDS(paste(path_final_res, 'cox_at_iptw_ipcw_pl_lag.rds', sep = '/'))

result_chart <- data.frame(matrix(nrow = 12, ncol = 3))
colnames(result_chart) <- c('estimate', 'lower_ci', 'upper_ci')
rownames(result_chart) <- c('ITT', 'ITT (IPTW)', 'AT', 'AT (IPTW)', 'AT (stratified lagged sIPCW)', 'AT (IPTW + stratified lagged IPCW)',
                            'AT (stratified non-lagged IPCW)', 'AT (IPTW + stratified non-lagged IPCW)', 'AT (pooled non-lagged IPCW)',
                            'AT (IPTW + pooled non-lagged IPCW)', 'AT (pooled lagged IPCW)', 'AT (IPTW + pooled lagged IPCW)')

result_chart[1, 'estimate'] <- exp(cox_itt$coef)
result_chart[1, 'lower_ci'] <- exp(confint(cox_itt))[1]
result_chart[1, 'upper_ci'] <- exp(confint(cox_itt))[2]

result_chart[2, 'estimate'] <- exp(cox_itt_iptw$coef)
result_chart[2, 'lower_ci'] <- exp(confint(cox_itt_iptw))[1]
result_chart[2, 'upper_ci'] <- exp(confint(cox_itt_iptw))[2]

result_chart[3, 'estimate'] <- exp(cox_at$coef)
result_chart[3, 'lower_ci'] <- exp(confint(cox_at))[1]
result_chart[3, 'upper_ci'] <- exp(confint(cox_at))[2]

result_chart[4, 'estimate'] <- exp(cox_at_siptw$coef)
result_chart[4, 'lower_ci'] <- exp(confint(cox_at_iptw))[1]
result_chart[4, 'upper_ci'] <- exp(confint(cox_at_iptw))[2]

result_chart[5, 'estimate'] <- exp(cox_at_ipcw_str_lag$coef)
result_chart[5, 'lower_ci'] <- exp(confint(cox_at_ipcw_str_lag))[1]
result_chart[5, 'upper_ci'] <- exp(confint(cox_at_ipcw_str_lag))[2]

result_chart[6, 'estimate'] <- exp(cox_at_iptw_ipcw_str_lag$coef)
result_chart[6, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_str_lag))[1]
result_chart[6, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_str_lag))[2]

result_chart[7, 'estimate'] <- exp(cox_at_ipcw_str_nonlag$coef)
result_chart[7, 'lower_ci'] <- exp(confint(cox_at_ipcw_str_nonlag))[1]
result_chart[7, 'upper_ci'] <- exp(confint(cox_at_ipcw_str_nonlag))[2]

result_chart[8, 'estimate'] <- exp(cox_at_iptw_ipcw_str_nonlag$coef)
result_chart[8, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_str_nonlag))[1]
result_chart[8, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_str_nonlag))[2]

result_chart[9, 'estimate'] <- exp(cox_at_ipcw_pl_nonlag$coef)
result_chart[9, 'lower_ci'] <- exp(confint(cox_at_ipcw_pl_nonlag))[1]
result_chart[9, 'upper_ci'] <- exp(confint(cox_at_ipcw_pl_nonlag))[2]

result_chart[10, 'estimate'] <- exp(cox_at_iptw_ipcw_pl_nonlag$coef)
result_chart[10, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_nonlag))[1]
result_chart[10, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_nonlag))[2]

result_chart[11, 'estimate'] <- exp(cox_at_ipcw_pl_lag$coef)
result_chart[11, 'lower_ci'] <- exp(confint(cox_at_ipcw_pl_lag))[1]
result_chart[11, 'upper_ci'] <- exp(confint(cox_at_ipcw_pl_lag))[2]

result_chart[12, 'estimate'] <- exp(cox_at_iptw_ipcw_pl_lag$coef)
result_chart[12, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_lag))[1]
result_chart[12, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_lag))[2]

result_chart <- cbind(model = row.names(result_chart), result_chart, row.names = NULL)
result_chart <- result_chart[c(2,4,6,8,10,12),]

saveRDS(result_chart, file = paste(path_vis, 'cox_results.xlsx', sep = '/'))

#### INCIDENCE RATE RATIO RESULTS ####

bootstrap_ci <- readRDS(paste(path_final_res, 'bootstrap_ci.rds', sep = '/'))

time <- data.frame(time = c(unique(cohort$month_year))) %>% 
  arrange(time)

# ITT IPTW: itt.2019-01.iptw.IR
ir_itt_iptw <- bootstrap_ci %>%
  filter(grepl("^ir\\.itt\\.iptw\\.2.*\\.IRR$", variable)) %>% 
  select(-variable)

ir_itt_iptw <- ir_itt_iptw[-1,]
ir_itt_iptw <- cbind(time, variable = 'ITT (IPTW)', ir_itt_iptw)

# AT IPTW: att.2019-01.iptw.IR
at_iptw_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_IR <- at_iptw_IR[-(1:4),]
at_iptw_IR <- cbind(time, variable = 'AT (IPTW)', at_iptw_IR)

# AT IPTW + stratified lagged IPCW: att.2019-01.iptw.ipcw.IR
at_iptw_ipcw_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_IR <- cbind(time, variable = 'AT (stratified lagged IPCW + IPTW)', at_iptw_ipcw_IR)

# AT IPTW + stratified non-lagged IPCW: att.2019-01.iptw.ipcw_nl.IR
at_iptw_ipcw_nl_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw_nl\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_nl_IR <- cbind(time, variable = 'AT (stratified non-lagged IPCW + IPTW)', at_iptw_ipcw_nl_IR)

# AT IPTW + pooled non-lagged IPCW: att.2019-01.iptw.ipcw_pl_nonlag.IR
at_iptw_ipcw_pl_nonlag_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw_pl_nonlag\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_pl_nonlag_IR <- cbind(time, variable = 'AT (pooled non-lagged IPCW + IPTW)', at_iptw_ipcw_pl_nonlag_IR)

# AT IPTW + pooled lagged IPCW: att.2019-01.iptw.ipcw_pl_lag.IR
at_iptw_ipcw_pl_lag_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw_pl_lag\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_pl_lag_IR <- cbind(time, variable = 'AT (pooled lagged IPCW + IPTW)', at_iptw_ipcw_pl_lag_IR)
IR_results <- rbind(itt_iptw_IR, at_iptw_IR, at_iptw_ipcw_IR, at_iptw_ipcw_nl_IR, at_iptw_ipcw_pl_nonlag_IR, at_iptw_ipcw_pl_lag_IR)

