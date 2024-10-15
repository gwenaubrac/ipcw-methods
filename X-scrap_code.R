## ---------------------------
##
## Program: X. Scrap Code
##
## Purpose: Scrap code that I may or may not need in the future. 
##
## Author: Gwen Aubrac
##
## Date Created: 2024/10/09
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

#### GET PERSON TIME AND EVENTS FOR EACH ANALYSIS TYPE ####

cohort_analytic_at <- readRDS(file = paste(path_main, 'cohort_analytic_at.rds', sep = '/'))

cohort_analytic_itt <- cohort %>% 
  mutate (
    person_time = as.numeric(itt_exit_date - entry_date),
    iptw_person_time = person_time * iptw,
    siptw_person_time = person_time * siptw,
    iptw_event = itt_event * iptw,
    siptw_event = itt_event * siptw,
  )

itt_person_time <- cohort_analytic_itt %>% 
  mutate(itt_pt = sum(person_time),
         itt_pt_iptw = sum(iptw_person_time)) %>% 
  select(id, itt_pt, itt_pt_iptw) 

itt_events <- cohort_analytic_itt %>% 
  mutate(total_events = sum(itt_event),
         total_events_iptw = sum(iptw_event)) %>% 
  select(year, total_events, total_events_iptw) 

at_person_time <- cohort_analytic_at %>% 
  group_by(id) %>%
  mutate(
    at_pt = sum(person_time),
    at_pt_iptw_ipcw_str_lag = sum(iptw_ipcw_str_lag_person_time),
    at_pt_iptw_ipcw_str_nonlag = sum(iptw_ipcw_str_nonlag_person_time),
    at_pt_iptw_ipcw_pl_nonlag = sum(iptw_ipcw_pl_nonlag_person_time),
    at_pt_iptw_ipcw_pl_lag = sum(iptw_ipcw_pl_lag_person_time)
  ) %>%
  select(
    id,
    at_pt,
    at_pt_iptw_ipcw_str_lag,
    at_pt_iptw_ipcw_str_nonlag,
    at_pt_iptw_ipcw_pl_nonlag,
    at_pt_iptw_ipcw_pl_lag
  ) %>%
  slice(1)

at_events <- cohort_analytic_at %>% 
  group_by(year) %>%
  mutate(
    total_events = sum(event_at_tstop),
    total_events_iptw_ipcw_str_lag = sum(iptw_ipcw_str_lag_event),
    total_events_iptw_ipcw_str_nonlag = sum(iptw_ipcw_str_nonlag_event),
    total_events_iptw_ipcw_pl_nonlag = sum(iptw_ipcw_pl_nonlag_event),
    total_events_iptw_ipcw_pl_lag = sum(iptw_ipcw_pl_lag_event)
  ) %>%
  select(
    year,
    total_events,
    total_events_iptw_ipcw_str_lag,
    total_events_iptw_ipcw_str_nonlag,
    total_events_iptw_ipcw_pl_nonlag,
    total_events_iptw_ipcw_pl_lag
  ) %>% 
  slice(1)

cohort_vis <- merge(cohort, itt_person_time, by.x = 'id', by.y = 'id', all.x = TRUE)
cohort_vis <- merge(cohort_vis, at_person_time, by.x= 'id', by.y = 'id', all.x = TRUE)

#### COVARIATE ASSOCIATION WITH OUTCOME ####

base_variables <- c(covariates, base_comorb) 
grouping_var <- 'region' 

trt_col_name <- 'trt_dummy'

covscoef <- data.frame(group = c(unique(cohort_vis[,grouping_var]))) %>% 
  arrange(group)

for (i in 1:length(base_variables)) {
  var <- base_variables[i]
  
  cols <- ((colnames(cohort_vis) == grouping_var) |
             (colnames(cohort_vis) == trt_col_name) |
             (colnames(cohort_vis) == var))
  
  # extract variable info from cohort_vis
  var_df <- cohort_vis[cols]
  var_i <- colnames(var_df) == var
  
  # check variable type
  if (is.character(var_df[,var_i])) {
    var_df[,var_i] <- as.factor(var_df[,var_i])
  }
  
  if (is.factor(var_df[,var_i]) & nlevels(var_df[,var_i])>2) { # if factor with >2 levels, create dummies
    var_df <- dummy_cols(var_df, select_columns = var, remove_first_dummy = TRUE)
    cols_keep <- colnames(var_df) != var
    var_df <- var_df[cols_keep]
  }
  
  names(var_df)[colnames(var_df) == trt_col_name] <- 'trt'
  names(var_df)[colnames(var_df) == grouping_var] <- 'group'
  
  # iterate through variables (may be multiple if var was a factor)
  for (j in 1:length(var_df)) {
    if (colnames(var_df)[j] == 'trt' | colnames(var_df)[j] == 'group') {
      next
    }
    
    name <- colnames(var_df[j])
    var_sum <- var_df 
    names(var_sum)[j] <- 'var'
    
    covscoef$var <- NA
    
    # calculate association with exposure through logistic regression
    for (i in 1:nrow(covscoef)) {
      group_i <- covscoef[i,1]
      covscoef[covscoef$'group' == group_i, 'var'] <- glm(trt ~ var, data = subset(var_sum, group == group_i), family = 'binomial')$coefficients[[2]]
    }
    
    names(covscoef)[names(covscoef) == 'var'] <- name
  }
}

rm(var_df, var_sum, cols, group_i, i, j, name, var, var_i)

covscoef <- covscoef[,-1] %>%
  pivot_longer(cols = everything(), names_to = "covariate", values_to = "exp(b)")

covscoef <- cbind(region = 'CPRD', covscoef, row.names = NULL)