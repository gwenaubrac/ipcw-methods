## ---------------------------
##
## Program: 10. Create plots
##
## Purpose: Create plots and tables to communicate findings. 
## 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-10-15
##
## ---------------------------
##
## Notes: 
##   
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

# for plots
library(lubridate)
library(dplyr)
library(magrittr)
library(fastDummies)
library(ggplot2)
library(tidyr)
library(cobalt)
library(scales)
library(survival)
library(cowplot)
library(table1)
library(writexl)
library(tidysmd)
library(tidyr)
library(gridExtra)

# for analyses

colors_trt = c('#0d0887', '#fdbd22')
colors_adj = c('#f56864', '#6a00a8')

'#6a00a8'
'#9a00b1'
'#d53e4f'
'#f56864'
'#f9a463'
'#fdbd22'
'#f9c43d'
'#f8e4a7'
'#f7f4f9'

#### DEFINE PATHS ####

path_intermediate_res_main <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, 'main', 'intermediate', sep = '/')
path_intermediate_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'intermediate', sep = '/')
path_final_res <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'final', sep = '/')
path_plots <- paste('Z:/EPI/Protocol 24_004042/Gwen - IPCW + vis/results', exposure, outcome, analysis, 'plots', sep = '/')

covariates <- readRDS(file = paste(path_intermediate_res_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_intermediate_res_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_intermediate_res_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_intermediate_res_main, 'dec_comorb.rds', sep = '/'))
times_dec <- readRDS(file = paste(path_intermediate_res, 'times_dec.rds', sep = '/'))

setwd(path_plots)

ref_trt <- case_when(
  exposure == 'antidepressant' ~ 'SNRI',
  exposure == 'antihypertensive' ~ 'ARB'
)

comp_trt <- case_when(
  exposure == 'antidepressant' ~ 'SSRI',
  exposure == 'antihypertensive' ~ 'ACEI'
)

#### TABLE 1 ####

# data formatting
cohort <- readRDS(paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))

cohort %<>%
  mutate(disc_one_year = as.factor(if_else(disc == 1 & disc_date < ymd(entry_date+365), 1, 0)))

tab_cohort <- cohort
tab_cohort %<>% mutate (censor = as.factor(censor))

label(tab_cohort$disc_one_year) <- 'Discontinued treatment within one year'
label(tab_cohort$sex) <- 'Sex'
label(tab_cohort$age_at_entry) <- 'Age'
label(tab_cohort$age_group) <- 'Age group'
label(tab_cohort$itt_follow_up) <- 'ITT follow-up'
label(tab_cohort$at_follow_up) <- 'AT follow-up'
label(tab_cohort$disc) <- 'Discontinued treatment'
label(tab_cohort$switch) <- 'Switched treatment'
label(tab_cohort$depression_base) <- 'Depression'
label(tab_cohort$hepatic_disease_base) <- 'Hepatic disease'
label(tab_cohort$hyperlipidemia_base) <- 'Hyperlipidemia'
label(tab_cohort$hypertension_base) <- 'Hypertension'
label(tab_cohort$chronic_renal_disease_base) <- 'Chronic renal disease'
label(tab_cohort$heart_failure_base) <- 'Heart failure'
label(tab_cohort$ischemic_heart_disease_base) <- 'Ischemic heart disease'
label(tab_cohort$lvh_base) <- 'Left ventricular hypertrophy'
label(tab_cohort$cerebrovascular_disease_base) <- 'Cerebrovascular disease'
label(tab_cohort$arrhythmia_base) <- 'Arrhythmia'
label(tab_cohort$pvd_base) <- 'Peripheral vascular disease'
label(tab_cohort$epilepsy_base) <- 'Epilepsy'
label(tab_cohort$valvular_heart_disease_base) <- 'Valvular heart disease'
label(tab_cohort$pacemaker_base) <- 'Pacemaker'
label(tab_cohort$cardiomyopathy_base) <- 'Cardiomyopathy'
label(tab_cohort$hypokalemia_base) <- 'Hypokalemia'
label(tab_cohort$hypomagnesemia_base) <- 'Hypomagnesemia'
label(tab_cohort$acute_renal_disease_base) <- 'Acute renal disease'
label(tab_cohort$stroke_base) <- 'Stroke'
label(tab_cohort$suicidal_ideation_self_harm_base) <- 'Suicidal ideation and self-harm'
label(tab_cohort$myocardial_infarction_base) <- 'Myocardial infarction'
label(tab_cohort$copd_base) <- 'COPD'
label(tab_cohort$hypocalcemia_base) <- 'Hypocalcemia'
label(tab_cohort$deprivation) <- 'Deprivation Index'
label(tab_cohort$ethnicity) <- 'Ethnicity'
label(tab_cohort$anxiety_base) <- 'Anxiety'

units(tab_cohort$age_at_entry) <- 'years'
units(tab_cohort$itt_follow_up) <- 'days'
units(tab_cohort$at_follow_up) <- 'days'

caption <- 'Baseline Characteristics of Cohort'

# order by most common baseline comorb
order_comorb <- cohort
order_comorb[base_comorb] <- lapply(order_comorb[base_comorb], function(x) as.numeric(as.character(x)))

comorb_freq <- order_comorb %>%
  summarise(across(all_of(base_comorb), ~ mean(.x, na.rm = TRUE) * 100))

comorb_freq_long <- comorb_freq %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Frequency") %>%
  arrange(desc(Frequency))

ordered_base_comorb <- comorb_freq_long$Variable

# create table
table1_formula <- as.formula(paste("~", paste(c("sex", "age_at_entry", "age_group", "deprivation", "ethnicity", "disc_one_year", "itt_follow_up", "at_follow_up", "disc", "switch", ordered_base_comorb), collapse = "+"), "|trt"))
table1 <- table1(
  table1_formula,
  data = tab_cohort,
  overall = c(right = 'Total'),
  caption = caption
)

table1
write.table (table1 , "table1.csv", col.names = T, row.names=F, append= F, sep=',')

# get SMDs
tidy_smd(cohort, c(age_at_entry), .group = trt, .wts = c(iptw))
tidy_smd(cohort, c(age_group), .group = trt, .wts = c(iptw))

#### COVARIATE ASSOCIATION OVER TIME ####

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_covariates.rds', sep = '/'))

# choose which variables to plot:

# most common comorbidities:
base_variables <- c('hypertension_base', 'depression_base', 'hyperlipidemia_base', 'suicidal_ideation_self_harm_base') 

# comorbidities most associated with trt:
base_variables <- c('age_group')

# all specified covariates
base_variables <- c(covariates) # specified covariates

grouping_var <- 'month_year' 
trt_col_name <- 'trt_dummy'

cov_desc <- data.frame(group = c(unique(cohort[,grouping_var]))) %>% 
  arrange(group)

for (i in 1:length(base_variables)) {
  var <- base_variables[i]
  
  cols <- ((colnames(cohort) == grouping_var) |
             (colnames(cohort) == trt_col_name) |
             (colnames(cohort) == var))
  
  # extract variable info from cohort
  var_df <- cohort[cols]
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
    
    cov_desc$var <- NA
    
    # calculate association with exposure through logistic regression
    for (i in 1:nrow(cov_desc)) {
      group_i <- cov_desc[i,1]
      cov_desc[cov_desc$'group' == group_i, 'var'] <- glm(trt ~ var, data = subset(var_sum, group == group_i), family = 'binomial')$coefficients[[2]]
    }
    
    names(cov_desc)[names(cov_desc) == 'var'] <- name
  }
}

cov_desc_long <- cov_desc %>% 
  pivot_longer(cols = -group, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

cov_desc_long <- cov_desc_long %>%
  mutate(
    variable = variable %>%
      sub("_base$", "", .) %>%
      gsub("_", " ", .)
  )

cov_desc_long$group <- as.character(cov_desc_long$group)
cov_desc_long$group <- as.Date(paste0(cov_desc_long$group, "-01"), format = "%Y-%m-%d")

p <- ggplot(cov_desc_long, aes(x = group, y = value, color = variable, group = variable)) +
  geom_point(size = 1) +
  geom_line(linewidth = 0.5) +
  labs(x = "calendar time", y = "coefficient of association") +
  scale_color_viridis_d(option = "C", begin = 0.1, end = 0.9, direction = -1, name = 'Covariate') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 10, face = 'bold'), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom", 
    legend.title = element_text(size = 8, face = 'bold'),
    legend.text = element_text(size = 8),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_date(
    breaks = date_breaks("6 months"), 
    labels = date_format("%b %Y")
  ) + 
  ggtitle("Association Between Covariates and Treatment over Time")

p

ggsave("covariate_association_plot.png", plot = p, width = 6, height = 4, units = "in", bg = 'white')

dev.off()

rm(cov_desc, cov_desc_long, p, var_df, var_sum)
rm(cols, cols_keep, group_i, grouping_var, i, j, name, trt_col_name, var, var_i)
rm(cohort)

#### PROPENSITY SCORE DENSITY ####

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_update_cov.rds', sep = '/'))

model <- readRDS(paste(path_final_res, 'iptw_model.rds', sep = '/'))

weights_out <- WeightIt::weightit(
  formula = model,
  data = cohort,
  method = 'glm')

weights_out$treat <- if_else(weights_out$treat == 0, ref_trt, comp_trt)

p <- bal.plot(
  weights_out,
  var.name = "prop.score",
  which = "both",
  type = "density",
  mirror = FALSE,
  colors = colors_trt, 
  grid = FALSE,
  position = 'bottom'
)

p <- p +
  ggtitle("Propensity Score Distributional Balance") +
  xlab("propensity score") +
  ylab("density")  

p

ggsave("propensity_score_plot.png", plot = p, width = 6, height = 4, units = "in", bg = 'white')

dev.off()

#### LOVE PLOT ####

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))
covs <- readRDS(file = paste(path_intermediate_res, 'iptw_covs.rds', sep = '/'))

bal_tab <- bal.tab(trt_dummy ~ covs, data = cohort,
                   weights = 'iptw',
                   binary = "std", continuous = "std", 
                   stats = c('mean.diffs'),
                   s.d.denom = 'pooled',
                   un = TRUE, 
                   thresholds = c(m = 0.1)
) 


iptw_balance_tab <- bal_tab$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (IPTW)' = Diff.Adj,
                'Balance' = M.Threshold)

iptw_balance_tab$Variable <- row.names(iptw_balance_tab)

write.table (iptw_balance_tab , "iptw_balance_tab.csv", col.names = T, row.names=T, append= F, sep=',')

var_names <- c(iptw_balance_tab$Variable)
loveplot_var_names <- sub("_base$", "", var_names)
loveplot_var_names <- gsub("_", " ", loveplot_var_names)
loveplot_var_names[loveplot_var_names=='sex Male'] <- 'male'
loveplot_var_names[loveplot_var_names=='pvd'] <- 'perivascular disease'
loveplot_var_names[loveplot_var_names=='lvh'] <- 'left ventricular heart disease'
loveplot_var_names[loveplot_var_names=='copd'] <- 'chronic obstructive pulmonary disease'
loveplot_var_names <- setNames(loveplot_var_names, var_names)
loveplot_var_names

p <- love.plot(
  bal_tab,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Covariate Balance with IPTW',
  colors = colors_adj,
  var.names = loveplot_var_names,
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2
)


p <- p + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title.x = element_text(size = 8, face = 'bold'),
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  
  legend.position = "bottom", 
  legend.title = element_text(size = 8, face = 'bold'),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(0.5, "lines")
) 

p

ggsave("loveplot.png", plot = p, width = 6, height = 6, units = "in", bg = 'white')

dev.off()

rm(covs, p, iptw_balance_tab, bal_tab)
rm(cohort)

#### CUMULATIVE INCIDENCE BY TRT ####

cohort_long <- readRDS(file = paste(path_intermediate_res, 'cohort_analytic_at.rds', sep = '/'))
cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))

## REFERENCE GROUP ##

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

cohort_ref <- cohort %>% 
  filter(trt_dummy == 0)

ref_itt <- survfit(Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ 1, data = cohort_ref)
ref_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ 1, data = cohort_ref, weights = iptw)
ref_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref)
ref_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = iptw)

ref_at_iptw_ipcw_str_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_str_lag*iptw)
ref_at_iptw_ipcw_str_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_str_nonlag*iptw)
ref_at_iptw_ipcw_pl_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_pl_nonlag*iptw)
ref_at_iptw_ipcw_pl_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_pl_lag*iptw)

# ITT vs AT

if (exposure == 'antidepressant') {
  y_limit <- case_when(
    analysis == 'main' ~ c(0, 0.08),
    analysis == '90_day_grace_period' ~ c(0, 0.08), # check
    analysis == 'flexible_grace_period' ~ c(0, 0.08), # check
    analysis == 'male' ~ c(0, 0.1), # check
    analysis == 'female' ~ c(0, 0.08), # check
    analysis == 'young' ~ c(0, 0.03), # check
    analysis == 'old' ~ c(0, 0.2), # check
    analysis == '2019' ~ c(0, 0.08), # check
    analysis == '2020' ~ c(0, 0.07), # check
    analysis == '2021' ~ c(0, 0.08), # check
    analysis == '2022' ~ c(0, 0.03), # check
    analysis == 'depressed' ~ c(0, 0.06), # check
    analysis == 'not_depressed' ~ c(0, 0.08), # check
    TRUE ~ c(NA, NA)
  )
}

if (exposure == 'antihypertensive') {
  y_limit <- case_when(
    analysis == 'main' ~ c(0, 0.035),
    analysis == 'flexible_grace_period' ~ c(0, 0.08), # check
    analysis == 'male' ~ c(0, 0.1), # check
    analysis == 'female' ~ c(0, 0.08), # check
    analysis == 'young' ~ c(0, 0.03), # check
    analysis == 'old' ~ c(0, 0.2), # check
    analysis == '2019' ~ c(0, 0.08), # check
    analysis == '2020' ~ c(0, 0.07), # check
    analysis == '2021' ~ c(0, 0.08), # check
    analysis == '2022' ~ c(0, 0.03), # check
    analysis == 'hypertensive' ~ c(0, 0.06), # check
    analysis == 'not_hypertensive' ~ c(0, 0.08), # check
    TRUE ~ c(NA, NA)
  )
}

png("ref_incidence_by_analysis.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", ref_trt, "Group"),
     col = '#f56864',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

lines(ref_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(ref_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)"),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8'),
       lty = c(2, 1, 2, 1), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

png(filename = "ref_incidence_by_model.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_at_iptw_ipcw_str_lag,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", ref_trt, "Group"),
     col = '#9a00b1',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw_ipcw_str_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_pl_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_pl_lag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + stratified lagged IPCW)", "AT (IPTW + stratified non-lagged IPCW)", 
                  "AT (IPTW + pooled non-lagged IPCW)", "AT (IPTW + pooled lagged IPCW"),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'),
       lty = c(2, 1, 1, 1, 1), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


## COMPARATOR GROUP ## 

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1)

cohort_comp <- cohort %>% 
  filter(trt_dummy == 1)

comp_itt <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp)
comp_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp, weights = iptw)
comp_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp)
comp_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw)

comp_at_iptw_ipcw_str_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_str_lag*iptw)
comp_at_iptw_ipcw_str_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_str_nonlag*iptw)
comp_at_iptw_ipcw_pl_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_pl_nonlag*iptw)
comp_at_iptw_ipcw_pl_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_pl_lag*iptw)

# ITT vs AT

png("comp_incidence_by_analysis.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", comp_trt, "Group"),
     col = '#f56864',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

lines(comp_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(comp_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)"),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8'),
       lty = c(2, 1, 2, 1), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

png(filename = "comp_incidence_by_model.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_at_iptw_ipcw_str_lag,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", comp_trt, "Group"),
     col = '#9a00b1',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw_ipcw_str_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_pl_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_pl_lag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + stratified lagged IPCW)", "AT (IPTW + stratified non-lagged IPCW)", 
                  "AT (IPTW + pooled non-lagged IPCW)", "AT (IPTW + pooled lagged IPCW"),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'),
       lty = c(2, 1, 1, 1, 1), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()

#### CUMULATIVE INCIDENCE ZOOM IN (365 DAYS) ####

cohort_long <- readRDS(file = paste(path_final_res, 'cohort_analytic_at.rds', sep = '/'))
cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))
times_dec <- readRDS(paste(path_intermediate_res, 'times_dec.rds', sep = '/'))
times_dec

## REFERENCE GROUP ## 

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

cohort_ref <- cohort %>% 
  filter(trt_dummy == 0)

ref_itt <- survfit(Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ 1, data = cohort_ref)
ref_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ 1, data = cohort_ref, weights = iptw)
ref_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref)
ref_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = iptw)

ref_at_iptw_ipcw_str_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_str_lag*iptw)
ref_at_iptw_ipcw_str_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_str_nonlag*iptw)
ref_at_iptw_ipcw_pl_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_pl_nonlag*iptw)
ref_at_iptw_ipcw_pl_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_pl_lag*iptw)

# ITT vs AT

x_limit <- c(0, 365)

if (exposure == 'antidepressant') {
  y_limit <- case_when(
    analysis == 'main' ~ c(0, 0.04),
    analysis == '90_day_grace_period' ~ c(0, 0.04), # check
    analysis == 'flexible_grace_period' ~ c(0, 0.04), # check
    analysis == 'male' ~ c(0, 0.05), # check
    analysis == 'female' ~ c(0, 0.04), # check
    analysis == 'young' ~ c(0, 0.1), # check
    analysis == 'old' ~ c(0, 0.2), # check
    analysis == '2019' ~ c(0, 0.08), # check
    analysis == '2020' ~ c(0, 0.03), # check
    analysis == '2021' ~ c(0, 0.04), # check
    analysis == '2022' ~ c(0, 0.015), # check
    analysis == 'depressed' ~ c(0, 0.03), # check
    analysis == 'not_depressed' ~ c(0, 0.04), # check
    TRUE ~ c(NA, NA)
  )
}

if (exposure == 'antihypertensive') {
  y_limit <- case_when(
    analysis == 'main' ~ c(0, 0.02),
    analysis == '90_day_grace_period' ~ c(0, 0.04), # check
    analysis == 'flexible_grace_period' ~ c(0, 0.04), # check
    analysis == 'male' ~ c(0, 0.05), # check
    analysis == 'female' ~ c(0, 0.04), # check
    analysis == 'young' ~ c(0, 0.1), # check
    analysis == 'old' ~ c(0, 0.2), # check
    analysis == '2019' ~ c(0, 0.08), # check
    analysis == '2020' ~ c(0, 0.03), # check
    analysis == '2021' ~ c(0, 0.04), # check
    analysis == '2022' ~ c(0, 0.015), # check
    analysis == 'depressed' ~ c(0, 0.03), # check
    analysis == 'not_depressed' ~ c(0, 0.04), # check
    TRUE ~ c(NA, NA)
  )
}

png("ref_incidence_by_analysis_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", ref_trt, "Group"),
     col = '#f56864',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(ref_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(ref_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)", "intervals"),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8', '#f8e4a7'),
       lty = c(2, 1, 2, 1, 3), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

png(filename = "ref_incidence_by_model_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_at_iptw_ipcw_str_lag,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", ref_trt, "Group"),
     col = '#9a00b1',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw_ipcw_str_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_pl_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_pl_lag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + stratified lagged IPCW)", "AT (IPTW + stratified non-lagged IPCW)", 
                  "AT (IPTW + pooled non-lagged IPCW)", "AT (IPTW + pooled lagged IPCW", 'intervals'),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d', '#f8e4a7'),
       lty = c(2, 1, 1, 1, 1, 3), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


## COMPARATOR GROUP ## 

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1) 

cohort_comp <- cohort %>% 
  filter(trt_dummy == 1)

comp_itt <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp)
comp_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp, weights = iptw)
comp_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp)
comp_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw)

comp_at_iptw_ipcw_str_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_str_lag*iptw)
comp_at_iptw_ipcw_str_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_str_nonlag*iptw)
comp_at_iptw_ipcw_pl_nonlag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_pl_nonlag*iptw)
comp_at_iptw_ipcw_pl_lag <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_pl_lag*iptw)

# ITT vs AT

png("comp_incidence_by_analysis_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", comp_trt, "Group"),
     col = '#f56864',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(comp_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(comp_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)", 'intervals'),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8', '#f8e4a7'),
       lty = c(2, 1, 2, 1, 3), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

png(filename = "comp_incidence_by_model_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_at_iptw_ipcw_str_lag,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = paste("Cumulative Incidence Rate in", comp_trt, "Group"),
     col = '#9a00b1',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw_ipcw_str_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_pl_nonlag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_pl_lag,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + stratified lagged IPCW)", "AT (IPTW + stratified non-lagged IPCW)", 
                  "AT (IPTW + pooled non-lagged IPCW)", "AT (IPTW + pooled lagged IPCW", 'intervals'),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d', '#f8e4a7'),
       lty = c(2, 1, 1, 1, 1, 3), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()

#### FOREST PLOT HR - STABILIZED ####

cox_itt <- readRDS(paste(path_final_res, 'cox_itt.rds', sep = '/'))
cox_itt_siptw <- readRDS(paste(path_final_res, 'cox_itt_siptw.rds', sep = '/'))
cox_at <- readRDS(paste(path_final_res, 'cox_at.rds', sep = '/'))
cox_at_siptw <- readRDS(paste(path_final_res, 'cox_at_siptw.rds', sep = '/'))
cox_at_sipcw_str_lag <- readRDS(paste(path_final_res, 'cox_at_sipcw_str_lag.rds', sep = '/'))
cox_at_siptw_sipcw_str_lag <- readRDS(paste(path_final_res, 'cox_at_siptw_sipcw_str_lag.rds', sep = '/'))
cox_at_sipcw_str_nonlag <- readRDS(paste(path_final_res, 'cox_at_sipcw_str_nonlag.rds', sep = '/'))
cox_at_siptw_sipcw_str_nonlag <- readRDS(paste(path_final_res, 'cox_at_siptw_sipcw_str_nonlag.rds', sep = '/'))
cox_at_sipcw_pl_nonlag <- readRDS(paste(path_final_res, 'cox_at_sipcw_pl_nonlag.rds', sep = '/'))
cox_at_siptw_sipcw_pl_nonlag <- readRDS(paste(path_final_res, 'cox_at_siptw_sipcw_pl_nonlag.rds', sep = '/'))
cox_at_sipcw_pl_lag <- readRDS(paste(path_final_res, 'cox_at_sipcw_pl_lag.rds', sep = '/'))
cox_at_siptw_sipcw_pl_lag <- readRDS(paste(path_final_res, 'cox_at_siptw_sipcw_pl_lag.rds', sep = '/'))

cox_result_chart <- data.frame(matrix(nrow = 12, ncol = 3))
colnames(cox_result_chart) <- c('estimate', 'lower_ci', 'upper_ci')
rownames(cox_result_chart) <- c('ITT', 'ITT (sIPTW)', 'AT', 'AT (sIPTW)', 'AT (stratified lagged sIPCW)', 'AT (sIPTW + stratified lagged sIPCW)',
                                 'AT (stratified non-lagged sIPCW)', 'AT (sIPTW + stratified non-lagged sIPCW)', 'AT (pooled non-lagged sIPCW)',
                                 'AT (sIPTW + pooled non-lagged sIPCW)', 'AT (pooled lagged sIPCW)', 'AT (sIPTW + pooled lagged sIPCW)')

cox_result_chart[1, 'estimate'] <- exp(cox_itt$coef)
cox_result_chart[1, 'lower_ci'] <- exp(confint(cox_itt))[1]
cox_result_chart[1, 'upper_ci'] <- exp(confint(cox_itt))[2]

cox_result_chart[2, 'estimate'] <- exp(cox_itt_siptw$coef)
cox_result_chart[2, 'lower_ci'] <- exp(confint(cox_itt_siptw))[1]
cox_result_chart[2, 'upper_ci'] <- exp(confint(cox_itt_siptw))[2]

cox_result_chart[3, 'estimate'] <- exp(cox_at$coef)
cox_result_chart[3, 'lower_ci'] <- exp(confint(cox_at))[1]
cox_result_chart[3, 'upper_ci'] <- exp(confint(cox_at))[2]

cox_result_chart[4, 'estimate'] <- exp(cox_at_siptw$coef)
cox_result_chart[4, 'lower_ci'] <- exp(confint(cox_at_siptw))[1]
cox_result_chart[4, 'upper_ci'] <- exp(confint(cox_at_siptw))[2]

cox_result_chart[5, 'estimate'] <- exp(cox_at_sipcw_str_lag$coef)
cox_result_chart[5, 'lower_ci'] <- exp(confint(cox_at_sipcw_str_lag))[1]
cox_result_chart[5, 'upper_ci'] <- exp(confint(cox_at_sipcw_str_lag))[2]

cox_result_chart[6, 'estimate'] <- exp(cox_at_siptw_sipcw_str_lag$coef)
cox_result_chart[6, 'lower_ci'] <- exp(confint(cox_at_siptw_sipcw_str_lag))[1]
cox_result_chart[6, 'upper_ci'] <- exp(confint(cox_at_siptw_sipcw_str_lag))[2]

cox_result_chart[7, 'estimate'] <- exp(cox_at_sipcw_str_nonlag$coef)
cox_result_chart[7, 'lower_ci'] <- exp(confint(cox_at_sipcw_str_nonlag))[1]
cox_result_chart[7, 'upper_ci'] <- exp(confint(cox_at_sipcw_str_nonlag))[2]

cox_result_chart[8, 'estimate'] <- exp(cox_at_siptw_sipcw_str_nonlag$coef)
cox_result_chart[8, 'lower_ci'] <- exp(confint(cox_at_siptw_sipcw_str_nonlag))[1]
cox_result_chart[8, 'upper_ci'] <- exp(confint(cox_at_siptw_sipcw_str_nonlag))[2]

cox_result_chart[9, 'estimate'] <- exp(cox_at_sipcw_pl_nonlag$coef)
cox_result_chart[9, 'lower_ci'] <- exp(confint(cox_at_sipcw_pl_nonlag))[1]
cox_result_chart[9, 'upper_ci'] <- exp(confint(cox_at_sipcw_pl_nonlag))[2]

cox_result_chart[10, 'estimate'] <- exp(cox_at_siptw_sipcw_pl_nonlag$coef)
cox_result_chart[10, 'lower_ci'] <- exp(confint(cox_at_siptw_sipcw_pl_nonlag))[1]
cox_result_chart[10, 'upper_ci'] <- exp(confint(cox_at_siptw_sipcw_pl_nonlag))[2]

cox_result_chart[11, 'estimate'] <- exp(cox_at_sipcw_pl_lag$coef)
cox_result_chart[11, 'lower_ci'] <- exp(confint(cox_at_sipcw_pl_lag))[1]
cox_result_chart[11, 'upper_ci'] <- exp(confint(cox_at_sipcw_pl_lag))[2]

cox_result_chart[12, 'estimate'] <- exp(cox_at_siptw_sipcw_pl_lag$coef)
cox_result_chart[12, 'lower_ci'] <- exp(confint(cox_at_siptw_sipcw_pl_lag))[1]
cox_result_chart[12, 'upper_ci'] <- exp(confint(cox_at_siptw_sipcw_pl_lag))[2]

cox_result_chart <- cbind(model = row.names(cox_result_chart), cox_result_chart, row.names = NULL)
write_xlsx(cox_result_chart, paste(path_final_res, 'cox_ratios_stab.xlsx', sep ='/'))

summary_cox_results <- cox_result_chart[c(2,4,6,8,10,12),]
model_order <- rev(c("ITT (sIPTW)", "AT (sIPTW)", "AT (sIPTW + stratified lagged sIPCW)", "AT (sIPTW + stratified non-lagged sIPCW)", 
                     "AT (sIPTW + pooled non-lagged sIPCW)", "AT (sIPTW + pooled lagged sIPCW)"))

summary_cox_results$model <- factor(summary_cox_results$model, levels = model_order)
model_colors <- rev(c('#f56864','#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'))

p <- ggplot(summary_cox_results, aes(x = estimate, y = model, color = model)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +  
  labs(x = "HR (95% CI)", y = "model", title = "Forest Plot of Hazard Ratios by Model") +  
  scale_color_manual(values = model_colors) +  
  theme(
    axis.text.y = element_text(size = 8), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(25, 0, 25, 0) 
  )

p

y_labels <- rev(paste(
  c(round(summary_cox_results$estimate, 2)),
  ' (',round(summary_cox_results$lower_ci, 2),', ',
  round(summary_cox_results$upper_ci, 2), ')',
  sep = ''
))

# in case some values are exactly the same for labeling:
y_labels <- rev(paste(
  c(round(summary_cox_results$estimate, 3)),
  ' (',round(summary_cox_results$lower_ci, 3),', ',
  round(summary_cox_results$upper_ci, 3), ')',
  sep = ''
))

y_labels

dummy_data <- data.frame(
  model = factor(y_labels, levels = y_labels),
  dummy = 1 
)

p2 <- ggplot(dummy_data, aes(x = dummy, y = model)) +
  scale_y_discrete(labels = y_labels) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(40, 0, 50, 10),
    legend.position = "none" 
  )

p2

combined_plot <- plot_grid(p, p2, ncol = 2, rel_widths = c(3, 1))
combined_plot

ggsave("forest_plot_HR_stab.png", plot = combined_plot, width = 6, height = 3, units = "in", bg = 'white')

dev.off()



#### FOREST PLOT HR - UNSTABILIZED ####

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

cox_result_chart <- data.frame(matrix(nrow = 12, ncol = 3))
colnames(cox_result_chart) <- c('estimate', 'lower_ci', 'upper_ci')
rownames(cox_result_chart) <- c('ITT', 'ITT (IPTW)', 'AT', 'AT (IPTW)', 'AT (stratified lagged IPCW)', 'AT (IPTW + stratified lagged IPCW)',
                            'AT (stratified non-lagged IPCW)', 'AT (IPTW + stratified non-lagged IPCW)', 'AT (pooled non-lagged IPCW)',
                            'AT (IPTW + pooled non-lagged IPCW)', 'AT (pooled lagged IPCW)', 'AT (IPTW + pooled lagged IPCW)')

cox_result_chart[1, 'estimate'] <- exp(cox_itt$coef)
cox_result_chart[1, 'lower_ci'] <- exp(confint(cox_itt))[1]
cox_result_chart[1, 'upper_ci'] <- exp(confint(cox_itt))[2]

cox_result_chart[2, 'estimate'] <- exp(cox_itt_iptw$coef)
cox_result_chart[2, 'lower_ci'] <- exp(confint(cox_itt_iptw))[1]
cox_result_chart[2, 'upper_ci'] <- exp(confint(cox_itt_iptw))[2]

cox_result_chart[3, 'estimate'] <- exp(cox_at$coef)
cox_result_chart[3, 'lower_ci'] <- exp(confint(cox_at))[1]
cox_result_chart[3, 'upper_ci'] <- exp(confint(cox_at))[2]

cox_result_chart[4, 'estimate'] <- exp(cox_at_iptw$coef)
cox_result_chart[4, 'lower_ci'] <- exp(confint(cox_at_iptw))[1]
cox_result_chart[4, 'upper_ci'] <- exp(confint(cox_at_iptw))[2]

cox_result_chart[5, 'estimate'] <- exp(cox_at_ipcw_str_lag$coef)
cox_result_chart[5, 'lower_ci'] <- exp(confint(cox_at_ipcw_str_lag))[1]
cox_result_chart[5, 'upper_ci'] <- exp(confint(cox_at_ipcw_str_lag))[2]

cox_result_chart[6, 'estimate'] <- exp(cox_at_iptw_ipcw_str_lag$coef)
cox_result_chart[6, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_str_lag))[1]
cox_result_chart[6, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_str_lag))[2]

cox_result_chart[7, 'estimate'] <- exp(cox_at_ipcw_str_nonlag$coef)
cox_result_chart[7, 'lower_ci'] <- exp(confint(cox_at_ipcw_str_nonlag))[1]
cox_result_chart[7, 'upper_ci'] <- exp(confint(cox_at_ipcw_str_nonlag))[2]

cox_result_chart[8, 'estimate'] <- exp(cox_at_iptw_ipcw_str_nonlag$coef)
cox_result_chart[8, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_str_nonlag))[1]
cox_result_chart[8, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_str_nonlag))[2]

cox_result_chart[9, 'estimate'] <- exp(cox_at_ipcw_pl_nonlag$coef)
cox_result_chart[9, 'lower_ci'] <- exp(confint(cox_at_ipcw_pl_nonlag))[1]
cox_result_chart[9, 'upper_ci'] <- exp(confint(cox_at_ipcw_pl_nonlag))[2]

cox_result_chart[10, 'estimate'] <- exp(cox_at_iptw_ipcw_pl_nonlag$coef)
cox_result_chart[10, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_nonlag))[1]
cox_result_chart[10, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_nonlag))[2]

cox_result_chart[11, 'estimate'] <- exp(cox_at_ipcw_pl_lag$coef)
cox_result_chart[11, 'lower_ci'] <- exp(confint(cox_at_ipcw_pl_lag))[1]
cox_result_chart[11, 'upper_ci'] <- exp(confint(cox_at_ipcw_pl_lag))[2]

cox_result_chart[12, 'estimate'] <- exp(cox_at_iptw_ipcw_pl_lag$coef)
cox_result_chart[12, 'lower_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_lag))[1]
cox_result_chart[12, 'upper_ci'] <- exp(confint(cox_at_iptw_ipcw_pl_lag))[2]

cox_result_chart <- cbind(model = row.names(cox_result_chart), cox_result_chart, row.names = NULL)
write_xlsx(cox_result_chart, paste(path_final_res, 'cox_ratios_unstab.xlsx', sep ='/'))

summary_cox_results <- cox_result_chart[c(2,4,6,8,10,12),]
model_order <- rev(c("ITT (IPTW)", "AT (IPTW)", "AT (IPTW + stratified lagged IPCW)", "AT (IPTW + stratified non-lagged IPCW)", 
                     "AT (IPTW + pooled non-lagged IPCW)", "AT (IPTW + pooled lagged IPCW)"))
summary_cox_results$model <- factor(summary_cox_results$model, levels = model_order)
model_colors <- rev(c('#f56864','#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'))

p <- ggplot(summary_cox_results, aes(x = estimate, y = model, color = model)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +  
  labs(x = "HR (95% CI)", y = "model", title = "Forest Plot of Hazard Ratios by Model") +  
  scale_color_manual(values = model_colors) +  
  theme(
    axis.text.y = element_text(size = 8), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(25, 0, 25, 0) 
  )

p

y_labels <- rev(paste(
  c(round(summary_cox_results$estimate, 2)),
  ' (',round(summary_cox_results$lower_ci, 2),', ',
  round(summary_cox_results$upper_ci, 2), ')',
  sep = ''
))

# in case some values are exactly the same for labeling:
y_labels <- rev(paste(
  c(round(summary_cox_results$estimate, 3)),
  ' (',round(summary_cox_results$lower_ci, 3),', ',
  round(summary_cox_results$upper_ci, 3), ')',
  sep = ''
))

y_labels

dummy_data <- data.frame(
  model = factor(y_labels, levels = y_labels),
  dummy = 1 
)

p2 <- ggplot(dummy_data, aes(x = dummy, y = model)) +
  scale_y_discrete(labels = y_labels) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(40, 0, 50, 10),
    legend.position = "none" 
  )

p2

combined_plot <- plot_grid(p, p2, ncol = 2, rel_widths = c(3, 1))
combined_plot

ggsave("forest_plot_HR_unstab.png", plot = combined_plot, width = 6, height = 3, units = "in", bg = 'white')

dev.off()


#### CENSORING DISTRIBUTION ####

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))
times_dec <- readRDS(paste(path_intermediate_res_main, 'times_dec.rds', sep = '/'))

cens_cohort <- cohort %>% 
  filter(censor == 1) %>% 
  mutate(censor_time = as.numeric(censor_date - entry_date)) %>% 
  select(id, trt_dummy, censor_time, censor) %>% 
  arrange(censor_time)

data_d1 <- subset(cens_cohort, censor_time < times_dec[[1]])
data_d2 <- subset(cens_cohort, censor_time >= times_dec[[1]] & censor_time < times_dec[[2]])
data_d3 <- subset(cens_cohort, censor_time >= times_dec[[2]] & censor_time < times_dec[[3]])
data_d4 <- subset(cens_cohort, censor_time >= times_dec[[3]] & censor_time < times_dec[[4]])
data_d5 <- subset(cens_cohort, censor_time >= times_dec[[4]] & censor_time < times_dec[[5]])
data_d6 <- subset(cens_cohort, censor_time >= times_dec[[5]] & censor_time < times_dec[[6]])
data_d7 <- subset(cens_cohort, censor_time >= times_dec[[6]] & censor_time < times_dec[[7]])
data_d8 <- subset(cens_cohort, censor_time >= times_dec[[7]] & censor_time < times_dec[[8]])
data_d9 <- subset(cens_cohort, censor_time >= times_dec[[8]] & censor_time < times_dec[[9]])
data_d10 <- subset(cens_cohort, censor_time >= times_dec[[9]])

cens_overall <- cens_cohort %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d1 <- data_d1 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d2 <- data_d2 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d3 <- data_d3 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d4 <- data_d4 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d5 <- data_d5 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d6 <- data_d6 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d7 <- data_d7 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d8 <- data_d8 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d9 <- data_d9 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

cens_d10 <- data_d10 %>% 
  group_by(censor_time) %>% 
  summarize (count = n())

plot_censoring_distribution <- function(data, x_var, y_var, title) {
  ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }})) +
    geom_line(color = "#6a00a8", linewidth = 0.5) +
    labs(x = "time (days)", y = "patients censored", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10, face = 'bold'),
      axis.title.y = element_text(size = 10, face = 'bold'),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 8, face = 'bold'),
      legend.text = element_text(size = 8),
      legend.box.margin = margin(0, 0, 0, 0)
    )
}

censor_d1 <- plot_censoring_distribution(cens_d1, censor_time, count, "Censoring in Interval 1")
ggsave("censor_dist_d1.png", plot = censor_d1, width = 4, height = 4, units = "in", bg = 'white')

censor_d2 <- plot_censoring_distribution(cens_d2, censor_time, count, "Censoring in Interval 2")
ggsave("censor_dist_d2.png", plot = censor_d2, width = 4, height = 4, units = "in", bg = 'white')

censor_d3 <- plot_censoring_distribution(cens_d3, censor_time, count, "Censoring in Interval 3")
ggsave("censor_dist_d3.png", plot = censor_d3, width = 4, height = 4, units = "in", bg = 'white')

censor_d4 <- plot_censoring_distribution(cens_d4, censor_time, count, "Censoring in Interval 4")
ggsave("censor_dist_d4.png", plot = censor_d4, width = 4, height = 4, units = "in", bg = 'white')

censor_d5 <- plot_censoring_distribution(cens_d5, censor_time, count, "Censoring in Interval 5")
ggsave("censor_dist_d5.png", plot = censor_d5, width = 4, height = 4, units = "in", bg = 'white')

censor_d6 <- plot_censoring_distribution(cens_d6, censor_time, count, "Censoring in Interval 6")
ggsave("censor_dist_d6.png", plot = censor_d6, width = 4, height = 4, units = "in", bg = 'white')

censor_d7 <- plot_censoring_distribution(cens_d7, censor_time, count, "Censoring in Interval 7")
ggsave("censor_dist_d7.png", plot = censor_d7, width = 4, height = 4, units = "in", bg = 'white')

censor_d8 <- plot_censoring_distribution(cens_d8, censor_time, count, "Censoring in Interval 8")
ggsave("censor_dist_d8.png", plot = censor_d8, width = 4, height = 4, units = "in", bg = 'white')

censor_d9 <- plot_censoring_distribution(cens_d9, censor_time, count, "Censoring in Interval 9")
ggsave("censor_dist_d9.png", plot = censor_d9, width = 4, height = 4, units = "in", bg = 'white')

censor_d10 <- plot_censoring_distribution(cens_d10, censor_time, count, "Censoring in Interval 10")
ggsave("censor_dist_d10.png", plot = censor_d10, width = 4, height = 4, units = "in", bg = 'white')

grid_cens_plots <- grid.arrange(censor_d1, censor_d2, censor_d3, censor_d4, censor_d5, 
             censor_d6, censor_d7, censor_d8, censor_d9, censor_d10, ncol=2)
ggsave("grid_cens_plots.png", plot = grid_cens_plots, width = 8, height = 10, units = "in", bg = 'white')

#### FOREST PLOT IRR - UNSTABILIZED ####

# choose one of saved CIs (from boot function or computed manually)
incidence_rates <- readRDS(paste(path_final_res, 'incidence_rates_ci_fx.rds', sep ='/'))
#incidence_rates <- readRDS(paste(path_final_res, 'incidence_rates_ci.rds', sep ='/'))

incidence_rates_chart <- incidence_rates[c(2,4,6,8,10,12),]
incidence_rates_chart %<>% select(-variable)

rownames(incidence_rates_chart) <- c(
  'ITT (IPTW)',
  'AT (IPTW)',
  'AT (IPTW + stratified lagged IPCW)',
  'AT (IPTW + stratified non-lagged IPCW)',
  'AT (IPTW + pooled non-lagged IPCW)',
  'AT (IPTW + pooled lagged IPCW)'
)

incidence_rates_chart <- cbind(model = row.names(incidence_rates_chart), incidence_rates_chart, row.names = NULL)
model_order <- rev(c("ITT (IPTW)", "AT (IPTW)", "AT (IPTW + stratified lagged IPCW)", "AT (IPTW + stratified non-lagged IPCW)", 
                     "AT (IPTW + pooled non-lagged IPCW)", "AT (IPTW + pooled lagged IPCW)"))

incidence_rates_chart$model <- factor(incidence_rates_chart$model, levels = model_order)
model_colors <- rev(c('#f56864','#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'))

p <- ggplot(incidence_rates_chart, aes(x = estimate, y = model, color = model)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +  
  labs(x = "IRR (95% CI)", y = "model", title = "Forest Plot of Incidence Rate Ratios by Model") +  
  scale_color_manual(values = model_colors) +  
  theme(
    axis.text.y = element_text(size = 8), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(25, 0, 25, 0) 
  )

p

y_labels <- rev(paste(
  c(round(incidence_rates_chart$estimate, 2)),
  ' (',round(incidence_rates_chart$lower_ci, 2),', ',
  round(incidence_rates_chart$upper_ci, 2), ')',
  sep = ''
))

# in case some values are exactly the same for labeling:
y_labels <- rev(paste(
  c(round(incidence_rates_chart$estimate, 3)),
  ' (',round(incidence_rates_chart$lower_ci, 3),', ',
  round(incidence_rates_chart$upper_ci, 3), ')',
  sep = ''
))

y_labels

dummy_data <- data.frame(
  model = factor(y_labels, levels = y_labels),
  dummy = 1 
)

p2 <- ggplot(dummy_data, aes(x = dummy, y = model)) +
  scale_y_discrete(labels = y_labels) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(40, 0, 50, 10),
    legend.position = "none" 
  )

p2

combined_plot <- plot_grid(p, p2, ncol = 2, rel_widths = c(3, 1))
combined_plot

ggsave("forest_plot_IRR_unstab.png", plot = combined_plot, width = 6, height = 3, units = "in", bg = 'white')

dev.off()

#### FOREST PLOT IRR - STABILIZED ####

incidence_rates <- readRDS(paste(path_final_res, 'incidence_rates_ci_fx.rds', sep ='/'))
#incidence_rates <- readRDS(paste(path_final_res, 'incidence_rates_ci.rds', sep ='/'))

incidence_rates_chart <- incidence_rates[c(25,26,28,30,32,34),]
incidence_rates_chart %<>% select(-variable)

rownames(incidence_rates_chart) <- c(
  'ITT (sIPTW)',
  'AT (sIPTW)',
  'AT (sIPTW + stratified lagged sIPCW)',
  'AT (sIPTW + stratified non-lagged sIPCW)',
  'AT (sIPTW + pooled non-lagged sIPCW)',
  'AT (sIPTW + pooled lagged sIPCW)'
)

incidence_rates_chart <- cbind(model = row.names(incidence_rates_chart), incidence_rates_chart, row.names = NULL)
model_order <- rev(c("ITT (sIPTW)", "AT (sIPTW)", "AT (sIPTW + stratified lagged sIPCW)", "AT (sIPTW + stratified non-lagged sIPCW)", 
                     "AT (sIPTW + pooled non-lagged sIPCW)", "AT (sIPTW + pooled lagged sIPCW)"))

incidence_rates_chart$model <- factor(incidence_rates_chart$model, levels = model_order)
model_colors <- rev(c('#f56864','#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'))

p <- ggplot(incidence_rates_chart, aes(x = estimate, y = model, color = model)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +  
  labs(x = "IRR (95% CI)", y = "model", title = "Forest Plot of Incidence Rate Ratios by Model") +  
  scale_color_manual(values = model_colors) +  
  theme(
    axis.text.y = element_text(size = 8), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(25, 0, 25, 0) 
  )

p

y_labels <- rev(paste(
  c(round(incidence_rates_chart$estimate, 2)),
  ' (',round(incidence_rates_chart$lower_ci, 2),', ',
  round(incidence_rates_chart$upper_ci, 2), ')',
  sep = ''
))

# in case some values are exactly the same for labeling:
y_labels <- rev(paste(
  c(round(incidence_rates_chart$estimate, 3)),
  ' (',round(incidence_rates_chart$lower_ci, 3),', ',
  round(incidence_rates_chart$upper_ci, 3), ')',
  sep = ''
))

y_labels

dummy_data <- data.frame(
  model = factor(y_labels, levels = y_labels),
  dummy = 1 
)

p2 <- ggplot(dummy_data, aes(x = dummy, y = model)) +
  scale_y_discrete(labels = y_labels) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(40, 0, 50, 10),
    legend.position = "none" 
  )

p2

combined_plot <- plot_grid(p, p2, ncol = 2, rel_widths = c(3, 1))
combined_plot

ggsave("forest_plot_IRR_stab.png", plot = combined_plot, width = 6, height = 3, units = "in", bg = 'white')

dev.off()


#### HISTOGRAM OF AGE OVER TIME BY TRT ####

# Proportion young and old by interval - percentile age groups
cohort_long <- readRDS(file = paste(path_intermediate_res, 'cohort_ipcw.rds', sep = '/'))

age_dist <- quantile(cohort$age_at_entry, probs = c(0.25, 0.5, 0.75))
age_dist

cohort_long %<>% mutate(combo_weight = iptw*ipcw_pl_nonlag,
                        age_25 = if_else(age_at_entry < age_dist[[1]], 1 , 0),
                        age_25_50 = if_else(age_at_entry >= age_dist[[1]] & age_at_entry < age_dist[[2]], 1, 0),
                        age_50_75 = if_else(age_at_entry >= age_dist[[2]] & age_at_entry < age_dist[[3]], 1, 0),
                        age_75 = if_else(age_at_entry >=age_dist[[3]], 1, 0))

hist_age_perc <- cohort_long %>% 
  group_by(trt_dummy, dec) %>% 
  summarize(prop_25 = sum(age_25 == 1) / n(), 
            prop_25_50 = sum(age_25_50 == 1) / n(),
            prop_50_75 = sum(age_50_75 == 1) / n(),
            prop_75 = sum(age_75 == 1) / n())

df_long <- pivot_longer(hist_age_perc, cols = c(prop_25, prop_25_50, prop_50_75, prop_75), 
                        names_to = "age_group", 
                        values_to = "proportion")
# unweighted
hist_age_perc <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9a00b1', '#f56864', '#f9a463', '#f8e4a7'),
                    labels = c('<25%', '25-50%', '50-75%', '>75%')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Age Groups Over Time by Treatment",
       x = "Interval",
       y = "Proportion",
       fill = "Age Percentile") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_age_perc
ggsave("hist_age_perc.png", plot = hist_age_perc, width = 6, height = 3, units = "in", bg = 'white')


# weighted with IPTW
hist_age_perc_weighted_iptw <- cohort_long %>% # weighted with IPTW
  group_by(trt_dummy, dec) %>%
  summarize(
    prop_25 = sum(iptw[age_25==1]) / sum(iptw),
    prop_25_50 = sum(iptw[age_25_50==1]) / sum(iptw),
    prop_50_75 = sum(iptw[age_50_75==1]) / sum(iptw),
    prop_75 = sum(iptw[age_75==1]) / sum(iptw),
  )

df_long <- pivot_longer(hist_age_perc_weighted_iptw, cols = c(prop_25, prop_25_50, prop_50_75, prop_75), 
                        names_to = "age_group", 
                        values_to = "proportion")

hist_age_perc_weighted_iptw <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9a00b1', '#f56864', '#f9a463', '#f8e4a7'),
                    labels = c('<25%', '25-50%', '50-75%', '>75%')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Age Groups Over Time by Treatment (IPTW + IPCW)",
       x = "Interval",
       y = "Proportion",
       fill = "Age Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_age_perc_weighted_iptw
ggsave("hist_age_perc_weighted_iptw.png", plot = hist_age_perc_weighted_iptw, width = 6, height = 3, units = "in", bg = 'white')

# weighted with IPTW + IPCW
hist_age_perc_weighted_combo <- cohort_long %>% 
  group_by(trt_dummy, dec) %>%
  summarize(
    prop_25 = sum(combo_weight[age_25==1]) / sum(combo_weight),
    prop_25_50 = sum(combo_weight[age_25_50==1]) / sum(combo_weight),
    prop_50_75 = sum(combo_weight[age_50_75==1]) / sum(combo_weight),
    prop_75 = sum(combo_weight[age_75==1]) / sum(combo_weight),
  )

df_long <- pivot_longer(hist_age_perc_weighted_combo, cols = c(prop_25, prop_25_50, prop_50_75, prop_75), 
                        names_to = "age_group", 
                        values_to = "proportion")

hist_age_perc_weighted_combo <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#9a00b1', '#f56864', '#f9a463', '#f8e4a7'),
                    labels = c('<25%', '25-50%', '50-75%', '>75%')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Age Groups Over Time by Treatment (IPTW + IPCW)",
       x = "Interval",
       y = "Proportion",
       fill = "Age Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_age_perc_weighted_combo
ggsave("hist_age_perc_weighted_combo.png", plot = hist_age_perc_weighted_combo, width = 6, height = 3, units = "in", bg = 'white')

# Proportion young and old by decile - age groups
# does not work for antihypertensive yet -- need to update new age groups
cohort_long <- readRDS(file = paste(path_intermediate_res, 'cohort_ipcw.rds', sep = '/'))
cohort_long %<>% mutate(combo_weight = iptw*ipcw_pl_nonlag)

hist_age <- cohort_long %>% 
  group_by(trt_dummy, dec) %>% 
  summarize(group_18_24 = sum(age_group == '18-24') / n(), 
            group_25_34 = sum(age_group == '25-34') / n(),
            group_35_44 = sum(age_group == '35-44') / n(),
            group_45_54 = sum(age_group == '45-54') / n(),
            group_55_64 = sum(age_group == '55-64') / n(),
            group_65_74 = sum(age_group == '65-74') / n(),
            group_75_84 = sum(age_group == '75-84') / n(),
            group_85_over = sum(age_group == '>85') / n()
            )

df_long <- pivot_longer(hist_age, cols = c(group_18_24, group_25_34, group_35_44, 
                                           group_45_54, group_55_64, group_65_74, group_75_84, group_85_over), 
                        names_to = "age_group", 
                        values_to = "proportion")

hist_age <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#6a00a8', '#9a00b1', '#d53e4f', '#f56864', '#f9a463', '#fdbd22','#f8e4a7', '#E0E0E0'),
                    labels = c('18-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '>85')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Age Groups Over Time by Treatment",
       x = "Interval",
       y = "Proportion",
       fill = "Age Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_age
ggsave("hist_age.png", plot = hist_age, width = 6, height = 3, units = "in", bg = 'white')


hist_age_weighted_iptw <- cohort_long %>% # weighted with IPTW only
  group_by(trt, dec) %>%
  summarize(
    group_18_24 = sum(siptw[age_group == '18-24']) / sum(siptw),
    group_25_34 = sum(siptw[age_group == '25-34']) / sum(siptw),
    group_35_44 = sum(siptw[age_group == '35-44']) / sum(siptw),
    group_45_54 = sum(siptw[age_group == '45-54']) / sum(siptw),
    group_55_64 = sum(siptw[age_group == '55-64']) / sum(siptw),
    group_65_74 = sum(siptw[age_group == '65-74']) / sum(siptw),
    group_75_84 = sum(siptw[age_group == '75-84']) / sum(siptw),
    group_85_over = sum(siptw[age_group == '>85']) / sum(siptw)
  )

df_long <- pivot_longer(hist_age_weighted_iptw, cols = c(group_18_24, group_25_34, group_35_44, 
                                                          group_45_54, group_55_64, group_65_74, group_75_84, group_85_over), 
                        names_to = "age_group", 
                        values_to = "proportion")

hist_age_weighted_iptw <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#6a00a8', '#9a00b1', '#d53e4f', '#f56864', '#f9a463', '#fdbd22','#f8e4a7', '#E0E0E0'),
                    labels = c('18-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '>85')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Age Groups Over Time by Treatment (IPTW)",
       x = "Interval",
       y = "Proportion",
       fill = "Age Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


hist_age_weighted_iptw
ggsave("hist_age_weighted_iptw.png", plot = hist_age_weighted_iptw, width = 6, height = 3, units = "in", bg = 'white')



hist_age_weighted_combo <- cohort_long %>% # weighted with IPTW + IPCW
  group_by(trt, dec) %>%
  summarize(
    group_18_24 = sum(combo_weight[age_group == '18-24']) / sum(combo_weight),
    group_25_34 = sum(combo_weight[age_group == '25-34']) / sum(combo_weight),
    group_35_44 = sum(combo_weight[age_group == '35-44']) / sum(combo_weight),
    group_45_54 = sum(combo_weight[age_group == '45-54']) / sum(combo_weight),
    group_55_64 = sum(combo_weight[age_group == '55-64']) / sum(combo_weight),
    group_65_74 = sum(combo_weight[age_group == '65-74']) / sum(combo_weight),
    group_75_84 = sum(combo_weight[age_group == '75-84']) / sum(combo_weight),
    group_85_over = sum(combo_weight[age_group == '>85']) / sum(combo_weight)
  )

df_long <- pivot_longer(hist_age_weighted_combo, cols = c(group_18_24, group_25_34, group_35_44, 
                                           group_45_54, group_55_64, group_65_74, group_75_84, group_85_over), 
                        names_to = "age_group", 
                        values_to = "proportion")

hist_age_weighted_combo <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = age_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#6a00a8', '#9a00b1', '#d53e4f', '#f56864', '#f9a463', '#fdbd22','#f8e4a7', '#E0E0E0'),
                    labels = c('18-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '>85')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Age Groups Over Time by Treatment (IPTW + IPCW)",
       x = "Interval",
       y = "Proportion",
       fill = "Age Group") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


hist_age_weighted_combo
ggsave("hist_age_weighted_combo.png", plot = hist_age_weighted_combo, width = 6, height = 3, units = "in", bg = 'white')


#### HISTOGRAM OF CENSORING OVER TIME BY TRT ####

# Proportion of censoring by treatment over time
hist_cens <- cohort_long %>% # unweighted
  group_by(trt_dummy, dec) %>% 
  summarize(prop_cens = sum(uncensored_at_tstop == 0) / n(), 
            prop_uncens = sum(uncensored_at_tstop == 1) / n())

df_long <- pivot_longer(hist_cens, cols = c(prop_cens, prop_uncens), 
                        names_to = "censoring", 
                        values_to = "proportion")

hist_cens <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = censoring)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6a00a8", "#f56864"),
                    labels = c('uncensored', 'censored')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Status Over Time by Treatment",
       x = "Interval",
       y = "Proportion",
       fill = "Status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_cens
ggsave("hist_cens.png", plot = hist_cens, width = 6, height = 3, units = "in", bg = 'white')

hist_cens_weighted <- cohort_long %>% # IPTW + IPCW weighted
  group_by(trt_dummy, dec) %>% 
  summarize(prop_cens = sum(combo_weight[uncensored_at_tstop == 0]) / sum(combo_weight), 
            prop_uncens = sum(combo_weight[uncensored_at_tstop == 1]) / sum(combo_weight))

df_long <- pivot_longer(hist_cens_weighted, cols = c(prop_cens, prop_uncens), 
                        names_to = "censoring", 
                        values_to = "proportion")

hist_cens_weighted <- ggplot(df_long, aes(x = factor(dec), y = proportion, fill = censoring)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6a00a8", "#f56864"),
                    labels = c('censored', 'uncensored')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Status Over Time by Treatment (IPTW + IPCW)",
       x = "Interval",
       y = "Proportion",
       fill = "Status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_cens_weighted
ggsave("hist_cens_weighted.png", plot = hist_cens_weighted, width = 6, height = 3, units = "in", bg = 'white')

#### HISTOGRAM OF CENSORING BY AGE AND TRT ####

# Proportion of discontinuation by age and treatment group
hist_cens_age <- cohort %>% 
  mutate(age_round = as.factor(round(age_at_entry / 10) * 10)) %>% 
  group_by(trt_dummy, age_round) %>% 
  summarize(prop_cens = sum(censor == 0) / n(), 
            prop_uncens = sum(censor == 1) / n())

df_long <- pivot_longer(hist_cens_age, cols = c(prop_cens, prop_uncens), 
                        names_to = "censoring", 
                        values_to = "proportion")

df_long$age_round <- factor(df_long$age_round, levels = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110))

hist_cens_age <- ggplot(df_long, aes(x = factor(age_round), y = proportion, fill = censoring)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6a00a8", "#f56864"),
                    labels = c('uncensored', 'censored')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Status by Age and Treatment",
       x = "Age",
       y = "Proportion",
       fill = "Status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7.5))

hist_cens_age
ggsave("hist_cens_age.png", plot = hist_cens_age, width = 6, height = 3, units = "in", bg = 'white')

hist_cens_age_weighted_iptw <- cohort_long %>% # weighted with IPTW
  mutate(age_round = as.factor(round(age_at_entry / 10) * 10)) %>% 
  group_by(trt_dummy, age_round) %>% 
  summarize(prop_cens = sum(iptw[censor == 0]) / sum(iptw), 
            prop_uncens = sum(iptw[censor == 1]) / sum(iptw))

df_long <- pivot_longer(hist_cens_age_weighted_iptw, cols = c(prop_cens, prop_uncens), 
                        names_to = "censoring", 
                        values_to = "proportion")

df_long$age_round <- factor(df_long$age_round, levels = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110))

hist_cens_age_weighted_iptw <- ggplot(df_long, aes(x = factor(age_round), y = proportion, fill = censoring)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6a00a8", "#f56864"),
                    labels = c('censored', 'uncensored')) +
  facet_wrap(~ trt_dummy, labeller = as_labeller(c('0' = ref_trt, '1' = comp_trt))) + 
  labs(title = "Status by Age and Treatment (IPTW)",
       x = "Age",
       y = "Proportion",
       fill = "Status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 7.5))

hist_cens_age_weighted_iptw
ggsave("hist_cens_age_weighted_iptw.png", plot = hist_cens_age_weighted_iptw, width = 6, height = 3, units = "in", bg = 'white')

#### HISTOGRAM OF BASELINE DEPRESSION BY AGE GROUP ####

cohort <- readRDS(file = paste(path_intermediate_res, 'cohort_iptw.rds', sep = '/'))
cohort %<>% mutate(binary_age_grp = if_else(age_at_entry >= 65, 'patients ≥65', 'patients <65'))

depression_by_age <- cohort %>% # unweighted
  group_by(trt, binary_age_grp) %>% 
  summarize(prop_dep = sum(depression_base == 1) / n(), 
            prop_no_dep = sum(depression_base == 0) / n())

depression_by_age_iptw <- cohort %>% # weighted (iptw only)
  group_by(trt, binary_age_grp) %>% 
  summarize(prop_dep = sum(iptw[depression_base==1]) / sum(iptw), 
            prop_no_dep = sum(iptw[depression_base==0]) / sum(iptw))

depression_by_age_long <- pivot_longer(depression_by_age, cols = c(prop_dep, prop_no_dep), 
                        names_to = "depression_status", 
                        values_to = "proportion")

depression_by_age_iptw_long <- pivot_longer(depression_by_age_iptw, cols = c(prop_dep, prop_no_dep), 
                                                names_to = "depression_status", 
                                                values_to = "proportion")

hist_depression_by_age <- ggplot(depression_by_age_long, aes(x = factor(binary_age_grp), y = proportion, fill = depression_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6a00a8", "#f9a463"),
                    labels = c('Depression', 'No depression')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Patients with Baseline Depression by Treatment",
       x = "Age group",
       y = "Proportion",
       fill = "Baseline Depression") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
hist_depression_by_age

ggsave("hist_depression_by_age.png", plot = hist_depression_by_age, width = 6, height = 3, units = "in", bg = 'white')


hist_depression_by_age_iptw <- ggplot(depression_by_age_iptw_long, aes(x = factor(binary_age_grp), y = proportion, fill = depression_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("#6a00a8", "#f9a463"),
                    labels = c('Depression', 'No depression')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Patients with Baseline Depression by Treatment (IPTW)",
       x = "Age group",
       y = "Proportion",
       fill = "Baseline Depression") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
hist_depression_by_age_iptw

ggsave("hist_depression_by_age_iptw.png", plot = hist_depression_by_age_iptw, width = 6, height = 3, units = "in", bg = 'white')

#### HISTOGRAM OF BASELINE DEPRESSION OVER TIME BY TREATMENT ####

cohort_long <- readRDS(file = paste(path_intermediate_res, 'cohort_ipcw.rds', sep = '/'))
cohort_long %<>% mutate(combo_weight = iptw*ipcw_pl_nonlag)

dep_over_time <- cohort_long %>% 
  group_by(trt, dec) %>% 
  summarize(prop_dep = sum(depression_base == 1) / n(), 
            prop_no_dep = sum(depression_base ==  0) / n()
  )

dep_over_time_long <- pivot_longer(dep_over_time, cols = c(prop_dep, prop_no_dep), 
                        names_to = "depression_status", 
                        values_to = "proportion")

hist_dep_over_time <- ggplot(dep_over_time_long, aes(x = factor(dec), y = proportion, fill = depression_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#6a00a8', '#f9a463'),
                    labels = c('Depression', 'No depression')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Baseline Depression Over Time by Treatment",
       x = "Interval",
       y = "Proportion",
       fill = "Depression status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_dep_over_time
ggsave("hist_dep_over_time.png", plot = hist_dep_over_time, width = 6, height = 3, units = "in", bg = 'white')



dep_over_time_iptw <- cohort_long %>% # weighted with IPTW only
  group_by(trt, dec) %>%
  summarize(
    prop_dep = sum(iptw[depression_base == 1]) / sum(iptw),
    prop_no_dep = sum(iptw[depression_base == 0]) / sum(iptw)
  )

dep_over_time_iptw_long <- pivot_longer(dep_over_time_iptw, cols = c(prop_dep, prop_no_dep), 
                        names_to = "depression_status", 
                        values_to = "proportion")

hist_dep_over_time_iptw <- ggplot(dep_over_time_iptw_long, aes(x = factor(dec), y = proportion, fill = depression_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#6a00a8', '#f9a463'),
                    labels = c('Depression', 'No depression')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Baseline Depression Over Time by Treatment (IPTW)",
       x = "Interval",
       y = "Proportion",
       fill = "Depression status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_dep_over_time_iptw
ggsave("hist_dep_over_time_iptw.png", plot = hist_dep_over_time_iptw, width = 6, height = 3, units = "in", bg = 'white')



dep_over_time_weighted_combo <- cohort_long %>% # weighted with IPTW + IPCW
  group_by(trt, dec) %>%
  summarize(
    prop_dep = sum(combo_weight[depression_base == 1]) / sum(combo_weight),
    prop_no_dep = sum(combo_weight[depression_base == 0]) / sum(combo_weight)
  )

dep_over_time_weighted_combo_long <- pivot_longer(dep_over_time_weighted_combo, cols = c(prop_dep, prop_no_dep), 
                                        names_to = "depression_status", 
                                        values_to = "proportion")

hist_dep_over_time_weighted_combo <- ggplot(dep_over_time_weighted_combo_long, aes(x = factor(dec), y = proportion, fill = depression_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c('#6a00a8', '#f9a463'),
                    labels = c('Depression', 'No depression')) +
  facet_wrap(~ trt, labeller = as_labeller(c('snri' = 'SNRI', 'ssri' = 'SSRI'))) + 
  labs(title = "Baseline Depression Over Time by Treatment (IPTW + IPCW)",
       x = "Interval",
       y = "Proportion",
       fill = "Depression status") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hist_dep_over_time_weighted_combo
ggsave("hist_dep_over_time_weighted_combo.png", plot = hist_dep_over_time_weighted_combo, width = 6, height = 3, units = "in", bg = 'white')
