**Compare IPCW models to Adjust for Informative Censoring using CPRD Data** <br />
Studies on treatment efficacy using observational data can be analyzed with different approaches to answer distinct questions: what is the effect of being prescribed the treatment (intention-to-treat approach)? And what is the effect of starting and adhering to a treatment (as-treated approach)? While the as-treated approach can generate useful insights, it is vulnerable to a type of selection bias called informative censoring if patients who are censored from the analysis for not adhering to their treatment differ systematically from those who remain in the analysis. Methods to address informative censoring are well-established, but the nuances of their implementation are less documented <br />

This project compares approaches to implementing inverse probability of censoring weights (IPCW), a weighting method that creates an artificial cohort in which there is no difference between patients who are censored from the analysis due to discontinuing their treatment and those who remain in the study. Logistic regression models are fitted in intervals of follow-up to estimate weights throughout the study, and results from a) lagged weights and b) non-lagged weights are compared. This is different from Cox model or Kaplan-Meier approaches, that estimate weights every time a patient is censored (which can be computationally intensive in large data, especially when bootstrapping is used). <br />

Access to CPRD data is supported by ISAC protocol 24_004042. This work is conducted at the Lady Davis Institute of the Jewish General Hospital in Montreal. <br />

| Program  | Description |
| ------------- | ------------- |
| 1-cohort_creation  | Create cohorts of first-time users of two treatments being compared within the CPRD  |
| 2-censoring  | Define cohort exit (censoring) and date of treatment discontinuation (separate programs for different grace periods)  |
| 3-covariates  | Define covariates and the first occurrence of comorbidities  |
| 4-outcome  | Define the occurrence and timing of all-cause mortality and follow-up  |
| 5-subgroups  | Split the data into subgroups of patients for subgroup analyses  |
| 6-iptw | Compute IPTW weights to adjust for confoudning |
| 7-ipcw  | Split data into follow-up intervals and compute different IPCW models to adjust for informative censoring  (original model: all intervals are included; no deaths model: intervals in which patients died are removed when fitting IPCW)|
| 8-analyses  | Analyze incidence rate ratios and hazard ratios of the outcome between two treatment groups  |
| 9-bootstrap  | Repeat steps 6-8 1,000 times to obtain robust confidence intervals for incidence rates  |
| 10-plot_results  | Create tables and plots to communicate findings  |

Packages used for this project:

R Packages Used for Analysis
1.
Canty A, Ripley B. boot: Bootstrap Functions (Originally by Angelo Canty for S). Published online April 8, 1999:1.3-31. doi:10.32614/CRAN.package.boot
2.
Greifer N. cobalt: Covariate Balance Tables and Plots. Published online April 2, 2024. Accessed September 16, 2024. https://cran.r-project.org/web/packages/cobalt/index.html
3.
Davison AC, Hinkley DV. Bootstrap Methods and Their Application. Cambridge University Press; 1997. doi:10.1017/CBO9780511802843
4.
Greifer N. WeightIt: Weighting for Covariate Balance in Observational Studies. Published online August 24, 2024. Accessed September 16, 2024. https://cran.r-project.org/web/packages/WeightIt/index.html
5.
Wickham H, Averick M, Bryan J, et al. Welcome to the Tidyverse. Journal of Open Source Software. 2019;4(43):1686. doi:10.21105/joss.01686
6.
Huling J, Bates D, Eddelbuettel D, Francois R, Qiu Y. fastglm: Fast and Stable Fitting of Generalized Linear Models using “RcppEigen.” Published online May 23, 2022. Accessed September 16, 2024. https://cran.r-project.org/web/packages/fastglm/index.html
7.
Wickham H. The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software. 2011;40:1-29. doi:10.18637/jss.v040.i01
8.
Kassambara A, Kosinski M, Biecek P, Fabian S. survminer: Drawing Survival Curves using “ggplot2.” Published online March 9, 2021. Accessed September 16, 2024. https://cran.r-project.org/web/packages/survminer/index.html
9.
Therneau TM, until 2009) TL (original S >R port and R maintainer, Elizabeth A, Cynthia C. survival: Survival Analysis. Published online June 5, 2024. Accessed September 16, 2024. https://cran.r-project.org/web/packages/survival/index.html
10.
Rich B. table1: Tables of Descriptive Statistics in HTML. Published online January 6, 2023. Accessed September 16, 2024. https://cran.r-project.org/web/packages/table1/index.html
11.
Barrett M. tidysmd: Tidy Standardized Mean Differences. Published online May 26, 2023. Accessed September 16, 2024. https://cran.r-project.org/web/packages/tidysmd/index.html
12.
Barrett T, Dowle M, Srinivasan A, et al. data.table: Extension of “data.frame.” Published online August 27, 2024. Accessed September 16, 2024. https://cran.r-project.org/web/packages/data.table/index.html
