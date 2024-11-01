**1. Compare IPCW models to Adjust for Informative Censoring using CPRD Data** <br />
Studies on treatment efficacy using observational data can be analyzed with different approaches to answer distinct questions: what is the effect of being prescribed the treatment (intention-to-treat approach)? And what is the effect of starting and adhering to a treatment (as-treated approach)? While the as-treated approach can generate useful insights, it is vulnerable to a type of selection bias called informative censoring if patients who are excluded from the analysis for not adhering to their treatment differ systematically from those who remain in the analysis. <br />

This project compares methods of implementing inverse probability of censoring weights (IPCW), which create an artificial cohort in which there is no difference between patients who are included or excluded from the analysis due to discontinuing their treatment. Logistic regression models are fitted in intervals of follow-up to estimate weights throughout the study. This is different from Cox model or Kaplan-Meier approaches, that estimate weights every time a patient is censored (which can be computationally intensive in large data, especially when bootstrapping is used). <br />

**2. Generate Back-end Data for Data Visualization in Multi-Database Studies** <br />
Distributed data networks and multi-database studies offer new ways of conducting health research that is rapid, rigorous, and reproducible. However, the study populations contained in each database may be  heterogeneous in terms of event rates and confounding structures. This project used CPRD data on different cohorts to generate data describing the rate of events, the prevalence of confounders, and the association between confounders and the event over time. This data will be combined with results from other data sources to test a visualization tool using R Shiny. The interactive tool will display features of the populations in each database to help researchers identify key differences between them at a glance. <br />

Access to CPRD data is supported by ISAC protocol 24_004042. This work is conducted at the Lady Davis Institute of the Jewish General Hospital in Montreal. <br />

| Program  | Description |
| ------------- | ------------- |
| 1-cohort_creation  | Create cohorts of first-time users of two treatments being compared within the CPRD  |
| 2-censoring  | Define cohort exit (censoring) and date of treatment discontinuation (separate programs for different grace periods)  |
| 3-covariates  | Define covariates and the first occurrence of comorbidities  |
| 4-outcome  | Define the occurrence and timing of all-cause mortality or a clinical event (separate programs)  |
| 5-subgroups  | Split the data into subgroups of patients for subgroup analyses  |
| 6-iptw | Compute IPTW weights to adjust for confoudning |
| 7-ipcw  | Split data into follow-up intervals and compute different IPCW models to adjust for informative censoring  |
| 8-analyses  | Analyze incidence rate ratios and hazard ratios of the outcome between two treatment groups  |
| 9-bootstrap  | Replicate steps 6-8 1,000 times to obtain robust confidence intervals for incidence rates  |
| 10-create_plots  | Create tables and plots to communicate findings  |
| 11-prepare_visualization  | Generate back-end data to describe event rates and confounding structures to test a visualization tool  |
