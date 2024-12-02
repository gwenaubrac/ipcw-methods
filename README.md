**Compare IPCW models to Adjust for Informative Censoring using CPRD Data** <br />
Studies on treatment efficacy using observational data can be analyzed with different approaches to answer distinct questions: what is the effect of being prescribed the treatment (intention-to-treat approach)? And what is the effect of starting and adhering to a treatment (as-treated approach)? While the as-treated approach can generate useful insights, it is vulnerable to a type of selection bias called informative censoring if patients who are censored from the analysis for not adhering to their treatment differ systematically from those who remain in the analysis. Methods to address informative censoring are well-established, but the nuances of their implementation are less documented <br />

This project compares approaches to implementing inverse probability of censoring weights (IPCW), a weighting method that creates an artificial cohort in which there is no difference between patients who are censored from the analysis due to discontinuing their treatment and those who remain in the study. Logistic regression models are fitted in intervals of follow-up to estimate weights throughout the study, and results from a) lagged weights and b) non-lagged weights are compared. This is different from Cox model or Kaplan-Meier approaches, that estimate weights every time a patient is censored (which can be computationally intensive in large data, especially when bootstrapping is used). <br />

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
