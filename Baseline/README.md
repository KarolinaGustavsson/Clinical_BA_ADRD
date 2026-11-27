# Baseline analysis: biological age and ADRD risk

Scripts for the **baseline (cross-sectional) analysis** of the association
between biological age and risk of Alzheimer's disease and related dementias
(ADRD) in 261,383 AMORIS individuals.

Every script begins with `source(here::here("config.R"))` and takes its paths,
output location and analysis parameters from there. Nothing hardcodes a
location. Outputs are written under `OUTPUT_DIR` (`./output` unless
`STUDY1_OUTPUT_DIR` says otherwise) and read back through the constants in
`config.R` §4.

---

## Workflow

Run in this order. Scripts 3–8 all need the objects script 2 leaves in the
session, and 4–7 need the outcome definitions from script 3.

### 1. `Selecting_biomarkers_B.Rmd` — biomarker screening

Screens the 17 extracted candidates on their correlation with chronological
age.

- Density plot for every candidate biomarker
- Spearman and Pearson correlation with chronological age, before and after
  outlier removal, overall and separately in men and women
- Candidates whose absolute correlation falls below 0.1 in both coefficients
  are dropped; 12 are retained

**Outputs:** correlation tables (`ADRDamoris_all_bm*.csv`,
`raw_ADRDamoris_all_bm*.csv`, `bm_all_new_training_set*.csv`), biomarker
density figure.

---

### 2. `Calculating_BA_B.Rmd` — cohort setup and biological age

- Load the baseline study population
- Define the 12 biomarkers and remove values beyond ±5 s.d. (`OUTLIER_SD`)
- Build the age-uniform training set (10,080 individuals, uniform age within
  each sex) and split off the test set
- Apply the Klemera–Doubal Method and PhenoAge from the `BioAge` package
- Derive the age gaps and their standardised versions
- Print the fitted KDM and PhenoAge model weights (Supplementary Table 12)

**Outputs:** `to1po2` (one row per individual) in the session, plus
`bioage_scale_df_table.csv`, `sampleIDs_baseline.csv` and
`training_set_baseline.csv`. **The last two are what the whole `Longitudinal/`
folder depends on** — see the note at the end of this file.

An optional section restricts the cohort to individuals who fasted at
baseline, for the fasting sensitivity analysis.

---

### 3. `CoxPH_B.Rmd` — outcome definitions and the main Cox models

This script defines the outcomes and covariates that scripts 4–7 all reuse, so
it cannot be skipped.

- Link SveDem and derive the ADRD indicator and the three subtypes from the
  Patient Register codes, anti-dementia dispensations and SveDem, assigning
  each case to exactly one subtype
- Derive `income_bins` and `NewUtbNiva`, and the Cox entry and exit ages
- Fit, for ADRD and each subtype, a sex-adjusted model and a multivariable
  model additionally adjusted for income quartile and education level, with
  chronological age as the underlying timescale
- Test the proportional hazards assumption with `cox.zph()`
- Sensitivity analyses: restricted to fasting individuals; with sampling date
  added as a covariate; and with a restricted AD definition (AD*, section
  2.2.1): early- and late-onset AD in SveDem and AD-specific ICD-10 codes
  outside SveDem (F00 and G30, excluding F00.2 and G30.8), with mixed dementia,
  medication-only and pre-ICD-10 cases censored at diagnosis
  (Supplementary Figs. 18–19)

**Outputs:** `combined_results_ADRD_main.csv` (main models;
`combined_results_ADRD_main_sens_provtagningsdatum.csv` when sampling date is
added for the sensitivity analysis),
`combined_results_Alzheimers_s_withoutG308.csv` (AD*),
`combined_results_SveDem_main.csv`.

---

### 4. `CoxPH_SveDem_censored_B.Rmd` — adjudicated SveDem subanalysis

The same models restricted to clinically adjudicated SveDem diagnoses. ADRD
cases ascertained **outside** SveDem are right-censored at their diagnosis
date rather than excluded, so they contribute follow-up up to the point at
which their outcome status becomes unknown.

**Outputs:** `combined_results_SveDem_rightcensored.csv`.

---

### 5. `CoxPH_subgroups_B.Rmd` — stratified models

- **By baseline age:** median split at 45 years
- **By sex:** men and women separately
- The Cox models from script 3 refitted within each stratum, with a test for
  difference between strata

**Outputs:** `combined_df_agegroup_45_main.csv`, `combined_sex_main.csv`.

---

### 6. `Table1_B.Rmd` — descriptive characteristics

- Table 1, stratified by the biological age groups, with means and standard
  deviations for continuous variables, N (%) for categorical variables,
  p-values and **standardised mean differences**
- The median split of the KDM age gap into "Low Age Gap" and "High Age Gap",
  written out per individual for reuse in the longitudinal tables and the
  cumulative incidence analysis
- Follow-up characteristics by ADRD status

**Outputs:** `table_kdm_agegap_stratified_two.csv`,
`followuptime_ADRD_b_smd.csv`, and
`table1_kdmgroups_baseline_perID.csv` — the last is **individual-level** and
must stay inside the secure environment.

---

### 7. `Fine_Gray/fg_project1_ADRD.R` — competing risks

Fine-Gray subdistribution hazard models with death as the competing event,
producing Supplementary Table 9. Run from a terminal rather than from
RStudio: the `finegray()` risk-set expansion is memory-hungry on a cohort of
this size, and results are appended to CSV one model at a time so a run can be
restarted after a crash.

```
Rscript fg_project1_ADRD.R inspect     # data diagnostics; run first
Rscript fg_project1_ADRD.R main        # ADRD and the three subtypes
Rscript fg_project1_ADRD.R svedem      # adjudicated SveDem outcomes
```

See [`Fine_Gray/README.md`](Fine_Gray/README.md) for the inputs, the
modelling choices and the fallbacks if memory runs out.

---

### 8. `CoxPH_lag_B.Rmd` — lag analysis

The Fig. 5a models with follow-up starting 5 or 10 years after baseline, to
reduce the influence of preclinical disease. Everyone alive, free of ADRD and
under follow-up at that point enters the risk set at baseline age plus the lag
(delayed entry, age as timescale), so no one is excluded on the basis of a
later event. Section 1 repeats the data preparation of `CoxPH_B.Rmd`; a lag of
0 is fitted as a check and must reproduce `combined_results_ADRD_main.csv`
exactly. The 10-year lag is Supplementary Figs. 16–17.

**Outputs:** `combined_results_ADRD_main_lag5y.csv`,
`combined_results_ADRD_main_lag10y.csv`, `lag_counts_N_events.csv`.

---

### 9. `Cumulative_incidence/` — absolute risk with death as a competing event

Aalen-Johansen cumulative incidence of ADRD and of death from age 65,
by biological age group (the median split from `Table1_B.Rmd`), on the same
analytic sample, event definitions and age timescale as the Cox and Fine-Gray
models. `cif_sexstd_project1.R` estimates the curves separately in men and
women and standardises them to the sex distribution of the sample
(Supplementary Figs. 26–27):

```
Rscript cif_sexstd_project1.R sex      # sex-standardised curves -> ./sexstd/
Rscript plot_cif_results.R sexstd      # figures from the saved curves
```

`fg_project1.R` in this folder is a later version of the script in
`Fine_Gray/`; its `cif` mode gives the unstandardised curves. The published
Supplementary Table 9 was produced by `Fine_Gray/fg_project1_ADRD.R`. See
[`Cumulative_incidence/README.md`](Cumulative_incidence/README.md).

---

### 10. `Forest_plots_B.Rmd` — forest plots

Draws the forest plots in Fig. 5 and Supplementary Figs. 14–25 and 29–30
from the results tables written by scripts 3–5 and 8 and by
`Longitudinal/BA_rate_L.Rmd`, so run it after that script.
A table at the top of the notebook maps each figure to its input file.

**Outputs:** PDF figures in `fig_path()`.

---

## Data and variables

**Inputs.** The baseline study population (`PATH_STUDYPOP`) and the SveDem
linkage (`PATH_SVEDEM`), both configured in `config.R`.

**The 12 retained biomarkers.** `S_Alb` (albumin), `S_Krea` (creatinine), `TC`
(total cholesterol), `TG` (triglycerides), `S_K` (potassium), `S_Urea` (urea),
`S_LD` (lactate dehydrogenase), `S_Ca` (calcium), `S_Urat` (urate),
`S_FAMN` (fructosamine), `S_Hapt` (haptoglobin), `fS_Gluk` (fasting glucose).

**Key derived variables.** `kdm`, `phenoage` (biological age, years);
`kdm_advance`, `phenoage_advance` (age gaps, years); `kdm_advance_sd`,
`phenoage_advance_sd` (standardised age gaps — the exposure in every model);
`ADRD`, `Alzheimers`, `Vascular`, `Other` (outcomes, 0/1); `status` (death,
the competing event).

Full definitions, including the coding of the covariates and the SveDem
indicators, are in [DATA_DICTIONARY.md](../DATA_DICTIONARY.md) §4.

---

## Notes

**Everything here is person-level.** `to1po2` holds one row per individual.
The longitudinal folder is where the measurement-level data live, and where
mixing the two grains becomes a hazard — see
[DATA_DICTIONARY.md](../DATA_DICTIONARY.md) §7.

**The handoff to the longitudinal analysis.** `Calculating_BA_B.Rmd` writes
`sampleIDs_baseline.csv` and `training_set_baseline.csv`. The training set
**must be reused, not regenerated**: KDM and PhenoAge are calibrated on it, so
a different training set puts baseline and longitudinal biological ages on
different scales.

**File encoding.** Output CSVs use semicolon delimiters, comma decimals and
Latin-1 encoding (`config.R` §6), and are read back with the same settings.
