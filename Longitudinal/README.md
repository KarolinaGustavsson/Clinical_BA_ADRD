# Longitudinal analysis: biological age trajectories and ADRD risk

Scripts for the **repeated-measures analysis**: 44,805 AMORIS participants
with three or more pre-diagnosis biomarker measurements, contributing 194,897
measurement occasions, with up to 26 years of register follow-up. The question
is whether **baseline biological age** and the **rate of change in biological
age** are independently associated with dementia risk.

Every script begins with `source(here::here("config.R"))`. Outputs are written
under `OUTPUT_DIR` and read back through the constants in `config.R` §4.

## Prerequisites

Complete the [baseline analysis](../Baseline) first. These scripts read
`sampleIDs_baseline.csv` and `training_set_baseline.csv` from it.

---

## Workflow

### 1. `Calculating_BA_L.Rmd` — biological age at every measurement

- Load the full AMORIS extract, **one row per measurement occasion**
- Restrict to the individuals selected at baseline (`PATH_SAMPLE_IDS`)
- Apply KDM and PhenoAge using the **baseline training set**
  (`PATH_TRAINING_SET`), not a newly derived one
- Keep individuals with three or more measurements

**Outputs:** `results_dt` in the session — one row per measurement occasion —
plus age-distribution and correlation figures for the subsample.

---

### 2. `Table1_L.Rmd` — cohort characteristics

- Baseline (first) measurement per individual
- Repeated-measurement pattern: number of measurements per participant, span
  between first and last, mean interval, and time from the last measurement to
  censoring
- Characteristics by biological age group, reusing the median split written by
  `Table1_B.Rmd` (`PATH_TABLE1_PERID`) so the groups match Table 1 exactly
- Follow-up characteristics by ADRD status

**Outputs:** `table_kdm_longitud_summary_part1*.csv`,
`table_kdm_longitud_baseline_summary_part2*.csv`,
`followuptime_byADRD_l.csv`.

---

### 3. `BA_rate_L.Rmd` — baseline level versus rate of change

**Preparation.** Link SveDem, harmonise the dementia codes across
ascertainment sources, derive entry and exit ages. Drop measurements recorded
**after** an ADRD diagnosis and require three or more pre-diagnosis
measurements; this excludes 61 of 44,866 individuals, leaving 44,805. The
script reports how many measurement rows are dropped.

**Two-stage estimation.** For each individual, fit
`lm(biological age ~ chronological age)` with age centred at that individual's
first measurement, so the intercept is estimated biological age at entry and
the slope is the change in biological age per year of chronological age.

**Survival models.** Fit Cox models with slope and intercept entered jointly,
adjusted for sex, with ADRD as the outcome and chronological age as the
timescale. **Section 3.3 collapses to one row per individual before fitting,
and that is what the manuscript reports** — see the note below.

**Sensitivity.** Randomly remove one measurement per individual, re-estimate
the slopes and intercepts, and refit; repeated twice, for both KDM and
PhenoAge (Supplementary Tables 10–11).

**Outputs:** `combined_results_ADRD_Baseline_forL.csv`,
`table9_random_drop_sensitivity_RMA.csv`, `bioage_scale_df_table_RMA.csv`.

---

### 4. `Trajectories_binned_L.Rmd` — 5-year age bins and sex differences

- Chronological age in 5-year bins from 17 up to the oldest observed age;
  results are reported for the 13 bins in 17–82
- `lm(biological age ~ age + age_bin + age:age_bin)` — bin main effects **and**
  the interaction, so each bin is fitted with its own intercept and its own
  slope, with 17–22 as the reference
- Cluster-robust standard errors (`sandwich`), clustered on `sampleID`, for
  the repeated measurements
- Global Wald test for whether the slopes differ across the 13 reported bins
  (17–82 years, 12 degrees of freedom)
- Sex differences from a joint interaction model including age, bin, sex and
  their interactions, giving sex-specific slopes and their differences with
  pointwise intervals, for the whole cohort and for men and women separately
- Spaghetti plots with the fitted bin segments and density shading

**Outputs:** `trajectory_table_*.csv` (slope per bin),
`trajectory_table_wald_*.csv` (Wald test), sex-difference tables, and the
trajectory figures.

---

### 5. `Trajectories_GAM_L.Rmd` — smooth trajectories

- `gam(biological age ~ s(chronological age))` with `mgcv`
- Predictions across 17–82 with **cluster-robust confidence bands**: the
  covariance matrix is built from the penalized normal equations
  (`A = X'X + S_lambda`) paired with a meat matrix clustered on `sampleID`.
  `sandwich::vcovCL()` is deliberately not used — it pairs the penalized bread
  with an unpenalized meat, which is inconsistent for a penalized GAM and was
  verified on simulated clustered data to understate the standard error by
  roughly twentyfold. The helper `gam_vcov_cluster()` in section 1.4 carries
  the full reasoning.
- Spaghetti plots on disjoint random subsamples, with two-dimensional kernel
  density shading

**Outputs:** `GAM_ci_per5year.csv` (fitted values and intervals at 5-year
landmark ages), spaghetti and density figures.

---

## Data and variables

**Input.** The full AMORIS extract (`PATH_AMORIS_ALL`), one row per
measurement occasion, plus the SveDem linkage (`PATH_SVEDEM`).

**Key derived variables.** `kdm_advance_sd`, `phenoage_advance_sd`
(standardised age gaps); `intercept` and `slope` (per-individual, from the
two-stage fit); `age_bin`; `n_measurements`; `ADRD`, `Alzheimers`, `Vascular`,
`Other`. Full definitions in
[DATA_DICTIONARY.md](../DATA_DICTIONARY.md) §4.

---

## Notes

**1. Analytic grain — the thing to get right here.**
`results_dt` holds one row per *measurement occasion*. The trajectory
`intercept` and `slope` are *person-level* quantities, but they are merged
back onto every measurement row for that person (`BA_rate_L.Rmd` §3.2). A Cox
model fitted on that merged frame enters each participant into the partial
likelihood once per measurement, which inflates the effective sample size and
understates the standard errors.

`BA_rate_L.Rmd` therefore collapses to one row per individual before fitting
the models it reports, in §3.3:

```r
results_kdm_1row <- results_kdm %>%
  distinct(sampleID, .keep_all = TRUE)
```

Any new model whose exposure is a
person-level summary must be fitted on the de-duplicated frame.

The trajectory models in scripts 4 and 5 *are* fitted on measurement-level
rows, which is correct: the exposure there is chronological age at the
measurement, not a person-level summary. Within-person correlation is handled
by cluster-robust standard errors, not by de-duplication.

See [DATA_DICTIONARY.md](../DATA_DICTIONARY.md) §7.

**2. The training set must come from the baseline analysis.** Biological age is
calibrated, not measured. Regenerating the training set here puts the
longitudinal biological ages on a different scale from the baseline ones and
the two are no longer comparable.

**3. Health check-up subsets.** Optional sections restrict to samples taken at
health check-ups (`debiteringskod %in% c("F", "H")`) rather than in clinical
care.

**4. File encoding.** Output CSVs use semicolon delimiters, comma decimals and
Latin-1 encoding (`config.R` §6).

---

## Approximate run times

```
Baseline analysis
    |
    +--> sampleIDs_baseline.csv, training_set_baseline.csv
          |
     Calculating_BA_L.Rmd        20-30 min
     Table1_L.Rmd                 5-10 min
     BA_rate_L.Rmd               30-45 min  (per-individual fits + Cox models)
     Trajectories_binned_L.Rmd   15-20 min
     Trajectories_GAM_L.Rmd      10-15 min
```
