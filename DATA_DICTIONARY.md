# Data dictionary and workflow

How the data move between scripts, and what every variable means.

The AMORIS data are restricted and cannot be redistributed, so the pipeline
cannot be executed from this repository alone (see [README.md](README.md)).
This document is written so that the analysis can be read, checked and
re-implemented without running it, and so that a researcher who *does* hold an
approved AMORIS extract knows exactly which inputs each script expects.

**Contents**

1. [Execution order](#1-execution-order)
2. [How the data move between scripts](#2-how-the-data-move-between-scripts)
3. [Source data files](#3-source-data-files)
4. [Variables](#4-variables)
5. [Intermediate files](#5-intermediate-files)
6. [Key in-memory objects](#6-key-in-memory-objects)
7. [Analytic grain — read before modifying](#7-analytic-grain--read-before-modifying)
8. [Which script produces which published result](#8-which-script-produces-which-published-result)
9. [Data access](#9-data-access)

---

## 1. Execution order

`config.R` is sourced by every script and is the only file that needs editing
for a new environment. Scripts then run in this order; each depends on objects
or files produced by the one before.

```
config.R                               edit once (or set AMORIS_DATA_ROOT)
   |
Baseline/
   1. Selecting_biomarkers_B.Rmd       biomarker screening
   2. Calculating_BA_B.Rmd             cohort setup, KDM + PhenoAge
   3. CoxPH_B.Rmd                      outcome definitions, main Cox models
   4. CoxPH_SveDem_censored_B.Rmd      SveDem subanalysis, non-SveDem censored
   5. CoxPH_subgroups_B.Rmd            age- and sex-stratified models
   6. Table1_B.Rmd                     Table 1, biological age groups
   7. Fine_Gray/fg_project1_ADRD.R     competing-risk models (Supp. Table 9)
   |
Longitudinal/
   8. Calculating_BA_L.Rmd             biological age at every measurement
   9. Table1_L.Rmd                     repeated-measures descriptives
  10. BA_rate_L.Rmd                    intercept vs. rate of change
  11. Trajectories_binned_L.Rmd        5-year age-bin trajectories, sex differences
  12. Trajectories_GAM_L.Rmd           smooth trajectories
   |
Review/                               analyses cited in the response to
                                      reviewers only — see Review/README.md
```

Steps 1–2 can be run independently of 3–7, but step 2 writes the sample and
training-set files that the whole `Longitudinal/` folder needs, so it cannot be
skipped.

---

## 2. How the data move between scripts

Two things are passed between scripts: **files on disk**, and **objects left in
the R session**. The second kind is easy to miss, because the notebooks were
written to be run in sequence in one session.

### 2.1 Files on disk

Written under `OUTPUT_DIR` (see `config.R` §3) and referenced everywhere by the
constants in `config.R` §4, never by a literal path.

| Produced by | File (config.R constant) | Consumed by | Why it must be reused, not regenerated |
|---|---|---|---|
| `Calculating_BA_B.Rmd` §2.7 | `sampleIDs_baseline.csv` (`PATH_SAMPLE_IDS`) | `Calculating_BA_L.Rmd` §1.2 | Fixes the analytic cohort, so baseline and longitudinal describe the same people |
| `Calculating_BA_B.Rmd` §2.8 | `training_set_baseline.csv` (`PATH_TRAINING_SET`) | `Calculating_BA_L.Rmd` §1.2 | **Critical.** KDM and PhenoAge are only comparable across the two analyses if both are calibrated on the *same* training set. Regenerating it puts the two sets of biological ages on different scales. |
| `Calculating_BA_B.Rmd` §2.5 | `bioage_scale_df_table.csv` (`PATH_BIOAGE_SCALE`) | reference; the s.d. used to convert age gaps to years | Reporting units: 1 s.d. of the KDM age gap = 6.5 years |
| `Table1_B.Rmd` §1.x | `table1_kdmgroups_baseline_perID.csv` (`PATH_TABLE1_PERID`) | `Table1_L.Rmd` §1.2, `Review/Baseline/Cumulative_incidence/` | Keeps "biologically younger/older" identical across Table 1, Table 1 (longitudinal) and the cumulative incidence curves |

`table1_kdmgroups_baseline_perID.csv` holds **one row per participant**
(`sampleID`, `KDM_group`). It is individual-level data: it stays inside the
secure environment and is excluded by `.gitignore`.

### 2.2 Objects carried in the R session

| Object | Created in | Grain | Used by |
|---|---|---|---|
| `amoris` | `Calculating_BA_B.Rmd` §1.2 | one row per individual | the rest of that notebook |
| `to1po2` | `Calculating_BA_B.Rmd` §2 | **one row per individual** | `CoxPH_B.Rmd`, `CoxPH_SveDem_censored_B.Rmd`, `CoxPH_subgroups_B.Rmd`, `Table1_B.Rmd` |
| `to1po2` (with outcomes) | `CoxPH_B.Rmd` §1 | one row per individual | `CoxPH_subgroups_B.Rmd`, `Table1_B.Rmd`, the Fine-Gray scripts |
| `n_amoris` | `Calculating_BA_L.Rmd` §1.2 | one row per **measurement occasion** | that notebook |
| `results_dt` | `Calculating_BA_L.Rmd` | **one row per measurement occasion** | `Table1_L.Rmd`, `BA_rate_L.Rmd`, `Trajectories_binned_L.Rmd`, `Trajectories_GAM_L.Rmd`, `Review/Longitudinal/*` |
| `results_kdm`, `results_phe` | `BA_rate_L.Rmd` §3 | **one row per individual** (slope + intercept) | the two-stage Cox models in that notebook |

The Fine-Gray scripts are the exception: they are run from a terminal rather
than inside a notebook, and load a saved workspace holding the prepared
`to1po2` and `to1po2_SveDem_D` instead of inheriting them from the session.
See `Baseline/Fine_Gray/README.md`.

### 2.3 Why the order matters

- Biological age is **calibrated**, not measured. The KDM and PhenoAge
  coefficients are estimated once, in `Calculating_BA_B.Rmd`, on an
  age-uniform training subset of 10,080 individuals, and then applied
  everywhere else. Any script that recalibrates produces numbers on a
  different scale.
- Outcome definitions live in `CoxPH_B.Rmd` §1, not in the source data. ADRD,
  the three subtypes and the SveDem indicators are all derived there, so
  nothing downstream can be run before it.
- `BA_rate_L.Rmd` collapses measurement-level rows to person-level slopes and
  intercepts. Everything after that point in that notebook is person-level —
  see §7.

---

## 3. Source data files

Three restricted SAS files. Their directory is set once in `config.R` §1, via
the `AMORIS_DATA_ROOT` environment variable or by editing `DATA_ROOT`.

### `StudyPop.sas7bdat` — baseline analytic population (`PATH_STUDYPOP`)

One row per individual, taken at the first time point at which all candidate
biomarkers were available. Carries all register linkages except SveDem.

### `AMORIS_20240417.sas7bdat` — full extract (`PATH_AMORIS_ALL`)

**One row per measurement occasion**, so an individual appears as many times as
they were sampled. The source of the repeated-measures subsample. The
distinction between this file and `StudyPop` is the single most important thing
to keep straight in this pipeline; see §7.

### `SveDem.sas7bdat` — dementia quality register linkage (`PATH_SVEDEM`)

Clinically adjudicated dementia diagnoses from SveDem, the Swedish registry for
cognitive/dementia disorders, with finer subtype resolution than the National
Patient Register. Used for the adjudicated subtype subanalysis.

---

## 4. Variables

Variable names are as they appear in the AMORIS extract; many are Swedish
abbreviations. Register-sourced variables carry missingness from the register
linkage.

### 4.1 Identifiers and dates

| Variable | Type | Description |
|---|---|---|
| `Id` / `sampleID` | num | Study-specific pseudonymised identifier. **Not** a personal identity number. Join key throughout. |
| `AterPNR` | num | Re-used identity number flag (data-quality marker from the register linkage) |
| `PersonNrAndr` | num | Changed identity number flag (data-quality marker) |
| `FODDATn` | date | Date of birth |
| `Provtagningsdatum` | date | Date of blood sampling for a given measurement occasion |
| `d1` | date | Sampling date as recorded in the register extract |
| `d2` | date | Date of last follow-up from the register extract |
| `firstDate` | date | First date of follow-up; equals the first sampling date |
| `lastDate` | date | Last date of follow-up: ADRD diagnosis, death, emigration, or 2020-12-31, whichever came first |
| `firstRTB` | date | First record in the population register |
| `lastIMMIGRATION` | date | Most recent immigration date |
| `DODSDATn` | date | Date of death |
| `ULORSAK` | chr | Underlying cause of death (ICD) |
| `INDATUM` | date | Date of first dementia diagnosis in the National Patient Register |
| `FDATUM` | date | Date of first anti-dementia drug dispensation |
| `Event` | num | Register code for how follow-up ended (diagnosis before sampling, diagnosis after sampling, death, emigration, lost to follow-up, population-register follow-up) |

### 4.2 Age

| Variable | Type | Description |
|---|---|---|
| `Alder` | chr | Chronological age as delivered in the extract |
| `age` | num | Chronological age in years at a given measurement occasion, derived from `FODDATn` and `Provtagningsdatum` |
| `age_w` | num | Age in whole years, used to build the age-uniform training set |
| `AgeAtFirstDate` | num | Age at entry: `(firstDate − FODDATn) / 365.25`. **Cox entry time.** |
| `AgeAtLastDate` | num | Age at exit: `(lastDate − FODDATn) / 365.25`. **Cox exit time.** |
| `age_of_diagnosis` | num | Age at ADRD diagnosis |
| `age_first_meas` | num | Age at an individual's first measurement; the centring point for the trajectory intercept |
| `birthcohort` | chr | Birth-year band, used descriptively in Table 1 |
| `cal_year`, `kalenderar` | num | Calendar year of measurement, derived from `Provtagningsdatum`. Used only in the calendar-year sensitivity analysis (`Review/Longitudinal/`) |

Chronological age is the underlying timescale in every Cox model
(`Surv(AgeAtFirstDate, AgeAtLastDate, event)`), which adjusts for chronological
age at baseline by construction.

### 4.3 Biomarkers

Seventeen candidates were extracted; twelve were retained after screening on
absolute correlation with chronological age within each sex
(`Selecting_biomarkers_B.Rmd`). Values more than 5 s.d. from the sample mean
are excluded (`OUTLIER_SD` in `config.R`).

**Retained in KDM and PhenoAge (12):**

| Variable | Description | Unit |
|---|---|---|
| `S_Alb` | Serum albumin | g/L |
| `S_Krea` | Serum creatinine | µmol/L |
| `TC` | Serum total cholesterol | mmol/L |
| `TG` | Fasting serum triglycerides | mmol/L |
| `S_K` | Serum potassium | mmol/L |
| `S_Urea` | Serum urea (blood urea nitrogen) | mmol/L |
| `S_LD` | Serum lactate dehydrogenase | µkat/L |
| `S_Ca` | Serum calcium | mmol/L |
| `S_Urat` | Serum urate (uric acid) | µmol/L |
| `S_FAMN` | Serum fructosamine | mmol/L |
| `S_Hapt` | Serum haptoglobin | g/L |
| `fS_Gluk` | Fasting serum glucose | mmol/L |

**Extracted but not retained:**

| Variable | Description | Unit |
|---|---|---|
| `S_Alp` | Serum alkaline phosphatase | µkat/L |
| `S_P` | Serum phosphate | mmol/L |
| `fS_Jaern` | Fasting serum iron | µmol/L |
| `fS_TIBC` | Total iron-binding capacity | µmol/L |
| `Fe_maet` | Iron saturation | % |

### 4.4 Sampling context

| Variable | Type | Description |
|---|---|---|
| `debiteringskod` | chr | Billing code for the sampling occasion. `"F"` and `"H"` indicate health check-ups; `"S"` indicates clinical care |
| `fastestatus` | chr | Fasting status as recorded: `"F"` overnight fast, `"N4"` 4-hour fast, `"N"`/`"S"` non-fasting |
| `fasting_status` | chr | Derived, human-readable version of `fastestatus` used in Table 1 |
| `BMI` | num | Body mass index (register-sourced; substantial missingness) |
| `DateBMI_extra` | date | Alternative date for the BMI record |

### 4.5 Demographic and socioeconomic covariates

All register-sourced and therefore subject to linkage missingness.

| Variable | Type | Description |
|---|---|---|
| `Kon` | chr | Sex as delivered, `"1"` or `"2"` |
| `gender` / `new_gender` | num | Recoded sex used in modelling: 1 = men, 2 = women |
| `Sex` | chr | Labelled sex (`"Men"`, `"Women"`) for tables |
| `FodGrEg4` | chr | Region or country of birth (`"00"` = Sweden) |
| `LISA` | date | Date of the LISA record supplying education |
| `Sun2000niva_old` | chr | Highest attained education level in LISA before diagnosis (SUN 2000 classification) |
| `UtbNiva` | chr | Education level in categories, from the population and housing census |
| `FoB_UtbNiva` | date | Date of the census education record |
| `NewUtbNiva` | num | **Derived.** Education level used as a covariate: `Sun2000niva_old` where numeric and not `"*"`, otherwise `UtbNiva`, otherwise `NA`. Kept numeric in all models. Missing for about 22% of the baseline cohort, concentrated in the youngest birth cohorts — characterised in the response to reviewers |
| `Education_level` | chr | Labelled version of `NewUtbNiva` for tables |
| `DispInk` | num | Disposable income |
| `DispInkFam` | num | Disposable family income |
| `DispInkKE` | num | Disposable family income per consumption unit |
| `ArbInk` | num | Last recorded employment income before diagnosis |
| `FoB_ArbInk` | date | Date of that employment-income record |
| `DispInk_quartile`, `ArbInk_quartile` | factor | **Derived.** Income quartiles |
| `income_bins` | num | **Derived.** Income quartile used as a covariate: `DispInk_quartile` where available, otherwise `ArbInk_quartile`, otherwise `NA`. Kept numeric in all models |
| `Income_quartile` | chr | Labelled version of `income_bins` for tables |

### 4.6 Biological age

Computed with `kdm_calc()` and `phenoage_calc()` from the **BioAge** package
(Kwon & Belsky, *GeroScience* 2021), calibrated on the age-uniform training
subset written by `Calculating_BA_B.Rmd`.

| Variable | Type | Description |
|---|---|---|
| `kdm` | num | Klemera–Doubal Method biological age, years |
| `phenoage` | num | PhenoAge biological age, years |
| `kdm_advance` | num | KDM age gap: `kdm − age`, years |
| `phenoage_advance` | num | PhenoAge age gap: `phenoage − age`, years |
| `kdm_advance_sd` | num | **The exposure in every Cox model.** KDM age gap standardised to s.d. units; 1 s.d. = 6.5 years |
| `phenoage_advance_sd` | num | PhenoAge age gap, standardised |
| `kdm_sd`, `phenoage_sd` | num | Standardised biological age (not the gap) |
| `kdm_a_group` | factor | **Derived.** Median split of `kdm_advance`: `"Low Age Gap"` / `"High Age Gap"`, i.e. biologically younger vs. older. Written out per individual as `KDM_group` |
| `KDM_group` | chr | `kdm_a_group` as carried in `table1_kdmgroups_baseline_perID.csv` |
| `kdm_plot`, `phenoage_plot` | num | Plotting copies, unused in models |

Biomarkers are **not** equally weighted. The fitted KDM and PhenoAge
coefficients are reported in Supplementary Table 12 and printed by
`Calculating_BA_B.Rmd`.

### 4.7 Outcomes

Derived in `CoxPH_B.Rmd` §1 from three ascertainment sources: the National
Patient Register (`Dementia`, ICD codes), anti-dementia drug dispensations
(`ATC`), and SveDem. Each case is assigned to exactly one subtype; ICD codes
indicating mixed Alzheimer's and vascular pathology are assigned to
Alzheimer's. The code lists are in Supplementary Tables 13–14.

| Variable | Type | Description |
|---|---|---|
| `Dementia` | chr | Dementia diagnosis in the National Patient Register, ICD code |
| `ATC` | chr | First anti-dementia drug dispensed |
| `ADRD` | 0/1 | Any ADRD: a Patient Register code, an anti-dementia dispensation, or a SveDem record |
| `subtype` | chr | **Derived.** `"Alzheimers"`, `"Vascular"` or `"Other"`; the single assignment from which the three indicators below are built |
| `Alzheimers` | 0/1 | Alzheimer's disease, including mixed AD/vascular |
| `Vascular` | 0/1 | Vascular dementia |
| `Other` | 0/1 | All other dementias |
| `status` | 0/1 | Death during follow-up (`1` if `DODSDATn` is present). The competing event in the Fine-Gray and cumulative incidence analyses |
| `Mortality` | chr | Labelled version of `status` for tables |
| `FollowUpTime` | num | `AgeAtLastDate − AgeAtFirstDate`, years |
| `first_diag_date` | date | Earliest diagnosis date across the ascertainment sources |

**SveDem-restricted indicators** (`S_` prefix), from `DIAGNOSIS_BR` in the
SveDem linkage:

| Variable | SveDem code | Description |
|---|---|---|
| `S_ADRD` | 1–9 | Any adjudicated dementia diagnosis |
| `S_e_alz` | 1 | Early-onset Alzheimer's disease |
| `S_l_alz` | 2 | Late-onset Alzheimer's disease |
| `S_alz_vas` | 3 | Mixed Alzheimer's and vascular dementia |
| `S_vas` | 4 | Vascular dementia |
| `S_lew` | 5 | Lewy body dementia |
| `S_fro` | 6 | Frontotemporal dementia |
| `S_par` | 7 | Parkinson's disease dementia |
| `S_uns` | 8 | Unspecified dementia |
| `S_oth` | 9 | Other dementia |
| `S_INDATUM` | — | `DIAGNOSIS_DATE_BR`, date of the adjudicated diagnosis |

In the SveDem subanalysis (`CoxPH_SveDem_censored_B.Rmd`) only
SveDem-recorded diagnoses count as events; ADRD cases ascertained outside
SveDem are **right-censored at their diagnosis date** rather than excluded.

### 4.8 SveDem clinical variables

The SveDem extract carries the full basic dementia work-up. Only
`DIAGNOSIS_BR` and `DIAGNOSIS_DATE_BR` enter the analyses; the remainder are
available in the linkage and listed here for completeness.

| Group | Variables |
|---|---|
| Registration | `CENTREID_BR`, `CENTRETYPEID_BR`, `BASAL_EXAMINATION_UNIT`, `BASAL_EXAMINATION_SVEDEM`, `LAN_CENTRE_BR`, `KOMMUN_CENTRE_BR`, `LAN_PATIENT`, `KOMMUN_PATIENT` |
| Work-up dates | `REFERRAL_DATE`, `INITIATE_DATE`, `D_BASAL_UTREDNING`, `D_ANTAL_UNDERSOKNINGAR_BR` |
| Anthropometry | `WEIGHT_BR`, `LENGTH_BR`, `D_BMI_BR` |
| Living situation | `RECIDENCY_BR`, `CORESIDENT_BR`, `OCCUPATION_STATUS_BR`, `DAYCARE_TODAY_BR`, `HOMECARE_TODAY_BR`, `HOMEHEALTHCARE_TODAY_BR` |
| Family history | `RELATIVES_COGNITION_BR`, `FAMILY_DEGREE1`, `FAMILY_DEGREE2` |
| Cognitive tests | `MMT_BR`, `MMT_VALUE_BR`, `MOCA_BR`, `MOCA_VALUE_BR`, `RUDAS_BR`, `RUDAS_VALUE_BR`, `COGNITION_BR`, `COGNITIVE_EXTENDED_BR` |
| Investigations | `BLOODCHEM_BR`, `LP_BR`, `CT_BR`, `MR_BR`, `SPECT_BR`, `EEG_BR` |
| Assessments | `FUNTIONAL_BR`, `MOTORIC_BR`, `LINGUAL_BR` |
| Diagnosis | `DIAGNOS_BR`, `DIAGNOSIS_ICD_BR`, `DIAGNOS_G_BR`, `DIAGNOSIS_DATE_BR` |
| Medication at diagnosis | `CHOLINESTERASE_INHIBITOR_BR`, `NMDA_ANTAGONIST_BR`, `VASCULAR_BR`, `ANTIDEPRESSION_BR`, `NEUROLEPTICA_BR`, `ANXIETY_SUPRESSOR_BR`, `BARBITUATES_BR`, `ANALAGETICS_BR`, `MEDICATION_N_BR` |
| Licences | `DRIVERS_LICENSE_BR`, `DRIVERS_AGREEMENT_BR`, `DRIVERS_REPORTED_BR`, `WEAPON_LICENSE_BR`, `WEAPON_UNSUITABLE_BR`, `WEAPON_REPORTED_BR` |
| Vital status | `D_DEAD`, `DODSDATUM_BEFREG`, `DATEOFBIRTH`, `GENDER` |

### 4.9 Trajectory and rate-of-change estimates

Person-level summaries produced in the repeated-measures subsample.

| Variable | Type | Description |
|---|---|---|
| `n_measurements` | num | Number of pre-diagnosis measurement occasions for an individual (minimum 3; mean 4.35) |
| `intercept` | num | Per-individual regression intercept: estimated biological age at that individual's **first** measurement, from `lm(biological age ~ centred chronological age)` |
| `slope` | num | Per-individual rate of change in biological age per year of chronological age, from the same fit |
| `age_bin` | factor | 5-year chronological-age band, 17–82, used in the binned trajectory model |
| `lower_CI`, `upper_CI` | num | Cluster-robust 95% confidence limits for a bin slope, clustered on `sampleID` |
| `diff_women_minus_men` | num | Sex difference in bin slope, from the joint interaction model |
| `wald_statistic`, `degrees_of_freedom`, `p_value` | num | Global Wald test for heterogeneity of slopes across bins |

The binned trajectory model is
`lm(biological age ~ age + age_bin + age:age_bin)` — bin main effects **and**
the interaction, so each bin has its own intercept and its own slope, with
17–22 as the reference. Standard errors are cluster-robust (`sandwich`),
clustered on the individual.

---

## 5. Intermediate files

All under `OUTPUT_DIR`. Written with `;` separators, `,` decimals and Latin-1
encoding (`config.R` §6), and read back with the same settings.

| File | Grain | Columns | Status |
|---|---|---|---|
| `sampleIDs_baseline.csv` | one row per individual | `sampleID` | **Individual-level.** Never leaves the secure environment |
| `training_set_baseline.csv` | one row per individual | `sampleID`, `age`, `gender`, the 12 biomarkers | **Individual-level.** Never leaves the secure environment |
| `table1_kdmgroups_baseline_perID.csv` | one row per individual | `sampleID`, `KDM_group` | **Individual-level.** Never leaves the secure environment |
| `bioage_scale_df_table.csv` | one row per measure | measure, mean, s.d. | Aggregate |
| `GAM_ci_per5year.csv` | one row per 5-year age point | age, fit, lower, upper | Aggregate (Supp. Tables 5, 8) |
| `trajectory_table_*.csv` | one row per age bin | bin, slope, CI, n | Aggregate (Supp. Tables 3–4, 6–7) |
| `combined_results_*.csv` | one row per model | outcome, measure, HR, CI, p | Aggregate (Fig. 5, Supp. Figs. 14–24) |
| `combined_df_agegroup_*.csv`, `combined_sex_*.csv` | one row per model | stratum, outcome, HR, CI, p | Aggregate (Fig. 5b–c) |
| `finegray_project1_*.csv` | one row per model | outcome, model, sHR, CI, p | Aggregate (Supp. Table 9) |
| `table9_random_drop_sensitivity_RMA.csv` | one row per model | run, exposure, HR, CI | Aggregate (Supp. Tables 10–11) |
| `cif_project1_curves.csv`, `cif_project1_summary.csv` | one row per group × time | group, state, time, estimate, CI, risk set | Aggregate (response to reviewers) |

The first three are inputs to later steps and are individual-level by
necessity. Everything else is aggregate. `.gitignore` excludes the whole
output directory regardless.

---

## 6. Key in-memory objects

| Object | Grain | Rows (as analysed) | Notes |
|---|---|---|---|
| `to1po2` | one row per individual | 261,383 | Baseline analytic cohort |
| `to1po2_SveDem_D` | one row per individual | SveDem-linked subset | SveDem subanalysis |
| `to1po2_5`, `to1po2_10` | one row per individual | — | Preclinical-exclusion sensitivity: ADRD cases diagnosed within 5 or 10 years of baseline removed. The script prints the number dropped at each cut-off |
| `results_dt` | one row per **measurement** | 194,897 | Repeated-measures subsample, 44,805 individuals |
| `results_kdm`, `results_phe` | one row per individual | 44,805 | Slope and intercept per person |

---

## 7. Analytic grain — read before modifying

Some objects hold one row per **individual** and others one row per
**measurement occasion**, and the two are joined repeatedly.

**Any model whose exposure is a person-level summary — a trajectory intercept
or slope — must be fitted on person-level rows.** Fitting it on
measurement-level rows counts participants with more measurements more than
once, inflates the effective sample size and understates the standard errors.
This is a real failure mode: it occurred in an earlier version of the two-stage
analysis and was corrected during peer review. The de-duplicated, person-level
models are what the manuscript reports.

Practical rules:

- `results_dt` is measurement-level. Anything merged onto it stays
  measurement-level.
- `results_kdm` and `results_phe` are person-level. The two-stage Cox models
  use these.
- Trajectory models *are* fitted on measurement-level rows, because the
  exposure there is chronological age at the measurement, not a person-level
  summary. Within-person correlation is handled by cluster-robust standard
  errors clustered on `sampleID`, not by de-duplication.
- Repeated measurements recorded **after** an ADRD diagnosis are dropped, and
  at least three pre-diagnosis measurements are required. This excludes 61 of
  44,866 individuals, leaving 44,805 individuals and 194,897 measurements.
  `BA_rate_L.Rmd` reports the number of rows dropped.

---

## 8. Which script produces which published result

What each script computes, and the manuscript items it feeds. Figure and table
numbers are given where the correspondence is one-to-one; otherwise the output
is described.

| Script | Computes | Feeds |
|---|---|---|
| `Baseline/Selecting_biomarkers_B.Rmd` | Biomarker densities; Spearman and Pearson correlations with chronological age, overall and by sex, before and after outlier removal | Fig. 2a; biomarker density supplementary figure |
| `Baseline/Calculating_BA_B.Rmd` | Cohort construction, age-uniform training set, KDM and PhenoAge calibration and the fitted model weights; biological age correlation matrix | Fig. 2b; Supp. Table 12 (model weights); age-distribution and scatter supplementary figures |
| `Baseline/CoxPH_B.Rmd` | Outcome and covariate definitions; baseline Cox models for ADRD and the three subtypes, sex-adjusted and multivariable; fasting and calendar-date sensitivity; preclinical exclusion at 5 and 10 years | Fig. 5a; sex-adjusted, PhenoAge and fasting supplementary figures; preclinical-exclusion figures in the response to reviewers |
| `Baseline/CoxPH_SveDem_censored_B.Rmd` | The same models restricted to adjudicated SveDem diagnoses, with ADRD cases ascertained outside SveDem right-censored at diagnosis | Fig. 5d and its PhenoAge counterpart |
| `Baseline/CoxPH_subgroups_B.Rmd` | Cox models stratified by baseline age (median split, 45 years) and by sex, with tests for difference between strata | Fig. 5b–c and their PhenoAge counterparts |
| `Baseline/Table1_B.Rmd` | Table 1, including standardised mean differences; the biologically younger/older median split; inclusion counts | Table 1; extended baseline characteristics; flowchart counts |
| `Baseline/Fine_Gray/fg_project1_ADRD.R` | Fine-Gray subdistribution hazard models with death as the competing event | **Supp. Table 9** |
| `Longitudinal/Calculating_BA_L.Rmd` | Biological age at every measurement occasion, calibrated on the baseline training set | Age-distribution and scatter supplementary figures for the longitudinal subset |
| `Longitudinal/Table1_L.Rmd` | Repeated-measures descriptives: measurements per person, spacing, follow-up, by biological age group and ADRD status | Descriptive supplementary tables for the longitudinal subset |
| `Longitudinal/BA_rate_L.Rmd` | Per-individual slope and intercept; two-stage Cox models; baseline models refitted in the longitudinal subset; random drop-one-measurement sensitivity | **Supp. Tables 10–11**; longitudinal-subset baseline supplementary figures |
| `Longitudinal/Trajectories_binned_L.Rmd` | 5-year age-bin trajectories with cluster-robust intervals; global Wald test; sex differences from the joint interaction model, with pointwise intervals | Fig. 4; **Supp. Tables 3–4, 6–7**; binned-trajectory supplementary figures |
| `Longitudinal/Trajectories_GAM_L.Rmd` | Smooth trajectories by GAM with cluster-robust confidence bands; density-shaded spaghetti plots | Fig. 3; **Supp. Tables 5, 8** |
| `Review/**` | Analyses cited in the response to reviewers only | See [Review/README.md](Review/README.md) |

## 9. Data access

AMORIS is governed by Swedish law on health data. Access is granted by the
cohort's data access committee subject to ethical review, and the data cannot
be redistributed by us. Researchers wishing to reproduce these analyses should
contact the AMORIS study group; the cohort profile
([Walldius et al., *Int J Epidemiol* 2017;46:1103](https://doi.org/10.1093/ije/dyw333))
describes the application process.

The analyses in this repository were approved by the Swedish Ethical Review
Authority (2022-06634-01).
