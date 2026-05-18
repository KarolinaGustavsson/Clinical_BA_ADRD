# Longitudinal Analysis: Biological Age and ADRD Risk

## Overview

This folder contains scripts for **longitudinal analysis** of 44,866 AMORIS participants with repeated biomarker measurements over 26 years of follow-up. The analysis examines whether **baseline biological age** and **rate of change in bioage** independently predict dementia risk.

## Prerequisites

Before running longitudinal scripts, complete the **Baseline analysis**.

---

## Workflow

Execute scripts in the following order:

### 1. **Calculating_BA_L.rmd** - Longitudinal Bioage Calculation
**Purpose:** Calculate biological ages for all repeated measurements using baseline training set

**Key Steps:**
- Load full AMORIS longitudinal cohort with all measurements
- Restrict to sampleIDs from baseline analysis (for consistency)
- Use **same training set as baseline** for KDM and PhenoAge algorithms
- Apply bioage algorithms to each measurement
- Merge with baseline data
- Retain individuals with ≥3 measurements for trajectory analysis

**Key Outputs:**
- `to1po2` - Longitudinal dataset with all measurements and calculated bioages
- Prepared for downstream trajectory and rate-of-change analyses

---

### 2. **Table1_L.rmd** - Longitudinal Cohort Characteristics
**Purpose:** Describe longitudinal cohort at baseline and repeated measurement patterns

**Key Steps:**
- Load preprocessed longitudinal cohort
- Extract baseline (first) measurement for each individual
- Characterize repeated measurement frequency:
  - Number of measurements per participant
  - Time span between measurements
  - Mean interval between measurements
- Generate **Table 1** in a similar way as the table for the baseline analysis
- Calculate p-values for group comparisons

**Key Outputs:**
- `table_kdm_longitud_summary_part1_p.csv` - Measurement frequency statistics
- `table_kdm_longitud_baseline_summary_part2_p.csv` - Baseline characteristics by bioage group
- `followuptime_byADRD_l.csv` - Follow-up characteristics by ADRD status

---

### 3. **BA_rate_L.rmd** - Time-to-Event Analysis
**Purpose:** Analyze whether **rate of change** predicts dementia risk beyond baseline biological age

**Key Steps:**

**Part 1: Data Preparation**
- Load longitudinal data with repeated measurements
- Link SveDem for dementia outcomes
- Harmonize dementia diagnostic codes across sources
- Calculate time-to-event and survival status

**Part 2: Trajectory Slope Calculation**
- For each individual: fit linear regression of biological age vs. chronological age
- Extract slope (rate of change in years of bioage per year of chronological age)
- Standardize slopes by sex and bioage group

**Part 3: Survival Analysis**
- Fit Cox models with ADRD as outcome
- Models stratified by sex 

**Key Outputs:**
- `combined_results_ADRD_RMA.csv` - Hazard ratios for baseline and rate of change

---

### 4. **Trajectories_binned_L.Rmd** - Trajectory Analysis with 5-Year Age Bins
**Purpose:** Visualize bioage trajectories across the lifespan and test for sex differences

**Key Steps:**

**Part 1: Data Preparation**
- Create 5-year age bins
- For each age bin: fit linear regression within bin
- Calculate slope (KDM or PhenoAge change per year of age)
- Compute 95% CI for each slope

**Part 2: Visualization**
- Create spaghetti plots (individual trajectories overlaid)
- Plot predicted trajectory line (fixed slope) with confidence interval
- Display trajectories for:
  - All individuals
  - Men only
  - Women only

**Part 3: Statistical Testing**
- Global Wald test: are slopes significantly different across age bins?
- Test for sex interaction in slopes
- Save slope estimates by age bin

**Key Outputs:**
- `trajectory_table_kdm.csv` / `trajectory_table_phe.csv` - Slope estimates by age bin
- `trajectory_table_wald_kdm_FH.csv` - Global Wald test results
- trajectory figures 

---

### 5. **Trajectories_GAM_L.rmd** - Smoothed Trajectory Modeling
**Purpose:** Fit smooth nonparametric trajectories using Generalized Additive Models

**Key Steps:**
- Apply Generalized Additive Models (GAM) to bioage trajectories
- Fit smooth spline: KDM (or PhenoAge) ~ s(chronological age)
- Generate predictions across age range with 95% CI
- Create spaghetti plots with GAM smooth overlay
- Subsample individuals for visualization (10 disjoint chunks)
- Add kernel density estimate shading to show data density

**Key Outputs:**
- `GAM_ci_per5year.csv` - Predicted bioages with CI at landmark ages
- PNG figures with spaghetti plots + GAM smooth + density shading:
  - `kdmspag_10grids1percent.png`
  - `kdmspag_10groups_density_NEW.png`

---

## Data & Variables

### Input Data
- Longitudinal AMORIS with **≥3 repeated measurements** per individual
- 44,866 participants

### Key Derived Variables
- **Baseline assessment:**
  - `kdm_advance_sd`, `phenoage_advance_sd` (standardized bioage difference)

- **Rate of change:**
  - `kdm_slope` - KDM change per year of chronological age
  - `phenoage_slope` - PhenoAge change per year

- **Outcomes:**
  - `ADRD`, `Alzheimers`, `Vascular`, `Other` (binary: 0/1)

---

## Important Notes

1. **Longitudinal Design:** Only individuals with ≥3 measurements are included in rate-of-change analyses (to ensure reliable slope estimates)

2. **Training Set Consistency:** Uses the **same bioage training set as baseline analysis** to ensure comparability

3. **Multiple Measurements per Individual:** Scripts correctly handle repeated measurements (person-time structure for survival analysis)

4. **Fasting Status:** Optional subgroup analyses available (restrict to health checkups: `debiteringskod %in% c("H")`)

5. **Output Format:** CSVs use semicolon delimiters and Latin1 encoding for Swedish system compatibility

---

## Workflow Summary

```
Baseline Analysis (1 hour)
    ↓
    └─→ sampleIDs_baseline.csv, training_set_baseline.csv
         ↓
    Longitudinal Analysis (2-3 hours)
         ├─→ Calculating_BA_L.rmd (20-30 min)
         ├─→ Table1_L.rmd (5-10 min)
         ├─→ BA_rate_L.rmd (30-45 min, Cox models)
         ├─→ Trajectories_binned_L.Rmd (15-20 min)
         └─→ Trajectories_GAM_L.rmd (10-15 min)
```


---

**Last updated:** May 2026
