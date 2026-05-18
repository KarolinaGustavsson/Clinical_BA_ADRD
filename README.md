# Longitudinal Biological Aging and Risk of Dementia in a Large Prospective Swedish Cohort

## Overview

This repository contains R scripts for analyzing the association between biological age and risk of Alzheimer's disease and related dementias (ADRD) in the AMORIS cohort. The analysis examines both **baseline bioage** and **rate of change in bioage** as exposures.

## Study Design

### Population
- **Baseline cohort (cross-sectional):** 261,383 AMORIS participants with complete biomarker data
- **Longitudinal cohort:** 44,866 individuals with repeated biomarker measurements over 26 years of follow-up

### Outcome
Alzheimer's Disease and Related Dementia (ADRD) diagnoses tracked through Swedish national health registers, categorized into:
- Alzheimer's disease (AD)
- Vascular dementia
- Other dementias

### Biological Age Assessment
- **Method:** Klemera-Doubál Method (KDM) and PhenoAge algorithms
- **Biomarkers (12):** S_Alb, S_Krea, TC, TG, S_K, S_Urea, S_LD, S_Ca, S_Urat, S_FAMN, S_Hapt, fS_Gluk
- **Key metrics:** 
  - Biological age (years)
  - Bioage advance (difference from chronological age, standardized)

## Repository Structure

```
Clinical_BA_ADRD/
├── README.md                          # This file
├── Baseline/                          # Cross-sectional baseline analysis
│   ├── README.md                      # Baseline analysis guide
│   ├── Calculating_BA_B.Rmd           # Main bioage calculation script
│   ├── Selecting_biomarkers_B.Rmd     # Biomarker validation
│   ├── Tabel1_B.rmd                   # Baseline Table 1
│   ├── CoxPH_B.Rmd                    # Main Cox PH analysis
│   └── CoxPH_subgroups_B.Rmd          # Stratified survival analysis
│
└── Longitudinal/                      # Repeated measurements analysis
    ├── README.md                      # Longitudinal analysis guide
    ├── Calculating_BA_L.rmd           # Bioage calculation for repeated measurements
    ├── Table1_L.rmd                   # Longitudinal cohort characteristics
    ├── BA_rate_L.rmd                  # Rate of change analysis
    ├── Trajectories_binned__L.Rmd     # Trajectory visualization (5-year bins)
    └── Trajectories_GAM_L.rmd         # Trajectory modeling with GAM
```

## Workflow

### Step 1: Baseline Analysis (Baseline folder)
Start with `Calculating_BA_B.Rmd` to set up the cohort, calculate biological ages, and generate outputs needed for the longitudinal analysis.

**See [Baseline/README.md](Baseline/README.md) for detailed instructions.**

### Step 2: Longitudinal Analysis (Longitudinal folder)
The longitudinal scripts use outputs from the baseline analysis.

**See [Longitudinal/README.md](Longitudinal/README.md) for detailed instructions.**

## Data Requirements

The scripts require access to the AMORIS data along with Swedish registers.

## Contact

Karolina Gustavsson  
Department of Medical Epidemiology and Biostatistics  
Karolinska Institute

## Citation

This repository accompanies an unpublished study. If you use the code or materials, please cite:

Gustavsson K, Sjölander A, Mostafaei S, Mak JKL, Religa D, Käll L, Modig K, Hägg S. *Longitudinal biological aging and risk of dementia in a large prospective Swedish cohort*. Manuscript in preparation, 2026.

A formal citation will be added once the manuscript is available as a preprint or publication.

---

**Last updated:** May 2026
