# Longitudinal Biological Aging and Risk of Dementia in a Large Prospective Swedish Cohort

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23420921.svg)](https://doi.org/10.5281/zenodo.23420921)

R code for the analyses in the manuscript. The analysis examines both
**baseline biological age** and the **rate of change in biological age** as
exposures for Alzheimer's disease and related dementias (ADRD) in the AMORIS
cohort.

The underlying data are restricted and cannot be redistributed, so the pipeline
cannot be executed from this repository alone. It is documented in enough
detail to be read and assessed without running it, and to be run by a
researcher who holds an approved AMORIS extract — see
[DATA_DICTIONARY.md](DATA_DICTIONARY.md).

---

## Study design

**Population.** 261,383 AMORIS participants with complete biomarker data at
baseline. Of these, 44,805 individuals with three or more pre-diagnosis
biomarker measurements form the repeated-measures subsample, contributing
194,897 measurement occasions.

**Follow-up.** 26 years, through Swedish national health registers, to
31 December 2020.

**Outcome.** ADRD, ascertained from the National Patient Register,
anti-dementia drug dispensations and SveDem (the Swedish registry for
cognitive/dementia disorders), and classified as Alzheimer's disease, vascular
dementia or other dementia. A subtype analysis counts only clinically
adjudicated SveDem diagnoses as events; cases ascertained elsewhere are
right-censored at diagnosis.

**Exposure.** Biological age from the Klemera–Doubal Method (KDM, primary) and
PhenoAge (secondary), computed from 12 clinical biomarkers with the
[BioAge](https://github.com/dayoonkwon/BioAge) package. The exposure in the
survival models is the **age gap** — biological minus chronological age —
standardised, where one standard deviation of the KDM age gap corresponds to
6.5 years. Biomarkers are not equally weighted; the fitted weights are in
Supplementary Table 12 and in [DATA_DICTIONARY.md](DATA_DICTIONARY.md) §4.6.

---

## Getting started

**1. Install the dependencies.** See [DEPENDENCIES.md](DEPENDENCIES.md).
Biological age estimation needs the `BioAge` package, which is not on CRAN:

```r
remotes::install_github("dayoonkwon/BioAge", ref = "380cd3e017ce2d3457d7fbec4f92e033fcaa49c7")
```

The commit is pinned deliberately: BioAge declares version `0.1.0` both at this
commit and at its current `HEAD`, so the version number alone does not identify
the code these results were produced with. See
[DEPENDENCIES.md](DEPENDENCIES.md).

**2. Point the code at your data.** Set one environment variable — for example
in `~/.Renviron`:

```sh
AMORIS_DATA_ROOT=/path/to/amoris
```

or edit `DATA_ROOT` at the top of [config.R](config.R). That is the only
location any script needs. Then check the inputs are reachable before starting
a long run:

```r
source("config.R")
check_config()
```

Outputs go to `./output` by default; set `STUDY1_OUTPUT_DIR` to send them
elsewhere.

**3. Run the scripts in order.** See the workflow below, and
[DATA_DICTIONARY.md](DATA_DICTIONARY.md) for the full execution order, the
variable definitions, and the files passed between steps.

**4. Record your environment.** After the pipeline completes:

```r
source("write_session_info.R")
```

---

## Repository structure

```
Clinical_BA_ADRD/
├── README.md                          # this file
├── config.R                           # paths and analysis parameters — edit this
├── DEPENDENCIES.md                    # packages and versions
├── DATA_DICTIONARY.md                 # variables, data flow, run order
├── write_session_info.R               # records software versions
├── .gitignore                         # keeps participant-level files out of git
│
├── Baseline/                          # cross-sectional baseline analysis
│   ├── README.md
│   ├── Selecting_biomarkers_B.Rmd     # biomarker screening
│   ├── Calculating_BA_B.Rmd           # cohort setup, KDM + PhenoAge calibration
│   ├── CoxPH_B.Rmd                    # outcome definitions, main Cox models
│   ├── CoxPH_SveDem_censored_B.Rmd    # adjudicated SveDem subanalysis
│   ├── CoxPH_subgroups_B.Rmd          # age- and sex-stratified models
│   ├── Table1_B.Rmd                   # Table 1, biological age groups
│   ├── Forest_plots_B.Rmd             # forest plots (Fig. 5, Supp. Figs. 14–21, 23–24)
│   └── Fine_Gray/                     # competing risks (Supp. Table 9)
│
├── Longitudinal/                      # repeated-measures analysis
│   ├── README.md
│   ├── Calculating_BA_L.Rmd           # biological age at every measurement
│   ├── Table1_L.Rmd                   # repeated-measures descriptives
│   ├── BA_rate_L.Rmd                  # intercept vs. rate of change
│   ├── Trajectories_binned_L.Rmd      # 5-year age-bin trajectories
│   └── Trajectories_GAM_L.Rmd         # smooth trajectories
│
└── Review/                            # analyses cited in the response to
    ├── README.md                      #   reviewers only, not in the paper
    ├── Baseline/Cumulative_incidence/  # absolute risk, competing events
    └── Longitudinal/                  # calendar-year adjustment, inclusion checks
```

---

## Workflow

### Step 1 — Baseline

Start with `Selecting_biomarkers_B.Rmd`, then `Calculating_BA_B.Rmd` to build
the cohort, calibrate KDM and PhenoAge, and write the sample and training-set
files the longitudinal analysis depends on. `CoxPH_B.Rmd` then defines the
outcomes and covariates that everything downstream uses.

See [Baseline/README.md](Baseline/README.md).

### Step 2 — Longitudinal

The longitudinal scripts read `sampleIDs_baseline.csv` and
`training_set_baseline.csv` from step 1. **The training set must be reused, not
regenerated** — biological age is calibrated on it, so a new training set puts
the two sets of biological ages on different scales.

See [Longitudinal/README.md](Longitudinal/README.md).

### Step 3 — Review analyses (optional)

Analyses run in response to peer review and cited in the response to
reviewers, which produce no figure or table in the manuscript or supplementary
information. See [Review/README.md](Review/README.md).

---

## Two things to know before modifying anything

**Analytic grain.** Some objects hold one row per *individual* and others one
row per *measurement occasion*, and the two are joined repeatedly. Any model
whose exposure is a person-level summary — a trajectory intercept or slope —
must be fitted on person-level rows, or participants with more measurements are
counted more than once. [DATA_DICTIONARY.md](DATA_DICTIONARY.md) §7 sets out
which object is which; it is worth reading before touching `BA_rate_L.Rmd`.

**Calibration order.** KDM and PhenoAge coefficients are estimated once, in
`Calculating_BA_B.Rmd`, and applied everywhere else. Scripts are not
independent of each other.

---

## Data availability

The scripts require the AMORIS cohort data with linkage to Swedish national
health registers.

**These data cannot be redistributed.** AMORIS is governed by Swedish law on
health data, and access is granted by the cohort's data access committee
subject to ethical review. Researchers wishing to reproduce these analyses
should contact the AMORIS study group; the cohort profile
([Walldius et al., *Int J Epidemiol* 2017;46:1103](https://doi.org/10.1093/ije/dyw333))
describes the application process.

Nothing derived from individual participants is included in this repository.
`.gitignore` excludes the source extracts, R workspaces, the output directory
and the individual-level intermediate files by name; see
[DATA_DICTIONARY.md](DATA_DICTIONARY.md) §5 for which intermediate files are
individual-level and which are aggregate.

The analyses were approved by the Swedish Ethical Review Authority
(2022-06634-01).

---

## License

The code is released under the [MIT License](LICENSE). The AMORIS data it
analyses are not covered by this licence and are not included; see Data
availability above.

---

## Contact

Karolina Gustavsson
Department of Medical Epidemiology and Biostatistics
Karolinska Institutet

## Citation

This repository accompanies a manuscript under review. If you use the code or
materials, please cite:

Gustavsson K, Sjölander A, Mostafaei S, Mak JKL, Religa D, Käll L, Modig K,
Hägg S. *Longitudinal biological aging and risk of dementia in a large
prospective Swedish cohort*. Manuscript under review, 2026.

The version of the code used for the revised manuscript is archived at Zenodo:

Gustavsson K, Sjölander A, Mostafaei S, Mak JKL, Religa D, Käll L, Modig K,
Hägg S. *Clinical_BA_ADRD* (v1.0.0). Zenodo, 2026.
https://doi.org/10.5281/zenodo.23420921

The concept DOI https://doi.org/10.5281/zenodo.23420920 always resolves to the
latest archived version. A formal citation of the article will be added on
publication.
