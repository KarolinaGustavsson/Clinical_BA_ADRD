# Dependencies

Software required to run the analysis pipeline.

> **The version column is not yet filled in.** Run
> [`write_session_info.R`](write_session_info.R) on the platform where the
> results were produced and paste its output into the tables below and into the
> appendix. Versions recorded on any other machine would not describe what
> produced the published numbers.

---

## R

| | |
|---|---|
| R version | 4.5.2 (as stated in the manuscript Methods) |
| Platform | `___` (from `sessionInfo()`) |
| OS | `___` |

---

## The BioAge package

Biological age is computed with `kdm_calc()` and `phenoage_calc()` from
**BioAge** (Kwon & Belsky, *GeroScience* 2021;43:2795–2808).

BioAge is not on CRAN. Install it from GitHub:

```r
# install.packages("remotes")
remotes::install_github("dayoonkwon/BioAge")
```

On a platform with no outbound internet — which is the normal situation for
restricted health data — download the source on a connected machine, transfer
it in, and install locally:

```r
install.packages("BioAge_<version>.tar.gz", repos = NULL, type = "source")
```

| | |
|---|---|
| BioAge version | `___` |
| Commit SHA | `___` |

The package also supplies the NHANES III reference dataset used in the
optional section of `Calculating_BA_B.Rmd` that sets the size of the training
set: `data("NHANES3", package = "BioAge")`.

### A note on the earlier `function.R`

Earlier versions of `Calculating_BA_B.Rmd` and `Calculating_BA_L.Rmd` began by
sourcing a local file called `function.R`, which held a copy of the BioAge
functions. It existed only to make the package available inside a secure
analysis environment with no CRAN or GitHub access, and it was never part of
this repository — which made those two scripts impossible to run anywhere
else. Both now depend on the published package instead, and the offline route
is the local source install shown above.

---

## CRAN packages

Every package the scripts attach or call with `::`. All are on CRAN.

| Package | Version | Used for |
|---|---|---|
| broom | `___` | Tidying model output |
| data.table | `___` | Fast I/O and grouped operations on the measurement-level data |
| dplyr | `___` | Data manipulation |
| fields | `___` | `interp.surface()`, for density shading in the trajectory figures |
| flexsurv | `___` | Parametric survival models |
| ggplot2 | `___` | Figures |
| ggtext | `___` | `geom_richtext()`, for marked-up annotations in the trajectory figures |
| gridExtra | `___` | Multi-panel figure layout |
| gtable | `___` | Panel alignment in composed figures |
| haven | `___` | Reading the SAS `.sas7bdat` extracts |
| here | `___` | Locating `config.R` from any script, without `setwd()` |
| htmlTable | `___` | Table rendering |
| lmtest | `___` | Coefficient tests against a robust covariance matrix |
| MASS | `___` | `kde2d()`, for the two-dimensional density shading |
| mgcv | `___` | Generalized additive models for the smooth trajectories |
| rmarkdown | `___` | Notebook rendering |
| sandwich | `___` | Cluster-robust standard errors for repeated measures |
| stringr | `___` | String handling |
| survey | `___` | Survey-weighted summaries |
| survival | `___` | Cox models, `finegray()`, and `survfit()` for cumulative incidence |
| table1 | `___` | Table 1 |
| tidyr | `___` | Reshaping |

Install them all at once:

```r
install.packages(c(
  "broom", "data.table", "dplyr", "fields", "flexsurv", "ggplot2",
  "ggtext", "gridExtra", "gtable", "haven", "here", "htmlTable",
  "lmtest", "MASS", "mgcv", "rmarkdown", "sandwich", "stringr",
  "survey", "survival", "table1", "tidyr"
))
```

No package beyond `survival` is needed for the competing-risk analyses:
Fine-Gray uses `survival::finegray()` and cumulative incidence uses
`survival::survfit()` on a multi-state fit, which requires **survival ≥ 3.0**.
`cmprsk::cuminc()` is deliberately not used — it takes only time and status, so
it cannot express the age-based entry/exit the rest of the analysis uses. The
reasoning is in
[`Review/Baseline/Cumulative_incidence/README.md`](Review/Baseline/Cumulative_incidence/README.md).

---

## Platform notes

The analyses were run on a secure platform with no direct CRAN access, where
packages were installed to a project library from an internal mirror. Both are
configured through the environment rather than in the scripts — see
[`config.R`](config.R) §7:

```sh
R_PROJECT_LIB=/path/to/project/library
R_CRAN_MIRROR=https://internal.mirror.example/cran
```

Leave both unset to install and load packages normally.

The Fine-Gray and cumulative incidence scripts are run from a terminal, not
from inside RStudio, because the `finegray()` risk-set expansion is
memory-hungry on a cohort of this size. `run_fg_project1.ps1` is a PowerShell
driver that gives each outcome its own R process so memory is released between
outcomes; the equivalent plain commands are in each folder's `README.md`.

---

## Recording the versions

```r
source("write_session_info.R")
```

Run it from the repository root, on the analysis platform, after the pipeline
completes. It loads every package listed above, writes `SESSION_INFO.txt`, and
formats the versions as a markdown table ready to paste into this file.

---

## Appendix: full `sessionInfo()`

```
___ paste the sessionInfo() section of SESSION_INFO.txt here ___
```
