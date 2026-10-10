# Dependencies

Software used to produce the results in the manuscript. Versions were recorded
on the analysis platform where the analyses were run, with
`write_session_info.R`.

---

## R

| | |
|---|---|
| R version | 4.5.3 (2026-03-11 ucrt) |
| Platform | `x86_64-w64-mingw32/x64` |
| OS | Windows 11 x64 (build 26100) |

Packages were installed to a project library on the analysis platform rather
than to the default system library, because the platform has no outbound
internet access. See [Platform notes](#platform-notes) below. The three
packages that ship with R — `MASS`, `mgcv` and `survival` — were loaded from
R's own library; everything else came from the project library.

---

## The BioAge package

Biological age is computed with `kdm_calc()` and `phenoage_calc()` from
**BioAge** (Kwon & Belsky, *GeroScience* 2021;43:2795–2808), which is not on
CRAN.

| | |
|---|---|
| Declared version | 0.1.0 |
| Source | https://github.com/dayoonkwon/BioAge |
| Commit used | `380cd3e017ce2d3457d7fbec4f92e033fcaa49c7` (2023-10-11) |

The functions were used from a local copy of the package source taken in 2024,
loaded directly rather than installed as a package, so R recorded no version
or commit for them. The commit above is identified from the upstream history:
it was the repository's `HEAD` for the whole of 2024 — there are no commits
between 2023-10-11 and 2026-04-02 — so it is the only state of the package the
analyses could have used.

### Install the version these results were produced with

**The declared version number does not identify the code.** It is `0.1.0` both
at the commit above and at the current `HEAD`, but the two differ. Pin the
commit:

```r
# install.packages("remotes")
remotes::install_github("dayoonkwon/BioAge",
                        ref = "380cd3e017ce2d3457d7fbec4f92e033fcaa49c7")
```

On a platform without outbound internet, download that ref on a connected
machine, transfer it in, and install locally:

```r
install.packages("BioAge-380cd3e.tar.gz", repos = NULL, type = "source")
```

### Why the commit has to be pinned

Installing from `HEAD` will **not** reproduce the PhenoAge values reported
here. On 2026-04-02, after these analyses were run, one constant in
`R/phenoage_calc.R` was changed:

```r
# at the commit used here
dat$phenoage0 = ((log(-.0055305 * (log(1 - m_orig))) / .09165)  + 141.50225)
# at HEAD today
dat$phenoage0 = ((log(-.0055305 * (log(1 - m_orig))) / .090165) + 141.50225)
```

`0.090165` is the constant published by Levine et al., and the upstream commit
message describes the change as restoring it, so the earlier `0.09165` was an
error in the package. It is the only difference in the package's R source
between the two commits; `kdm_calc()` is unaffected, so the KDM results — the
primary measure throughout the manuscript — are identical under either version.

The change shifts absolute PhenoAge by roughly 1.5 years and rescales it by
about 1.6%. The exposure used in the models is the **standardised PhenoAge age
gap**, which is insensitive to a shift and largely insensitive to a rescaling,
so the reported hazard ratios per standard deviation are not materially
affected; absolute PhenoAge values are.

### A note on the earlier `function.R`

Earlier versions of `Calculating_BA_B.Rmd` and `Calculating_BA_L.Rmd` began by
sourcing a local file called `function.R`, which carried the BioAge functions.
It existed only to make them available on a platform with no CRAN or GitHub
access, and it was never part of this repository — which made those two scripts
impossible to run anywhere else. They now depend on the package, installed at
the pinned commit above.

---

## CRAN packages

Every package the scripts attach or call with `::`.

| Package | Version | Used for |
|---|---|---|
| broom | 1.0.5 | Tidying model output |
| data.table | 1.14.8 | Fast I/O and grouped operations on the measurement-level data |
| dplyr | 1.1.3 | Data manipulation |
| fields | 17.1 | `interp.surface()`, for density shading in the trajectory figures |
| flexsurv | 2.2.2 | Parametric survival models |
| ggplot2 | 4.0.1 | Figures |
| ggtext | 0.1.2 | `geom_richtext()`, for marked-up annotations in the trajectory figures |
| gridExtra | 2.3 | Multi-panel figure layout |
| gtable | 0.3.6 | Panel alignment in composed figures |
| haven | 2.5.3 | Reading the SAS `.sas7bdat` extracts |
| here | 1.0.1 | Locating `config.R` from any script, without `setwd()` |
| htmlTable | 2.4.1 | Table rendering |
| lmtest | 0.9-40 | Coefficient tests against a robust covariance matrix |
| MASS | 7.3-65 | `kde2d()`, for the two-dimensional density shading |
| mgcv | 1.9-4 | Generalized additive models for the smooth trajectories |
| rmarkdown | 2.30 | Notebook rendering |
| sandwich | 3.1-1 | Cluster-robust standard errors for repeated measures |
| stringr | 1.5.0 | String handling |
| survey | 4.2-1 | Survey-weighted summaries |
| survival | 3.8-6 | Cox models, `finegray()`, and `survfit()` for cumulative incidence |
| table1 | 1.4.3 | Table 1 |
| tidyr | 1.3.0 | Reshaping |

To install the same versions, rather than whatever is current:

```r
# install.packages("remotes")
remotes::install_version("broom",      "1.0.5")
remotes::install_version("data.table", "1.14.8")
remotes::install_version("dplyr",      "1.1.3")
# ... and so on for the table above
```

Or, to install current versions and accept the difference:

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
it cannot express the age-based entry and exit the rest of the analysis uses.
The reasoning is in
[`Review/Baseline/Cumulative_incidence/README.md`](Review/Baseline/Cumulative_incidence/README.md).

---

## Platform notes

The analyses ran on a secure Windows analysis platform with no outbound
internet access, where packages were installed to a project library from an
internal mirror. Both are configured through the environment rather than in the
scripts — see [`config.R`](config.R) §7:

```sh
R_PROJECT_LIB=/path/to/project/library
R_CRAN_MIRROR=https://internal.mirror.example/cran
```

Leave both unset to install and load packages normally.

The Fine-Gray and cumulative incidence scripts are run from a terminal rather
than from inside RStudio, because the `finegray()` risk-set expansion is
memory-hungry on a cohort of this size. `run_fg_project1.ps1` is a PowerShell
driver that gives each outcome its own R process so memory is released between
outcomes; the equivalent plain commands are in each folder's `README.md`.

---

## Re-recording these versions

```r
source("write_session_info.R")
```

Run it from the repository root on the analysis platform. It records the
versions of every package listed above and writes `SESSION_INFO.txt`, with the
versions formatted as a markdown table ready to paste into this file.

---

## Appendix: full `sessionInfo()`

```
R version 4.5.3 (2026-03-11 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Swedish_Sweden.utf8  LC_CTYPE=Swedish_Sweden.utf8
[3] LC_MONETARY=Swedish_Sweden.utf8 LC_NUMERIC=C
[5] LC_TIME=Swedish_Sweden.utf8

time zone: Europe/Stockholm
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets
[7] methods   base

other attached packages:
 [1] tidyr_1.3.0        table1_1.4.3       survey_4.2-1
 [4] Matrix_1.7-4       stringr_1.5.0      sandwich_3.1-1
 [7] rmarkdown_2.30     mgcv_1.9-4         nlme_3.1-168
[10] MASS_7.3-65        lmtest_0.9-40      zoo_1.8-12
[13] htmlTable_2.4.1    here_1.0.1         haven_2.5.3
[16] gtable_0.3.6       gridExtra_2.3      ggtext_0.1.2
[19] ggplot2_4.0.1      flexsurv_2.2.2     survival_3.8-6
[22] fields_17.1        RColorBrewer_1.1-3 viridisLite_0.4.2
[25] spam_2.11-3        dplyr_1.1.3        data.table_1.14.8
[28] broom_1.0.5

loaded via a namespace (and not attached):
 [1] dotCall64_1.2       xfun_0.40           htmlwidgets_1.6.2
 [4] lattice_0.22-9      numDeriv_2016.8-1.1 quadprog_1.5-8
 [7] vctrs_0.6.3         tools_4.5.3         generics_0.1.3
[10] tibble_3.2.1        fansi_1.0.4         pkgconfig_2.0.3
[13] checkmate_2.2.0     S7_0.2.1            lifecycle_1.0.3
[16] compiler_4.5.3      farver_2.1.1        statmod_1.5.0
[19] mitools_2.4         mstate_0.3.2        htmltools_0.5.6
[22] maps_3.4.3          Formula_1.2-5       pillar_1.9.0
[25] muhaz_1.2.6.4       digest_0.6.33       tidyselect_1.2.0
[28] stringi_1.7.12      mvtnorm_1.2-3       purrr_1.0.2
[31] forcats_1.0.0       splines_4.5.3       fastmap_1.1.1
[34] rprojroot_2.0.4     cli_3.6.1           magrittr_2.0.3
[37] utf8_1.2.3          withr_3.0.2         scales_1.4.0
[40] backports_1.4.1     deSolve_1.38        hms_1.1.3
[43] evaluate_1.0.3      knitr_1.44          rlang_1.1.1
[46] gridtext_0.1.5      Rcpp_1.0.11         DBI_1.1.3
[49] glue_1.6.2          xml2_1.3.5          rstudioapi_0.15.0
[52] R6_2.5.1
```
