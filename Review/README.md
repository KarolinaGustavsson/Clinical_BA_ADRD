# Review analyses

Code for analyses that were run **in response to peer review** and are cited in
the response to reviewers, but that do not produce a figure or table in the
manuscript or the supplementary information.

Everything behind a published result lives in [`Baseline/`](../Baseline) and
[`Longitudinal/`](../Longitudinal) instead. The split is by where the result
appears, not by how the analysis was done — a sensitivity analysis reported in
the supplementary information (for example the Fine-Gray models in
Supplementary Table 9) stays in the main folders. Analyses that share code
with the main models also stay there: the preclinical-exclusion models
(Response Figs. 2–3) are in `CoxPH_B.Rmd`.

The subfolders mirror the main structure: `Baseline/` for the cross-sectional
cohort, `Longitudinal/` for the repeated-measures subsample. Scripts source
`config.R` from the repository root and follow the same conventions as the
rest of the repository.

---

## What is here

| Script | Question it answers | Reported as |
|---|---|---|
| [`Baseline/Cumulative_incidence/`](Baseline/Cumulative_incidence) | Absolute risk rather than relative: cumulative incidence of ADRD and of death by biological age group, with death as a competing event | Figure in the response to reviewers |
| [`Longitudinal/Trajectories_binned_calyear_L.Rmd`](Longitudinal/Trajectories_binned_calyear_L.Rmd) | Could secular assay drift between 1985 and 1996 be mistaken for within-person aging? Trajectories and sex differences re-estimated with calendar year of measurement added to every model | Narrative in the response to reviewers; estimates essentially unchanged |
| [`Longitudinal/Inspect_ADRD_before_last_measurement_L.Rmd`](Longitudinal/Inspect_ADRD_before_last_measurement_L.Rmd) | Where does the ADRD diagnosis fall within each individual's series of measurements? The basis for requiring three pre-diagnosis measurements | Inclusion counts in the Methods and the response |

---

## Baseline/Cumulative_incidence

Cumulative incidence functions estimated with the Aalen-Johansen estimator, on
the same analytic sample, event definitions and age timescale as the Cox and
Fine-Gray models.

`fg_project1.R` here is a later version of the script in
[`../Baseline/Fine_Gray`](../Baseline/Fine_Gray): it adds a `cif` mode and
keeps the earlier `main`, `svedem` and `inspect` modes. **The published
Fine-Gray subdistribution hazard ratios in Supplementary Table 9 were produced
by the earlier script**, `../Baseline/Fine_Gray/fg_project1_ADRD.R`, and the
`cif` mode never calls `finegray()`, so running it cannot move those numbers.
Both versions are kept so that each published result can be traced to the code
that actually produced it.

```
Rscript fg_project1.R inspect          # data diagnostics; run first
Rscript fg_project1.R cif              # cumulative incidence, all 4 outcomes
Rscript fg_project1.R cif ADRD         # one outcome
Rscript plot_cif_results.R             # figures from the saved curves
```

Run it from a terminal rather than inside RStudio, and read
[`README.md`](Baseline/Cumulative_incidence/README.md) first — it documents the
inputs, the derivation of the PhenoAge groups from the KDM split, the
checkpointing, and the memory behaviour. `plot_cif_results.R` writes the
figure captions, axis definitions and the estimates behind them to
`README_cif_figures.md` in the run directory; that file holds results and so
is not kept in this repository.

This analysis needs `table1_kdmgroups_baseline_perID.csv`, which is
individual-level. It stays inside the secure environment; see
[`../DATA_DICTIONARY.md`](../DATA_DICTIONARY.md) §5.

---

## Longitudinal/Trajectories_binned_calyear_L.Rmd

`Trajectories_binned_L.Rmd` with one calendar-year term added to every model:

```r
lm(phenoage ~ age * age_bin)               # published
lm(phenoage ~ age * age_bin + cal_year_c)  # here (section 5.0)
```

The `CAL_ADJ` switch at the top selects a linear year term (`"linear"`, the
default, and what the response reports), year fixed effects (`"factor"`, which
also absorbs a reagent or instrument change that shifts values in a single
year) or `"none"`, which reproduces the model without calendar year as a check.
The calendar term does not interact with age, so the bin slopes are recovered
from the coefficient vector exactly as in the model without it. Every output file carries the `TAG`
suffix so a run cannot overwrite the published results.

Note: the slope models in sections 2.1 and 2.3 still use the earlier
shared-intercept specification (`age + age:age_bin`), and `CAL_ADJ` defaults
to `"linear"`. Section 5.0 repeats the comparison with the current
specification (`age * age_bin`, as in Supplementary Tables 4-5) and a linear
calendar-year term, fitted with and without calendar year on identical rows.
That is the comparison quoted in the response to reviewers: the largest change
in any slope over the 13 reported bins (17-82 years) was 0.009. The
sex-difference models in section 3 already use the current specification.

---

## Not included here

Two analyses discussed in the response to reviewers produce no code in this
folder, by design:

- **The landmark design** at the third biomarker measurement was considered and
  argued against on interpretive grounds rather than estimated, so there is no
  script to publish.
- **Regression to the mean** is answered algebraically in the response — the
  estimated slope is unbiased for the slope in the target regression — and the
  argument rests on no new estimation.
