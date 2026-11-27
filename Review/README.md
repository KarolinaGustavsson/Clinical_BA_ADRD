# Review analyses

Code for analyses that were run **in response to peer review** and are cited in
the response to reviewers, but that do not produce a figure or table in the
manuscript or the supplementary information.

Everything behind a published result lives in [`Baseline/`](../Baseline) and
[`Longitudinal/`](../Longitudinal) instead. The split is by where the result
appears, not by how the analysis was done — a sensitivity analysis reported in
the supplementary information (for example the Fine-Gray models in
Supplementary Table 9, the cumulative incidence curves in Supplementary
Figs. 26–27 and the calendar-year adjustment in Supplementary Table 15) stays
in the main folders. Analyses that share code with the main models also stay
there: the 5-year lag is in `Baseline/CoxPH_lag_B.Rmd`, next to the 10-year
lag, and the slope-scaling check for the drop-one models is section 3.9 of
`Longitudinal/BA_rate_L.Rmd`.

Scripts source `config.R` from the repository root and follow the same
conventions as the rest of the repository.

---

## What is here

| Script | Question it answers | Reported as |
|---|---|---|
| [`Longitudinal/Inspect_ADRD_before_last_measurement_L.Rmd`](Longitudinal/Inspect_ADRD_before_last_measurement_L.Rmd) | Where does the ADRD diagnosis fall within each individual's series of measurements? The basis for requiring three pre-diagnosis measurements | Inclusion counts in the Methods and the response |

---

## Not included here

Two analyses discussed in the response to reviewers produce no code in this
folder, by design:

- **The landmark design** at the third biomarker measurement was not
  estimated. The longitudinal Cox models instead start follow-up at the last
  measurement used to estimate the slope and intercept (delayed entry), which
  separates the exposure period from the risk period while keeping all
  pre-diagnosis measurements (`Longitudinal/BA_rate_L.Rmd`, section 3.6).
- **Regression to the mean** is answered algebraically in the response — the
  estimated slope is unbiased for the slope in the target regression — and the
  argument rests on no new estimation.
