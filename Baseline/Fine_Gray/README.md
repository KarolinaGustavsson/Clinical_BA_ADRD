# Fine-Gray competing risk analysis — project 1 (dementia)

Fine-Gray subdistribution hazard models for dementia, with **death as the
competing event** and **attained age as the timescale**, mirroring the Cox
models in `CoxPH_B.Rmd` from the `Clinical_BA_ADRD` repo.

Adapted from the cancer side project's Fine-Gray scripts.

**This script produced the subdistribution hazard ratios in Supplementary
Table 9.** A later version of it, adding cumulative incidence functions
(Aalen-Johansen) for the response to reviewers, is in
[`../../Review/Baseline/Cumulative_incidence`](../../Review/Baseline/Cumulative_incidence).
That version keeps the `main`, `svedem` and `inspect` modes unchanged, but the
published Supplementary Table 9 numbers come from the script in this folder.

---

## Running it on the secure analysis platform

Put these three files in a folder together with `cox_main_2025_06_06.RData`:

```
fg_project1_ADRD.R              config, all logic, and the three run modes below
test_fg_project1_ADRD.R         80 synthetic checks, no real data needed
run_fg_project1.ps1        optional: one process per outcome (Windows)
```

`fg_project1_ADRD.R` takes a mode as its first argument:

```powershell
Rscript test_fg_project1_ADRD.R           # ~30 s, confirms nothing broke in transit
Rscript fg_project1_ADRD.R inspect        # data diagnostics, run once before trusting anything
Rscript fg_project1_ADRD.R main           # 4 outcomes
Rscript fg_project1_ADRD.R svedem         # 10 outcomes
```

**On the test suite:** 78 of the 80 checks pass. The two failures are both in
section 10, "Time rounding", and are a stale expectation in the test rather
than a fault in the analysis: the check compares counts of distinct exit times
and expansion size, which is not what the rounding actually changes. The
rounding check was rewritten for the later version of this script in
`../../Review/Baseline/Cumulative_incidence`, whose suite passes in full, and
the rewritten check shows that rounding *excludes* a small number of people —
those whose exit age is rounded onto their entry age — rather than moving
anyone between event categories. Supplementary Table 9 was produced without
the rounding escape hatch, so neither the failures nor the rounding affect it.

No editing required. The script finds the RData file itself. No data
preparation step is needed either — the four timeline ages are derived from the
existing `cox_main_2025_06_06.RData` on the fly (see *The four timeline ages*
below).

`fg_project1_ADRD.R inspect` prints only counts and summaries, no individual-level
data, so its output is safe to copy out of the secure analysis platform. Check in particular:

- section 1: how many `ADRD` cases have **no** diagnosis date (should be 0)
- section 2: whether `lastDate` equals `DODSDATn` for the dead
- section 3: whether follow-up stops at diagnosis
- section 5: deaths or diagnoses recorded *before* baseline

**If a run dies, just start it again.** Every model is written to CSV the
moment it finishes, and completed models are skipped on restart. Nothing is
recomputed and nothing is duplicated.

### If it still runs out of memory

Run one outcome per process, so memory is fully released between outcomes:

```powershell
.\run_fg_project1.ps1              # or: .\run_fg_project1.ps1 -Which svedem
```

Or drive individual outcomes by hand:

```powershell
Rscript fg_project1_ADRD.R main ADRD
Rscript fg_project1_ADRD.R svedem S_vas
```

Last resort: set `fg_config$time_round_digits <- 2` near the top of
`fg_project1_ADRD.R`. This rounds exit ages to ~4 days, collapsing tied event times
and shrinking the expansion substantially. Report anything produced this way as
a sensitivity analysis, not as the primary result.

### Output

| File | Contents |
|---|---|
| `finegray_project1_main.csv` | 4 outcomes × 2 models × 2 BA = 16 rows |
| `finegray_project1_svedem.csv` | 10 outcomes × 2 models × 2 BA = 40 rows |
| `finegray_project1_*.log` | timestamped progress, timings, memory |

Semicolon-separated with comma decimals and Latin1 encoding, same as the Cox
output. Read back with `read.csv2()`. Section 10.0 of the Rmd merges them.

A row with a non-empty `error` column means that model failed; everything else
still ran. Check `error` before interpreting results.

---

## The analysis

**Outcomes.** Main: `ADRD`, `Alzheimers`, `Vascular`, `Other`.
SveDem: `S_ADRD`, `S_e_alz`, `S_l_alz`, `S_alz_vas`, `S_vas`, `S_lew`,
`S_fro`, `S_par`, `S_uns`, `S_oth`.

**BA measures.** `kdm_advance_sd` and `phenoage_advance_sd`, each in its own
model.

**Adjustment sets**, matching `CoxPH_B.Rmd`:

- Sex-adjusted: `gender`
- Multivariable: `gender + income_bins + NewUtbNiva`

**Timescale.** `Surv(AgeAtFirstDate, AgeAtLastDate, fg_status)` — attained age
with delayed entry, as in the Cox models.

---

## Decisions worth knowing about

### The four timeline ages

Cox only ever needed one exit column, so `CoxPH_B.Rmd` overwrites
`AgeAtLastDate` with the age at dementia diagnosis for anyone who has one.
That is fine for Cox — everyone who is not a case is simply censored — but it
destroys exactly what Fine-Gray needs: **when the competing death happened**.

`derive_fg_timeline()` splits the timeline back apart:

| Column | Source |
|---|---|
| `age_entry` | `AgeAtFirstDate`, unchanged |
| `age_dementia` | `AgeAtLastDate`, **reused as-is** where it really is a diagnosis age |
| `age_death` | recomputed from `DODSDATn` and `FODDATn` |
| `age_censor` | `lastDate` / `d2`, which was never overwritten |

`age_dementia` deliberately **reuses** `AgeAtLastDate` rather than recomputing
it, so the dementia timing is identical to the Cox models by construction. The
ATC/`FDATUM` rule for the meds-only group is inherited automatically and cannot
drift out of sync with `CoxPH_B.Rmd`.

The one thing that *is* recomputed is **which rows that reuse is valid for**.
`CoxPH_B.Rmd` only overwrites `AgeAtLastDate` when `INDATUM`/`S_INDATUM` exists,
or `FDATUM` for the meds-only group. Someone with `ADRD == 1` but none of those
keeps the `lastDate` age, which is not a diagnosis age at all — so the date
source is re-derived with the same precedence rather than trusting `ADRD == 1`.
Such rows are logged as a warning and treated as non-cases;
`fg_project1_ADRD.R inspect` section 1 counts them before you run anything.

Exit time and status are then composed per outcome on a **first-event-wins**
basis: whichever of diagnosis, death or end of follow-up comes first. Ties go
to dementia, then death. This is correct even when a death is recorded before
the diagnosis date.

### Other dementia subtypes are competing events

Because the registers record only one dementia type per person, a vascular case
can never become an Alzheimer case — their hazard for that outcome is zero.
Leaving them in the risk set would pad it with people who cannot have the
event, which is exactly what competing risks exists to handle.

This is also the only option that works if follow-up stops at diagnosis: if
`d2` is the diagnosis date for dementia cases, there is no post-diagnosis
follow-up to censor them into. `fg_project1_ADRD.R inspect` section 3 checks this.

`ADRD` and `S_ADRD` are unaffected — there is no "other subtype" there, so the
rule reduces to: event, else death, else censored. The test suite asserts that
the two settings give bit-identical `ADRD` coding.

As a check, the counts should satisfy:
`competing_deaths(subtype) = competing_deaths(ADRD) + (ADRD events − subtype events)`,
and `censored` should be the same for every outcome.

**Sensitivity analysis:** set `fg_config$competing_other_dementia <- FALSE` to
censor other subtypes at their diagnosis age instead. That is informative
censoring, so treat it as a robustness check rather than a primary result.

### Covariate coding matches CoxPH_B.Rmd

`income_bins` and `NewUtbNiva` are **numeric**, exactly as in the Cox models, so
the hazard ratios are directly comparable. This also happens to help memory:
as factors, `NewUtbNiva` alone would add ~5 dummy columns to every row of the
expansion.

`gender` is left as stored; `coxph` and `finegray` coerce character to factor
themselves.

### One expansion per outcome

Each outcome is expanded **once** and all four models are fitted on that single
expansion. This is valid because `finegray()` derives its weights from
`(entry, exit, status, id)` only — the right-hand-side variables are merely
carried along.

This is what makes the run feasible: 14 expansions instead of 28, and the
expansion dominates both runtime and memory. The claim is verified by
`test_fg_project1_ADRD.R` on synthetic data and by section 6.0 of the Rmd on the
real cohort. **If section 6.0 ever returns `FALSE`, stop and expand per model.**

### Shared complete-case sample

Both adjustment sets and both BA measures use the same rows: complete cases on
`AgeAtFirstDate`, `AgeAtLastDate`, both BA measures, and all covariates.

So `N` is identical across all four models for a given outcome, which makes
them directly comparable — but the sex-adjusted `N` will be **smaller** than
the corresponding sex-adjusted Cox model in `CoxPH_B.Rmd`, which only needed
`gender` to be non-missing. This is a deliberate consequence of sharing one
expansion.

### Robust standard errors

`cluster(fg_id)` with robust standard errors is the default, which is correct
for weighted Fine-Gray data. If the robust sandwich estimator runs out of
memory, the code falls back to the model-based variance and records
`robust = FALSE` in that row. Those standard errors are anti-conservative —
check the `robust` column before reporting.

### Weights above 1 are expected

With age as the timescale there is delayed entry, and `finegray()` weights
under left truncation are ratios of censoring survival that routinely exceed 1.
Only positivity is guaranteed. Not a bug.

---

## Files

| File | Purpose |
|---|---|
| `fg_project1_ADRD.R` | Everything: config, timeline derivation, cleaning, expansion, fitting, resume logic, plus the `main` / `svedem` / `inspect` run modes |
| `test_fg_project1_ADRD.R` | 80 checks on synthetic data; `source()`s `fg_project1_ADRD.R` without triggering a real run |
| `run_fg_project1.ps1` | Optional: one process per outcome, maximum memory safety (Windows) |
| `fine_gray_competing_risk_project1_attained_age_with_tests.Rmd` | Optional interactive notebook: same checks plus an expansion-size probe. Not needed to run the analysis - the terminal is the preferred path either way. |

`fg_project1_ADRD.R` used to be four separate files (`fg_project1_helpers.R`,
`fg_project1_main.R`, `fg_project1_svedem.R`, `inspect_data_for_fg.R`). They were
merged for a smaller footprint; nothing about the logic changed, only how it's
invoked - see the "Command-line entry point" section at the bottom of the file.

## Requirements

R with `survival` (tested against 3.8.3). No other packages.
