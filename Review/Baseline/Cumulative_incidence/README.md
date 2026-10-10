# Fine-Gray competing risk analysis — project 1 (dementia)

Fine-Gray subdistribution hazard models for dementia, with **death as the
competing event** and **attained age as the timescale**, mirroring the Cox
models in `CoxPH_B.Rmd` from the `Clinical_BA_ADRD` repo.

Adapted from the cancer side project's Fine-Gray scripts.

---

## Running it on the secure analysis platform

Put these three files in a folder together with `coxph_B_correctsubtypes.RData`:

```
fg_project1.R              config, all logic, and the three run modes below
test_fg_project1.R         80 synthetic checks, no real data needed
run_fg_project1.ps1        optional: one process per outcome (Windows)
```

`fg_project1.R` takes a mode as its first argument:

```powershell
Rscript test_fg_project1.R           # ~30 s, confirms nothing broke in transit
Rscript fg_project1.R inspect        # data diagnostics, run once before trusting anything
Rscript fg_project1.R main           # 4 outcomes
Rscript fg_project1.R svedem         # 10 outcomes
Rscript fg_project1.R cif            # absolute risk; minutes, no expansion
```

No editing required. The script finds the RData file itself. No data
preparation step is needed either — the four timeline ages are derived from the
existing `coxph_B_correctsubtypes.RData` on the fly (see *The four timeline ages*
below).

`fg_project1.R inspect` prints only counts and summaries, no individual-level
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
Rscript fg_project1.R main ADRD
Rscript fg_project1.R svedem S_vas
```

**`fg_config$time_round_digits` is already set to 1** — this was not a default,
it is what the published run used (`finegray_project1_main.log`, 2026-09-07),
because the exact expansion would not fit in memory. It rounds exit ages to
~5 weeks, collapsing tied event times. 3 ≈ one day, 2 ≈ four days, 1 ≈ five
weeks, `NULL` = exact.

Leave it alone. The `cif` mode inherits it so the cumulative incidence is
computed on exactly the timeline behind Supplementary Table 9; changing it would
silently desynchronise the two.

Rounding can in principle create ties, and ties resolve to dementia before
death, so it could move someone between event categories. `compare_rounding_effect()`
checks whether it actually did — the `cif` mode runs it automatically and writes
`cif_project1_rounding_check.csv`. If the counts are identical under rounded and
exact timelines, the rounding only collapsed tied times and reclassified nobody,
and Supplementary Table 9's counts are unaffected by it.

### Output

| File | Contents |
|---|---|
| `finegray_project1_main.csv` | 4 outcomes × 2 models × 2 BA = 16 rows |
| `finegray_project1_svedem.csv` | 10 outcomes × 2 models × 2 BA = 40 rows |
| `finegray_project1_*.log` | timestamped progress, timings, memory |
| `cif_project1_curves.csv` | cumulative incidence on a time grid, long format |
| `cif_project1_summary.csv` | absolute risk differences at the reporting times |
| `cif_project1_rounding_check.csv` | exit-age rounding: excluded vs reclassified |
| `cif_project1.log` | group merge diagnostics, rule check, timings |

Figures are a separate step. `plot_cif_results.R` reads
`cif_project1_curves.csv` — aggregate, so it runs inside or outside the secure analysis platform —
and writes `figure_cif_ADRD_*.pdf` plus `README_cif_figures.md`, which carries
the caption, the axis units and the key numbers. Keeping it separate means the
figures can be reworked without re-running the analysis.

```powershell
Rscript plot_cif_results.R .
```

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
`fg_project1.R inspect` section 1 counts them before you run anything.

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
follow-up to censor them into. `fg_project1.R inspect` section 3 checks this.

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
`test_fg_project1.R` on synthetic data and by section 6.0 of the Rmd on the
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

## Cumulative incidence (the `cif` mode)

Added for **reviewer comment 4**, which asks for absolute risk measures to
illustrate the competing risk of death. Supplementary Table 9 reports
subdistribution hazard ratios, which are still a *relative* measure — so it does
not answer the comment on its own.

**This mode cannot change the Fine-Gray results.** It is a different estimator,
not a different Fine-Gray run: Aalen-Johansen cumulative incidence computed with
`survfit()` on the output of `clean_fg_data()`. The sample, the
first-event-wins coding and the "other subtypes are competing events" rule are
shared with the Fine-Gray models by construction, but `finegray()` is never
called. No expansion, no memory wall, separate output files, and the published
sHRs are untouched. The test suite asserts this mechanically rather than
trusting it.

```powershell
Rscript fg_project1.R cif            # all 4 main outcomes
Rscript fg_project1.R cif ADRD       # just one
```

### What it needs

`table1_kdmgroups_baseline_perID.csv` — the per-ID biologically-younger /
biologically-older split already used for Table 1. Put it next to the script, or
set `STUDY1_OUTPUT_DIR` so `cif_config$group_file_candidates` finds it where
`Table1_B.Rmd` wrote it. Nothing is hard-coded:
the ID column is whichever name the file and the cohort share, matched
case-insensitively, and the group column is found by its values.

### The PhenoAge groups are derived, not invented

Table 1's split is for KDM. Rather than guess whether it is on the clock value or
the age gap, every plausible rule is scored against the supplied KDM groups —
`advance > 0`, `advance_sd > 0`, `clock > chronological age`, and the two median
splits. Only a rule that reproduces them **exactly** is accepted, and that same
rule is then applied to PhenoAge. If nothing matches exactly the run stops and
says so, including when it is a near miss, because a near miss usually means the
same rule computed on a slightly different sample.

Several rules can be exactly equivalent on this data — `advance > 0` and
`clock > chronological age` are the same statement when the advance is defined as
the difference. That is fine, as long as they still agree once transferred to
PhenoAge; if they diverge there they were only coincidentally equal, and the run
stops.

### Both event types, and the naive comparison

Three curves per group:

- **dementia**, Aalen-Johansen
- **death**, Aalen-Johansen — the half that actually answers the comment. The
  biologically older group dies faster, which depletes the at-risk pool and is
  *why* the subdistribution HRs attenuate.
- **naive 1 − KM**, treating death as censoring. This is exactly the mistake the
  reviewer warns about, so reporting it quantifies how much it matters. It is
  always ≥ the Aalen-Johansen curve.

### Two time axes

| Axis | What it is | Notes |
|---|---|---|
| `age` | attained age, anchored at `cif_config$index_age` (default 65) | The models' own timescale. A CIF "from time zero" is undefined under delayed entry, so people still under observation at the index age start their clock there and the curve reads as "cumulative incidence between 65 and *t*". |
| `years` | time since baseline | No left truncation. Matches the convention of the UK Biobank omics paper, but it averages absolute risks over whatever baseline ages are in the cohort. |

Both are produced every run. Which one leads the supplementary figure is
`cif_config$primary_axis`; decide it from the **baseline age by group** block the
log prints before the curves. If the two groups enter at similar ages the
years axis is defensible; if not, the age axis has to lead.

This is also why `survfit()` is used rather than `cmprsk::cuminc()`: `cuminc()`
takes only `(ftime, fstatus)` and has no delayed-entry interface, so it cannot do
the attained-age axis at all. `survfit()` also keeps the dependency list at
`survival` alone.

### Output is aggregate only

Same discipline as `inspect`: grid rows with `n_risk` below
`cif_config$min_risk_set` (default 10) are suppressed, so both CSVs can be copied
out of the secure analysis platform and the figure rebuilt outside it from the numbers alone. The PNGs
are produced in base graphics for the same reason the rest of this uses no extra
packages.

### Check before trusting it

The log prints, in this order: how many cohort rows matched the group file (the
unmatched are excluded and counted, never silently dropped); which level is being
treated as biologically older and why; the agreement of every candidate rule;
the group sizes; and baseline age by group. Read all of it. Getting the older
level backwards would reverse the figure without erroring.

---

## Files

| File | Purpose |
|---|---|
| `fg_project1.R` | Everything: config, timeline derivation, cleaning, expansion, fitting, resume logic, plus the `main` / `svedem` / `inspect` / `cif` run modes |
| `test_fg_project1.R` | 136 checks on synthetic data; `source()`s `fg_project1.R` without triggering a real run |
| `table1_kdmgroups_baseline_perID.csv` | Input for the `cif` mode only — Table 1's biological age groups, per ID |
| `run_fg_project1.ps1` | Optional: one process per outcome, maximum memory safety (Windows) |
| `fine_gray_competing_risk_project1_attained_age_with_tests.Rmd` | Optional interactive notebook: same checks plus an expansion-size probe. Not needed to run the analysis - the terminal is the preferred path either way. |

`fg_project1.R` used to be four separate files (`fg_project1_helpers.R`,
`fg_project1_main.R`, `fg_project1_svedem.R`, `inspect_data_for_fg.R`). They were
merged for a smaller footprint; nothing about the logic changed, only how it's
invoked - see the "Command-line entry point" section at the bottom of the file.

## Requirements

R with `survival` (tested against 3.8.3). No other packages.
