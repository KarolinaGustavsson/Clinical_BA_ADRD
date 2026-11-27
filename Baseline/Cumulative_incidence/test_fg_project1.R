#!/usr/bin/env Rscript
# =====================================================================
# Tests for the project 1 Fine-Gray helpers
# =====================================================================
#
# Runs entirely on small synthetic data - no real data, no VDI needed.
# Run before shipping to the VDI, and after any edit to the helpers:
#
#   Rscript test_fg_project1.R
#
# Exits non-zero if anything fails.

# Source the implementation file belonging to THIS test file, by stripping the
# leading "test_" from our own name. So test_fg_project1_cif.R sources
# fg_project1_cif.R even when an older fg_project1.R sits in the same folder -
# otherwise a renamed pair could silently test the wrong copy, which during a
# transfer to the VDI is exactly when it would matter most.
.tf_cli    <- commandArgs(trailingOnly = FALSE)
.tf_file   <- sub("^--file=", "", .tf_cli[grepl("^--file=", .tf_cli)])
.tf_target <- if (length(.tf_file) > 0) {
  sub("^test_", "", sub(".*[\\\\/]", "", .tf_file))
} else "fg_project1.R"
if (!file.exists(.tf_target)) {
  if (!file.exists("fg_project1.R")) {
    stop("Cannot find the script to test. Looked for '", .tf_target,
         "' and 'fg_project1.R' in ", getwd())
  }
  .tf_target <- "fg_project1.R"
}
cat("Testing:", .tf_target, "\n")
source(.tf_target)

failures <- 0L
checks   <- 0L

check <- function(label, condition, detail = "") {
  checks <<- checks + 1L
  ok <- isTRUE(condition)
  cat(if (ok) "  PASS  " else "  FAIL  ", label, "\n", sep = "")
  if (!ok) {
    failures <<- failures + 1L
    if (nzchar(detail)) cat("        ", detail, "\n", sep = "")
  }
  invisible(ok)
}

# Synthetic cohort already carrying the four timeline ages.
make_test_data <- function(n = 800, seed = 42) {
  set.seed(seed)
  d <- data.frame(
    age_entry = runif(n, 45, 70),
    kdm_advance_sd      = rnorm(n),
    phenoage_advance_sd = rnorm(n),
    gender      = sample(c("1", "2"), n, replace = TRUE),
    income_bins = sample(1:4, n, replace = TRUE),
    NewUtbNiva  = sample(1:6, n, replace = TRUE)
  )
  d$age_censor <- d$age_entry + rexp(n, 1 / 12)

  fate <- sample(
    c("alz", "vas", "oth", "dead", "censor"), n, replace = TRUE,
    prob = c(0.08, 0.05, 0.05, 0.35, 0.47)
  )
  d$Alzheimers <- as.integer(fate == "alz")
  d$Vascular   <- as.integer(fate == "vas")
  d$Other      <- as.integer(fate == "oth")
  d$ADRD       <- as.integer(fate %in% c("alz", "vas", "oth"))

  # Dementia cases are diagnosed somewhere inside their follow-up.
  d$age_dementia <- ifelse(d$ADRD == 1,
                           d$age_entry + (d$age_censor - d$age_entry) * runif(n, 0.1, 0.9),
                           NA_real_)
  # Deaths: everyone in the "dead" group, plus 60% of dementia cases later on.
  dies <- fate == "dead" | (d$ADRD == 1 & runif(n) < 0.6)
  d$age_death <- ifelse(dies, d$age_censor, NA_real_)

  # Missingness the cleaner has to drop.
  d$income_bins[sample(n, 20)]    <- NA
  d$kdm_advance_sd[sample(n, 15)] <- NA
  # An impossible exit age that must also be dropped.
  d$age_censor[1] <- d$age_entry[1] - 1
  d$age_dementia[1] <- NA
  d$age_death[1] <- NA
  d
}

test_data <- make_test_data()

# ---------------------------------------------------------------------
cat("\n1. Event coding: first event wins\n")
# ---------------------------------------------------------------------

# Each row is a case the coding rules must get right. Outcome = Alzheimers.
coding_cases <- data.frame(
  age_entry    = rep(60, 7),
  age_dementia = c(70, 70, 70, NA, NA, 70, 70),
  age_death    = c(NA, 75, 75, 75, NA, NA, 65),
  age_censor   = c(80, 80, 80, 80, 80, 80, 80),
  Alzheimers   = c(1L, 1L, 0L, 0L, 0L, 0L, 1L),
  Vascular     = c(0L, 0L, 1L, 0L, 0L, 1L, 0L),
  ADRD         = c(1L, 1L, 1L, 0L, 0L, 1L, 1L),
  kdm_advance_sd = 0, phenoage_advance_sd = 0,
  gender = "1", income_bins = 1, NewUtbNiva = 1
)
coded <- clean_fg_data(coding_cases, "Alzheimers", "test")

check("own subtype, alive -> event",
      coded$fg_status[1] == "event")
check("own subtype then death -> event at diagnosis, not death",
      coded$fg_status[2] == "event" && coded$AgeAtLastDate[2] == 70)
check("other subtype then death -> competing at diagnosis age",
      coded$fg_status[3] == "death" && coded$AgeAtLastDate[3] == 70)
check("no dementia, died -> competing death at the death age",
      coded$fg_status[4] == "death" && coded$AgeAtLastDate[4] == 75)
check("no dementia, alive -> censored at end of follow-up",
      coded$fg_status[5] == "censor" && coded$AgeAtLastDate[5] == 80)
check("other subtype, alive -> competing at diagnosis age",
      coded$fg_status[6] == "death" && coded$AgeAtLastDate[6] == 70)
check("death recorded before diagnosis -> death wins",
      coded$fg_status[7] == "death" && coded$AgeAtLastDate[7] == 65)

# ADRD has no "other subtype", so the rule must reduce to the plain one.
adrd_coded <- clean_fg_data(coding_cases, "ADRD", "test")
check("ADRD: every dementia case is an event",
      all(adrd_coded$fg_status[c(1, 2, 3, 6)] == "event"))
check("ADRD: death without dementia is competing",
      adrd_coded$fg_status[4] == "death")
check("ADRD: death before diagnosis is still a competing death",
      adrd_coded$fg_status[7] == "death")

# The sensitivity setting: other subtypes censored instead of competing.
censored_variant <- clean_fg_data(coding_cases, "Alzheimers", "test",
                                  competing_other_dementia = FALSE)
check("competing_other_dementia=FALSE censors other subtypes",
      all(censored_variant$fg_status[c(3, 6)] == "censor"))
check("competing_other_dementia=FALSE leaves real deaths competing",
      censored_variant$fg_status[4] == "death")
check("competing_other_dementia=FALSE leaves own-subtype events alone",
      all(censored_variant$fg_status[c(1, 2)] == "event"))
check("ADRD is identical under both competing-event settings",
      identical(as.character(adrd_coded$fg_status),
                as.character(clean_fg_data(coding_cases, "ADRD", "test",
                                           competing_other_dementia = FALSE)$fg_status)))

# ---------------------------------------------------------------------
cat("\n2. Timeline derivation from CoxPH_B-style columns\n")
# ---------------------------------------------------------------------

# Raw layout as it comes out of CoxPH_B.Rmd, where AgeAtLastDate has already
# been overwritten with the diagnosis age for anyone who has a diagnosis date.
born <- as.Date("1930-01-01")
raw <- data.frame(
  FODDATn   = rep(born, 5),
  firstDate = rep(as.Date("1990-01-01"), 5),                      # age 60
  lastDate  = as.Date(c("2010-01-01", "2005-01-01", "2010-01-01",
                        "2000-01-01", "2010-01-01")),
  DODSDATn  = as.Date(c(NA, "2005-01-01", NA, "2000-01-01", NA)),
  INDATUM   = as.Date(c("2000-01-01", NA, NA, NA, NA)),
  S_INDATUM = as.Date(c(NA, "1998-01-01", NA, NA, NA)),
  FDATUM    = as.Date(c(NA, NA, "1995-01-01", NA, NA)),
  ATC       = c("", "", "N06DA", "", ""),
  Dementia  = c("G30", "", "", "", ""),
  S_ADRD    = c(0L, 1L, 0L, 0L, 0L),
  ADRD      = c(1L, 1L, 1L, 0L, 1L),   # row 5: ADRD=1 but NO date anywhere
  AgeAtFirstDate = rep(60, 5),
  # What CoxPH_B would have produced: diagnosis age where a date existed,
  # otherwise the lastDate age.
  AgeAtLastDate = c(70, 68, 65, 70, 80)
)

tl <- derive_fg_timeline(raw, quiet = TRUE)

check("timeline returns exactly the four age columns",
      identical(names(tl), timeline_columns))
check("age_entry is taken from AgeAtFirstDate",
      all(tl$age_entry == 60))
check("age_dementia reuses AgeAtLastDate for an ICD diagnosis",
      tl$age_dementia[1] == 70)
check("age_dementia reuses AgeAtLastDate for a SveDem diagnosis",
      tl$age_dementia[2] == 68)
check("age_dementia reuses AgeAtLastDate for the meds-only group (FDATUM)",
      tl$age_dementia[3] == 65,
      paste("got", tl$age_dementia[3]))
check("age_dementia is NA when there is no dementia at all",
      is.na(tl$age_dementia[4]))
check("ADRD==1 with no date anywhere does NOT get a dementia age",
      is.na(tl$age_dementia[5]),
      "otherwise the lastDate age would be mistaken for a diagnosis age")

check("age_death is recovered from DODSDATn",
      isTRUE(all.equal(tl$age_death[2], 75, tolerance = 1e-3)),
      paste("got", tl$age_death[2]))
check("age_death is NA for the living",
      is.na(tl$age_death[1]) && is.na(tl$age_death[3]))
check("age_censor comes from lastDate, which was never overwritten",
      isTRUE(all.equal(tl$age_censor[1], 80, tolerance = 1e-3)),
      paste("got", tl$age_censor[1]))
check("age_censor is unaffected by the diagnosis overwrite",
      isTRUE(all.equal(tl$age_censor[3], 80, tolerance = 1e-3)))

check("already-derived timeline columns are passed through untouched", {
  pre <- data.frame(age_entry = 1, age_dementia = 2, age_death = 3, age_censor = 4)
  identical(derive_fg_timeline(pre, quiet = TRUE), pre)
})
check("a missing date column raises a clear error", {
  err <- tryCatch(derive_fg_timeline(raw[, setdiff(names(raw), "DODSDATn")], quiet = TRUE),
                  error = conditionMessage)
  is.character(err) && grepl("death", err)
})

# End to end from raw columns: row 2 was diagnosed at 68 and died at 75, so
# the diagnosis must win and the death must not be silently used as the exit.
raw_full <- cbind(raw, kdm_advance_sd = 0, phenoage_advance_sd = 0,
                  gender = "1", income_bins = 1, NewUtbNiva = 1,
                  Alzheimers = c(1L, 0L, 1L, 0L, 0L))
raw_coded <- clean_fg_data(raw_full, "Alzheimers", "test")
check("raw pipeline: diagnosis beats a later death",
      raw_coded$fg_status[2] == "death" && raw_coded$AgeAtLastDate[2] == 68)
check("raw pipeline: own subtype is an event at the diagnosis age",
      raw_coded$fg_status[1] == "event" && raw_coded$AgeAtLastDate[1] == 70)

# ---------------------------------------------------------------------
cat("\n3. Cleaning\n")
# ---------------------------------------------------------------------

cleaned <- clean_fg_data(test_data, "Alzheimers", "test")

check("rows with missing covariates or BA are dropped",
      nrow(cleaned) < nrow(test_data))
check("no missing values remain in model variables",
      !anyNA(cleaned[, c("AgeAtFirstDate", "AgeAtLastDate", fg_config$ba_measures,
                         all_model_covariates)]))
check("exit age always after entry age",
      all(cleaned$AgeAtLastDate > cleaned$AgeAtFirstDate))
check("fg_id is 1..N with no gaps after row dropping",
      identical(cleaned$fg_id, seq_len(nrow(cleaned))))
check("income_bins stays numeric (matches CoxPH_B)",
      is.numeric(cleaned$income_bins))
check("NewUtbNiva stays numeric (matches CoxPH_B)",
      is.numeric(cleaned$NewUtbNiva))
check("gender is coerced to a factor",
      is.factor(cleaned$gender))
check("all three status levels are present",
      all(c("censor", "event", "death") %in% levels(cleaned$fg_status)))

counts <- fg_counts(cleaned)
check("counts partition N exactly",
      counts$events + counts$competing_deaths + counts$censored == counts$N,
      paste(counts$events, "+", counts$competing_deaths, "+",
            counts$censored, "!=", counts$N))

check("missing outcome indicator is treated as not-an-event, not an error", {
  na_data <- test_data
  na_data$Alzheimers[1:10] <- NA
  nrow(clean_fg_data(na_data, "Alzheimers", "test")) > 0
})

check("a missing required column raises a clear error", {
  bad <- test_data[, setdiff(names(test_data), "NewUtbNiva")]
  err <- tryCatch(clean_fg_data(bad, "Alzheimers", "test"), error = conditionMessage)
  is.character(err) && grepl("NewUtbNiva", err)
})

# ---------------------------------------------------------------------
cat("\n4. Expansion invariance (the one-expansion-per-outcome claim)\n")
# ---------------------------------------------------------------------

# The whole memory strategy rests on finegray() weights depending only on
# (entry, exit, status, id) and not on the carried covariates. Verify it
# rather than trusting it.
full    <- expand_fg_data(cleaned)
minimal <- expand_fg_data(cleaned, ba_measures = "kdm_advance_sd",
                          covariates = "gender")

check("expansion row count is independent of carried covariates",
      nrow(full) == nrow(minimal),
      paste(nrow(full), "vs", nrow(minimal)))
check("expansion weights are independent of carried covariates",
      isTRUE(all.equal(full$fgwt, minimal$fgwt)))
check("expansion risk intervals are independent of carried covariates",
      isTRUE(all.equal(full$fgstart, minimal$fgstart)) &&
      isTRUE(all.equal(full$fgstop,  minimal$fgstop)))

check("expanded data keeps only the columns coxph needs",
      setequal(names(full),
               c("fgstart", "fgstop", "fgstatus", "fgwt", "fg_id",
                 fg_config$ba_measures, all_model_covariates)))
check("expanded event count matches the cleaned event count",
      sum(full$fgstatus == 1) == counts$events,
      paste(sum(full$fgstatus == 1), "vs", counts$events))
# Note: weights are NOT bounded by 1 here. With age as the timescale there is
# delayed entry, and finegray() weights under left truncation are ratios of
# censoring survival that routinely exceed 1. Only positivity is guaranteed.
check("all weights are positive and finite",
      all(full$fgwt > 0) && all(is.finite(full$fgwt)))
check("expansion is at least as large as the input",
      nrow(full) >= nrow(cleaned))

# ---------------------------------------------------------------------
cat("\n5. Model fitting\n")
# ---------------------------------------------------------------------

fit_mv <- fit_with_fallback(full, "kdm_advance_sd",
                            adjustment_sets[["Multivariable model"]], "test")

check("multivariable fit returns a finite beta",
      is.finite(fit_mv$beta), paste("error:", fit_mv$error))
check("multivariable fit returns a positive se",
      is.finite(fit_mv$se) && fit_mv$se > 0)
check("robust estimator was used", isTRUE(fit_mv$robust))
check("no error recorded on a clean fit", is.na(fit_mv$error))

# Compare against a directly written model to catch formula-assembly bugs.
reference <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
    kdm_advance_sd + gender + income_bins + NewUtbNiva + cluster(fg_id),
  weights = fgwt, data = full, robust = TRUE
)
check("beta matches a hand-written reference model",
      isTRUE(all.equal(fit_mv$beta, unname(coef(reference)["kdm_advance_sd"]))),
      paste(fit_mv$beta, "vs", coef(reference)["kdm_advance_sd"]))
check("se matches the robust sandwich se of the reference model",
      isTRUE(all.equal(fit_mv$se,
                       unname(sqrt(diag(vcov(reference)))["kdm_advance_sd"]))),
      paste(fit_mv$se, "vs", sqrt(diag(vcov(reference)))["kdm_advance_sd"]))

fit_sex <- fit_with_fallback(full, "phenoage_advance_sd",
                             adjustment_sets[["Sex-adjusted model"]], "test")
check("sex-adjusted fit works for the other BA measure",
      is.finite(fit_sex$beta) && is.finite(fit_sex$se),
      paste("error:", fit_sex$error))

check("a non-existent BA variable is reported, not thrown", {
  bad <- fit_with_fallback(full, "not_a_variable", "gender", "test")
  is.na(bad$beta) && !is.na(bad$error)
})

# ---------------------------------------------------------------------
cat("\n6. Result rows\n")
# ---------------------------------------------------------------------

row <- make_result_row("Main", "ADRD", "Multivariable model", "kdm_advance_sd",
                       beta = log(1.5), se = 0.1, N = 100, events = 20)
check("result row has exactly the expected columns",
      identical(names(row), result_columns))
check("HR is exp(beta)", isTRUE(all.equal(row$HR, 1.5)))
check("confidence interval brackets the HR",
      row$CI_lower < row$HR && row$HR < row$CI_upper)
check("p-value is the two-sided normal test",
      isTRUE(all.equal(row$p, 2 * pnorm(-abs(log(1.5) / 0.1)))))

blank <- make_result_row("Main", "ADRD", "Multivariable model", "kdm_advance_sd",
                         error = "expansion failed")
check("a failed row carries NA estimates rather than garbage",
      is.na(blank$HR) && is.na(blank$p) && !is.na(blank$error))

# ---------------------------------------------------------------------
cat("\n7. Incremental writing and resume\n")
# ---------------------------------------------------------------------

tmp_dir  <- tempfile("fgtest"); dir.create(tmp_dir)
tmp_file <- file.path(tmp_dir, "results.csv")
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

check("no completed keys before the file exists",
      length(completed_keys(tmp_file)) == 0)

append_result(row, tmp_file)
check("first append creates the file with a header",
      file.exists(tmp_file) && length(readLines(tmp_file)) == 2)

append_result(make_result_row("Main", "ADRD", "Sex-adjusted model",
                              "kdm_advance_sd", beta = 0.1, se = 0.05),
              tmp_file)
check("second append adds a row without repeating the header",
      length(readLines(tmp_file)) == 3)

done <- completed_keys(tmp_file)
check("both written results are detected as complete",
      length(done) == 2 &&
      result_key("Main", "ADRD", "Multivariable model", "kdm_advance_sd") %in% done)
check("an unwritten result is not reported as complete",
      !(result_key("Main", "Vascular", "Multivariable model", "kdm_advance_sd") %in% done))

reread <- read.table(tmp_file, sep = ";", dec = ",", header = TRUE,
                     quote = "", fileEncoding = "Latin1")
check("numbers survive the semicolon/comma-decimal round trip",
      isTRUE(all.equal(reread$HR[1], 1.5, tolerance = 1e-6)),
      paste("got", reread$HR[1]))

# ---------------------------------------------------------------------
cat("\n8. End-to-end run with resume\n")
# ---------------------------------------------------------------------

e2e_file <- file.path(tmp_dir, "e2e.csv")
small    <- make_test_data(n = 300, seed = 7)
specs    <- list(list(outcome_label = "ADRD", event_var = "ADRD"),
                 list(outcome_label = "Vascular", event_var = "Vascular"))

sink(nullfile())
run_outcome_block(small, specs, "Main", e2e_file)
sink()

res <- read.table(e2e_file, sep = ";", dec = ",", header = TRUE,
                  quote = "", fileEncoding = "Latin1")
check("end-to-end produces 2 outcomes x 2 models x 2 BA = 8 rows",
      nrow(res) == 8, paste("got", nrow(res)))
check("every end-to-end model produced an estimate",
      all(is.finite(res$HR)),
      paste("failed:", paste(res$error[!is.finite(res$HR)], collapse = "; ")))
check("both BA measures appear in the results",
      setequal(res$BA, fg_config$ba_measures))
check("both adjustment sets appear in the results",
      setequal(res$model, names(adjustment_sets)))
check("N is identical across models within an outcome (shared complete-case sample)",
      length(unique(res$N[res$outcome == "ADRD"])) == 1)

# Simulate a crash-and-restart: drop the last 3 rows, rerun, expect them back.
writeLines(readLines(e2e_file)[1:6], e2e_file)
sink(nullfile())
run_outcome_block(small, specs, "Main", e2e_file)
sink()

res2 <- read.table(e2e_file, sep = ";", dec = ",", header = TRUE,
                   quote = "", fileEncoding = "Latin1")
check("restart refills exactly the missing rows",
      nrow(res2) == 8, paste("got", nrow(res2)))
check("restart does not duplicate completed work",
      !any(duplicated(result_key(res2$dataset, res2$outcome, res2$model, res2$BA))))
check("restarted results match the originals",
      isTRUE(all.equal(sort(res$HR), sort(res2$HR), tolerance = 1e-8)))

# A fully complete file should do no work at all.
before <- readLines(e2e_file)
sink(nullfile())
run_outcome_block(small, specs, "Main", e2e_file)
sink()
check("a complete file is left untouched on rerun",
      identical(before, readLines(e2e_file)))

# ---------------------------------------------------------------------
cat("\n9. Degenerate outcomes are recorded, not fatal\n")
# ---------------------------------------------------------------------

zero_file <- file.path(tmp_dir, "zero.csv")
no_events <- small
no_events$S_lew <- 0L

sink(nullfile())
run_outcome_block(no_events, list(list(outcome_label = "S_lew", event_var = "S_lew")),
                  "SveDem", zero_file)
sink()

zero <- read.table(zero_file, sep = ";", dec = ",", header = TRUE, quote = "",
                   fileEncoding = "Latin1")
check("an outcome with zero events writes rows instead of crashing",
      nrow(zero) == 4, paste("got", nrow(zero)))
check("zero-event rows carry an explanatory error",
      all(!is.na(zero$error)))

# ---------------------------------------------------------------------
cat("\n10. Time rounding (memory escape hatch)\n")
# ---------------------------------------------------------------------

# The unrounded reference has to be built explicitly: fg_config$time_round_digits
# is now 1, matching the published run, so `cleaned` and `full` above are already
# rounded and comparing against them would compare a thing with itself.
unrounded      <- clean_fg_data(test_data, "Alzheimers", "test",
                                time_round_digits = NULL)
unrounded_full <- expand_fg_data(unrounded)
rounded        <- clean_fg_data(test_data, "Alzheimers", "test",
                                time_round_digits = 1)

check("rounding reduces the number of distinct exit times",
      length(unique(rounded$AgeAtLastDate)) < length(unique(unrounded$AgeAtLastDate)))
check("rounding still yields a usable expansion",
      nrow(expand_fg_data(rounded)) > 0)
check("rounding shrinks the expansion",
      nrow(expand_fg_data(rounded)) < nrow(unrounded_full))
check("the default cleaning now applies the published rounding",
      identical(cleaned$AgeAtLastDate, rounded$AgeAtLastDate))

# ---------------------------------------------------------------------
cat("\n11. Cumulative incidence (the cif mode)\n")
# ---------------------------------------------------------------------

quietly <- function(expr) {
  sink(nullfile())
  on.exit(sink(), add = TRUE)
  force(expr)
}

# --- 11a. extra_vars cannot move the Fine-Gray results ----------------
# The cif mode needs the group labels carried through clean_fg_data(). That must
# be provably inert for the Fine-Gray path, because those results are already
# written up.

check("extra_vars = character(0) is identical to not passing it at all",
      identical(clean_fg_data(test_data, "Alzheimers", "test"),
                clean_fg_data(test_data, "Alzheimers", "test",
                              extra_vars = character(0))))

with_extra <- test_data
with_extra$some_group <- factor(
  ifelse(with_extra$kdm_advance_sd > 0, "older", "younger"))
with_extra$all_missing <- NA_character_

carried      <- clean_fg_data(with_extra, "Alzheimers", "test",
                              extra_vars = "some_group")
carried_na   <- clean_fg_data(with_extra, "Alzheimers", "test",
                              extra_vars = "all_missing")
baseline_out <- clean_fg_data(with_extra, "Alzheimers", "test")

fg_cols <- c("AgeAtFirstDate", "AgeAtLastDate", "fg_status", "fg_id",
             fg_config$ba_measures, all_model_covariates)

check("carrying an extra variable leaves every Fine-Gray column untouched",
      identical(baseline_out[, fg_cols], carried[, fg_cols]))
check("an all-missing extra variable does NOT shrink N",
      nrow(carried_na) == nrow(baseline_out),
      paste(nrow(carried_na), "vs", nrow(baseline_out)))
check("an all-missing extra variable leaves the Fine-Gray columns untouched",
      identical(baseline_out[, fg_cols], carried_na[, fg_cols]))
check("the extra variable actually arrives in the output",
      "some_group" %in% names(carried))
check("extra_vars overlapping covariates does not duplicate a column",
      !any(duplicated(names(clean_fg_data(with_extra, "Alzheimers", "test",
                                          extra_vars = "gender")))))

# --- 11b. Aalen-Johansen against a hand-computed example --------------
# Five people, one group, no ties. Worked through by hand:
#   t=1 event : dCIF_e = 1.000 * 1/5 = 0.2            S = 0.8
#   t=2 death : dCIF_d = 0.800 * 1/4 = 0.2            S = 0.6
#   t=3 event : dCIF_e = 0.600 * 1/3 = 0.2 -> 0.4     S = 0.4
#   t=4 censor: nothing changes                       S = 0.4
#   t=5 death : dCIF_d = 0.400 * 1/1 = 0.4 -> 0.6     S = 0.0
# so CIF_dementia = .2 .2 .4 .4 .4 and CIF_death = 0 .2 .2 .2 .6,
# and 0.4 + 0.6 + 0.0 = 1 exactly.

hand <- data.frame(
  AgeAtFirstDate = rep(0, 5),
  AgeAtLastDate  = c(1, 2, 3, 4, 5),
  fg_status = factor(c("event", "death", "event", "censor", "death"),
                     levels = c("censor", "event", "death")),
  fg_id = 1:5,
  grp = factor(rep("Biologically younger", 5),
               levels = c("Biologically younger", "Biologically older"))
)

hand_rows <- quietly(cif_for_outcome(hand, "grp", "years"))
at <- function(rows, st, t) {
  v <- rows$estimate[rows$state == st & rows$time == t]
  if (length(v) == 0) NA_real_ else v[1]
}

check("AJ dementia CIF matches the hand calculation",
      isTRUE(all.equal(vapply(1:5, function(t) at(hand_rows, "dementia", t), numeric(1)),
                       c(0.2, 0.2, 0.4, 0.4, 0.4), tolerance = 1e-8)),
      paste(round(vapply(1:5, function(t) at(hand_rows, "dementia", t), numeric(1)), 4),
            collapse = " "))
check("AJ death CIF matches the hand calculation",
      isTRUE(all.equal(vapply(1:5, function(t) at(hand_rows, "death", t), numeric(1)),
                       c(0.0, 0.2, 0.2, 0.2, 0.6), tolerance = 1e-8)),
      paste(round(vapply(1:5, function(t) at(hand_rows, "death", t), numeric(1)), 4),
            collapse = " "))
check("the two CIFs sum to 1 once everyone has left the risk set",
      isTRUE(all.equal(at(hand_rows, "dementia", 5) + at(hand_rows, "death", 5), 1,
                       tolerance = 1e-8)))
check("before the first event the CIF is 0 and the risk set is everyone",
      at(hand_rows, "dementia", 0.5) == 0 &&
      hand_rows$n_risk[hand_rows$state == "dementia" & hand_rows$time == 0.5][1] == 5)
check("a single-level grouping still produces one labelled stratum",
      identical(unique(hand_rows$group), "Biologically younger"))
check("per-group denominators are carried on every row",
      all(hand_rows$n_group == 5) &&
      all(hand_rows$events_group == 2) && all(hand_rows$deaths_group == 2))

# --- 11c. The naive 1 - KM comparison --------------------------------
# Treating death as censoring must overstate absolute risk, never understate it.

check("naive 1 - KM is >= the AJ CIF at every grid time", {
  aj <- hand_rows[hand_rows$state == "dementia", ]
  km <- hand_rows[hand_rows$state == "naive_km", ]
  m  <- merge(aj[, c("time", "estimate")], km[, c("time", "estimate")],
              by = "time", suffixes = c("_aj", "_km"))
  nrow(m) > 0 && all(m$estimate_km >= m$estimate_aj - 1e-10)
})
check("naive 1 - KM is strictly larger once a competing death has occurred",
      at(hand_rows, "naive_km", 3) > at(hand_rows, "dementia", 3) + 1e-8,
      paste(at(hand_rows, "naive_km", 3), "vs", at(hand_rows, "dementia", 3)))

# --- 11d. Shape guarantees on a larger synthetic cohort --------------

cif_input <- clean_fg_data(test_data, "Alzheimers", "test")
cif_input$grp <- factor(
  ifelse(cif_input$kdm_advance_sd > 0, "Biologically older", "Biologically younger"),
  levels = c("Biologically younger", "Biologically older"))

big_rows <- quietly(cif_for_outcome(cif_input, "grp", "years"))

check("every estimate is a probability",
      all(big_rows$estimate >= 0 & big_rows$estimate <= 1, na.rm = TRUE))
check("every CIF is monotone non-decreasing in time", {
  ok <- TRUE
  for (g in unique(big_rows$group)) {
    for (st in unique(big_rows$state)) {
      s <- big_rows[big_rows$group == g & big_rows$state == st, ]
      s <- s[order(s$time), ]
      if (nrow(s) > 1 && any(diff(s$estimate) < -1e-10)) ok <- FALSE
    }
  }
  ok
})
check("dementia and death CIFs never sum above 1",
      {
        d <- big_rows[big_rows$state == "dementia", c("group", "time", "estimate")]
        k <- big_rows[big_rows$state == "death",    c("group", "time", "estimate")]
        m <- merge(d, k, by = c("group", "time"), suffixes = c("_d", "_k"))
        all(m$estimate_d + m$estimate_k <= 1 + 1e-10)
      })
check("confidence limits bracket the estimate where they exist", {
  s <- big_rows[big_rows$state %in% c("dementia", "death") &
                is.finite(big_rows$CI_lower) & is.finite(big_rows$CI_upper), ]
  nrow(s) > 0 && all(s$CI_lower <= s$estimate + 1e-10) &&
    all(s$CI_upper >= s$estimate - 1e-10)
})
check("both states and the naive curve are all present",
      setequal(unique(big_rows$state), c("dementia", "death", "naive_km")))
check("the cumulative event count never decreases", {
  s <- big_rows[big_rows$state == "dementia" &
                big_rows$group == "Biologically older", ]
  s <- s[order(s$time), ]
  all(diff(s$n_event_cum) >= 0)
})

# --- 11e. The attained-age axis and its index age --------------------

age_cases <- data.frame(
  AgeAtFirstDate = c(50, 70, 60, 60),
  AgeAtLastDate  = c(62, 80, 90, 64),   # row 1 leaves before 65, row 4 after it
  fg_status = factor(c("event", "event", "death", "censor"),
                     levels = c("censor", "event", "death")),
  fg_id = 1:4,
  grp = factor(rep("Biologically younger", 4),
               levels = c("Biologically younger", "Biologically older"))
)
aged <- prepare_cif_axis(age_cases, "age", index_age = 65)

check("someone who already left before the index age is dropped",
      !(1 %in% aged$fg_id) && nrow(aged) == 2,
      paste("kept ids:", paste(aged$fg_id, collapse = ",")))
check("someone who entered before the index age starts there instead",
      aged$cif_start[aged$fg_id == 3] == 65)
check("someone who entered after the index age keeps their own entry age",
      aged$cif_start[aged$fg_id == 2] == 70)
check("the exit age is never altered by the index age",
      all(aged$cif_stop == age_cases$AgeAtLastDate[match(aged$fg_id, age_cases$fg_id)]))
check("the years axis starts everyone at zero", {
  y <- prepare_cif_axis(age_cases, "years")
  all(y$cif_start == 0) &&
    isTRUE(all.equal(y$cif_stop, age_cases$AgeAtLastDate - age_cases$AgeAtFirstDate))
})
check("an unknown axis is rejected",
      is.character(tryCatch(prepare_cif_axis(age_cases, "calendar"),
                            error = conditionMessage)))
check("the age grid starts at the index age and includes the reporting ages",
      min(cif_grid("age", 65)) == 65 &&
      all(cif_config$report_ages %in% cif_grid("age", 65)))
check("the years grid includes the reporting years",
      all(cif_config$report_years %in% cif_grid("years")))

# --- 11f. Reusing Table 1's grouping rather than inventing one --------

make_cif_test_data <- function(n = 600, seed = 11) {
  d <- make_test_data(n = n, seed = seed)
  d$lopnr <- seq_len(nrow(d))
  # The real cohort and the real Table 1 file both key on sampleID.
  d$sampleID <- d$lopnr
  d$AgeAtFirstDate <- d$age_entry
  # Mean deliberately away from 0 so "advance > 0" and "advance_sd > 0" are
  # genuinely different rules - otherwise the rule check proves nothing.
  d$kdm_advance      <- d$kdm_advance_sd * 3 + 1
  d$kdm              <- d$AgeAtFirstDate + d$kdm_advance
  d$phenoage_advance <- d$phenoage_advance_sd * 2 - 0.5
  d$phenoage         <- d$AgeAtFirstDate + d$phenoage_advance
  d
}

cif_data <- make_cif_test_data()
truth    <- cif_data$kdm_advance > 0          # the rule Table 1 "used"

# Only the first 550 IDs with a usable KDM value, so unmatched rows are
# exercised too. Column names deliberately cased differently from the cohort
# data, since the real group file and the real cohort spell the ID differently.
in_file <- cif_data$lopnr <= 550 & !is.na(truth)
groups_file <- data.frame(
  LopNr = cif_data$lopnr[in_file],
  kdm_group = ifelse(truth[in_file], "Biologically older", "Biologically younger"),
  stringsAsFactors = FALSE
)

cn <- quietly(identify_group_columns(groups_file, cif_data))
check("the ID column is found despite differing capitalisation",
      cn$id_in_groups == "LopNr" && cn$id_in_data == "lopnr")
check("the group column is identified",
      cn$group_col == "kdm_group")

check("the biologically older level is read from the label",
      quietly(group_older_level(groups_file$kdm_group,
                               cif_data$kdm_advance_sd[in_file])) ==
        "Biologically older")
check("the older level is still found when the labels are uninformative", {
  coded <- ifelse(truth[in_file], "B", "A")
  quietly(group_older_level(coded, cif_data$kdm_advance_sd[in_file])) == "B"
})

verdict <- quietly(verify_group_rule(cif_data, ifelse(in_file, truth, NA),
                                     prefix = "kdm"))
check("the true grouping rule is reproduced exactly",
      "advance > 0" %in% verdict$exact,
      paste("exact:", paste(verdict$exact, collapse = ", ")))
check("a rule that is not the truth is not reported as exact",
      !("advance_sd > 0" %in% verdict$exact))
check("every candidate rule is scored and reported",
      nrow(verdict$scores) >= 4 && all(verdict$scores$pct <= 100))

check("a random grouping matches no rule exactly", {
  set.seed(99)
  noise <- sample(c(TRUE, FALSE), nrow(cif_data), replace = TRUE)
  v <- quietly(verify_group_rule(cif_data, noise, prefix = "kdm"))
  length(v$exact) == 0
})

attached <- quietly(attach_cif_groups(cif_data, groups_file))

check("attaching the groups adds both grouping columns",
      all(c("kdm_group", "phenoage_group") %in% names(attached)))
check("the KDM groups come from the file verbatim",
      all((attached$kdm_group == "Biologically older") == truth, na.rm = TRUE))
check("people absent from the group file get NA, not a guessed group",
      all(is.na(attached$kdm_group[!in_file])) &&
      !any(is.na(attached$kdm_group[in_file])))
check("the PhenoAge grouping uses the same rule, on PhenoAge values",
      all((attached$phenoage_group == "Biologically older") ==
            (attached$phenoage_advance > 0), na.rm = TRUE))
check("PhenoAge groups are only assigned where a KDM group exists",
      identical(is.na(attached$kdm_group), is.na(attached$phenoage_group)))
check("the grouping levels put 'younger' first, so it is the reference",
      identical(levels(attached$kdm_group),
                c("Biologically younger", "Biologically older")))
check("the rule that was used is recorded on the result",
      !is.null(attr(attached, "cif_rule")))

# The real Table 1 split turned out to be a MEDIAN split of the age gap, keyed
# on sampleID and labelled "High Age Gap" / "Low Age Gap". Both median rules
# reproduce it exactly, so the transfer has to cope with more than one match.
med_truth <- cif_data$kdm_advance > median(cif_data$kdm_advance, na.rm = TRUE)
med_file  <- data.frame(
  sampleID  = cif_data$sampleID[!is.na(med_truth)],
  KDM_group = ifelse(med_truth[!is.na(med_truth)], "High Age Gap", "Low Age Gap"),
  stringsAsFactors = FALSE
)

check("the sampleID key links without falling back to a guess", {
  cn2 <- quietly(identify_group_columns(med_file, cif_data))
  cn2$id_in_groups == "sampleID" && cn2$id_in_data == "sampleID" &&
    cn2$group_col == "KDM_group"
})
check("'High Age Gap' is recognised as the biologically older level",
      quietly(group_older_level(med_file$KDM_group,
                                cif_data$kdm_advance_sd[!is.na(med_truth)])) ==
        "High Age Gap")

med_attached <- quietly(attach_cif_groups(cif_data, med_file))

check("a median-split grouping is reproduced exactly",
      all((med_attached$kdm_group == "Biologically older") == med_truth,
          na.rm = TRUE))
check("a median split marks about half the cohort as older", {
  p <- mean(med_attached$kdm_group == "Biologically older", na.rm = TRUE)
  p > 0.45 && p < 0.55
})
check("two equally exact median rules do not block the transfer",
      !all(is.na(med_attached$phenoage_group)))

# Prevalence and agreement are different quantities and must both be reported:
# a rule can agree with Table 1 for ~98% of people while calling only ~50% of
# them older. Conflating the two is what makes a zero-split look implausible.
med_verdict <- quietly(verify_group_rule(
  cif_data, ifelse(is.na(med_truth), NA, med_truth), prefix = "kdm"))

check("the rule table reports prevalence as well as agreement",
      "marks_older_pct" %in% names(med_verdict$scores))
check("a median-split rule marks ~50% older however well it agrees", {
  r <- med_verdict$scores[med_verdict$scores$rule == "advance > median", ]
  r$pct == 100 && r$marks_older_pct > 45 && r$marks_older_pct < 55
})
check("a zero-split rule also marks ~50% older, not ~98%", {
  r <- med_verdict$scores[med_verdict$scores$rule == "advance > 0", ]
  r$marks_older_pct > 20 && r$marks_older_pct < 80
})

check("a small disagreement between exact rules is tolerated and reported", {
  d <- cif_data
  # Break clock > chronological age for 3 people, so it diverges from
  # advance > 0 by 0.5% of the cohort - under the 1% tolerance.
  flip <- which(d$phenoage_advance > 0)[1:3]
  d$phenoage[flip] <- d$AgeAtFirstDate[flip] - 1
  a <- tryCatch(quietly(attach_cif_groups(d, groups_file)),
                error = function(e) e)
  !inherits(a, "error") && !all(is.na(a$phenoage_group))
})
check("a large disagreement between exact rules is still refused", {
  d <- cif_data
  flip <- which(d$phenoage_advance > 0)
  d$phenoage[flip] <- d$AgeAtFirstDate[flip] - 1   # about half the cohort
  err <- tryCatch(quietly(attach_cif_groups(d, groups_file)),
                  error = conditionMessage)
  is.character(err) && grepl("tolerance", err)
})

check("a group file that cannot be linked raises a clear error", {
  bad <- data.frame(unrelated_key = 1:10, kdm_group = "Biologically older",
                    stringsAsFactors = FALSE)
  err <- tryCatch(quietly(attach_cif_groups(cif_data, bad)),
                  error = conditionMessage)
  is.character(err) && grepl("shared|link", err)
})
check("a grouping no rule explains is refused rather than approximated", {
  set.seed(123)
  scrambled <- data.frame(
    LopNr = cif_data$lopnr,
    kdm_group = sample(c("Biologically older", "Biologically younger"),
                       nrow(cif_data), replace = TRUE),
    stringsAsFactors = FALSE)
  err <- tryCatch(quietly(attach_cif_groups(cif_data, scrambled)),
                  error = conditionMessage)
  is.character(err) && grepl("exactly", err)
})

# --- 11g. Risk differences -------------------------------------------

diff_rows <- data.frame(
  group = c("Biologically younger", "Biologically older"),
  state = "dementia", time = 80,
  estimate = c(0.10, 0.16), se = c(0.01, 0.02),
  CI_lower = NA_real_, CI_upper = NA_real_,
  n_risk = c(500, 500), n_event_cum = c(50, 80),
  n_group = 500, events_group = 50, deaths_group = 100,
  stringsAsFactors = FALSE
)
dd <- cif_differences(diff_rows, "age", 80)

check("the risk difference is older minus younger",
      isTRUE(all.equal(dd$difference, 0.06)))
check("the difference CI uses the pooled standard error",
      isTRUE(all.equal(dd$diff_upper - dd$difference,
                       1.96 * sqrt(0.01^2 + 0.02^2))))
check("the difference interval brackets the difference",
      dd$diff_lower < dd$difference && dd$difference < dd$diff_upper)
check("a reporting time with only one group present is skipped",
      is.null(cif_differences(diff_rows[1, ], "age", 80)))

# --- 11h. End to end, with suppression and resume ---------------------

cif_dir <- tempfile("ciftest"); dir.create(cif_dir)
on.exit(unlink(cif_dir, recursive = TRUE), add = TRUE)
curves_file  <- file.path(cif_dir, "curves.csv")
summary_file <- file.path(cif_dir, "summary.csv")

old_plots <- cif_config$make_plots
cif_config$make_plots <- FALSE

quietly(run_cif_block(
  data          = attached,
  specs         = list(list(outcome_label = "ADRD", event_var = "ADRD")),
  dataset_label = "Main",
  curves_file   = curves_file,
  summary_file  = summary_file,
  axes          = "years"
))

curves <- read.table(curves_file, sep = ";", dec = ",", header = TRUE,
                     quote = "", fileEncoding = "Latin1")

check("the curves file is written with exactly the declared columns",
      setequal(names(curves), cif_result_columns),
      paste(setdiff(names(curves), cif_result_columns), collapse = ", "))
check("both groupings are present in the output",
      setequal(unique(curves$grouping), c("kdm_group", "phenoage_group")))
check("thin risk sets are suppressed from the written output",
      all(curves$n_risk >= cif_config$min_risk_set))
check("numbers survive the semicolon / comma-decimal round trip",
      all(curves$estimate >= 0 & curves$estimate <= 1))
check("the summary file carries the risk differences",
      file.exists(summary_file) &&
      nrow(read.table(summary_file, sep = ";", dec = ",", header = TRUE,
                      quote = "", fileEncoding = "Latin1")) > 0)

before_rerun <- readLines(curves_file)
quietly(run_cif_block(
  data          = attached,
  specs         = list(list(outcome_label = "ADRD", event_var = "ADRD")),
  dataset_label = "Main",
  curves_file   = curves_file,
  summary_file  = summary_file,
  axes          = "years"
))
check("a completed curve set is not recomputed or duplicated",
      identical(before_rerun, readLines(curves_file)))

cif_config$make_plots <- old_plots

check("the cif mode writes to its own files, never the Fine-Gray ones",
      !any(grepl("^finegray", list.files(cif_dir))))

# --- 11i. The exit-age rounding check ---------------------------------
# The published Fine-Gray run used time_round_digits = 1. The cif mode inherits
# it so both are computed on the same timeline, and reports whether it mattered.

check("the shipped rounding matches the published Fine-Gray run",
      identical(fg_config$time_round_digits, 1))

rc <- quietly(compare_rounding_effect(
  attached, list(list(outcome_label = "ADRD", event_var = "ADRD")), "Main"))

check("the rounding check reports both settings side by side",
      all(c("N_rounded", "events_rounded", "N_exact", "events_exact",
            "excluded_by_rounding", "reclassified", "counts_identical")
          %in% names(rc)))
check("the rounding check flags a real difference rather than hiding it", {
  # Entry/exit ages here are spread over decades, so rounding to a whole year
  # must collapse something and change the counts.
  coarse <- quietly(compare_rounding_effect(
    attached, list(list(outcome_label = "ADRD", event_var = "ADRD")), "Main",
    rounded = 0))
  !coarse$counts_identical[1]
})
check("comparing a setting against itself reports identical",
      isTRUE(quietly(compare_rounding_effect(
        attached, list(list(outcome_label = "ADRD", event_var = "ADRD")),
        "Main", rounded = 1, exact = 1))$counts_identical[1]))

# Excluding someone and reclassifying someone have completely different
# implications for the Methods, so the check must tell them apart rather than
# reporting "the counts differ" and leaving it there.
check("the person ID is found case-insensitively",
      identical(find_id_column(attached), "sampleID"))

excl <- quietly(compare_rounding_effect(
  attached, list(list(outcome_label = "ADRD", event_var = "ADRD")), "Main",
  rounded = 0))

check("people lost to rounding are counted as excluded, not as reclassified",
      excl$excluded_by_rounding[1] > 0)
check("the excluded count explains the whole change in N",
      excl$N_exact[1] - excl$N_rounded[1] ==
        excl$excluded_by_rounding[1] - excl$gained_by_rounding[1])
check("reclassification is reported as its own number",
      is.finite(excl$reclassified[1]))
check("without an ID column the check says so instead of guessing", {
  no_id <- attached
  no_id$sampleID <- NULL; no_id$lopnr <- NULL
  r <- quietly(compare_rounding_effect(
    no_id, list(list(outcome_label = "ADRD", event_var = "ADRD")), "Main",
    rounded = 0))
  is.na(r$reclassified[1]) && is.na(r$excluded_by_rounding[1])
})

# --- 11j. The figure actually renders ---------------------------------
# run_cif_block() wraps plot_cif() in tryCatch and only logs a failure, so a
# broken plot call produces a truncated PNG and an easily missed log line
# rather than an error. Render one for real instead of trusting that.

plot_file <- file.path(cif_dir, "plot_check.png")
plot_err <- tryCatch({
  plot_cif(big_rows, plot_file, "ADRD - test", "years")
  NA_character_
}, error = conditionMessage)

check("plot_cif runs without error", is.na(plot_err), paste("error:", plot_err))
check("plot_cif writes a complete-looking PNG",
      file.exists(plot_file) && file.info(plot_file)$size > 10000,
      paste("size:", if (file.exists(plot_file)) file.info(plot_file)$size else "missing"))

check("plot_cif copes with a single group present", {
  one <- big_rows[big_rows$group == big_rows$group[1], ]
  f <- file.path(cif_dir, "plot_one_group.png")
  err <- tryCatch({ plot_cif(one, f, "one group", "years"); NA_character_ },
                  error = conditionMessage)
  is.na(err) && file.exists(f)
})
check("a cohort with no competing deaths warns but still returns a curve", {
  no_death <- hand
  no_death$fg_status <- factor(
    c("event", "censor", "event", "censor", "censor"),
    levels = c("censor", "event", "death"))
  r <- quietly(cif_for_outcome(no_death, "grp", "years"))
  !is.null(r) && "dementia" %in% r$state
})
check("an outcome with no events of its own is refused, not silently plotted", {
  no_event <- hand
  no_event$fg_status <- factor(
    c("death", "death", "censor", "censor", "death"),
    levels = c("censor", "event", "death"))
  err <- tryCatch(quietly(cif_for_outcome(no_event, "grp", "years")),
                  error = conditionMessage)
  is.character(err) && grepl("event", err)
})

check("plot_cif copes with an all-missing curve", {
  blanked <- big_rows
  blanked$estimate[blanked$state == "death"] <- NA_real_
  f <- file.path(cif_dir, "plot_na.png")
  err <- tryCatch({ plot_cif(blanked, f, "na curve", "age"); NA_character_ },
                  error = conditionMessage)
  is.na(err) && file.exists(f)
})

# ---------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(checks - failures, "/", checks, " checks passed\n", sep = "")
if (failures > 0) {
  cat(failures, " FAILED\n", sep = "")
  quit(save = "no", status = 1)
}
cat("All good.\n")
