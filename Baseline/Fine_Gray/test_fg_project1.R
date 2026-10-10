#!/usr/bin/env Rscript
# =====================================================================
# Tests for the project 1 Fine-Gray helpers
# =====================================================================
#
# Runs entirely on small synthetic data - no real data, no VDI needed.
# Run before shipping to the VDI, and after any edit to the helpers:
#
#   Rscript test_fg_project1_ADRD.R
#
# Exits non-zero if anything fails.

source("fg_project1_ADRD.R")

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

rounded <- clean_fg_data(test_data, "Alzheimers", "test", time_round_digits = 1)
check("rounding reduces the number of distinct exit times",
      length(unique(rounded$AgeAtLastDate)) < length(unique(cleaned$AgeAtLastDate)))
check("rounding still yields a usable expansion",
      nrow(expand_fg_data(rounded)) > 0)
check("rounding shrinks the expansion",
      nrow(expand_fg_data(rounded)) < nrow(full))

# ---------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(checks - failures, "/", checks, " checks passed\n", sep = "")
if (failures > 0) {
  cat(failures, " FAILED\n", sep = "")
  quit(save = "no", status = 1)
}
cat("All good.\n")
