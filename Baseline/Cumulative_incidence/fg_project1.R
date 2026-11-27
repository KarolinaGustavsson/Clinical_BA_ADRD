#!/usr/bin/env Rscript
# =====================================================================
# Fine-Gray competing risk analysis for project 1 (dementia)
# =====================================================================
#
# One file, four modes. Config, timeline derivation, cleaning, expansion,
# fitting and resume logic all live here; the bottom of the file just picks
# which dataset/diagnostics to run based on the first command-line argument.
#
#   Rscript fg_project1.R main   [outcome ...]   # ADRD Alzheimers Vascular Other
#   Rscript fg_project1.R svedem [outcome ...]   # S_ADRD S_e_alz S_l_alz S_alz_vas
#                                                 # S_vas S_lew S_fro S_par S_uns S_oth
#   Rscript fg_project1.R inspect                # data diagnostics; run this first
#   Rscript fg_project1.R cif    [outcome ...]   # absolute risk: Aalen-Johansen
#                                                 # cumulative incidence by BA group
#
# Put this next to cox_main_2025_06_06.RData and run from a terminal with
# everything else closed - the finegray() expansion is what runs out of
# memory inside RStudio/Rmd. If a run dies, just start it again: every model
# is appended to CSV the moment it finishes, and completed models are
# skipped on restart.
#
# Data expectations (prepared upstream by Calculating_BA_B.Rmd + CoxPH_B.Rmd):
#   to1po2           main analysis dataset
#   to1po2_SveDem_D  SveDem-restricted dataset
# Both must already contain AgeAtFirstDate, AgeAtLastDate, status, ADRD,
# the outcome indicators, the BA measures and the covariates.
#
# Modelling choices, all deliberate - see README_fg_project1.md:
#   * income_bins and NewUtbNiva stay NUMERIC, matching coxph in CoxPH_B.Rmd.
#   * Other dementia subtypes are competing events, alongside death.
#   * One finegray() expansion per outcome serves all 4 models.
#   * The cif mode never calls finegray(), so it cannot move the published
#     subdistribution hazard ratios. It reuses clean_fg_data() only.
#
# This file is also source()'d by test_fg_project1.R to reuse the functions
# below without triggering a real run - see "Command-line entry point" at
# the bottom for how that is guarded.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(survival)
})

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

fg_config <- list(
  # Searched in order in the working directory. coxph_B_correctsubtypes is what
  # the published Fine-Gray run actually loaded (see finegray_project1_main.log,
  # 2026-09-07) - it is the corrected-subtypes version of the prepared data. The
  # cox_main names are kept as a fallback for older copies.
  rdata_candidates = c(
    "coxph_B_correctsubtypes.RData",
    "coxph_B_correctsubtypes.rdata",
    "coxph_B_correctsubtypes.RDATA",
    "coxph_B_correctsubtypes",
    "cox_main_2025_06_06.RData",
    "cox_main_2025_06_06.rdata",
    "cox_main_2025_06_06.RDATA",
    "cox_main_2025_06_06"
  ),

  ba_measures = c("kdm_advance_sd", "phenoage_advance_sd"),

  # gender is left exactly as stored; coxph and finegray coerce character to
  # factor themselves, so this reproduces CoxPH_B.Rmd behaviour either way.
  # income_bins and NewUtbNiva are numeric there, so they stay numeric here.
  factor_covariates  = "gender",
  numeric_covariates = c("income_bins", "NewUtbNiva"),

  # Marks people whose AgeAtLastDate is a dementia diagnosis date rather than
  # a death/censoring date. Only used in legacy mode (see derive_fg_timeline).
  dementia_flag_var = "ADRD",

  # Should a different dementia subtype count as a competing event?
  # TRUE  - correct when the registers only ever record one dementia type, so a
  #         vascular case can never become an Alzheimer case. Also the only
  #         option that works when follow-up stops at diagnosis.
  # FALSE - other subtypes are censored at their diagnosis age instead. This is
  #         informative censoring; keep it only as a sensitivity analysis.
  # Either way ADRD / S_ADRD are unaffected: there is no "other subtype" there.
  competing_other_dementia = TRUE,

  # Column names in the prepared data, used to build the four age columns.
  # Vectors are tried in order and the first one present is used.
  timeline_columns = list(
    birth        = "FODDATn",
    death        = "DODSDATn",
    end          = c("lastDate", "d2"),
    entry_age    = "AgeAtFirstDate",
    cox_exit_age = "AgeAtLastDate",
    dx           = c("INDATUM", "S_INDATUM"),
    meds_date    = "FDATUM",
    atc          = "ATC",
    icd          = "Dementia",
    svedem       = "S_ADRD"
  ),

  # Rounding the exit age collapses tied event times and shrinks the finegray
  # expansion, which is where the memory goes. NULL = no rounding (exact).
  # 3 ~ one day, 2 ~ four days, 1 ~ five weeks.
  #
  # 1 is NOT the default for convenience - it is what the published run used
  # (finegray_project1_main.log, 2026-09-07), because the exact expansion would
  # not fit in memory. The cif mode inherits it so the cumulative incidence is
  # computed on exactly the timeline behind Supplementary Table 9; changing it
  # here would silently desynchronise the two. The cif mode reports the event
  # counts under both settings so the effect of the rounding is visible rather
  # than assumed - see compare_rounding_effect().
  time_round_digits = 1
)

adjustment_sets <- list(
  "Sex-adjusted model"  = "gender",
  "Multivariable model" = c("gender", "income_bins", "NewUtbNiva")
)

# All covariates any model needs. Used for the shared complete-case sample.
all_model_covariates <- unique(unlist(adjustment_sets, use.names = FALSE))

main_specs <- list(
  list(outcome_label = "ADRD",       event_var = "ADRD"),
  list(outcome_label = "Alzheimers", event_var = "Alzheimers"),
  list(outcome_label = "Vascular",   event_var = "Vascular"),
  list(outcome_label = "Other",      event_var = "Other")
)

svedem_specs <- list(
  list(outcome_label = "S_ADRD",    event_var = "S_ADRD"),
  list(outcome_label = "S_e_alz",   event_var = "S_e_alz"),
  list(outcome_label = "S_l_alz",   event_var = "S_l_alz"),
  list(outcome_label = "S_alz_vas", event_var = "S_alz_vas"),
  list(outcome_label = "S_vas",     event_var = "S_vas"),
  list(outcome_label = "S_lew",     event_var = "S_lew"),
  list(outcome_label = "S_fro",     event_var = "S_fro"),
  list(outcome_label = "S_par",     event_var = "S_par"),
  list(outcome_label = "S_uns",     event_var = "S_uns"),
  list(outcome_label = "S_oth",     event_var = "S_oth")
)

# ---------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------

fg_log_file <- NULL

log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n")
  cat(line)
  flush.console()
  if (!is.null(fg_log_file)) {
    cat(line, file = fg_log_file, append = TRUE)
  }
  invisible(NULL)
}

print_gc <- function(label = "gc") {
  used <- sum(gc()[, "used"] * c(8, 8) / 1e6)
  log_msg("Memory (", label, "): ", round(used, 1), " Mb in use")
  invisible(NULL)
}

elapsed_min <- function(start_time) {
  as.numeric(difftime(Sys.time(), start_time, units = "mins"))
}

#' Install an error handler that logs a traceback before dying.
#' Deliberately does NOT save.image() - the workspace here is several Gb.
install_crash_handler <- function(label) {
  options(error = function() {
    cat("\n==============================\n")
    cat("ERROR at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " in ", label, "\n", sep = "")
    traceback(2)
    cat("Results completed before the crash are already in the CSV.\n")
    cat("Re-running the script resumes from there.\n")
    cat("==============================\n")
    quit(save = "no", status = 1)
  })
}

# ---------------------------------------------------------------------
# Loading the prepared data
# ---------------------------------------------------------------------

find_rdata_file <- function(candidates = fg_config$rdata_candidates) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "Could not find the prepared RData file in: ", getwd(),
      "\nTried: ", paste(candidates, collapse = ", ")
    )
  }
  existing[1]
}

#' Load the prepared workspace and return only the two datasets we need.
#'
#' Everything else from the Cox workspace is dropped immediately, which is the
#' single biggest memory saving available before the expansion starts.
load_prepared_data <- function(required = c("to1po2", "to1po2_SveDem_D")) {
  file <- find_rdata_file()
  log_msg("Loading ", file)
  start <- Sys.time()

  env <- new.env(parent = emptyenv())
  loaded <- load(file, envir = env)
  log_msg("Loaded ", length(loaded), " objects in ", round(elapsed_min(start), 2), " min")

  missing_objs <- setdiff(required, loaded)
  if (length(missing_objs) > 0) {
    stop(
      "Missing objects in ", file, ": ", paste(missing_objs, collapse = ", "),
      "\nFound: ", paste(loaded, collapse = ", ")
    )
  }

  out <- mget(required, envir = env)
  rm(env)
  gc()
  for (nm in required) {
    log_msg("  ", nm, ": ", nrow(out[[nm]]), " rows x ", ncol(out[[nm]]), " cols")
  }
  print_gc("after load")
  out
}

# ---------------------------------------------------------------------
# Result storage with resume support
# ---------------------------------------------------------------------

result_columns <- c(
  "dataset", "outcome", "model", "BA",
  "beta", "se", "HR", "CI_lower", "CI_upper", "p",
  "N", "events", "competing_deaths", "censored",
  "fg_expanded_N", "expansion_elapsed_min", "model_elapsed_min",
  "robust", "error"
)

result_key <- function(dataset, outcome, model, BA) {
  paste(dataset, outcome, model, BA, sep = "|")
}

#' Keys already present in the results file, so a restart can skip them.
completed_keys <- function(results_file) {
  if (!file.exists(results_file)) return(character(0))
  done <- tryCatch(
    read.table(results_file, sep = ";", dec = ",", header = TRUE,
               quote = "", fileEncoding = "Latin1", colClasses = "character"),
    error = function(e) {
      warning("Could not read ", results_file, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(done) || nrow(done) == 0) return(character(0))
  result_key(done$dataset, done$outcome, done$model, done$BA)
}

#' Append one result row, writing the header only when creating the file.
append_result <- function(row, results_file) {
  row <- as.data.frame(row, stringsAsFactors = FALSE)[, result_columns, drop = FALSE]
  new_file <- !file.exists(results_file)
  write.table(
    row, results_file,
    sep = ";", dec = ",",
    row.names = FALSE,
    col.names = new_file,
    append = !new_file,
    quote = FALSE,
    fileEncoding = "Latin1"
  )
  invisible(row)
}

make_result_row <- function(dataset, outcome, model, BA,
                            beta = NA_real_, se = NA_real_,
                            N = NA_integer_, events = NA_integer_,
                            competing_deaths = NA_integer_, censored = NA_integer_,
                            fg_expanded_N = NA_integer_,
                            expansion_elapsed_min = NA_real_,
                            model_elapsed_min = NA_real_,
                            robust = NA, error = NA_character_) {
  if (is.na(beta) || is.na(se)) {
    HR <- CI_lower <- CI_upper <- p <- NA_real_
  } else {
    HR       <- exp(beta)
    CI_lower <- exp(beta - 1.96 * se)
    CI_upper <- exp(beta + 1.96 * se)
    p        <- 2 * pnorm(-abs(beta / se))
  }

  data.frame(
    dataset = dataset, outcome = outcome, model = model, BA = BA,
    beta = beta, se = se, HR = HR, CI_lower = CI_lower, CI_upper = CI_upper, p = p,
    N = N, events = events, competing_deaths = competing_deaths, censored = censored,
    fg_expanded_N = fg_expanded_N,
    expansion_elapsed_min = expansion_elapsed_min,
    model_elapsed_min = model_elapsed_min,
    robust = robust, error = error,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------
# Data cleaning
# ---------------------------------------------------------------------

require_columns <- function(data, columns, context) {
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", context, ": ", paste(missing_cols, collapse = ", "))
  }
}

# ---------------------------------------------------------------------
# Timeline: the four age columns Fine-Gray needs
# ---------------------------------------------------------------------

timeline_columns <- c("age_entry", "age_dementia", "age_death", "age_censor")

first_present <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) NULL else hit[1]
}

age_between <- function(date, birth) {
  as.numeric(difftime(as.Date(date), as.Date(birth), units = "days")) / 365.25
}

#' Build the four age columns Fine-Gray needs from the prepared Cox data.
#'
#' Cox only ever needed one exit column, so CoxPH_B.Rmd overwrites
#' AgeAtLastDate with the age at dementia diagnosis for anyone who has one.
#' That is fine for Cox - everyone who is not a case is simply censored - but
#' it destroys the information Fine-Gray needs, namely *when* the competing
#' death happened. This splits the timeline back apart:
#'
#'   age_entry     age at baseline                        (AgeAtFirstDate)
#'   age_dementia  age at dementia diagnosis, else NA     (AgeAtLastDate, reused)
#'   age_death     age at DODSDATn, else NA               (recovered here)
#'   age_censor    age at end of follow-up                (lastDate / d2)
#'
#' age_dementia deliberately REUSES AgeAtLastDate rather than recomputing it,
#' so the dementia timing is identical to the Cox models by construction. The
#' ATC/FDATUM rule for the meds-only group is therefore inherited automatically
#' and cannot drift out of sync.
#'
#' The one thing that must be recomputed is which rows that reuse is valid for.
#' CoxPH_B.Rmd only overwrites AgeAtLastDate when INDATUM/S_INDATUM exists, or
#' FDATUM for the meds-only group. An ADRD case with none of those keeps the
#' lastDate age, which is not a diagnosis age at all - so the date source is
#' recomputed here with the same precedence rather than trusting ADRD == 1.
derive_fg_timeline <- function(data, cols = fg_config$timeline_columns,
                               quiet = FALSE) {

  if (all(timeline_columns %in% names(data))) {
    return(data[, timeline_columns, drop = FALSE])
  }

  need_col <- function(key) {
    hit <- first_present(data, cols[[key]])
    if (is.null(hit)) {
      stop("derive_fg_timeline: none of these columns found for '", key, "': ",
           paste(cols[[key]], collapse = ", "))
    }
    hit
  }

  entry_col <- need_col("entry_age")
  cox_col   <- need_col("cox_exit_age")
  birth_col <- need_col("birth")
  death_col <- need_col("death")
  end_col   <- need_col("end")

  n <- nrow(data)

  # Same precedence as CoxPH_B.Rmd section 1.4: FDATUM for the meds-only group,
  # then INDATUM/S_INDATUM overrides it where present.
  dx_cols <- intersect(cols$dx, names(data))
  if (length(dx_cols) == 0) {
    stop("derive_fg_timeline: no dementia date column found (looked for ",
         paste(cols$dx, collapse = ", "), ")")
  }
  dx_dates <- lapply(dx_cols, function(cn) as.Date(data[[cn]]))
  earliest_dx <- Reduce(function(a, b) pmin(a, b, na.rm = TRUE), dx_dates)

  meds_col   <- first_present(data, cols$meds_date)
  atc_col    <- first_present(data, cols$atc)
  icd_col    <- first_present(data, cols$icd)
  svedem_col <- first_present(data, cols$svedem)

  meds_only <- rep(FALSE, n)
  if (!is.null(atc_col) && !is.null(icd_col) && !is.null(svedem_col)) {
    meds_only <- data[[atc_col]] != "" & !is.na(data[[atc_col]]) &
      (data[[icd_col]] == "" | is.na(data[[icd_col]])) &
      !is.na(data[[svedem_col]]) & data[[svedem_col]] == 0
  }
  has_meds_date <- if (is.null(meds_col)) rep(FALSE, n) else !is.na(data[[meds_col]])

  has_dementia_date <- !is.na(earliest_dx) | (meds_only & has_meds_date)

  out <- data.frame(
    age_entry    = as.numeric(data[[entry_col]]),
    # Reuse the Cox exit age wherever it really is a dementia date.
    age_dementia = ifelse(has_dementia_date, as.numeric(data[[cox_col]]), NA_real_),
    age_death    = age_between(data[[death_col]], data[[birth_col]]),
    age_censor   = age_between(data[[end_col]],   data[[birth_col]]),
    stringsAsFactors = FALSE
  )

  if (!quiet) {
    log_msg("Timeline derived from ", entry_col, " / ", cox_col, " / ",
            death_col, " / ", end_col)
    log_msg("  with a dementia date: ", sum(has_dementia_date),
            "  |  with a death date: ", sum(!is.na(out$age_death)),
            "  |  censor age missing: ", sum(is.na(out$age_censor)))
    flag_col <- first_present(data, fg_config$dementia_flag_var)
    if (!is.null(flag_col)) {
      orphan <- !is.na(data[[flag_col]]) & data[[flag_col]] == 1 & !has_dementia_date
      if (any(orphan)) {
        log_msg("  WARNING: ", sum(orphan), " rows have ", flag_col,
                " == 1 but no diagnosis date; they cannot be timed and are",
                " treated as non-cases. Check the 'inspect' mode section 1.")
      }
    }
  }

  out
}

#' Build the Fine-Gray input for one outcome.
#'
#' Exit time and status come from the four timeline ages on a first-event-wins
#' basis: whichever of dementia diagnosis, death or end of follow-up happens
#' first is what the person exits on. Ties go to dementia, then death.
#'
#' Event coding:
#'   1 = the dementia outcome of interest
#'   2 = competing event (death, and other dementia subtypes by default)
#'   0 = censored
#'
#' A different dementia subtype is a competing event because the registers
#' record only one dementia type per person - a vascular case can never become
#' an Alzheimer case, so leaving them in the risk set would pad it with people
#' whose hazard is zero. Set fg_config$competing_other_dementia = FALSE to
#' censor them at diagnosis instead (informative censoring; sensitivity only).
#' ADRD / S_ADRD are unaffected either way - there is no "other subtype" there.
#'
#' extra_vars carries additional columns through to the output without letting
#' them influence the sample. They are NOT part of the complete-case filter, so
#' a missing value in one of them can never silently shrink N. Used by the cif
#' mode to carry the biological-age group labels. Default character(0) means
#' every Fine-Gray call site behaves exactly as it did before this argument
#' existed - asserted by a snapshot check in test_fg_project1.R section 11.
clean_fg_data <- function(data, event_var, dataset_label,
                          ba_measures = fg_config$ba_measures,
                          covariates = all_model_covariates,
                          competing_other_dementia = fg_config$competing_other_dementia,
                          time_round_digits = fg_config$time_round_digits,
                          quiet = TRUE,
                          extra_vars = character(0)) {

  require_columns(data, unique(c(event_var, ba_measures, covariates, extra_vars)),
                  context = paste0(dataset_label, " / ", event_var))

  timeline <- derive_fg_timeline(data, quiet = quiet)

  df <- data.frame(
    timeline,
    fg_outcome = data[[event_var]],
    data[, unique(c(ba_measures, covariates, extra_vars)), drop = FALSE],
    stringsAsFactors = FALSE
  )

  if (!is.null(time_round_digits)) {
    for (v in timeline_columns) df[[v]] <- round(df[[v]], time_round_digits)
  }

  # First event wins. na.rm handles "never diagnosed" and "still alive".
  df$exit_age <- pmin(df$age_dementia, df$age_death, df$age_censor, na.rm = TRUE)

  # Ties: dementia first, then death. eq() treats NA as "did not happen".
  eq <- function(x) !is.na(x) & x == df$exit_age
  first_dementia <- eq(df$age_dementia)
  first_death    <- !first_dementia & eq(df$age_death)

  outcome_yes <- !is.na(df$fg_outcome) & df$fg_outcome == 1

  df$fg_event_num <- 0L
  df$fg_event_num[first_dementia & outcome_yes] <- 1L
  df$fg_event_num[first_death] <- 2L
  if (competing_other_dementia) {
    df$fg_event_num[first_dementia & !outcome_yes] <- 2L
  }

  df$fg_status <- factor(
    df$fg_event_num,
    levels = c(0, 1, 2),
    labels = c("censor", "event", "death")
  )

  # Compact integer id for delayed-entry expansion and robust clustering.
  df$fg_id <- seq_len(nrow(df))

  # Match coxph coercion: character/factor becomes a factor, numeric stays numeric.
  for (cv in intersect(fg_config$factor_covariates, names(df))) {
    if (!is.numeric(df[[cv]])) df[[cv]] <- droplevels(factor(df[[cv]]))
  }
  for (cv in intersect(fg_config$numeric_covariates, names(df))) {
    df[[cv]] <- as.numeric(df[[cv]])
  }

  # extra_vars is deliberately absent here - see the note in the header.
  complete_vars <- unique(c("age_entry", "exit_age", ba_measures, covariates))
  keep <- complete.cases(df[, complete_vars, drop = FALSE]) &
    !is.na(df$fg_status) &
    is.finite(df$exit_age) &
    df$exit_age > df$age_entry

  df <- df[keep, , drop = FALSE]
  df <- droplevels(df)
  # Renumber so ids are 1..N again after dropping rows.
  df$fg_id <- seq_len(nrow(df))

  # Keep the names the rest of the pipeline expects.
  names(df)[names(df) == "age_entry"] <- "AgeAtFirstDate"
  names(df)[names(df) == "exit_age"]  <- "AgeAtLastDate"

  df[, unique(c("AgeAtFirstDate", "AgeAtLastDate", "fg_status", "fg_id",
                ba_measures, covariates, extra_vars)), drop = FALSE]
}

fg_counts <- function(df) {
  list(
    N                = nrow(df),
    events           = sum(df$fg_status == "event"),
    competing_deaths = sum(df$fg_status == "death"),
    censored         = sum(df$fg_status == "censor")
  )
}

# ---------------------------------------------------------------------
# Fine-Gray expansion
# ---------------------------------------------------------------------

#' Expand once per outcome and reuse for every model.
#'
#' finegray() derives its weights from the censoring distribution of
#' (entry, exit, status, id) only - the right-hand side variables are just
#' carried along into the output. So carrying every BA measure and every
#' covariate at once gives an expansion valid for all four models, instead of
#' repeating the expensive expansion per model. verify_expansion_invariance()
#' in the Rmd checks this rather than taking it on faith.
expand_fg_data <- function(df,
                           ba_measures = fg_config$ba_measures,
                           covariates = all_model_covariates) {
  rhs <- unique(c(ba_measures, covariates, "fg_id"))
  fg_formula <- as.formula(paste0(
    "Surv(AgeAtFirstDate, AgeAtLastDate, fg_status) ~ ",
    paste(rhs, collapse = " + ")
  ))

  fg <- survival::finegray(fg_formula, data = df, etype = "event", id = fg_id)

  # finegray keeps extra bookkeeping columns; coxph only needs these.
  fg[, c("fgstart", "fgstop", "fgstatus", "fgwt", "fg_id", ba_measures, covariates),
     drop = FALSE]
}

# ---------------------------------------------------------------------
# Weighted Cox model on the expanded data
# ---------------------------------------------------------------------

fit_fg_model <- function(fg_expanded, ba_var, covariates, use_robust = TRUE) {
  terms <- c(ba_var, covariates)
  if (use_robust) terms <- c(terms, "cluster(fg_id)")

  fit_formula <- as.formula(
    paste0("Surv(fgstart, fgstop, fgstatus) ~ ", paste(terms, collapse = " + "))
  )

  survival::coxph(
    fit_formula,
    weights = fgwt,
    data = fg_expanded,
    robust = use_robust,
    x = FALSE, y = FALSE, model = FALSE,
    singular.ok = TRUE
  )
}

#' Fit one model, falling back to the non-robust variant if the robust
#' sandwich estimator runs out of memory. Returns beta/se, never the fit.
fit_with_fallback <- function(fg_expanded, ba_var, covariates, label) {
  start <- Sys.time()

  fit <- tryCatch(
    fit_fg_model(fg_expanded, ba_var, covariates, use_robust = TRUE),
    error = function(e) e
  )
  used_robust <- TRUE

  if (inherits(fit, "error")) {
    log_msg("  robust fit failed (", conditionMessage(fit), "); trying non-robust")
    gc()
    robust_msg <- conditionMessage(fit)
    fit <- tryCatch(
      fit_fg_model(fg_expanded, ba_var, covariates, use_robust = FALSE),
      error = function(e) e
    )
    used_robust <- FALSE
    if (inherits(fit, "error")) {
      return(list(
        beta = NA_real_, se = NA_real_, robust = NA,
        elapsed = elapsed_min(start),
        error = paste0("robust: ", robust_msg, " | non-robust: ", conditionMessage(fit))
      ))
    }
  }

  coefs <- coef(fit)
  if (!(ba_var %in% names(coefs)) || is.na(coefs[[ba_var]])) {
    return(list(
      beta = NA_real_, se = NA_real_, robust = used_robust,
      elapsed = elapsed_min(start),
      error = paste0("No estimable coefficient for ", ba_var, " in ", label)
    ))
  }

  # vcov() is the robust sandwich when robust = TRUE, the model-based
  # variance otherwise, so this picks up the right se either way.
  se <- sqrt(diag(vcov(fit)))[[ba_var]]

  list(
    beta = unname(coefs[[ba_var]]), se = se, robust = used_robust,
    elapsed = elapsed_min(start), error = NA_character_
  )
}

# ---------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------

#' Run every model for one outcome, appending each result as it completes.
#'
#' Skips work already present in results_file, so a crashed run resumes by
#' simply being started again.
run_outcome <- function(data, spec, dataset_label, results_file,
                        ba_measures = fg_config$ba_measures,
                        sets = adjustment_sets,
                        done_keys = character(0)) {

  wanted <- expand.grid(
    model = names(sets), BA = ba_measures,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  wanted$key <- result_key(dataset_label, spec$outcome_label, wanted$model, wanted$BA)
  todo <- wanted[!(wanted$key %in% done_keys), , drop = FALSE]

  if (nrow(todo) == 0) {
    log_msg("Skipping ", dataset_label, " / ", spec$outcome_label, " (already complete)")
    return(invisible(character(0)))
  }

  log_msg("=== ", dataset_label, " / ", spec$outcome_label,
          " (", nrow(todo), " of ", nrow(wanted), " models to run) ===")

  df <- clean_fg_data(data, spec$event_var, dataset_label)
  counts <- fg_counts(df)
  log_msg("N=", counts$N, "  events=", counts$events,
          "  competing deaths=", counts$competing_deaths,
          "  censored=", counts$censored)

  fail_all <- function(msg) {
    for (i in seq_len(nrow(todo))) {
      append_result(make_result_row(
        dataset = dataset_label, outcome = spec$outcome_label,
        model = todo$model[i], BA = todo$BA[i],
        N = counts$N, events = counts$events,
        competing_deaths = counts$competing_deaths, censored = counts$censored,
        error = msg
      ), results_file)
    }
    log_msg("  ", msg)
    invisible(todo$key)
  }

  if (counts$N == 0 || counts$events == 0) {
    return(fail_all("No usable observations or no events after cleaning"))
  }

  expansion_start <- Sys.time()
  fg_expanded <- tryCatch(expand_fg_data(df), error = function(e) e)
  rm(df); gc()

  if (inherits(fg_expanded, "error")) {
    return(fail_all(paste("Fine-Gray expansion failed:", conditionMessage(fg_expanded))))
  }

  expansion_elapsed <- elapsed_min(expansion_start)
  expanded_n <- nrow(fg_expanded)
  log_msg("  expanded to ", expanded_n, " rows in ",
          round(expansion_elapsed, 2), " min")
  print_gc("after expansion")

  for (i in seq_len(nrow(todo))) {
    model_name <- todo$model[i]
    ba_var     <- todo$BA[i]
    label      <- paste(dataset_label, spec$outcome_label, model_name, ba_var, sep = " / ")
    log_msg("  fitting ", label)

    fit <- fit_with_fallback(fg_expanded, ba_var, sets[[model_name]], label)

    append_result(make_result_row(
      dataset = dataset_label, outcome = spec$outcome_label,
      model = model_name, BA = ba_var,
      beta = fit$beta, se = fit$se,
      N = counts$N, events = counts$events,
      competing_deaths = counts$competing_deaths, censored = counts$censored,
      fg_expanded_N = expanded_n,
      expansion_elapsed_min = expansion_elapsed,
      model_elapsed_min = fit$elapsed,
      robust = fit$robust, error = fit$error
    ), results_file)

    log_msg("    HR=", if (is.na(fit$beta)) "NA" else round(exp(fit$beta), 4),
            "  (", round(fit$elapsed, 2), " min)")
    gc()
  }

  rm(fg_expanded); gc()
  print_gc("after outcome")
  invisible(todo$key)
}

#' Run a list of outcomes against one dataset, resuming from results_file.
run_outcome_block <- function(data, specs, dataset_label, results_file,
                              only_outcomes = NULL) {
  if (!is.null(only_outcomes)) {
    specs <- Filter(function(s) s$outcome_label %in% only_outcomes, specs)
    if (length(specs) == 0) {
      stop("No outcomes matched: ", paste(only_outcomes, collapse = ", "))
    }
  }

  done <- completed_keys(results_file)
  if (length(done) > 0) log_msg("Resuming: ", length(done), " results already in ", results_file)

  for (spec in specs) {
    done <- c(done, run_outcome(data, spec, dataset_label, results_file, done_keys = done))
  }
  invisible(done)
}

# ---------------------------------------------------------------------
# Data diagnostics (formerly inspect_data_for_fg.R)
# ---------------------------------------------------------------------
#
# Prints only counts, dates and summaries - no individual-level data - so the
# output is safe to copy out of the VDI and paste into a discussion.
#
# Answers the questions the prep code has to make assumptions about:
#   1. Does AgeAtLastDate hold a dementia date for every dementia case?
#   2. Is d2 / lastDate the death date for the dead, or something else?
#   3. For dementia cases, is d2 the diagnosis date? (If so, we cannot follow
#      them past diagnosis, which decides how subtypes must be handled.)
#   4. Which dementia date source applies to whom (INDATUM / S_INDATUM / FDATUM)?
#   5. Are there deaths recorded before baseline, or diagnoses before baseline?
run_inspection <- function(d) {
  hdr <- function(x) cat("\n\n========== ", x, " ==========\n", sep = "")

  hdr("0. Which columns exist")
  wanted <- c("FODDATn", "DODSDATn", "Provtagningsdatum", "d1", "d2",
              "firstDate", "lastDate", "INDATUM", "S_INDATUM",
              "FDATUM", "ATC", "Dementia", "S_ADRD", "ADRD",
              "Alzheimers", "Vascular", "Other",
              "AgeAtFirstDate", "AgeAtLastDate", "status")
  present <- wanted[wanted %in% names(d)]
  cat("Present :", paste(present, collapse = ", "), "\n")
  cat("Missing :", paste(setdiff(wanted, names(d)), collapse = ", "), "\n")
  cat("N rows  :", nrow(d), "\n")

  has <- function(x) x %in% names(d)
  as_date <- function(x) as.Date(x)

  hdr("1. Does AgeAtLastDate hold a dementia date for every dementia case?")
  # CoxPH_B.Rmd only overwrites AgeAtLastDate when INDATUM/S_INDATUM (or FDATUM
  # for the meds-only group) exists. An ADRD case with none of those keeps the
  # lastDate age, which is NOT a dementia date - it would be silently mis-timed.
  if (has("ADRD") && has("INDATUM")) {
    dx_date <- if (has("S_INDATUM")) pmin(as_date(d$INDATUM), as_date(d$S_INDATUM), na.rm = TRUE)
               else as_date(d$INDATUM)
    meds_only <- if (has("ATC") && has("Dementia") && has("S_ADRD")) {
      d$ATC != "" & !is.na(d$ATC) & (d$Dementia == "" | is.na(d$Dementia)) & d$S_ADRD == 0
    } else rep(FALSE, nrow(d))
    fd <- if (has("FDATUM")) as_date(d$FDATUM) else as.Date(rep(NA, nrow(d)))

    # Same precedence as CoxPH_B: FDATUM for meds-only, then dx date overrides.
    has_dem_date <- !is.na(dx_date) | (meds_only & !is.na(fd))
    dem <- d$ADRD == 1 & !is.na(d$ADRD)

    cat("ADRD cases                          :", sum(dem), "\n")
    cat("  with a usable dementia date       :", sum(dem & has_dem_date), "\n")
    cat("  WITHOUT any dementia date         :", sum(dem & !has_dem_date), "\n")
    cat("    -> these have AgeAtLastDate = lastDate age, not a diagnosis age.\n")
    cat("    -> they must be dropped or given an explicit rule.\n")
    cat("Non-ADRD rows that do have a date   :", sum(!dem & has_dem_date), "\n")
    cat("\nDate source used (among those with a date):\n")
    print(table(source = ifelse(!is.na(dx_date), "INDATUM/S_INDATUM",
                         ifelse(meds_only & !is.na(fd), "FDATUM (meds-only)", "none")))[
          c("INDATUM/S_INDATUM", "FDATUM (meds-only)")])
  } else {
    cat("Cannot check: need ADRD and INDATUM.\n")
  }

  hdr("2. What is d2 / lastDate for the dead?")
  end_col <- if (has("lastDate")) "lastDate" else if (has("d2")) "d2" else NA
  if (!is.na(end_col) && has("DODSDATn")) {
    dead <- !is.na(d$DODSDATn)
    diff_days <- as.numeric(as_date(d[[end_col]][dead]) - as_date(d$DODSDATn[dead]))
    cat("Using end-of-followup column:", end_col, "\n")
    cat("Dead individuals:", sum(dead), "\n")
    cat("Of those,", end_col, "equals DODSDATn exactly:",
        sum(diff_days == 0, na.rm = TRUE), "\n")
    cat("Differs (days), summary:\n")
    print(summary(diff_days[diff_days != 0]))
    cat("\nIf these are all 0, ", end_col, " is the death date for the dead\n",
        "and age_death can come from either source.\n", sep = "")
  } else {
    cat("Cannot check: need DODSDATn and lastDate/d2.\n")
  }

  hdr("3. Is d2 the diagnosis date for dementia cases?")
  if (!is.na(end_col) && has("INDATUM") && has("ADRD")) {
    dem <- d$ADRD == 1 & !is.na(d$ADRD)
    dx <- if (has("S_INDATUM")) pmin(as_date(d$INDATUM), as_date(d$S_INDATUM), na.rm = TRUE)
          else as_date(d$INDATUM)
    cmp <- as.numeric(as_date(d[[end_col]][dem]) - dx[dem])
    cat("Dementia cases:", sum(dem), "\n")
    cat(end_col, "equals the diagnosis date:", sum(cmp == 0, na.rm = TRUE), "\n")
    cat("Summary of", end_col, "minus diagnosis date, in days:\n")
    print(summary(cmp))
    cat("\nAll zeros means follow-up STOPS at diagnosis, so people with one\n",
        "dementia subtype cannot be followed onward for another subtype.\n",
        "That settles it: other subtypes must be competing events.\n", sep = "")

    cat("\nOf dementia cases, how many also have a death date:",
        sum(dem & !is.na(d$DODSDATn)), "\n")
    if (has("DODSDATn")) {
      after <- as.numeric(as_date(d$DODSDATn[dem]) - dx[dem])
      cat("Days from diagnosis to death (dementia cases who died):\n")
      print(summary(after))
      cat("Negative values = death recorded BEFORE diagnosis:",
          sum(after < 0, na.rm = TRUE), "\n")
    }
  }

  hdr("4. Which dementia date source applies")
  if (has("INDATUM")) {
    src <- data.frame(
      has_INDATUM   = !is.na(d$INDATUM),
      has_S_INDATUM = if (has("S_INDATUM")) !is.na(d$S_INDATUM) else FALSE,
      has_FDATUM    = if (has("FDATUM")) !is.na(d$FDATUM) else FALSE
    )
    print(colSums(src))
    cat("\nCombinations (INDATUM / S_INDATUM / FDATUM):\n")
    print(table(paste(src$has_INDATUM, src$has_S_INDATUM, src$has_FDATUM)))

    if (has("ADRD")) {
      no_date <- d$ADRD == 1 & !src$has_INDATUM & !src$has_S_INDATUM & !src$has_FDATUM
      cat("\nADRD cases with NO usable diagnosis date:", sum(no_date, na.rm = TRUE),
          "\n(these cannot be given an event time and will be dropped)\n")
    }

    if (has("ATC") && has("Dementia") && has("S_ADRD")) {
      meds_only <- d$ATC != "" & !is.na(d$ATC) &
                   (d$Dementia == "" | is.na(d$Dementia)) & d$S_ADRD == 0
      cat("\n'Meds-only' group (ATC but no ICD and not in SveDem):",
          sum(meds_only, na.rm = TRUE), "\n")
      cat("Of those, how many have FDATUM:",
          if (has("FDATUM")) sum(meds_only & !is.na(d$FDATUM), na.rm = TRUE) else NA, "\n")
    }
  }

  hdr("5. Timing sanity around baseline")
  start_col <- if (has("firstDate")) "firstDate" else if (has("d1")) "d1" else
               if (has("Provtagningsdatum")) "Provtagningsdatum" else NA
  if (!is.na(start_col)) {
    cat("Using baseline column:", start_col, "\n")
    b <- as_date(d[[start_col]])
    cat("Baseline date range:", format(range(b, na.rm = TRUE)), "\n")
    if (has("DODSDATn")) {
      cat("Deaths recorded BEFORE baseline:",
          sum(as_date(d$DODSDATn) < b, na.rm = TRUE), "\n")
    }
    if (has("INDATUM")) {
      cat("Dementia diagnoses BEFORE baseline (prevalent, should be excluded):",
          sum(as_date(d$INDATUM) < b, na.rm = TRUE), "\n")
    }
    if (!is.na(end_col)) {
      cat("End-of-followup before baseline:",
          sum(as_date(d[[end_col]]) < b, na.rm = TRUE), "\n")
      cat("End-of-followup range:", format(range(as_date(d[[end_col]]), na.rm = TRUE)), "\n")
    }
  }

  hdr("6. Existing age columns, for comparison")
  for (v in c("AgeAtFirstDate", "AgeAtLastDate")) {
    if (has(v)) { cat(v, ":\n"); print(summary(d[[v]])) }
  }
  if (has("status")) { cat("\nstatus (death indicator):\n"); print(table(d$status, useNA = "ifany")) }

  hdr("7. Outcome counts")
  for (v in c("ADRD", "Alzheimers", "Vascular", "Other", "S_ADRD")) {
    if (has(v)) cat(sprintf("%-12s %d\n", v, sum(d[[v]] == 1, na.rm = TRUE)))
  }
  if (has("Alzheimers") && has("Vascular") && has("Other")) {
    cat("\nPeople counted in more than one subtype (should be 0):",
        sum((d$Alzheimers == 1) + (d$Vascular == 1) + (d$Other == 1) > 1, na.rm = TRUE), "\n")
  }

  cat("\n\nDone. This output contains no individual-level data.\n")
  invisible(NULL)
}

# =====================================================================
# Cumulative incidence functions (the 'cif' mode)
# =====================================================================
#
# Answers reviewer comment 4, which asks for ABSOLUTE risk. Supplementary
# Table 9 reports subdistribution hazard ratios, which are still relative.
#
# This is a different estimator, not a different Fine-Gray run:
# Aalen-Johansen cumulative incidence, computed with survfit() on exactly the
# output of clean_fg_data(). So the sample, the first-event-wins coding and the
# "other subtypes are competing events" rule are shared with the Fine-Gray
# models by construction - but finegray() is never called, so there is no
# expansion, no memory wall, and no way for the published sHRs to move.
#
# Why survfit() and not cmprsk::cuminc(): cuminc() takes only (ftime, fstatus)
# and has no delayed-entry interface, so it cannot do the attained-age timescale
# that the Cox and Fine-Gray models use. survfit() accepts
# Surv(start, stop, status_factor) and needs no package beyond survival.
#
# Both event types are reported. The dementia CIF alone does not answer the
# comment - the point is that the biologically older group is dying faster, so
# the mortality CIF is half the answer.

cif_config <- list(
  # Table 1's biologically-younger / biologically-older split, per ID.
  # Local filename first so the file can simply be dropped next to the script.
  group_file_candidates = c(
    # Next to the script is the normal case: copy the file in alongside it.
    "table1_kdmgroups_baseline_perID.csv",
    # Otherwise look where Table1_B.Rmd wrote it, i.e. under the output
    # directory named by STUDY1_OUTPUT_DIR (see config.R). Dropped from the
    # vector when the variable is unset.
    if (nzchar(Sys.getenv("STUDY1_OUTPUT_DIR"))) {
      file.path(Sys.getenv("STUDY1_OUTPUT_DIR"),
                "table1_kdmgroups_baseline_perID.csv")
    }
  ),

  # Tried case-insensitively against both the group file and the cohort data.
  # sampleID is what table1_kdmgroups_baseline_perID.csv actually uses.
  id_candidates = c("sampleid", "lopnr", "id", "personid", "pid", "serial"),

  # When several rules reproduce the KDM groups exactly but disagree once
  # transferred to another clock, prefer the one built on the same variable the
  # published Cox and Fine-Gray models use as the exposure, i.e. *_advance_sd.
  prefer_rule_order = c("advance_sd > median", "advance > median",
                        "advance_sd > 0", "advance > 0",
                        "clock > chronological age"),

  # How much disagreement between equally-exact rules is tolerable before the
  # transferred grouping is considered ambiguous. A handful of rows, usually
  # from differing missingness, is not a reason to abandon the analysis; a
  # substantial split means the rules were only coincidentally equivalent.
  rule_disagreement_tol = 0.01,

  chron_age_col = "AgeAtFirstDate",

  # Where the clock variables live, first present wins.
  clock_columns = list(
    kdm = list(clock      = c("kdm", "kdm0", "KDM", "kdm_age"),
               advance    = c("kdm_advance", "kdm_adv"),
               advance_sd = "kdm_advance_sd"),
    phenoage = list(clock      = c("phenoage", "PhenoAge", "phenoage_age"),
                    advance    = c("phenoage_advance", "phenoage_adv"),
                    advance_sd = "phenoage_advance_sd")
  ),

  # Attained-age axis: the models' own timescale. A CIF "from time zero" is
  # undefined under delayed entry, so it is anchored at index_age and reads as
  # "cumulative incidence between index_age and t".
  index_age   = 65,
  report_ages = c(70, 75, 80, 85, 90),
  max_age     = 100,
  age_grid_by = 1,

  # Years-since-baseline axis: Shayan's convention. His 5 and 9 year horizons
  # are far short of this cohort's ~22.5 years of follow-up.
  report_years  = c(5, 10, 15, 20),
  max_years     = 25,
  years_grid_by = 0.5,

  # Which axis leads the supplementary figure. Decide after looking at the
  # baseline-age balance between groups that run_cif_block() prints.
  primary_axis = "age",

  # Grid rows thinner than this are suppressed, so the CSVs stay aggregate
  # enough to copy out of the VDI - same discipline as the inspect mode.
  min_risk_set = 10,

  # Rules are required to reproduce the supplied KDM groups exactly. Anything
  # between this and exact is reported as a near miss rather than accepted.
  near_match_pct = 99.5,

  # The publication figures come from plot_cif_results.R, which works off
  # cif_project1_curves.csv and can be iterated on without re-running anything.
  # Set TRUE for quick-look PNGs while debugging a run.
  make_plots = FALSE
)

cif_result_columns <- c(
  "dataset", "outcome", "grouping", "axis", "group", "state",
  "time", "estimate", "se", "CI_lower", "CI_upper",
  "n_risk", "n_event_cum", "n_group", "events_group", "deaths_group"
)

cif_key <- function(dataset, outcome, grouping, axis) {
  paste(dataset, outcome, grouping, axis, sep = "|")
}

completed_cif_keys <- function(results_file) {
  if (!file.exists(results_file)) return(character(0))
  done <- tryCatch(
    read.table(results_file, sep = ";", dec = ",", header = TRUE,
               quote = "", fileEncoding = "Latin1", colClasses = "character"),
    error = function(e) {
      warning("Could not read ", results_file, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(done) || nrow(done) == 0) return(character(0))
  unique(cif_key(done$dataset, done$outcome, done$grouping, done$axis))
}

#' Append CIF rows, writing the header only when creating the file.
append_cif_rows <- function(rows, results_file) {
  if (is.null(rows) || nrow(rows) == 0) return(invisible(NULL))
  rows <- as.data.frame(rows, stringsAsFactors = FALSE)[, cif_result_columns, drop = FALSE]
  new_file <- !file.exists(results_file)
  write.table(
    rows, results_file,
    sep = ";", dec = ",",
    row.names = FALSE,
    col.names = new_file,
    append = !new_file,
    quote = FALSE,
    fileEncoding = "Latin1"
  )
  invisible(rows)
}

# ---------------------------------------------------------------------
# The exposure grouping: reuse Table 1's split, do not invent one
# ---------------------------------------------------------------------

#' Read the per-ID biological age group file written for Table 1.
#'
#' Same read arguments Karolina uses for it interactively, so a file that reads
#' correctly there reads correctly here.
load_cif_groups <- function(candidates = cif_config$group_file_candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Could not find the Table 1 group file. Tried:\n  ",
         paste(candidates, collapse = "\n  "),
         "\nPut table1_kdmgroups_baseline_perID.csv next to this script, or ",
         "edit cif_config$group_file_candidates.")
  }
  log_msg("Reading biological age groups from ", hit[1])
  g <- read.table(hit[1], sep = ";", dec = ",", header = TRUE,
                  stringsAsFactors = FALSE, fileEncoding = "Latin1")
  log_msg("  ", nrow(g), " rows x ", ncol(g), " cols: ",
          paste(names(g), collapse = ", "))
  g
}

#' Work out which column is the ID and which is the group label.
#'
#' Nothing is hard-coded: the ID column is whichever candidate appears in both
#' the group file and the cohort data (matched case-insensitively, since the two
#' sources spell it differently), and the group column is the remaining column
#' with a handful of distinct values.
identify_group_columns <- function(groups, data,
                                   id_candidates = cif_config$id_candidates) {
  g_lower <- tolower(names(groups))
  d_lower <- tolower(names(data))

  shared <- intersect(intersect(g_lower, d_lower), tolower(id_candidates))
  if (length(shared) == 0) {
    # Fall back to any column name shared by both, preferring narrow ones.
    shared <- intersect(g_lower, d_lower)
    if (length(shared) == 0) {
      stop("No column is shared between the group file (",
           paste(names(groups), collapse = ", "),
           ") and the cohort data, so the two cannot be linked. ",
           "Add the ID column name to cif_config$id_candidates.")
    }
    log_msg("  WARNING: no recognised ID name; linking on '", shared[1],
            "' because it is the only shared column. Check this is the person ID.")
  }
  id_in_groups <- names(groups)[match(shared[1], g_lower)]
  id_in_data   <- names(data)[match(shared[1], d_lower)]

  others <- setdiff(names(groups), id_in_groups)
  n_levels <- vapply(others, function(cn) length(unique(groups[[cn]][!is.na(groups[[cn]])])),
                     integer(1))
  binary_ish <- others[n_levels >= 2 & n_levels <= 4]

  if (length(binary_ish) == 0) {
    stop("No column in the group file looks like a 2-4 level group. ",
         "Columns and their distinct-value counts: ",
         paste(sprintf("%s=%d", others, n_levels), collapse = ", "))
  }
  named <- binary_ish[grepl("group|grp|kdm|ba", binary_ish, ignore.case = TRUE)]
  group_col <- if (length(named) > 0) named[1] else binary_ish[1]

  if (length(binary_ish) > 1) {
    log_msg("  Several columns could be the group (",
            paste(binary_ish, collapse = ", "), "); using '", group_col, "'")
  }

  list(id_in_groups = id_in_groups, id_in_data = id_in_data, group_col = group_col)
}

#' Which value of the group column means "biologically older"?
#'
#' Decided by name where the labels say so, otherwise by which level actually
#' has the higher mean age advance - and the choice is always logged, because
#' getting it backwards would silently reverse the figure.
group_older_level <- function(values, advance) {
  levs <- sort(unique(values[!is.na(values)]))
  if (length(levs) != 2) {
    stop("Expected exactly 2 biological age groups, found ", length(levs), ": ",
         paste(levs, collapse = ", "))
  }

  mean_adv <- vapply(levs, function(l) mean(advance[!is.na(values) & values == l],
                                            na.rm = TRUE), numeric(1))
  by_advance <- levs[which.max(mean_adv)]

  named <- levs[grepl("old|high|accel", as.character(levs), ignore.case = TRUE)]
  older <- if (length(named) == 1) named else by_advance

  log_msg("  Group levels: ",
          paste(sprintf("%s (mean advance %.3f)", levs, mean_adv), collapse = ", "))
  log_msg("  Treating '", older, "' as biologically older",
          if (length(named) == 1) " (from the label)" else " (from the higher mean advance)")

  if (!identical(as.character(older), as.character(by_advance))) {
    log_msg("  WARNING: the label and the mean age advance disagree about which ",
            "group is older. Check the group file before using the figure.")
  }
  older
}

# ---------------------------------------------------------------------
# Verifying the grouping rule, then transferring it to PhenoAge
# ---------------------------------------------------------------------
#
# Karolina's question was whether the Table 1 split is on the clock value or on
# the age gap. Rather than guess, every plausible rule is tested against the
# supplied KDM groups; only a rule that reproduces them exactly is accepted, and
# that same rule is then applied to PhenoAge. So the PhenoAge split is the same
# construction as the KDM one rather than a new choice.

cif_group_rules <- list(
  "advance > 0" = list(
    needs = "advance",
    f = function(d, cols) d[[cols$advance]] > 0
  ),
  "advance_sd > 0" = list(
    needs = "advance_sd",
    f = function(d, cols) d[[cols$advance_sd]] > 0
  ),
  "clock > chronological age" = list(
    needs = "clock",
    f = function(d, cols) d[[cols$clock]] > d[[cif_config$chron_age_col]]
  ),
  "advance > median" = list(
    needs = "advance",
    f = function(d, cols) {
      v <- d[[cols$advance]]
      v > median(v, na.rm = TRUE)
    }
  ),
  "advance_sd > median" = list(
    needs = "advance_sd",
    f = function(d, cols) {
      v <- d[[cols$advance_sd]]
      v > median(v, na.rm = TRUE)
    }
  )
)

#' Locate the clock / advance / advance_sd columns for one BA measure.
resolve_clock_columns <- function(data, prefix,
                                 spec = cif_config$clock_columns) {
  if (is.null(spec[[prefix]])) stop("No clock_columns entry for '", prefix, "'")
  out <- lapply(spec[[prefix]], function(cands) first_present(data, cands))
  names(out) <- names(spec[[prefix]])
  out
}

#' Which rules are computable for this prefix?
applicable_rules <- function(data, cols, rules = cif_group_rules) {
  ok <- vapply(rules, function(r) {
    if (is.null(cols[[r$needs]])) return(FALSE)
    if (r$needs == "clock" && !(cif_config$chron_age_col %in% names(data))) return(FALSE)
    TRUE
  }, logical(1))
  rules[ok]
}

#' Score every applicable rule against the supplied grouping.
#'
#' is_older is TRUE / FALSE / NA per row of data, NA meaning the person is not
#' in the group file. Returns the agreement table plus the rules that matched
#' exactly, and never picks a rule for you when none matched.
verify_group_rule <- function(data, is_older, prefix = "kdm",
                              rules = cif_group_rules) {
  cols  <- resolve_clock_columns(data, prefix)
  rules <- applicable_rules(data, cols, rules)

  if (length(rules) == 0) {
    stop("None of the candidate grouping rules can be computed for '", prefix,
         "'. Resolved columns: ",
         paste(sprintf("%s=%s", names(cols),
                       vapply(cols, function(x) if (is.null(x)) "<absent>" else x,
                              character(1))),
               collapse = ", "))
  }

  log_msg("  Candidate rules for ", prefix, ": ", paste(names(rules), collapse = ", "))

  scores <- do.call(rbind, lapply(names(rules), function(nm) {
    pred <- tryCatch(rules[[nm]]$f(data, cols), error = function(e) NULL)
    if (is.null(pred)) {
      return(data.frame(rule = nm, n = 0L, agree = 0L, pct = NA_real_,
                        marks_older_pct = NA_real_, stringsAsFactors = FALSE))
    }
    usable <- !is.na(is_older) & !is.na(pred)
    data.frame(
      rule  = nm,
      n     = sum(usable),
      agree = sum(pred[usable] == is_older[usable]),
      pct   = 100 * sum(pred[usable] == is_older[usable]) / max(sum(usable), 1),
      # Prevalence: what share of people this rule calls biologically older,
      # which is a completely different quantity from the agreement and is
      # reported alongside it so the two cannot be confused. A median split is
      # 50% by construction; a split at zero is near 50% too, because the age
      # gap is centred near zero by the way the clocks are calibrated.
      marks_older_pct = 100 * sum(pred, na.rm = TRUE) / max(sum(!is.na(pred)), 1),
      stringsAsFactors = FALSE
    )
  }))
  scores <- scores[order(-scores$pct), , drop = FALSE]

  log_msg("  Each candidate rule for ", prefix, ", scored two ways:")
  log_msg("    'agrees' = share of people it puts in the SAME group as Table 1")
  log_msg("    'older'  = share of people the rule itself calls older ",
          "(a prevalence, not an agreement)")
  log_msg(sprintf("    %-26s %8s %8s", "rule", "agrees", "older"))
  log_msg(sprintf("    %-26s %8s %7.1f%%   <- Table 1's own grouping",
                  "(the supplied groups)", "-",
                  100 * sum(is_older, na.rm = TRUE) / max(sum(!is.na(is_older)), 1)))
  for (i in seq_len(nrow(scores))) {
    log_msg(sprintf("    %-26s %7.2f%% %7.1f%%   (%d / %d agree)",
                    scores$rule[i], scores$pct[i], scores$marks_older_pct[i],
                    scores$agree[i], scores$n[i]))
  }

  exact <- scores$rule[!is.na(scores$pct) & scores$n > 0 & scores$agree == scores$n]
  list(scores = scores, exact = exact, cols = cols, rules = rules)
}

#' Apply the verified rule(s) to another clock, measuring any disagreement.
#'
#' Several rules can reproduce the KDM groups exactly - a median split is
#' invariant under any monotone rescaling, so "advance > median" and
#' "advance_sd > median" are the same statement whenever the standardisation is
#' a single global linear transform. If they then diverge on another clock, the
#' usual causes are a standardisation applied within strata (so it is not a
#' global transform) or differing missingness between the two variables (so the
#' two medians are taken over different samples).
#'
#' Neither is automatically a reason to abandon the grouping, so this quantifies
#' the disagreement rather than refusing on sight: it picks the rule built on the
#' variable the published models use, reports exactly how many people the other
#' rules would classify differently, and stops only if that is more than
#' cif_config$rule_disagreement_tol of the cohort.
transfer_group_rule <- function(data, rule_names, prefix, rules = cif_group_rules,
                                tol = cif_config$rule_disagreement_tol) {
  cols <- resolve_clock_columns(data, prefix)
  usable <- rule_names[vapply(rule_names, function(nm) {
    r <- rules[[nm]]
    !is.null(r) && !is.null(cols[[r$needs]]) &&
      (r$needs != "clock" || cif_config$chron_age_col %in% names(data))
  }, logical(1))]

  if (length(usable) == 0) {
    stop("The rule(s) that matched the KDM groups (",
         paste(rule_names, collapse = ", "),
         ") cannot be computed for '", prefix,
         "' - resolved columns: ",
         paste(sprintf("%s=%s", names(cols),
                       vapply(cols, function(x) if (is.null(x)) "<absent>" else x,
                              character(1))),
               collapse = ", "))
  }

  preds <- lapply(usable, function(nm) rules[[nm]]$f(data, cols))
  names(preds) <- usable

  # Preference order, so the choice is a documented decision rather than
  # whichever rule happened to be listed first.
  chosen <- usable[order(match(usable, cif_config$prefer_rule_order,
                               nomatch = length(cif_config$prefer_rule_order) + 1L))][1]

  if (length(usable) > 1) {
    log_msg("  Several rules reproduced the KDM groups exactly; comparing them ",
            "on '", prefix, "':")
    for (k in names(cols)) {
      if (!is.null(cols[[k]])) {
        log_msg(sprintf("    %-22s %d missing of %d",
                        cols[[k]], sum(is.na(data[[cols[[k]]]])), nrow(data)))
      }
    }

    base <- preds[[chosen]]
    worst <- 0
    for (nm in setdiff(usable, chosen)) {
      p    <- preds[[nm]]
      both <- !is.na(base) & !is.na(p)
      differ    <- sum(base[both] != p[both])
      na_differ <- sum(is.na(base) != is.na(p))
      pct <- 100 * (differ + na_differ) / nrow(data)
      worst <- max(worst, pct)
      log_msg(sprintf(
        "    '%s' vs '%s': %d of %d jointly classified differ, plus %d with differing missingness = %.3f%% of the cohort",
        chosen, nm, differ, sum(both), na_differ, pct))
    }

    if (worst > 100 * tol) {
      stop("The rules that reproduced the KDM groups exactly (",
           paste(usable, collapse = ", "),
           ") disagree for '", prefix, "' on ", sprintf("%.3f%%", worst),
           " of the cohort, which is above the ", sprintf("%.3f%%", 100 * tol),
           " tolerance. They were only coincidentally equivalent, so the ",
           prefix, " grouping cannot be derived safely. Supply the Table 1 ",
           "script, or raise cif_config$rule_disagreement_tol deliberately.")
    }
    log_msg(sprintf("    Largest disagreement %.3f%% is within the %.3f%% ",
                    worst, 100 * tol),
            "tolerance, so the choice of rule is immaterial here.")
  }

  log_msg("  ", prefix, " groups derived with: ", chosen)
  preds[[chosen]]
}

#' Merge the Table 1 groups onto the cohort and build both grouping columns.
#'
#' Returns the data with kdm_group and phenoage_group added, plus the merge
#' diagnostics. The KDM groups come from the file verbatim - the rule is only
#' used to carry the same definition across to PhenoAge.
attach_cif_groups <- function(data, groups = NULL) {
  if (is.null(groups)) groups <- load_cif_groups()

  cn <- identify_group_columns(groups, data)
  log_msg("  Linking on ", cn$id_in_data, " (data) <-> ", cn$id_in_groups,
          " (group file); group column '", cn$group_col, "'")

  keys_data   <- as.character(data[[cn$id_in_data]])
  keys_groups <- as.character(groups[[cn$id_in_groups]])

  dup <- sum(duplicated(keys_groups))
  if (dup > 0) {
    log_msg("  WARNING: ", dup, " duplicated IDs in the group file; keeping the first")
    groups <- groups[!duplicated(keys_groups), , drop = FALSE]
    keys_groups <- as.character(groups[[cn$id_in_groups]])
  }

  raw_group <- groups[[cn$group_col]][match(keys_data, keys_groups)]

  matched   <- sum(!is.na(raw_group))
  unmatched <- sum(is.na(raw_group))
  log_msg("  Matched ", matched, " of ", nrow(data), " cohort rows (",
          unmatched, " not in the group file; they are excluded from the ",
          "curves and counted, never silently dropped)")

  kdm_cols <- resolve_clock_columns(data, "kdm")
  adv_for_direction <- if (!is.null(kdm_cols$advance_sd)) {
    data[[kdm_cols$advance_sd]]
  } else if (!is.null(kdm_cols$advance)) {
    data[[kdm_cols$advance]]
  } else {
    stop("Need a KDM advance column to check which group is biologically older")
  }

  older_level <- group_older_level(raw_group, adv_for_direction)
  is_older <- ifelse(is.na(raw_group), NA, raw_group == older_level)

  verdict <- verify_group_rule(data, is_older, prefix = "kdm")

  if (length(verdict$exact) == 0) {
    best <- verdict$scores[1, ]
    msg <- paste0(
      "No candidate rule reproduces the Table 1 KDM groups exactly. ",
      "Best was '", best$rule, "' at ", sprintf("%.2f%%", best$pct), ". "
    )
    if (!is.na(best$pct) && best$pct >= cif_config$near_match_pct) {
      msg <- paste0(msg,
        "That is a near miss, which usually means the same rule computed on a ",
        "different sample (e.g. a median taken within the Table 1 subset). ",
        "The PhenoAge grouping cannot be derived safely from a near miss - ",
        "supply the Table 1 script so the exact definition can be reused.")
    } else {
      msg <- paste0(msg,
        "Supply the Table 1 script so the exact definition can be reused, or ",
        "add the rule to cif_group_rules.")
    }
    stop(msg)
  }

  log_msg("  KDM groups reproduced exactly by: ", paste(verdict$exact, collapse = ", "))

  pheno_older <- transfer_group_rule(data, verdict$exact, prefix = "phenoage")

  lab <- function(x) factor(
    ifelse(is.na(x), NA_character_,
           ifelse(x, "Biologically older", "Biologically younger")),
    levels = c("Biologically younger", "Biologically older")
  )

  data$kdm_group      <- lab(is_older)
  # PhenoAge uses the KDM rule, but only for people who have a Table 1 group,
  # so both groupings describe the same set of people.
  data$phenoage_group <- lab(ifelse(is.na(is_older), NA, pheno_older))

  log_msg("  kdm_group      : ", paste(sprintf("%s=%d", levels(data$kdm_group),
          tabulate(data$kdm_group, 2)), collapse = ", "))
  log_msg("  phenoage_group : ", paste(sprintf("%s=%d", levels(data$phenoage_group),
          tabulate(data$phenoage_group, 2)), collapse = ", "))
  log_msg("  Groups agree for ",
          sum(data$kdm_group == data$phenoage_group, na.rm = TRUE), " of ",
          sum(!is.na(data$kdm_group)), " classified people")

  # The two groupings should describe the same people. They will not if PhenoAge
  # is missing for someone who does have a Table 1 KDM group, which would make
  # the two panels of the figure rest on different samples.
  pheno_gap <- sum(!is.na(data$kdm_group) & is.na(data$phenoage_group))
  if (pheno_gap > 0) {
    log_msg("  WARNING: ", pheno_gap, " people have a KDM group but no PhenoAge ",
            "value, so the phenoage_group curves cover a slightly smaller ",
            "sample than the kdm_group ones. Report both Ns.")
  }

  attr(data, "cif_rule")  <- verdict$exact[1]
  attr(data, "cif_merge") <- list(matched = matched, unmatched = unmatched,
                                  older_level = older_level)
  data
}

# ---------------------------------------------------------------------
# Aalen-Johansen cumulative incidence
# ---------------------------------------------------------------------

#' Step a survfit / survfitms object onto a fixed time grid.
#'
#' Deliberately avoids summary.survfit(): its output shape has shifted between
#' survival versions, and walking the fit by hand is both version-proof and
#' explicit about what happens before the first event (estimate 0, risk set as
#' at entry).
read_fit_onto_grid <- function(fit, grid, default_stratum = "all") {
  multistate <- !is.null(fit$states)

  # survfit omits $strata entirely when the grouping has a single level, so the
  # caller has to say what that one stratum is called. Without this the label
  # would silently become "all" and nothing would join back to it.
  lens <- fit$strata
  if (is.null(lens)) lens <- stats::setNames(length(fit$time), default_stratum)
  ends   <- cumsum(lens)
  starts <- ends - lens + 1L

  mat <- function(x, n) {
    if (is.null(x)) return(NULL)
    if (is.matrix(x)) x else matrix(x, nrow = n)
  }

  n_all <- length(fit$time)
  if (multistate) {
    est_m <- mat(fit$pstate,  n_all)
    lo_m  <- mat(fit$lower,   n_all)
    hi_m  <- mat(fit$upper,   n_all)
    se_m  <- mat(fit$std.err, n_all)
    states <- fit$states
  } else {
    # Plain Kaplan-Meier: report 1 - S as the naive cumulative incidence, which
    # flips the confidence limits.
    est_m <- mat(1 - fit$surv, n_all)
    lo_m  <- if (is.null(fit$upper)) NULL else mat(1 - fit$upper, n_all)
    hi_m  <- if (is.null(fit$lower)) NULL else mat(1 - fit$lower, n_all)
    se_m  <- NULL
    states <- "naive_km"
  }
  nr_m <- mat(fit$n.risk,  n_all)
  ne_m <- mat(fit$n.event, n_all)

  risk_col <- if (multistate && "(s0)" %in% colnames(nr_m)) {
    match("(s0)", colnames(nr_m))
  } else 1L

  out <- list()
  for (k in seq_along(lens)) {
    rows <- starts[k]:ends[k]
    tt   <- fit$time[rows]
    j    <- findInterval(grid, tt)
    jj   <- pmax(j, 1L)
    before <- j == 0L

    nr <- ifelse(before, nr_m[rows[1], risk_col], nr_m[rows, risk_col][jj])

    for (ci in seq_along(states)) {
      st <- states[ci]
      # Skip the initial ("still event-free") state. It is always first for a
      # survfitms built from a factor status, and is usually but not always
      # spelled "(s0)" - so do not rely on the name alone.
      if (multistate && (ci == 1L || identical(st, "(s0)"))) next

      est <- ifelse(before, 0, est_m[rows, ci, drop = TRUE][jj])
      lo  <- if (is.null(lo_m)) NA_real_ else ifelse(before, NA_real_, lo_m[rows, ci, drop = TRUE][jj])
      hi  <- if (is.null(hi_m)) NA_real_ else ifelse(before, NA_real_, hi_m[rows, ci, drop = TRUE][jj])
      se  <- if (is.null(se_m)) NA_real_ else ifelse(before, NA_real_, se_m[rows, ci, drop = TRUE][jj])

      ev_col <- if (!is.null(colnames(ne_m)) && st %in% colnames(ne_m)) {
        match(st, colnames(ne_m))
      } else if (ncol(ne_m) >= ci) ci else 1L
      cum_ev <- cumsum(ne_m[rows, ev_col, drop = TRUE])
      nev <- ifelse(before, 0, cum_ev[jj])

      out[[length(out) + 1L]] <- data.frame(
        stratum = names(lens)[k],
        state   = st,
        time    = grid,
        estimate = as.numeric(est),
        se       = as.numeric(se),
        CI_lower = as.numeric(lo),
        CI_upper = as.numeric(hi),
        n_risk      = as.numeric(nr),
        n_event_cum = as.numeric(nev),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

#' Put one outcome on one time axis, ready for survfit.
#'
#' age   - attained age with delayed entry, anchored at index_age. People who
#'         already left before index_age are dropped; people who entered before
#'         it have their clock started there.
#' years - time since baseline, which has no left truncation at all. Shayan's
#'         axis; kept as the secondary panel.
prepare_cif_axis <- function(df, axis, index_age = cif_config$index_age) {
  if (identical(axis, "age")) {
    d <- df[df$AgeAtLastDate > index_age, , drop = FALSE]
    d$cif_start <- pmax(d$AgeAtFirstDate, index_age)
    d$cif_stop  <- d$AgeAtLastDate
  } else if (identical(axis, "years")) {
    d <- df
    d$cif_start <- 0
    d$cif_stop  <- d$AgeAtLastDate - d$AgeAtFirstDate
  } else {
    stop("Unknown axis '", axis, "' - expected 'age' or 'years'")
  }
  d[d$cif_stop > d$cif_start, , drop = FALSE]
}

cif_grid <- function(axis, index_age = cif_config$index_age) {
  if (identical(axis, "age")) {
    sort(unique(c(seq(index_age, cif_config$max_age, by = cif_config$age_grid_by),
                  cif_config$report_ages)))
  } else {
    sort(unique(c(seq(0, cif_config$max_years, by = cif_config$years_grid_by),
                  cif_config$report_years)))
  }
}

#' Aalen-Johansen CIF for both event types, plus the naive 1 - KM comparison.
#'
#' The naive curve treats death as censoring, which is exactly the mistake the
#' reviewer is warning about, so reporting both quantifies how much it matters.
cif_for_outcome <- function(df, group_var, axis,
                            index_age = cif_config$index_age) {
  d <- prepare_cif_axis(df, axis, index_age)
  d$cif_grp <- d[[group_var]]
  d <- d[!is.na(d$cif_grp), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$cif_grp <- droplevels(factor(d$cif_grp))

  grid <- cif_grid(axis, index_age)
  lone <- paste0("cif_grp=", levels(d$cif_grp)[1])

  # Count what actually happened, not what the factor levels allow: survfit
  # derives its state names from the levels, so an outcome with zero events
  # still reports an "event" state and would yield a flat-zero curve.
  n_event <- sum(d$fg_status == "event")
  n_death <- sum(d$fg_status == "death")

  if (n_event == 0) {
    stop("No events of this outcome among the ", nrow(d),
         " people left after preparing the '", axis,
         "' axis, so there is no cumulative incidence to estimate.")
  }
  if (n_death == 0) {
    log_msg("    WARNING: no competing deaths here, so the CIF reduces to ",
            "1 - KM. Check the event coding before using this curve.")
  }

  aj <- suppressMessages(survival::survfit(
    Surv(cif_start, cif_stop, fg_status) ~ cif_grp, data = d, id = fg_id
  ))

  # The multi-state interface (a factor status giving the Aalen-Johansen
  # estimator) needs survival >= 3.0. Without it survfit returns an ordinary
  # Kaplan-Meier fit, which read_fit_onto_grid would report as 1 - S - silently
  # handing back the naive curve as if it were the CIF. Refuse that outright.
  if (is.null(aj$states)) {
    stop("survfit() returned an ordinary survival fit, not a multi-state one, ",
         "so this would be 1 - KM rather than the Aalen-Johansen CIF. ",
         "The multi-state interface needs survival >= 3.0; this session has ",
         as.character(utils::packageVersion("survival")), ".")
  }

  aj_rows <- read_fit_onto_grid(aj, grid, default_stratum = lone)

  naive <- suppressMessages(survival::survfit(
    Surv(cif_start, cif_stop, fg_status == "event") ~ cif_grp, data = d
  ))
  naive_rows <- read_fit_onto_grid(naive, grid, default_stratum = lone)

  rows <- rbind(aj_rows, naive_rows)
  rows$group <- sub("^cif_grp=", "", rows$stratum)
  rows$stratum <- NULL

  # Per-group denominators, so every curve can be read without the raw data.
  per_group <- do.call(rbind, lapply(levels(d$cif_grp), function(g) {
    sub <- d[d$cif_grp == g, , drop = FALSE]
    data.frame(group = g, n_group = nrow(sub),
               events_group = sum(sub$fg_status == "event"),
               deaths_group = sum(sub$fg_status == "death"),
               stringsAsFactors = FALSE)
  }))
  rows <- merge(rows, per_group, by = "group", all.x = TRUE)

  # Keep the state names readable in the CSV.
  rows$state[rows$state == "event"] <- "dementia"
  rows$state[rows$state == "death"] <- "death"

  attr(rows, "n_used") <- nrow(d)
  rows
}

#' Absolute risk differences between the two groups at the reporting times.
#'
#' se comes from survfit's std.err, which for a multi-state fit is on the
#' probability-in-state scale, so a Wald interval on the difference is the
#' natural thing. The naive curve is excluded: for a plain Kaplan-Meier fit
#' std.err is on the cumulative-hazard scale, not this one.
cif_differences <- function(rows, axis, report_times,
                            min_risk_set = cif_config$min_risk_set) {
  keep <- rows$time %in% report_times & rows$state %in% c("dementia", "death") &
    rows$n_risk >= min_risk_set
  r <- rows[keep, , drop = FALSE]
  if (nrow(r) == 0) return(NULL)

  out <- list()
  for (st in unique(r$state)) {
    for (tt in sort(unique(r$time[r$state == st]))) {
      sub <- r[r$state == st & r$time == tt, , drop = FALSE]
      young <- sub[sub$group == "Biologically younger", , drop = FALSE]
      old   <- sub[sub$group == "Biologically older", , drop = FALSE]
      if (nrow(young) != 1 || nrow(old) != 1) next

      diff <- old$estimate - young$estimate
      se   <- sqrt(old$se^2 + young$se^2)
      out[[length(out) + 1L]] <- data.frame(
        axis = axis, state = st, time = tt,
        est_younger = young$estimate, est_older = old$estimate,
        difference = diff,
        diff_lower = diff - 1.96 * se,
        diff_upper = diff + 1.96 * se,
        n_risk_younger = young$n_risk, n_risk_older = old$n_risk,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

#' Dementia and death CIFs side by side with the naive comparison.
#'
#' Base graphics on purpose - the README promises survival and nothing else.
plot_cif <- function(rows, file, title, axis) {
  # Reference group first, so "younger" gets the cool colour and reads as the
  # baseline. Alphabetical order would put "older" first and invert that.
  preferred <- c("Biologically younger", "Biologically older")
  present   <- unique(rows$group)
  groups    <- c(intersect(preferred, present), sort(setdiff(present, preferred)))
  cols      <- stats::setNames(c("#2166ac", "#b2182b")[seq_along(groups)], groups)
  xlab   <- if (identical(axis, "age")) "Attained age, years" else
            "Years since baseline"

  ok <- is.finite(rows$estimate)
  ymax <- max(rows$estimate[ok & rows$state %in% c("dementia", "death")], 0.05, na.rm = TRUE)

  grDevices::png(file, width = 1800, height = 800, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  # lend goes on par(), not on the individual calls: legend() takes no "..." and
  # would reject it, which would leave the key drawn with square dot ends while
  # the curves had round ones. Set here, everything inherits it - lines and
  # legend swatches alike.
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1), lend = "round")

  # Second curve is a dot pattern: "12" is one unit on, two off, both scaled by
  # lwd, so the gap is twice the dot diameter - tight enough to read as a
  # continuous curve, open enough never to be mistaken for the solid one. With
  # par(lend = "round") the dots draw as circles, and legend() accepts the same
  # pattern string, so the key matches the curves exactly.
  lty_for <- c("solid", "12")
  lwd_for <- c(2, 2.8)

  draw <- function(states, main, ymax_local) {
    plot(NA, xlim = range(rows$time), ylim = c(0, ymax_local),
         xlab = xlab, ylab = "Cumulative incidence", main = main)
    for (g in groups) {
      for (i in seq_along(states)) {
        sub <- rows[rows$group == g & rows$state == states[i] & is.finite(rows$estimate), ]
        sub <- sub[order(sub$time), ]
        if (nrow(sub) > 0) {
          graphics::lines(sub$time, sub$estimate, col = cols[[g]],
                          lty = lty_for[i], lwd = lwd_for[i], type = "s")
        }
      }
    }
  }

  # Group colour and curve type are separate legends, so neither has to encode
  # both and the reader is never asked to decode a four-way key.
  # seg.len widens the line sample so a spaced dot pattern is recognisable in
  # the key; at the default length it would show only a dot or two.
  add_legends <- function(labels) {
    graphics::legend("topleft", legend = groups, col = cols[groups],
                     lwd = 2, bty = "n")
    graphics::legend("bottomright", legend = labels, col = "grey30",
                     lty = lty_for, lwd = lwd_for, seg.len = 3.5, bty = "n")
  }

  draw(c("dementia", "death"), title, ymax)
  add_legends(c("Dementia", "Death"))

  dem_max <- max(rows$estimate[ok & rows$state %in% c("dementia", "naive_km")],
                 0.05, na.rm = TRUE)
  draw(c("dementia", "naive_km"), "Dementia only", dem_max)
  add_legends(c("Aalen-Johansen", "Naive 1 - KM"))

  invisible(file)
}

# ---------------------------------------------------------------------
# Orchestration for the cif mode
# ---------------------------------------------------------------------

#' Baseline age by group - the check that decides which axis leads the figure.
#'
#' A years-since-baseline CIF is only interpretable if the two groups enter at
#' similar ages; otherwise it averages absolute risks over different age bands.
describe_baseline_age_by_group <- function(df, group_var) {
  g <- droplevels(factor(df[[group_var]]))
  levs <- levels(g)
  log_msg("  Baseline age by ", group_var, ":")
  means <- stats::setNames(rep(NA_real_, length(levs)), levs)
  for (lev in levs) {
    a <- df$AgeAtFirstDate[!is.na(g) & g == lev]
    means[[lev]] <- mean(a)
    log_msg(sprintf("    %-22s n=%6d  mean=%5.1f  sd=%4.1f  range=%.1f-%.1f",
                    lev, length(a), mean(a), stats::sd(a), min(a), max(a)))
  }
  if (length(levs) == 2) {
    log_msg(sprintf("    difference in mean baseline age: %.2f years",
                    means[[2]] - means[[1]]))
    log_msg("    A large difference here means the years-since-baseline curves ",
            "average absolute risks over different age bands, so the ",
            "attained-age axis should lead the figure.")
  }
  invisible(NULL)
}

#' Did the exit-age rounding change any event classification?
#'
#' The published Fine-Gray run rounded exit ages to ~5 weeks so the expansion
#' would fit in memory. Rounding can create ties, and ties resolve to dementia
#' before death, so in principle it can move someone between event categories.
#' Whether it actually did is answerable for free, because clean_fg_data() is
#' cheap without the expansion. Identical counts mean the rounding only
#' collapsed tied times and reclassified nobody - which is the reassurance
#' Supplementary Table 9 needs, and which was never checked when it was run.
#' Find the person ID column, case-insensitively.
#'
#' Walks the candidates in their configured order rather than the data's column
#' order, so when a cohort carries more than one identifier the preferred one
#' wins - sampleID, which is what the Table 1 group file keys on.
find_id_column <- function(data, candidates = cif_config$id_candidates) {
  lower <- tolower(names(data))
  for (cand in tolower(candidates)) {
    hit <- which(lower == cand)
    if (length(hit) > 0) return(names(data)[hit[1]])
  }
  NULL
}

compare_rounding_effect <- function(data, specs, dataset_label,
                                    rounded = fg_config$time_round_digits,
                                    exact = NULL,
                                    id_col = find_id_column(data)) {
  log_msg("Effect of the exit-age rounding (published run used: ",
          if (is.null(rounded)) "none" else rounded, "):")

  # Counts alone cannot separate the two things rounding can do: EXCLUDE
  # someone (follow-up rounds to zero, so exit_age > age_entry fails) or
  # RECLASSIFY someone (rounding creates a tie, and ties resolve to dementia
  # before death). Those have completely different implications for the
  # Methods, so compare per person rather than in aggregate.
  if (is.null(id_col)) {
    log_msg("  WARNING: no person ID column found (looked for ",
            paste(cif_config$id_candidates, collapse = ", "),
            "), so excluded and reclassified people cannot be told apart. ",
            "Reporting counts only.")
  }

  out <- do.call(rbind, lapply(specs, function(spec) {
    ev <- if (is.null(id_col)) character(0) else id_col
    a  <- clean_fg_data(data, spec$event_var, dataset_label,
                        time_round_digits = rounded, extra_vars = ev)
    b  <- clean_fg_data(data, spec$event_var, dataset_label,
                        time_round_digits = exact,   extra_vars = ev)
    ca <- fg_counts(a); cb <- fg_counts(b)

    dropped <- gained <- reclassified <- NA_integer_
    if (!is.null(id_col)) {
      ka <- as.character(a[[id_col]])
      kb <- as.character(b[[id_col]])
      dropped <- length(setdiff(kb, ka))   # in the exact run, lost by rounding
      gained  <- length(setdiff(ka, kb))
      both    <- intersect(ka, kb)
      reclassified <- sum(
        as.character(a$fg_status)[match(both, ka)] !=
        as.character(b$fg_status)[match(both, kb)]
      )
    }
    rm(a, b); gc()

    log_msg(sprintf(
      "  %-12s rounded: N=%d ev=%d death=%d cens=%d | exact: N=%d ev=%d death=%d cens=%d",
      spec$outcome_label,
      ca$N, ca$events, ca$competing_deaths, ca$censored,
      cb$N, cb$events, cb$competing_deaths, cb$censored))
    if (!is.null(id_col)) {
      log_msg(sprintf(
        "  %-12s   -> %d excluded by rounding, %d gained, %d of %d shared people reclassified",
        "", dropped, gained, reclassified, length(both)))
    }

    data.frame(
      outcome = spec$outcome_label,
      N_rounded = ca$N, events_rounded = ca$events,
      deaths_rounded = ca$competing_deaths, censored_rounded = ca$censored,
      N_exact = cb$N, events_exact = cb$events,
      deaths_exact = cb$competing_deaths, censored_exact = cb$censored,
      excluded_by_rounding = dropped, gained_by_rounding = gained,
      reclassified = reclassified,
      counts_identical = identical(ca, cb),
      stringsAsFactors = FALSE
    )
  }))

  if (all(out$counts_identical)) {
    log_msg("  The rounding changed nothing at all.")
  } else if (!is.null(id_col) && all(out$reclassified == 0, na.rm = TRUE)) {
    log_msg("  The rounding reclassified NOBODY. Its only effect was to ",
            "exclude ", max(out$excluded_by_rounding, na.rm = TRUE),
            " people whose follow-up rounded to zero, i.e. under ~5 weeks. ",
            "That is an exclusion to state in the Methods, not a ",
            "misclassification, and the event definitions behind ",
            "Supplementary Table 9 are unaffected.")
  } else if (!is.null(id_col)) {
    log_msg("  WARNING: the rounding reclassified ",
            sum(out$reclassified, na.rm = TRUE),
            " person-outcomes between event categories. That belongs in the ",
            "Methods explicitly.")
  } else {
    log_msg("  Counts differ, but without an ID column the cause cannot be ",
            "attributed. Re-run with the person ID present.")
  }
  out
}

#' Run every CIF for one dataset, appending as it goes.
run_cif_block <- function(data, specs, dataset_label,
                          curves_file = "cif_project1_curves.csv",
                          summary_file = "cif_project1_summary.csv",
                          groupings = c("kdm_group", "phenoage_group"),
                          axes = c("age", "years"),
                          only_outcomes = NULL,
                          index_age = cif_config$index_age) {

  if (!is.null(only_outcomes)) {
    specs <- Filter(function(s) s$outcome_label %in% only_outcomes, specs)
    if (length(specs) == 0) {
      stop("No outcomes matched: ", paste(only_outcomes, collapse = ", "))
    }
  }

  done <- completed_cif_keys(curves_file)
  if (length(done) > 0) {
    log_msg("Resuming: ", length(done), " curve sets already in ", curves_file)
  }

  summaries <- list()

  for (spec in specs) {
    todo <- expand.grid(grouping = groupings, axis = axes,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    todo$key <- cif_key(dataset_label, spec$outcome_label, todo$grouping, todo$axis)
    todo <- todo[!(todo$key %in% done), , drop = FALSE]
    if (nrow(todo) == 0) {
      log_msg("Skipping ", dataset_label, " / ", spec$outcome_label, " (already complete)")
      next
    }

    log_msg("=== ", dataset_label, " / ", spec$outcome_label, " ===")
    df <- clean_fg_data(data, spec$event_var, dataset_label,
                        extra_vars = groupings)
    counts <- fg_counts(df)
    log_msg("N=", counts$N, "  events=", counts$events,
            "  competing deaths=", counts$competing_deaths,
            "  censored=", counts$censored)
    log_msg("  without a Table 1 group: ", sum(is.na(df[[groupings[1]]])),
            " (excluded from the curves)")
    describe_baseline_age_by_group(df, groupings[1])

    for (i in seq_len(nrow(todo))) {
      grouping <- todo$grouping[i]
      axis     <- todo$axis[i]
      label <- paste(dataset_label, spec$outcome_label, grouping, axis, sep = " / ")
      log_msg("  ", label)

      rows <- tryCatch(cif_for_outcome(df, grouping, axis, index_age),
                       error = function(e) e)
      if (inherits(rows, "error")) {
        log_msg("    FAILED: ", conditionMessage(rows))
        next
      }
      if (is.null(rows) || nrow(rows) == 0) {
        log_msg("    no usable rows")
        next
      }

      report_times <- if (identical(axis, "age")) cif_config$report_ages else
                      cif_config$report_years
      diffs <- cif_differences(rows, axis, report_times)
      if (!is.null(diffs)) {
        diffs$dataset  <- dataset_label
        diffs$outcome  <- spec$outcome_label
        diffs$grouping <- grouping
        summaries[[length(summaries) + 1L]] <- diffs
      }

      # Suppress thin risk sets so the CSV stays safe to copy out of the VDI.
      thin <- rows$n_risk < cif_config$min_risk_set
      if (any(thin)) {
        log_msg("    suppressing ", sum(thin), " grid rows with n_risk < ",
                cif_config$min_risk_set)
        rows <- rows[!thin, , drop = FALSE]
      }
      if (nrow(rows) == 0) {
        log_msg("    everything suppressed; nothing written")
        next
      }

      rows$dataset  <- dataset_label
      rows$outcome  <- spec$outcome_label
      rows$grouping <- grouping
      rows$axis     <- axis
      append_cif_rows(rows, curves_file)

      if (isTRUE(cif_config$make_plots)) {
        png_file <- paste0("cif_project1_", dataset_label, "_",
                           spec$outcome_label, "_", grouping, "_", axis, ".png")
        ok <- tryCatch({
          plot_cif(rows, png_file,
                   paste(spec$outcome_label, "-", grouping), axis)
          TRUE
        }, error = function(e) { log_msg("    plot failed: ", conditionMessage(e)); FALSE })
        if (ok) log_msg("    wrote ", png_file)
      }

      done <- c(done, todo$key[i])
      gc()
    }
    rm(df); gc()
  }

  if (length(summaries) > 0) {
    all_diffs <- do.call(rbind, summaries)
    cols <- c("dataset", "outcome", "grouping", "axis", "state", "time",
              "est_younger", "est_older", "difference", "diff_lower", "diff_upper",
              "n_risk_younger", "n_risk_older")
    write.table(all_diffs[, cols, drop = FALSE], summary_file,
                sep = ";", dec = ",", row.names = FALSE,
                col.names = !file.exists(summary_file),
                append = file.exists(summary_file),
                quote = FALSE, fileEncoding = "Latin1")
    log_msg("Wrote ", nrow(all_diffs), " summary rows to ", summary_file)
  }

  invisible(done)
}

# =======================================================================
# Command-line entry point
# =======================================================================
#
# Only runs when this file is executed directly via Rscript (checked by
# looking for "fg_project1.R" in the --file= argument Rscript passes to
# itself). When test_fg_project1.R does source("fg_project1.R") to reuse the
# functions above, --file= points at test_fg_project1.R instead, so this
# block is skipped and no analysis runs as a side effect of testing.

.fg_cli_raw  <- commandArgs(trailingOnly = FALSE)
.fg_file_arg <- sub("^--file=", "", .fg_cli_raw[grepl("^--file=", .fg_cli_raw)])
# Any fg_project1*.R counts as a direct run, so a copy carrying a suffix - e.g.
# fg_project1_cif.R while transferring a new version to the VDI - still works.
# Matching only the exact name would make a renamed copy exit silently with no
# output and no error, which is the worst possible failure here.
# test_fg_project1*.R does NOT match, because of the leading anchor, so
# source()ing this file from the tests still runs no analysis.
.fg_invoked_directly <- length(.fg_file_arg) > 0 &&
  grepl("^fg_project1.*\\.R$", sub(".*[\\\\/]", "", .fg_file_arg))

if (.fg_invoked_directly) {
  options(warn = 1)

  .fg_args <- commandArgs(trailingOnly = TRUE)
  .fg_mode <- if (length(.fg_args) >= 1) .fg_args[1] else ""
  .fg_only_outcomes <- if (length(.fg_args) >= 2) .fg_args[-1] else NULL

  .fg_usage <- function() {
    cat("Usage:\n")
    cat("  Rscript fg_project1.R main   [outcome ...]\n")
    cat("      outcomes: ADRD Alzheimers Vascular Other\n")
    cat("  Rscript fg_project1.R svedem [outcome ...]\n")
    cat("      outcomes: S_ADRD S_e_alz S_l_alz S_alz_vas S_vas S_lew S_fro S_par S_uns S_oth\n")
    cat("  Rscript fg_project1.R inspect\n")
    cat("      data diagnostics - run this once before trusting main/svedem output\n")
    cat("  Rscript fg_project1.R cif    [outcome ...]\n")
    cat("      Aalen-Johansen cumulative incidence by biological age group;\n")
    cat("      absolute risk for reviewer comment 4. Needs\n")
    cat("      table1_kdmgroups_baseline_perID.csv. No finegray() expansion,\n")
    cat("      so this does not touch the Fine-Gray results.\n")
  }

  if (.fg_mode == "main") {
    # Same behaviour as the old fg_project1_main.R: results appended to
    # finegray_project1_main.csv one row at a time; completed models are
    # skipped if this is a restart after a crash.
    fg_log_file <- "finegray_project1_main.log"
    install_crash_handler("fg_project1.R main")

    log_msg("Starting main Fine-Gray analysis in ", getwd())
    if (!is.null(.fg_only_outcomes)) {
      log_msg("Restricted to: ", paste(.fg_only_outcomes, collapse = ", "))
    }

    .fg_start <- Sys.time()
    .fg_datasets <- load_prepared_data(required = "to1po2")

    run_outcome_block(
      data          = .fg_datasets$to1po2,
      specs         = main_specs,
      dataset_label = "Main",
      results_file  = "finegray_project1_main.csv",
      only_outcomes = .fg_only_outcomes
    )

    log_msg("Done in ", round(elapsed_min(.fg_start), 2),
            " min. Results: finegray_project1_main.csv")

  } else if (.fg_mode == "svedem") {
    # Same behaviour as the old fg_project1_svedem.R: results appended to
    # finegray_project1_svedem.csv one row at a time.
    fg_log_file <- "finegray_project1_svedem.log"
    install_crash_handler("fg_project1.R svedem")

    log_msg("Starting SveDem Fine-Gray analysis in ", getwd())
    if (!is.null(.fg_only_outcomes)) {
      log_msg("Restricted to: ", paste(.fg_only_outcomes, collapse = ", "))
    }

    .fg_start <- Sys.time()
    .fg_datasets <- load_prepared_data(required = "to1po2_SveDem_D")

    run_outcome_block(
      data          = .fg_datasets$to1po2_SveDem_D,
      specs         = svedem_specs,
      dataset_label = "SveDem",
      results_file  = "finegray_project1_svedem.csv",
      only_outcomes = .fg_only_outcomes
    )

    log_msg("Done in ", round(elapsed_min(.fg_start), 2),
            " min. Results: finegray_project1_svedem.csv")

  } else if (.fg_mode == "inspect") {
    .fg_datasets <- load_prepared_data(required = "to1po2")
    run_inspection(.fg_datasets$to1po2)

  } else if (.fg_mode == "cif") {
    # Absolute risk for reviewer comment 4. Separate estimator, separate output
    # files; nothing here can change finegray_project1_*.csv.
    fg_log_file <- "cif_project1.log"
    install_crash_handler("fg_project1.R cif")

    log_msg("Starting cumulative incidence analysis in ", getwd())
    if (!is.null(.fg_only_outcomes)) {
      log_msg("Restricted to: ", paste(.fg_only_outcomes, collapse = ", "))
    }

    .fg_start <- Sys.time()
    .fg_datasets <- load_prepared_data(required = "to1po2")
    .cif_data <- attach_cif_groups(.fg_datasets$to1po2)
    rm(.fg_datasets); gc()

    log_msg("Index age for the attained-age axis: ", cif_config$index_age)
    log_msg("Exit-age rounding: ",
            if (is.null(fg_config$time_round_digits)) "none"
            else fg_config$time_round_digits,
            " (matching the published Fine-Gray run)")

    # Reconciliation: these counts must match finegray_project1_final.csv, and
    # the comparison says whether the rounding mattered at all.
    .cif_rounding <- compare_rounding_effect(
      .cif_data,
      if (is.null(.fg_only_outcomes)) main_specs else
        Filter(function(s) s$outcome_label %in% .fg_only_outcomes, main_specs),
      "Main")
    write.table(.cif_rounding, "cif_project1_rounding_check.csv",
                sep = ";", dec = ",", row.names = FALSE, quote = FALSE,
                fileEncoding = "Latin1")

    run_cif_block(
      data          = .cif_data,
      specs         = main_specs,
      dataset_label = "Main",
      only_outcomes = .fg_only_outcomes
    )

    log_msg("Done in ", round(elapsed_min(.fg_start), 2),
            " min. Results: cif_project1_curves.csv, cif_project1_summary.csv")

  } else {
    .fg_usage()
    quit(save = "no", status = 1)
  }
}
