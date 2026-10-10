#!/usr/bin/env Rscript
# =====================================================================
# Fine-Gray competing risk analysis for project 1 (dementia)
# =====================================================================
#
# One file, three modes. Config, timeline derivation, cleaning, expansion,
# fitting and resume logic all live here; the bottom of the file just picks
# which dataset/diagnostics to run based on the first command-line argument.
#
#   Rscript fg_project1_ADRD.R main   [outcome ...]   # ADRD Alzheimers Vascular Other
#   Rscript fg_project1_ADRD.R svedem [outcome ...]   # S_ADRD S_e_alz S_l_alz S_alz_vas
#                                                 # S_vas S_lew S_fro S_par S_uns S_oth
#   Rscript fg_project1_ADRD.R inspect                # data diagnostics; run this first
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
# Modelling choices, all deliberate - see README_fg_project1_ADRD.md:
#   * income_bins and NewUtbNiva stay NUMERIC, matching coxph in CoxPH_B.Rmd.
#   * Other dementia subtypes are competing events, alongside death.
#   * One finegray() expansion per outcome serves all 4 models.
#
# This file is also source()'d by test_fg_project1_ADRD.R to reuse the functions
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
  # Searched in order in the working directory.
  rdata_candidates = c(
    "coxph_B_correctsubtypes.RData",
    "coxph_B_correctsubtypes.rdata",
    "coxph_B_correctsubtypes.RDATA",
    "coxph_B_correctsubtypes"
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
  # 3 ~ one day, 2 ~ four days, 1 ~ five weeks. Only touch this if the exact
  # run runs out of memory, and report it as a sensitivity analysis.
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
clean_fg_data <- function(data, event_var, dataset_label,
                          ba_measures = fg_config$ba_measures,
                          covariates = all_model_covariates,
                          competing_other_dementia = fg_config$competing_other_dementia,
                          time_round_digits = fg_config$time_round_digits,
                          quiet = TRUE) {

  require_columns(data, unique(c(event_var, ba_measures, covariates)),
                  context = paste0(dataset_label, " / ", event_var))

  timeline <- derive_fg_timeline(data, quiet = quiet)

  df <- data.frame(
    timeline,
    fg_outcome = data[[event_var]],
    data[, unique(c(ba_measures, covariates)), drop = FALSE],
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

  df[, c("AgeAtFirstDate", "AgeAtLastDate", "fg_status", "fg_id",
         ba_measures, covariates), drop = FALSE]
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

# =======================================================================
# Command-line entry point
# =======================================================================
#
# Only runs when this file is executed directly via Rscript (checked by
# looking for "fg_project1_ADRD.R" in the --file= argument Rscript passes to
# itself). When test_fg_project1_ADRD.R does source("fg_project1_ADRD.R") to reuse the
# functions above, --file= points at test_fg_project1_ADRD.R instead, so this
# block is skipped and no analysis runs as a side effect of testing.

.fg_cli_raw  <- commandArgs(trailingOnly = FALSE)
.fg_file_arg <- sub("^--file=", "", .fg_cli_raw[grepl("^--file=", .fg_cli_raw)])
.fg_invoked_directly <- length(.fg_file_arg) > 0 &&
  identical(sub(".*[\\\\/]", "", .fg_file_arg), "fg_project1_ADRD.R")

if (.fg_invoked_directly) {
  options(warn = 1)

  .fg_args <- commandArgs(trailingOnly = TRUE)
  .fg_mode <- if (length(.fg_args) >= 1) .fg_args[1] else ""
  .fg_only_outcomes <- if (length(.fg_args) >= 2) .fg_args[-1] else NULL

  .fg_usage <- function() {
    cat("Usage:\n")
    cat("  Rscript fg_project1_ADRD.R main   [outcome ...]\n")
    cat("      outcomes: ADRD Alzheimers Vascular Other\n")
    cat("  Rscript fg_project1_ADRD.R svedem [outcome ...]\n")
    cat("      outcomes: S_ADRD S_e_alz S_l_alz S_alz_vas S_vas S_lew S_fro S_par S_uns S_oth\n")
    cat("  Rscript fg_project1_ADRD.R inspect\n")
    cat("      data diagnostics - run this once before trusting main/svedem output\n")
  }

  if (.fg_mode == "main") {
    # Same behaviour as the old fg_project1_ADRD_main.R: results appended to
    # finegray_project1_main.csv one row at a time; completed models are
    # skipped if this is a restart after a crash.
    fg_log_file <- "finegray_project1_main.log"
    install_crash_handler("fg_project1_ADRD.R main")

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
    # Same behaviour as the old fg_project1_ADRD_svedem.R: results appended to
    # finegray_project1_svedem.csv one row at a time.
    fg_log_file <- "finegray_project1_svedem.log"
    install_crash_handler("fg_project1_ADRD.R svedem")

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

  } else {
    .fg_usage()
    quit(save = "no", status = 1)
  }
}
