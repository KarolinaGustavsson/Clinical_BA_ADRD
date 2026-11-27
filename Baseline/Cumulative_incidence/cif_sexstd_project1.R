#!/usr/bin/env Rscript
# =====================================================================
# Sex-standardised cumulative incidence (Reviewer 3, round 2, comment 4)
# =====================================================================
#
# The round-1 curves (fg_project1.R cif) compare the two biological age groups
# crudely. The biologically older group has more men (56.0% vs 48.1%), so part
# of the mortality difference may be sex composition. This script directly
# standardises the Aalen-Johansen curves to the sex distribution of the
# analysed population:
#
#   CIF_std(group, t) = sum over sex s of  w_s * CIF(group, sex s, t)
#
# where w_s is the share of sex s among everyone on the axis (both groups
# combined). Each sex-specific curve is the same estimator as before
# (cif_for_outcome(), delayed entry, death competing), so with one sex the
# result equals the crude curve exactly. Variance: sexes are independent
# samples, so se_std = sqrt(sum w_s^2 se_s^2).
#
# A second mode, "mv", standardises to sex, income quartile and education with
# a model-based (g-formula) approach instead - see mv_standardised_cif().
#
# Run on the VDI, in the folder holding fg_project1.R and the RData:
#
#   Rscript cif_sexstd_project1.R sex        # sex-standardised -> ./sexstd/
#   Rscript plot_cif_results.R sexstd        # figures + README into ./sexstd/
#
#   Rscript cif_sexstd_project1.R mv 0       # multivariable, no bootstrap -> ./mvstd/
#   Rscript cif_sexstd_project1.R mv 100     # ... with 100 bootstrap replicates (slow)
#   Rscript plot_cif_results.R mvstd
#
# Writes only into ./sexstd/ or ./mvstd/, so cif_project1_* and finegray_project1_* are
# untouched. Output is aggregate (risk sets < cif_config$min_risk_set in either
# sex are suppressed), so ./sexstd/ can be copied out of the VDI.
#
# The file name deliberately does not start with "fg_project1", so sourcing
# fg_project1.R below does not trigger its command-line entry point.

source("fg_project1.R")

#' Sex-standardised CIF for one outcome, grouping and axis.
#' Returns rows in the same shape as cif_for_outcome() (dementia and death only).
sex_standardised_cif <- function(df, group_var, axis,
                                 sex_var = "gender",
                                 index_age = cif_config$index_age,
                                 min_risk_set = cif_config$min_risk_set) {
  d_axis <- prepare_cif_axis(df, axis, index_age)
  d_axis <- d_axis[!is.na(d_axis[[group_var]]) & !is.na(d_axis[[sex_var]]), , drop = FALSE]
  w <- prop.table(table(as.character(d_axis[[sex_var]])))
  log_msg("    standard population (", axis, "): ",
          paste0(sex_var, " ", names(w), " = ", sprintf("%.4f", w), collapse = ", "))

  per_sex <- lapply(names(w), function(s) {
    sub <- df[!is.na(df[[sex_var]]) & as.character(df[[sex_var]]) == s, , drop = FALSE]
    r <- cif_for_outcome(sub, group_var, axis, index_age)
    r <- r[r$state %in% c("dementia", "death"), , drop = FALSE]
    r$w <- as.numeric(w[[s]])
    r$sex <- s
    r
  })
  r <- do.call(rbind, per_sex)

  parts <- split(r, list(r$group, r$state, r$time), drop = TRUE)
  out <- lapply(parts, function(x) {
    if (nrow(x) != length(w)) return(NULL)          # a sex missing at this time
    se2 <- ifelse(is.na(x$se), 0, x$se^2)            # NA only before the first event
    est <- sum(x$w * x$estimate)
    se  <- sqrt(sum(x$w^2 * se2))
    data.frame(
      group = x$group[1], state = x$state[1], time = x$time[1],
      estimate = est, se = se,
      CI_lower = max(0, est - 1.96 * se), CI_upper = min(1, est + 1.96 * se),
      n_risk = sum(x$n_risk), min_sex_risk = min(x$n_risk),
      n_event_cum = sum(x$n_event_cum),
      n_group = sum(x$n_group[!duplicated(x$sex)]),
      events_group = sum(x$events_group[!duplicated(x$sex)]),
      deaths_group = sum(x$deaths_group[!duplicated(x$sex)]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(out$group, out$state, out$time), , drop = FALSE]

  # Aggregate-only output: drop grid rows where either sex is thin.
  thin <- out$min_sex_risk < min_risk_set
  if (any(thin)) log_msg("    suppressing ", sum(thin), " grid rows with a sex-specific risk set < ",
                         min_risk_set)
  out <- out[!thin, , drop = FALSE]
  out$min_sex_risk <- NULL
  rownames(out) <- NULL
  attr(out, "weights") <- w
  out
}

#' Multivariable-standardised CIF (model-based, g-formula).
#'
#' Within each biological age group, cause-specific Cox models for dementia and
#' death (a multi-state coxph, same delayed entry and timescale as the crude
#' curves) are fitted with the multivariable covariates. Each group has its own
#' baseline hazards, so the curve shapes are not forced to be proportional
#' between groups; only the covariate effects are assumed proportional. The
#' predicted CIF is averaged over the covariate distribution of everyone on the
#' axis (both groups combined), so both groups are compared at the same
#' sex, income and education mix. CIs are bootstrap percentiles (resampling
#' individuals within group); n_boot = 0 gives point estimates only.
mv_standardised_cif <- function(df, group_var, axis,
                                covars = adjustment_sets[["Multivariable model"]],
                                index_age = cif_config$index_age,
                                min_risk_set = cif_config$min_risk_set,
                                n_boot = 0, seed = 2026) {
  d <- prepare_cif_axis(df, axis, index_age)
  d <- d[!is.na(d[[group_var]]) & complete.cases(d[, covars, drop = FALSE]), , drop = FALSE]
  if (nrow(d) == 0) stop("No one at risk on the '", axis, "' axis after index age ", index_age)
  d$cif_grp <- droplevels(factor(d[[group_var]]))
  grid <- cif_grid(axis, index_age)

  # Standard population: observed covariate patterns, both groups combined.
  pat <- aggregate(list(n = rep(1L, nrow(d))), by = d[, covars, drop = FALSE], FUN = sum)
  w   <- pat$n / sum(pat$n)
  newdata <- pat[, covars, drop = FALSE]
  log_msg("    standard population (", axis, "): ", nrow(pat), " covariate patterns, N = ", sum(pat$n))

  form <- stats::as.formula(paste("Surv(cif_start, cif_stop, fg_status) ~",
                                  paste(covars, collapse = " + ")))

  # Standardised CIF for one group's data -> matrix [grid, c(dementia, death)]
  std_curve <- function(sub) {
    fit <- survival::coxph(form, data = sub, id = fg_id, model = TRUE)
    sf  <- survival::survfit(fit, newdata = newdata)
    j   <- findInterval(grid, sf$time)
    out <- matrix(0, nrow = length(grid), ncol = 2, dimnames = list(NULL, c("dementia", "death")))
    for (k in 1:2) {
      st  <- c("event", "death")[k]
      ps  <- sf$pstate[, , match(st, sf$states), drop = FALSE][, , 1]
      if (is.null(dim(ps))) ps <- matrix(ps, ncol = nrow(newdata))
      avg <- as.vector(ps %*% w)                     # average over patterns
      out[, k] <- ifelse(j == 0, 0, avg[pmax(j, 1L)])
    }
    out
  }

  groups <- levels(d$cif_grp)
  est <- lapply(groups, function(g) std_curve(d[d$cif_grp == g, , drop = FALSE]))
  names(est) <- groups

  boot <- NULL
  if (n_boot > 0) {
    set.seed(seed)
    log_msg("    bootstrap: ", n_boot, " replicates")
    boot <- lapply(groups, function(g) {
      sub <- d[d$cif_grp == g, , drop = FALSE]
      reps <- vapply(seq_len(n_boot), function(b) {
        bs <- sub[sample.int(nrow(sub), replace = TRUE), , drop = FALSE]
        bs$fg_id <- seq_len(nrow(bs))
        as.vector(std_curve(bs))
      }, numeric(length(grid) * 2))
      reps
    })
    names(boot) <- groups
  }

  rows <- list()
  for (g in groups) {
    sub <- d[d$cif_grp == g, , drop = FALSE]
    n_risk <- vapply(grid, function(t) sum(sub$cif_start < t & sub$cif_stop >= t), numeric(1))
    for (k in 1:2) {
      st <- c("dementia", "death")[k]
      ev <- if (st == "dementia") "event" else "death"
      idx <- (k - 1) * length(grid) + seq_along(grid)
      se <- lo <- hi <- rep(NA_real_, length(grid))
      if (!is.null(boot)) {
        b  <- boot[[g]][idx, , drop = FALSE]
        se <- apply(b, 1, stats::sd)
        lo <- apply(b, 1, stats::quantile, 0.025)
        hi <- apply(b, 1, stats::quantile, 0.975)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        group = g, state = st, time = grid,
        estimate = est[[g]][, k], se = se, CI_lower = lo, CI_upper = hi,
        n_risk = n_risk,
        n_event_cum = vapply(grid, function(t) sum(sub$fg_status == ev & sub$cif_stop <= t), numeric(1)),
        n_group = nrow(sub), events_group = sum(sub$fg_status == "event"),
        deaths_group = sum(sub$fg_status == "death"),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  thin <- out$n_risk < min_risk_set
  if (any(thin)) log_msg("    suppressing ", sum(thin), " grid rows with risk set < ", min_risk_set)
  out[!thin, , drop = FALSE]
}

run_sexstd_block <- function(data, specs, dataset_label,
                             out_dir = "sexstd",
                             groupings = c("kdm_group", "phenoage_group"),
                             axes = c("age", "years"),
                             method = c("sex", "mv"), n_boot = 0) {
  method <- match.arg(method)
  dir.create(out_dir, showWarnings = FALSE)
  curves_file  <- file.path(out_dir, "cif_project1_curves.csv")
  summary_file <- file.path(out_dir, "cif_project1_summary.csv")
  if (file.exists(curves_file) || file.exists(summary_file)) {
    stop("Output already exists in ", out_dir, " - delete it before re-running.")
  }
  summaries <- list()
  for (spec in specs) {
    log_msg("=== ", dataset_label, " / ", spec$outcome_label, " (sex-standardised) ===")
    df <- clean_fg_data(data, spec$event_var, dataset_label, extra_vars = groupings)
    counts <- fg_counts(df)
    log_msg("N=", counts$N, "  events=", counts$events,
            "  competing deaths=", counts$competing_deaths, "  censored=", counts$censored)
    for (grouping in groupings) for (axis in axes) {
      log_msg("  ", grouping, " / ", axis)
      rows <- tryCatch(
        if (method == "sex") sex_standardised_cif(df, grouping, axis)
        else mv_standardised_cif(df, grouping, axis, n_boot = n_boot),
        error = function(e) e)
      if (inherits(rows, "error")) { log_msg("    FAILED: ", conditionMessage(rows)); next }

      report_times <- if (identical(axis, "age")) cif_config$report_ages else cif_config$report_years
      diffs <- cif_differences(rows, axis, report_times)
      if (!is.null(diffs)) {
        diffs$dataset <- dataset_label; diffs$outcome <- spec$outcome_label
        diffs$grouping <- grouping
        summaries[[length(summaries) + 1L]] <- diffs
      }
      rows$dataset <- dataset_label; rows$outcome <- spec$outcome_label
      rows$grouping <- grouping; rows$axis <- axis
      append_cif_rows(rows, curves_file)
    }
    rm(df); gc()
  }
  if (length(summaries) > 0) {
    all_diffs <- do.call(rbind, summaries)
    cols <- c("dataset", "outcome", "grouping", "axis", "state", "time",
              "est_younger", "est_older", "difference", "diff_lower", "diff_upper",
              "n_risk_younger", "n_risk_older")
    write.table(all_diffs[, cols, drop = FALSE], summary_file, sep = ";", dec = ",",
                row.names = FALSE, quote = FALSE, fileEncoding = "Latin1")
    log_msg("Wrote ", nrow(all_diffs), " summary rows to ", summary_file)
  }
  invisible(curves_file)
}

# ---------------------------------------------------------------------
# Command-line entry point (only when run as a script, not when sourced)
# ---------------------------------------------------------------------
.sx_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(.sx_file) > 0 && grepl("^cif_sexstd_project1\\.R$", basename(.sx_file))) {
  .sx_args   <- commandArgs(trailingOnly = TRUE)
  .sx_method <- if (length(.sx_args) >= 1) .sx_args[1] else "sex"
  .sx_boot   <- if (length(.sx_args) >= 2) as.integer(.sx_args[2]) else 0L
  if (!.sx_method %in% c("sex", "mv")) stop("First argument must be 'sex' or 'mv'")
  .sx_dir    <- if (.sx_method == "sex") "sexstd" else "mvstd"

  fg_log_file <- paste0("cif_", .sx_dir, "_project1.log")
  install_crash_handler("cif_sexstd_project1.R")
  .sx_start <- Sys.time()
  log_msg(if (.sx_method == "sex") "Sex-standardised" else "Multivariable-standardised",
          " cumulative incidence in ", getwd(),
          if (.sx_method == "mv") paste0(" (bootstrap replicates: ", .sx_boot, ")") else "")
  .sx_data <- load_prepared_data(required = "to1po2")
  .sx_cif  <- attach_cif_groups(.sx_data$to1po2)
  rm(.sx_data); gc()
  # ADRD only, as in the published figure. Add outcomes here if needed.
  run_sexstd_block(.sx_cif, Filter(function(s) s$outcome_label == "ADRD", main_specs), "Main",
                   out_dir = .sx_dir, method = .sx_method, n_boot = .sx_boot)
  log_msg("Done in ", round(elapsed_min(.sx_start), 2),
          " min. Now run: Rscript plot_cif_results.R ", .sx_dir)
}
