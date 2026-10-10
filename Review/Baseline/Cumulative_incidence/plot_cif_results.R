#!/usr/bin/env Rscript
# =====================================================================
# Publication figures from the exported CIF results
# =====================================================================
#
# Reads cif_project1_curves.csv only. That file is aggregate, so this runs
# anywhere - inside the VDI or outside it - and needs no data access.
#
#   Rscript plot_cif_results.R [directory containing cif_project1_curves.csv]
#
# Writes, into the same directory:
#   figure_cif_ADRD_<grouping>_<axis>.pdf   the figures
#   README_cif_figures.md                   what each one shows, with the
#                                           headings, units and numbers
#
# Figures carry no title and no numbers at risk: everything a caption needs is
# generated into the README instead, so the figure stays clean and the wording
# lives somewhere editable. Legends sit under the x axis.
#
# Subtypes are deliberately not plotted. Their competing event is death OR
# another dementia subtype, which needs its own explanation and is more detail
# than this reviewer answer calls for. Change cif_plot_config$outcomes to
# reinstate them.
#
# base graphics only, same as the rest of the pipeline.

cif_plot_config <- list(
  outcomes = "ADRD",

  # Stop the axis once the smaller group drops below this many at risk. Beyond
  # that the curve is driven by a handful of people and reads as structure that
  # is not there.
  min_risk = 1000,

  outcome_labels  = c(ADRD = "ADRD"),
  grouping_labels = c(kdm_group      = "KDM biological age",
                      phenoage_group = "PhenoAge biological age"),
  axis_labels     = c(age = "Age", years = "Time since baseline"),
  axis_units      = c(age = "years", years = "years"),

  group_colours = c("Biologically younger" = "#2166ac",
                    "Biologically older"   = "#b2182b"),
  state_lty     = c(dementia = "solid", death = "12"),
  state_lwd     = c(dementia = 2.2,     death = 2.8),
  state_labels  = c(dementia = "ADRD",  death = "Death"),

  # Inches - a single-column supplementary figure.
  width = 7, height = 5, pointsize = 11
)

read_curves <- function(dir = ".") {
  f <- file.path(dir, "cif_project1_curves.csv")
  if (!file.exists(f)) stop("Cannot find ", f)
  read.table(f, sep = ";", dec = ",", header = TRUE, quote = "",
             fileEncoding = "Latin1", stringsAsFactors = FALSE)
}

read_summary <- function(dir = ".") {
  f <- file.path(dir, "cif_project1_summary.csv")
  if (!file.exists(f)) return(NULL)
  read.table(f, sep = ";", dec = ",", header = TRUE, quote = "",
             fileEncoding = "Latin1", stringsAsFactors = FALSE)
}

#' Turn (x, y) into the vertices of a step function.
step_xy <- function(x, y) {
  n <- length(x)
  if (n < 2) return(list(x = x, y = y))
  list(x = rep(x, each = 2)[-1],
       y = rep(y, each = 2)[-(2 * n)])
}

#' A confidence band that follows the steps rather than cutting the corners.
step_band <- function(x, lo, hi, col) {
  ok <- is.finite(lo) & is.finite(hi)
  if (sum(ok) < 2) return(invisible(NULL))
  up   <- step_xy(x[ok], hi[ok])
  down <- step_xy(x[ok], lo[ok])
  polygon(c(up$x, rev(down$x)), c(up$y, rev(down$y)), col = col, border = NA)
}

#' Trim the sparse tail: keep times where BOTH groups still have min_risk.
trim_to_risk_set <- function(d, min_risk) {
  risk_by_time <- aggregate(n_risk ~ time, data = d, FUN = min)
  keep <- risk_by_time$time[risk_by_time$n_risk >= min_risk]
  if (length(keep) < 2) {
    warning("Fewer than two time points survive min_risk = ", min_risk,
            "; keeping everything")
    keep <- sort(unique(d$time))
  }
  d[d$time %in% keep, ]
}

#' One figure: cumulative incidence of ADRD and of death, by group.
#'
#' No title and no numbers at risk by design - see the README this script
#' writes. Both legends sit below the x axis, so nothing overlaps the curves.
plot_cif_figure <- function(curves, outcome, grouping, axis,
                            file = NULL, ci = TRUE,
                            states = c("dementia", "death"),
                            min_risk = cif_plot_config$min_risk,
                            cfg = cif_plot_config) {

  d <- curves[curves$outcome == outcome & curves$grouping == grouping &
              curves$axis == axis & curves$state %in% states, ]
  if (nrow(d) == 0) stop("No rows for ", outcome, " / ", grouping, " / ", axis)

  groups <- intersect(names(cfg$group_colours), unique(d$group))
  d      <- trim_to_risk_set(d, min_risk)
  xlim   <- range(d$time)
  ymax   <- max(d$estimate[is.finite(d$estimate)], na.rm = TRUE)

  if (!is.null(file)) {
    grDevices::pdf(file, width = cfg$width, height = cfg$height,
                   pointsize = cfg$pointsize)
    # Only the device close is registered here. Restoring par() afterwards
    # would run AFTER dev.off() and force R to open a fresh device just to
    # receive the settings - which is where a stray Rplots.pdf comes from.
    on.exit(grDevices::dev.off(), add = TRUE)
  } else {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
  }

  # No top margin to speak of: there is no title. The bottom carries the x
  # label and two legend rows.
  graphics::par(mar = c(8.2, 4.4, 1.0, 1.0), lend = "round",
                mgp = c(2.3, 0.7, 0))

  plot(NA, xlim = xlim, ylim = c(0, ymax * 1.04),
       xlab = "", ylab = "Cumulative incidence", yaxt = "n", bty = "n",
       xaxs = "i", yaxs = "i")
  at <- pretty(c(0, ymax))
  graphics::axis(2, at = at, labels = paste0(round(100 * at), "%"), las = 1)
  graphics::abline(h = at, col = "grey92", lwd = 1)
  title(xlab = cfg$axis_labels[[axis]], line = 2.3)

  for (st in states) {
    for (g in groups) {
      s <- d[d$group == g & d$state == st, ]
      s <- s[order(s$time), ]
      if (nrow(s) == 0) next
      if (ci) {
        step_band(s$time, s$CI_lower, s$CI_upper,
                  col = grDevices::adjustcolor(cfg$group_colours[[g]],
                                               alpha.f = 0.16))
      }
      xy <- step_xy(s$time, s$estimate)
      graphics::lines(xy$x, xy$y, col = cfg$group_colours[[g]],
                      lty = cfg$state_lty[[st]], lwd = cfg$state_lwd[[st]])
    }
  }

  # Two horizontal legend rows under the axis. xpd = NA lets them sit in the
  # margin; inset is negative to push them below the plotting region.
  graphics::legend("bottom", inset = c(0, -0.30), xpd = NA, horiz = TRUE,
                   legend = groups, col = cfg$group_colours[groups],
                   lwd = 2.4, bty = "n", cex = 0.92, seg.len = 2.0,
                   text.width = NA)
  if (length(states) > 1) {
    graphics::legend("bottom", inset = c(0, -0.39), xpd = NA, horiz = TRUE,
                     legend = unname(cfg$state_labels[states]), col = "grey35",
                     lty = unname(cfg$state_lty[states]),
                     lwd = unname(cfg$state_lwd[states]),
                     bty = "n", cex = 0.92, seg.len = 3.0, text.width = NA)
  }

  invisible(file)
}

#' How much would ignoring the competing risk overstate absolute risk?
#'
#' The naive 1 - KM curve treats death as censoring. It is not plotted - it is
#' a methods point, not a result - but the number belongs in the response.
naive_overstatement <- function(curves, outcome, grouping, axis, at_time) {
  d  <- curves[curves$outcome == outcome & curves$grouping == grouping &
               curves$axis == axis & curves$time == at_time, ]
  aj <- d[d$state == "dementia", c("group", "estimate")]
  km <- d[d$state == "naive_km", c("group", "estimate")]
  m  <- merge(aj, km, by = "group", suffixes = c("_aj", "_naive"))
  m$absolute_overstatement <- m$estimate_naive - m$estimate_aj
  m$relative_overstatement <- m$estimate_naive / m$estimate_aj - 1
  m$time <- at_time
  m
}

# ---------------------------------------------------------------------
# The README that carries what the figures deliberately do not
# ---------------------------------------------------------------------

pct <- function(x, dp = 1) sprintf(paste0("%.", dp, "f%%"), 100 * x)

write_figure_readme <- function(dir, curves, made, cfg = cif_plot_config) {
  smry <- read_summary(dir)
  L <- c(
    "# Cumulative incidence figures",
    "",
    sprintf("Generated %s by `plot_cif_results.R` from `cif_project1_curves.csv`.",
            format(Sys.Date())),
    "",
    "## What these show",
    "",
    "Absolute risk of dementia, and of death, estimated with the",
    "**Aalen-Johansen estimator** with death as a competing event. This is the",
    "absolute-risk counterpart of the Fine-Gray subdistribution hazard ratios in",
    "Supplementary Table 9, which are relative measures and so do not answer the",
    "reviewer's request on their own.",
    "",
    "Same analytic sample and same event coding as the Fine-Gray models: they",
    "share `clean_fg_data()`, including the exit-age rounding, so the curves sit",
    "on exactly the timeline behind Supplementary Table 9. `finegray()` is never",
    "called, so none of this changes the published hazard ratios.",
    "",
    "## Unit of the vertical axis",
    "",
    "**Cumulative incidence is a probability**, plotted as a percentage. At a",
    "given age it is the probability that a person has been diagnosed with ADRD",
    "(or has died) by that age, among those alive and dementia-free at age 65.",
    "It is not a rate and has no time unit.",
    "",
    "Because death is treated as a competing event rather than as censoring, the",
    "two curves for a group are on a common scale and can be added: together with",
    "the probability of remaining alive and dementia-free they sum to 1. A",
    "Kaplan-Meier style `1 - S` curve does not have this property and overstates",
    "absolute risk.",
    "",
    "## Unit of the horizontal axis",
    "",
    "**Age, in years** - attained age, the same timescale as the Cox and",
    "Fine-Gray models. Curves start at age 65: with delayed entry a cumulative",
    "incidence 'from time zero' is undefined, so it is anchored at an index age",
    "and reads as the incidence accumulating between age 65 and the age shown.",
    "The `_years` figures use time since baseline in years instead.",
    "",
    "The axis is cut where the smaller group falls below",
    sprintf("%s people still at risk, so no part of a curve rests on a handful",
            format(cfg$min_risk, big.mark = " ")),
    "of individuals.",
    "",
    "Shaded bands are 95% confidence intervals.",
    "",
    "## Suggested caption",
    "",
    "> Cumulative incidence of ADRD and of death by biological age group, with",
    "> death as a competing event (Aalen-Johansen estimator). Age in years;",
    "> curves are anchored at age 65 and shaded bands are 95% confidence",
    "> intervals. Biological age groups are the median split of the KDM age gap",
    "> used in Table 1.",
    "",
    "## Figures",
    ""
  )

  for (i in seq_len(nrow(made))) {
    L <- c(L, sprintf("- `%s` - %s, by %s, on the %s axis.",
                      basename(made$file[i]),
                      cfg$outcome_labels[[made$outcome[i]]],
                      cfg$grouping_labels[[made$grouping[i]]],
                      if (made$axis[i] == "age") "attained age" else
                        "time since baseline"))
  }

  if (!is.null(smry)) {
    s <- smry[smry$outcome == "ADRD" & smry$grouping == "kdm_group" &
              smry$axis == "age", ]
    if (nrow(s) > 0) {
      L <- c(L, "", "## Key numbers (ADRD, KDM groups, attained age)", "",
             "Younger and older are the Table 1 median split of the KDM age gap.",
             "",
             "| Event | Age | Younger | Older | Difference (95% CI) |",
             "|---|---|---|---|---|")
      s <- s[order(s$state, s$time), ]
      for (i in seq_len(nrow(s))) {
        L <- c(L, sprintf("| %s | %d | %s | %s | %s (%s to %s) |",
                          ifelse(s$state[i] == "dementia", "ADRD", "Death"),
                          round(s$time[i]),
                          pct(s$est_younger[i]), pct(s$est_older[i]),
                          pct(s$difference[i], 2),
                          pct(s$diff_lower[i], 2), pct(s$diff_upper[i], 2)))
      }
    }
  }

  ov <- tryCatch(naive_overstatement(curves, "ADRD", "kdm_group", "age", 90),
                 error = function(e) NULL)
  if (!is.null(ov) && nrow(ov) > 0) {
    L <- c(L, "",
           "## Why the competing risk has to be handled", "",
           "Treating death as censoring (a naive `1 - Kaplan-Meier` curve)",
           "overstates the absolute risk of ADRD at age 90 by:", "")
    for (i in seq_len(nrow(ov))) {
      L <- c(L, sprintf("- %s: %s naive vs %s correct, overstated by %s (%s relative)",
                        ov$group[i], pct(ov$estimate_naive[i]),
                        pct(ov$estimate_aj[i]),
                        pct(ov$absolute_overstatement[i]),
                        sprintf("%.0f%%", 100 * ov$relative_overstatement[i])))
    }
    L <- c(L, "",
           "This is a methods point rather than a result, so it is reported here",
           "and in the response letter rather than drawn as a panel.")
  }

  L <- c(L, "",
         "## Not shown", "",
         "Dementia subtypes. Their competing event is death **or another",
         "dementia subtype**, because the registers record only one type per",
         "person - a different quantity needing its own explanation, and more",
         "detail than this answer calls for. Set `cif_plot_config$outcomes` to",
         "reinstate them.", "")

  f <- file.path(dir, "README_cif_figures.md")
  writeLines(L, f)
  f
}

# ---------------------------------------------------------------------

if (!interactive()) {
  args   <- commandArgs(trailingOnly = TRUE)
  dir    <- if (length(args) >= 1) args[1] else "."
  curves <- read_curves(dir)

  combos <- unique(curves[curves$outcome %in% cif_plot_config$outcomes,
                          c("outcome", "grouping", "axis")])
  if (nrow(combos) == 0) {
    stop("No rows for outcome(s): ",
         paste(cif_plot_config$outcomes, collapse = ", "))
  }

  made <- NULL
  for (i in seq_len(nrow(combos))) {
    f <- file.path(dir, sprintf("figure_cif_%s_%s_%s.pdf",
                                combos$outcome[i], combos$grouping[i],
                                combos$axis[i]))
    ok <- tryCatch({
      plot_cif_figure(curves, combos$outcome[i], combos$grouping[i],
                      combos$axis[i], file = f)
      TRUE
    }, error = function(e) {
      cat("FAILED", f, ":", conditionMessage(e), "\n"); FALSE
    })
    if (ok) {
      cat("wrote", f, "\n")
      made <- rbind(made, data.frame(file = f, outcome = combos$outcome[i],
                                     grouping = combos$grouping[i],
                                     axis = combos$axis[i],
                                     stringsAsFactors = FALSE))
    }
  }

  if (!is.null(made)) cat("wrote", write_figure_readme(dir, curves, made), "\n")
}
