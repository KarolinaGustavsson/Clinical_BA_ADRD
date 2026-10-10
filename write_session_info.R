# =============================================================================
# write_session_info.R
# =============================================================================
# Records the exact software versions used to produce the results.
# Run on the analysis platform after the pipeline completes, then paste the
# output into DEPENDENCIES.md.
#
#     source("write_session_info.R")
#
# Writes SESSION_INFO.txt to the repository root.
# =============================================================================

# Every package the pipeline attaches or calls with ::, taken from the scripts
# themselves. Keep this list in step with DEPENDENCIES.md.
pkgs <- c(
  # Biological age estimation (not on CRAN)
  "BioAge",
  # Data handling
  "data.table", "dplyr", "haven", "here", "stringr", "tidyr",
  # Modelling
  "broom", "flexsurv", "lmtest", "MASS", "mgcv", "sandwich",
  "survey", "survival",
  # Tables and figures
  "fields", "ggplot2", "ggtext", "gridExtra", "gtable", "htmlTable",
  "rmarkdown", "table1"
)

installed <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
missing   <- setdiff(pkgs, installed)

invisible(lapply(installed, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

versions <- data.frame(
  package = installed,
  version = vapply(installed, function(p) as.character(utils::packageVersion(p)), character(1)),
  row.names = NULL,
  stringsAsFactors = FALSE
)
versions <- versions[order(tolower(versions$package)), ]

con <- file("SESSION_INFO.txt", open = "wt")
on.exit(close(con))

writeLines(c(
  "=============================================================",
  " Clinical_BA_ADRD -- session information",
  paste0(" Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "=============================================================",
  ""
), con)

writeLines("--- Package versions (markdown table for DEPENDENCIES.md) ---", con)
writeLines("", con)
writeLines("| Package | Version |", con)
writeLines("|---|---|", con)
writeLines(sprintf("| %s | %s |", versions$package, versions$version), con)
writeLines("", con)

if (length(missing)) {
  writeLines("--- NOT INSTALLED (versions cannot be recorded) ---", con)
  writeLines(paste0("  ", missing), con)
  writeLines("", con)
  warning(
    "These packages are not installed and were skipped: ",
    paste(missing, collapse = ", ")
  )
}

writeLines("--- Full sessionInfo() ---", con)
writeLines("", con)
writeLines(capture.output(sessionInfo()), con)

message("Wrote SESSION_INFO.txt (", nrow(versions), " packages recorded).")
if (length(missing)) {
  message("Missing: ", paste(missing, collapse = ", "))
}
