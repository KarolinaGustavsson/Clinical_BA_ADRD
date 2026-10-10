# =============================================================================
# config.R  --  central configuration for the Clinical_BA_ADRD pipeline
# =============================================================================
#
# This is the only file that needs editing to run the pipeline in a new
# environment. Every analysis script sources it and takes paths, output
# locations and analysis parameters from it instead of hardcoding them.
#
# Usage, at the top of every analysis script:
#
#     source(here::here("config.R"))
#     amoris <- haven::read_sas(PATH_STUDYPOP)
#
# Access to the underlying AMORIS data is restricted; see README.md.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Where the restricted data live
# -----------------------------------------------------------------------------
# Point this at the directory holding the three SAS extracts. Either set the
# environment variable AMORIS_DATA_ROOT (preferred -- it keeps the local path
# out of the scripts and out of version control), for example in ~/.Renviron:
#
#     AMORIS_DATA_ROOT=/path/to/amoris
#
# or edit the fallback below. The analyses for this paper were run on a secure
# analysis platform at Karolinska Institutet where this was a mapped network
# drive; the path itself is site-specific and carries no meaning elsewhere.

DATA_ROOT <- path.expand(
  Sys.getenv("AMORIS_DATA_ROOT", unset = "~/amoris_data")
)


# -----------------------------------------------------------------------------
# 2. Source data files (SAS format, restricted access)
# -----------------------------------------------------------------------------

# Baseline study population: one row per individual, at the first time point
# where all candidate biomarkers were available, with register linkages
# except SveDem.
PATH_STUDYPOP <- file.path(DATA_ROOT, "StudyPop.sas7bdat")

# Full AMORIS extract with every measurement occasion. One row per
# measurement; the source of the repeated-measures subsample.
PATH_AMORIS_ALL <- file.path(DATA_ROOT, "AMORIS_20240417.sas7bdat")

# SveDem (Swedish registry for cognitive/dementia disorders) linkage, used for
# the clinically adjudicated subtype analysis.
PATH_SVEDEM <- file.path(DATA_ROOT, "SveDem.sas7bdat")


# -----------------------------------------------------------------------------
# 3. Output locations
# -----------------------------------------------------------------------------
# Tables, intermediate CSVs and figures are written here. Override with
# STUDY1_OUTPUT_DIR if the data platform keeps outputs on a separate volume.
# Nothing is created when this file is sourced; out_path() and fig_path()
# create the directory the first time something is actually written.

# Anchored to the repository root, NOT to the working directory. RStudio runs
# a notebook with the working directory set to that notebook's own folder, so
# a working-directory default would send Baseline/ output to Baseline/output
# and Longitudinal/ output to Longitudinal/output -- and the files the
# longitudinal scripts read back from the baseline run would not be found.

PROJECT_ROOT <- local({
  # here::here() locates the root from any subdirectory, using the .here
  # marker file shipped with this repository.
  if (requireNamespace("here", quietly = TRUE)) {
    root <- tryCatch(here::here(), error = function(e) NULL)
    if (!is.null(root) && file.exists(file.path(root, "config.R"))) return(root)
  }
  # Fallback if `here` is not installed: walk up from the working directory
  # looking for the directory that holds config.R.
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "config.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  getwd()
})

OUTPUT_DIR <- path.expand(
  Sys.getenv("STUDY1_OUTPUT_DIR", unset = file.path(PROJECT_ROOT, "output"))
)
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

# Build a path under OUTPUT_DIR, creating the directory if it does not exist.
out_path <- function(...) {
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  file.path(OUTPUT_DIR, ...)
}

# Same, for figures.
fig_path <- function(...) {
  if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)
  file.path(FIGURE_DIR, ...)
}


# -----------------------------------------------------------------------------
# 4. Intermediate files passed between scripts
# -----------------------------------------------------------------------------
# Written by the Baseline pipeline and read back by the Longitudinal pipeline.
# Column definitions and which script writes which file are in
# DATA_DICTIONARY.md.
#
# These are plain paths; the directory is created on first write.

PATH_SAMPLE_IDS   <- file.path(OUTPUT_DIR, "sampleIDs_baseline.csv")
PATH_TRAINING_SET <- file.path(OUTPUT_DIR, "training_set_baseline.csv")
PATH_BIOAGE_SCALE <- file.path(OUTPUT_DIR, "bioage_scale_df_table.csv")

# Biologically younger / older grouping, one row per individual. Written by
# Table1_B.Rmd, read by Table1_L.Rmd and by the cumulative incidence analysis.
# Individual-level: keep it inside the secure environment, never in the
# repository.
PATH_TABLE1_PERID <- file.path(OUTPUT_DIR, "table1_kdmgroups_baseline_perID.csv")


# -----------------------------------------------------------------------------
# 5. Analysis parameters
# -----------------------------------------------------------------------------
# Values that define the analysis, stated once here so a reader does not have
# to recover them from the scripts.

# The twelve biomarkers retained after screening on absolute correlation with
# chronological age, computed separately in men and women. Selection is shown
# in Selecting_biomarkers_B.Rmd, which also lists the candidates dropped.
BIOMARKERS <- c(
  "S_Alb",    # Serum albumin
  "S_Krea",   # Serum creatinine
  "TC",       # Total cholesterol
  "TG",       # Triglycerides
  "S_K",      # Serum potassium
  "S_Urea",   # Serum urea
  "S_LD",     # Lactate dehydrogenase
  "S_Ca",     # Serum calcium
  "S_Urat",   # Uric acid (urate)
  "S_FAMN",   # Fructosamine
  "S_Hapt",   # Haptoglobin
  "fS_Gluk"   # Fasting serum glucose
)

BIOMARKER_R_THRESHOLD <- 0.1    # Screening threshold on |r| with chronological age
OUTLIER_SD            <- 5      # Biomarker values beyond +/- 5 SD are excluded
TRAINING_SET_SIZE     <- 10080  # Matches the NHANES III training set used by
                                # BioAge, while keeping a uniform age
                                # distribution within each sex
MIN_MEASUREMENTS      <- 3      # Minimum pre-diagnosis measurement occasions
                                # for the repeated-measures subsample
AGE_BIN_WIDTH         <- 5      # Width of the trajectory age bins, years
AGE_BIN_RANGE         <- c(17, 82)        # Age range covered by the bins
FOLLOWUP_END          <- as.Date("2020-12-31")  # Last date of register linkage
SEED                  <- 20250417         # Training/test split and subsampling


# -----------------------------------------------------------------------------
# 6. Output formatting
# -----------------------------------------------------------------------------
# CSVs are written with semicolon delimiters, comma decimals and Latin-1
# encoding, so they open directly in a Swedish-locale spreadsheet on the
# analysis platform. Anything read back by a later script uses the same
# settings.

CSV_SEP      <- ";"
CSV_DEC      <- ","
CSV_ENCODING <- "latin1"


# -----------------------------------------------------------------------------
# 7. Optional: platform without direct CRAN access
# -----------------------------------------------------------------------------
# A secure analysis platform may have no outbound internet, with packages
# installed to a project library from an internal mirror. Set these through
# the environment if that applies; both are unset by default and ordinary
# installation is used.
#
#     R_PROJECT_LIB=/path/to/project/library
#     R_CRAN_MIRROR=https://internal.mirror.example/cran

LOCAL_LIB   <- Sys.getenv("R_PROJECT_LIB",  unset = "")
CRAN_MIRROR <- Sys.getenv("R_CRAN_MIRROR",  unset = "")

if (nzchar(LOCAL_LIB) && dir.exists(LOCAL_LIB)) .libPaths(LOCAL_LIB)
if (nzchar(CRAN_MIRROR)) options(repos = CRAN_MIRROR)


# -----------------------------------------------------------------------------
# 8. Sanity check
# -----------------------------------------------------------------------------
# Call check_config() after sourcing, to confirm the inputs are reachable
# before a long run starts.

check_config <- function() {
  paths <- c(
    StudyPop   = PATH_STUDYPOP,
    AMORIS_all = PATH_AMORIS_ALL,
    SveDem     = PATH_SVEDEM
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    warning(
      "Configured input files not found:\n",
      paste0("  ", names(missing), ": ", missing, collapse = "\n"),
      "\nSet AMORIS_DATA_ROOT, or edit DATA_ROOT in config.R.",
      call. = FALSE
    )
  } else {
    message("config.R: all source data files found under ", DATA_ROOT)
  }
  invisible(length(missing) == 0)
}
