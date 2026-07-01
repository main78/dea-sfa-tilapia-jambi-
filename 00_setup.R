# ============================================================
# 00_setup.R  —  GLOBAL CONFIGURATION FOR J2 REPLICATION SCRIPTS
# ============================================================
# Replication code for:
#   "When and Why Do DEA and SFA Diverge? Evidence from
#    Spatially Heterogeneous Tilapia Aquaculture"
#   Journal of Productivity Analysis (submitted)
#
# Author : Ma'in (NIM P3C123010), FEB UNJA
#
# USAGE:
#   1. Clone or download this repository.
#   2. Place the survey data file in a subfolder called  data/
#      (data are not included in this repository — see README.md
#      and the manuscript's Data Availability statement).
#   3. Open R / RStudio and set the working directory to the
#      repository root (the folder containing this file):
#         setwd("path/to/dea-sfa-tilapia-jambi")
#   4. Run scripts in order: J2_01 → J2_02 → ... → J2_08.
#      Each script begins with source("00_setup.R"), so you do
#      not need to source this file manually.
#
# NOTE: J2_08_MonteCarlo_SmallSample.R is fully self-contained
#   (no survey data required) and can be run independently.
# ============================================================


# ============================================================
# 1. REQUIRED PACKAGES
# ============================================================
# Install missing packages before first run:
#   install.packages(c("frontier","ggplot2","dplyr","tidyr",
#                      "openxlsx","officer","stargazer","lmtest"))

suppressPackageStartupMessages({

  # --- Stochastic Frontier Analysis ---
  library(frontier)      # sfa(), efficiencies()

  # --- Data manipulation ---
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)

  # --- Visualisation ---
  library(ggplot2)
  library(patchwork)     # multi-panel plots
  library(scales)

  # --- Output export ---
  library(openxlsx)      # Excel workbooks
  library(officer)       # Word documents
  library(stargazer)     # LaTeX tables (optional)

  # --- Statistical tests ---
  library(lmtest)        # coeftest()
  library(car)           # vif()
})

cat("OK: All packages loaded.\n\n")


# ============================================================
# 2. REPRODUCIBILITY
# ============================================================
set.seed(2026)


# ============================================================
# 3. DATA PATH  (relative to repository root)
# ============================================================
# Place the survey data file here before running J2_01 to J2_07.
# J2_08 (Monte Carlo) does not require survey data.
PATH_DATA  <- file.path("data", "Data_Master_Disertasi2.xlsx")
SHEET_NAME <- "Data_Master"
SKIP_ROWS  <- 2   # Row 1 = title; Row 2 = group header; Row 3 = column names


# ============================================================
# 4. OUTPUT DIRECTORIES  (created automatically on first run)
# ============================================================
DIR_ROOT   <- "Journal_2_rev"
DIR_RDS    <- file.path(DIR_ROOT, "rds")
DIR_XLSX   <- file.path(DIR_ROOT, "xlsx")
DIR_TEX    <- file.path(DIR_ROOT, "tex")
DIR_WORD   <- file.path(DIR_ROOT, "word")
DIR_PLOTS  <- file.path(DIR_ROOT, "plots")

for (d in c(DIR_RDS, DIR_XLSX, DIR_TEX, DIR_WORD, DIR_PLOTS)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("  Created:", d, "\n")
  }
}
cat("OK: Output directories ready under", DIR_ROOT, "\n\n")


# ============================================================
# 5. VARIABLE DEFINITIONS
# ============================================================

# DEA inputs (4 variables; X4_FERT and X5_MEDS excluded due to
# structural zeros — see companion paper for justification)
INPUT_VARS  <- c("X1_SIZE", "X2_SEED", "X3_FEED", "X6_LABR")
OUTPUT_VAR  <- "Y_OUT"        # Tilapia harvest (kg/cycle)
COST_VARS   <- c("C1_SIZE", "C2_SEED", "C3_FEED", "C6_LABR")

# Cluster identifier (1 = K1 lowland, 2 = K2 mid-elevation, 3 = K3 highland)
CLUSTER_VAR <- "klaster"

# Elevation clusters (for labelling in plots and tables)
CLUSTER_LABELS <- c(
  "K1" = "K1 (0\u2013100 m a.s.l.)",
  "K2" = "K2 (>100\u2013500 m a.s.l.)",
  "K3" = "K3 (>500 m a.s.l.)"
)

# Production-process indicators used in RQ3 and RQ4
INDIKATOR_H <- c("Density", "FCR", "SR", "Laborint")

# SFA formula (Cobb-Douglas, log-transformed inputs and output)
SFA_FORMULA <- ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6


# ============================================================
# 6. ANALYSIS CONSTANTS
# ============================================================
B_BOOTSTRAP <- 2000    # DEA bootstrap replications (Simar & Wilson 2007)
B_SFA_BOOT  <- 300     # SFA bootstrap replications (J2_03)
B_PAIRED    <- 2000    # Paired bootstrap replications for rho comparison (J2_06)
CV_K        <- 10      # Cross-validation folds (J2_07)
CV_R        <- 50      # Cross-validation repeats (J2_07)
MC_S        <- 1000    # Monte Carlo replications per sample size (J2_08)
ALPHA       <- 0.05    # Significance level
P_TRIM      <- 0.99    # Outlier trimming threshold (Winsorization at P99)


# ============================================================
# 7. HELPER: SIGNIFICANCE LABELS
# ============================================================
sig_label <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.10  ~ ".",
    TRUE      ~ "ns"
  )
}


# ============================================================
# 8. CHECK DATA FILE EXISTS
# ============================================================
if (file.exists(PATH_DATA)) {
  cat("OK: Data file found ->", PATH_DATA, "\n\n")
} else {
  cat("WARNING: Data file NOT found at:", PATH_DATA, "\n")
  cat("  Place the survey data file at:\n")
  cat("    <repo-root>/data/Data_Master_Disertasi2.xlsx\n")
  cat("  See README.md for data availability information.\n")
  cat("  Note: J2_08_MonteCarlo_SmallSample.R can still run\n")
  cat("  without survey data.\n\n")
}

cat("00_setup.R complete. Proceed with J2_00_Helpers.R, then J2_01 onward.\n\n")
