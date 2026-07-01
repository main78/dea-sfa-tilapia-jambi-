# ============================================================
# J2_05_Check_Divergence_By_Cluster.R
# QUICK CHECK: WHICH CLUSTER SHOWS THE LARGEST DEA-SFA DIVERGENCE?
# ============================================================
# Journal-2 pipeline (rebuilt) - Ma'in (NIM P3C123010), FEB UNJA
# Run AFTER J2_04. Not part of the core RQ1-RQ4 pipeline - this is a
# small follow-up check to resolve one open question flagged in the
# manuscript draft: does the cluster with the largest absolute
# DEA-SFA divergence (this paper) coincide with the cluster the
# companion dissertation's Bland-Altman analysis flagged as the locus
# of disagreement (K3), or does it instead point to K2 (the cluster
# flagged here as numerically fragile in RQ1-RQ2)?
#
# OUTPUT: printed to console only; also saved to
#   Journal_2_rev/xlsx/J2_05_Divergence_By_Cluster.xlsx
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

cat("\n============================================================\n")
cat(" J2_05: DIVERGENCE MAGNITUDE BY ELEVATION CLUSTER\n")
cat("============================================================\n")

j2_04 <- readRDS(file.path(DIR_RDS, "J2_04_Results.rds"))
data_sfa <- j2_04$data_sfa

divergence_by_cluster <- data_sfa %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    Mean_Abs_Divergence   = round(mean(Abs_Divergence), 4),
    Median_Abs_Divergence = round(median(Abs_Divergence), 4),
    SD_Abs_Divergence     = round(sd(Abs_Divergence), 4),
    .groups = "drop"
  )

cat("\n  Absolute divergence |TE_DEA - TE_SFA(Model B)| by cluster:\n")
print(divergence_by_cluster)

largest_cluster <- divergence_by_cluster$Cluster[which.max(divergence_by_cluster$Mean_Abs_Divergence)]
cat(sprintf("\n  Cluster with the LARGEST mean absolute divergence: %s\n", largest_cluster))
cat("  Compare this to:\n")
cat("    - RQ1/RQ2 fragility locus (parameter instability): K2\n")
cat("    - Companion dissertation's Bland-Altman locus (score disagreement): K3\n")

all_tables <- list(
  Divergence_By_Cluster = list(df = as.data.frame(divergence_by_cluster),
    caption = "Table. Absolute DEA-SFA divergence by elevation cluster.")
)
write_xlsx_multi(all_tables, file.path(DIR_XLSX, "J2_05_Divergence_By_Cluster.xlsx"))

cat("\nOK J2_05_Check_Divergence_By_Cluster.R complete.\n")
