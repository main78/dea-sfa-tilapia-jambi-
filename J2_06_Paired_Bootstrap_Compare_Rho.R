# ============================================================
# J2_06_Paired_Bootstrap_Compare_Rho.R
# FORMAL TEST: ARE DIFFERENCES IN SPEARMAN RHO ACROSS SFA SPECIFICATIONS
# STATISTICALLY MEANINGFUL? (addresses reviewer's "Masalah 2")
# ============================================================
# Journal-2 pipeline (rebuilt) - Ma'in (NIM P3C123010), FEB UNJA
# Run AFTER J2_02 (reads J2_02_Results.rds; does NOT refit any model).
#
# MOTIVATION:
#   Table 4 in the manuscript reports four Spearman correlations between
#   SFA technical efficiency and the DEA meta-frontier:
#     per-cluster SFA = 0.7598 | pooled SFA = 0.8130
#     Model A (slope interactions) = 0.8051 | Model B (intercept shifters) = 0.8147
#   A reviewer could reasonably ask whether 0.813 vs 0.815 (pooled vs Model B)
#   is a statistically meaningful difference, or just noise.
#
# WHY NOT A FISHER r-to-z TEST:
#   Fisher's r-to-z comparison assumes the two correlations come from
#   INDEPENDENT samples. Here, all four correlations are computed on the
#   SAME 380 farm-level pairs - only the SFA estimator differs. This is a
#   paired-sample situation, so the appropriate test is a paired bootstrap
#   of the DIFFERENCE in rho, resampling farms (not refitting any model),
#   which automatically preserves the pairing across specifications.
#
# METHOD:
#   For each of the 6 pairwise comparisons among {per-cluster, pooled,
#   Model A, Model B}, draw B = 2000 bootstrap resamples of the 380 DMUs
#   (with replacement). In each replicate, recompute Spearman rho(TE_spec, DEA)
#   for BOTH specifications on the SAME resampled indices, then take the
#   difference. The resulting bootstrap distribution of differences gives a
#   percentile 95% CI and a two-sided bootstrap p-value for each comparison.
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/word/tex: J2_06_RQ_PairedBootstrap_RhoComparison
#   plots/J2_06_*.png
#   rds/J2_06_Results.rds
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

set.seed(2026)

cat("\n============================================================\n")
cat(" J2_06: PAIRED BOOTSTRAP COMPARISON OF SPEARMAN RHO\n")
cat("============================================================\n")

# ============================================================
# STEP 1: LOAD DATA (no refitting - all TE columns already saved)
# ============================================================
cat("\n--- STEP 1: Load Data from J2_02 ---\n")

j2_02 <- readRDS(file.path(DIR_RDS, "J2_02_Results.rds"))
data_sfa <- j2_02$data_sfa

stopifnot(all(c("TE_DEA", "TE_SFA_cluster", "TE_SFA_pooled", "TE_SFA_modelA", "TE_SFA_modelB") %in% names(data_sfa)))

specs <- list(
  PerCluster = data_sfa$TE_SFA_cluster,
  Pooled     = data_sfa$TE_SFA_pooled,
  ModelA     = data_sfa$TE_SFA_modelA,
  ModelB     = data_sfa$TE_SFA_modelB
)
te_dea <- data_sfa$TE_DEA
n <- nrow(data_sfa)
cat(sprintf("  OK. n=%d. Specifications loaded: %s\n", n, paste(names(specs), collapse = ", ")))

# ============================================================
# STEP 2: POINT-ESTIMATE RHO (cross-check against Table 4)
# ============================================================
cat("\n--- STEP 2: Point-Estimate Rho (cross-check vs. manuscript Table 4) ---\n")

point_rho <- sapply(specs, function(x) {
  ok <- complete.cases(x, te_dea)
  cor(x[ok], te_dea[ok], method = "spearman")
})
cat("  Point-estimate Spearman rho per specification:\n")
print(round(point_rho, 4))
cat("\n  (Compare against manuscript Table 4: PerCluster=0.7598, Pooled=0.8130,\n")
cat("   ModelA=0.8051, ModelB=0.8147. Small differences <0.001 may arise from\n")
cat("   rounding or NA handling and are not substantively meaningful.)\n")

# ============================================================
# STEP 3: PAIRED BOOTSTRAP OF RHO DIFFERENCES (B = 2000)
# ============================================================
cat("\n--- STEP 3: Paired Bootstrap (B = 2000) ---\n")

B <- 2000
spec_names <- names(specs)
pairs <- combn(spec_names, 2, simplify = FALSE)  # all 6 pairwise comparisons

boot_diff <- matrix(NA_real_, nrow = B, ncol = length(pairs))
colnames(boot_diff) <- sapply(pairs, function(p) paste(p[1], "minus", p[2]))

t_start <- Sys.time()
for (b in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  te_dea_b <- te_dea[idx]

  rho_b <- sapply(specs, function(x) {
    x_b <- x[idx]
    ok <- complete.cases(x_b, te_dea_b)
    if (sum(ok) < 10) return(NA_real_)
    cor(x_b[ok], te_dea_b[ok], method = "spearman")
  })

  for (j in seq_along(pairs)) {
    p <- pairs[[j]]
    boot_diff[b, j] <- rho_b[p[1]] - rho_b[p[2]]
  }

  if (b %% 500 == 0) {
    elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 2)
    cat(sprintf("  ... replication %d/%d done | %.2f min elapsed\n", b, B, elapsed))
  }
}

cat(sprintf("\n  Bootstrap complete. %d/%d replications fully valid.\n",
            sum(complete.cases(boot_diff)), B))

# ============================================================
# STEP 4: SUMMARY TABLE - CI AND BOOTSTRAP p-VALUE PER COMPARISON
# ============================================================
cat("\n--- STEP 4: Summary of Pairwise Differences ---\n")

bootstrap_p_two_sided <- function(boot_vec, point_diff) {
  ok <- !is.na(boot_vec)
  boot_vec <- boot_vec[ok]
  if (point_diff >= 0) {
    p1 <- mean(boot_vec < 0)
  } else {
    p1 <- mean(boot_vec > 0)
  }
  min(1, 2 * p1)
}

comparison_tbl <- do.call(rbind, lapply(seq_along(pairs), function(j) {
  p <- pairs[[j]]
  point_diff <- point_rho[p[1]] - point_rho[p[2]]
  boot_vec <- boot_diff[, j]
  data.frame(
    Comparison        = paste0(p[1], " minus ", p[2]),
    Rho_Spec1         = round(point_rho[p[1]], 4),
    Rho_Spec2         = round(point_rho[p[2]], 4),
    Point_Diff        = round(point_diff, 4),
    Bootstrap_SE      = round(sd(boot_vec, na.rm = TRUE), 4),
    CI_Lower_95       = round(quantile(boot_vec, 0.025, na.rm = TRUE), 4),
    CI_Upper_95       = round(quantile(boot_vec, 0.975, na.rm = TRUE), 4),
    Bootstrap_p_value = round(bootstrap_p_two_sided(boot_vec, point_diff), 4)
  )
}))
comparison_tbl$Significant <- sapply(comparison_tbl$Bootstrap_p_value, sig_label)
rownames(comparison_tbl) <- NULL

cat("\n  Pairwise comparison of Spearman rho (paired bootstrap, B = 2000):\n")
print(comparison_tbl)

n_sig <- sum(comparison_tbl$Significant != "ns")
cat(sprintf("\n  %d of %d pairwise comparisons are statistically significant (p<0.10).\n",
            n_sig, nrow(comparison_tbl)))

# Highlight the comparison most relevant to the manuscript's claim (Model B best)
modelb_rows <- comparison_tbl[grepl("ModelB", comparison_tbl$Comparison), ]
cat("\n  Comparisons involving Model B specifically (the specification favoured in the manuscript):\n")
print(modelb_rows)

# ============================================================
# STEP 5: PLOTS (English, no titles, black footnote - house style)
# ============================================================
cat("\n--- STEP 5: Plots ---\n")

# Focus on the two comparisons most central to the manuscript's narrative:
# Model B vs Pooled (the "is 0.815 really better than 0.813?" question), and
# Model B vs Per-cluster (the architectural-symmetry headline result).
# NOTE: built dynamically from actual column names produced by combn() above,
# rather than hardcoded strings, so the direction always matches regardless
# of how spec_names happens to be ordered.
key_pairs <- colnames(boot_diff)[
  grepl("ModelB", colnames(boot_diff)) &
    (grepl("Pooled", colnames(boot_diff)) | grepl("PerCluster", colnames(boot_diff)))
]
cat(sprintf("\n  Key pairs selected for histogram: %s\n", paste(key_pairs, collapse = " | ")))

boot_long <- do.call(rbind, lapply(key_pairs, function(nm) {
  data.frame(Comparison = nm, Diff = boot_diff[, colnames(boot_diff) == nm])
}))
point_df <- comparison_tbl[comparison_tbl$Comparison %in% key_pairs, c("Comparison", "Point_Diff")]

p1 <- ggplot(boot_long, aes(x = Diff)) +
  geom_histogram(fill = "#2E75B6", alpha = 0.75, bins = 40, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  geom_vline(data = point_df, aes(xintercept = Point_Diff), color = "#C0504D", linewidth = 0.8) +
  facet_wrap(~Comparison, scales = "free", nrow = 1) +
  labs(x = "Bootstrap difference in Spearman rho", y = "Frequency",
       caption = wrap_caption("B = 2000 paired bootstrap replications (farms resampled, not models refitted); vertical red line marks the full-sample point estimate.")) +
  theme_j2_en()
save_png_j2(p1, "J2_06_Bootstrap_Diff_Histograms", w = 10, h = 4.5)

p2 <- ggplot(comparison_tbl, aes(x = reorder(Comparison, Point_Diff), y = Point_Diff)) +
  geom_hline(yintercept = 0, color = "#000000", linewidth = 0.4) +
  geom_pointrange(aes(ymin = CI_Lower_95, ymax = CI_Upper_95), color = "#2E75B6", linewidth = 0.6, size = 0.4) +
  coord_flip() +
  labs(x = NULL, y = "Difference in Spearman rho (95% bootstrap CI)",
       caption = wrap_caption("Paired bootstrap, B = 2000; n = 380.")) +
  theme_j2_en()
save_png_j2(p2, "J2_06_Coefplot_All_Pairwise_Differences", w = 8, h = 5)

# ============================================================
# STEP 6: EXPORT
# ============================================================
all_tables <- list(
  Point_Estimate_Rho = list(
    df = data.frame(Specification = names(point_rho), Rho = round(point_rho, 4)),
    caption = "Table A1. Point-estimate Spearman rho per SFA specification (cross-check against Table 4)."
  ),
  Pairwise_Rho_Comparison = list(
    df = comparison_tbl,
    caption = "Table A2. Paired bootstrap comparison of Spearman rho across SFA specifications (B = 2000)."
  )
)
finalize_outputs(all_tables, "J2_06_RQ_PairedBootstrap_RhoComparison",
                  "Formal Comparison of Spearman Rho Across SFA Specifications (Paired Bootstrap)")

# ============================================================
# STEP 7: SAVE RDS
# ============================================================
saveRDS(list(point_rho = point_rho, boot_diff = boot_diff, comparison_tbl = comparison_tbl),
        file.path(DIR_RDS, "J2_06_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_06_Results.rds")))

cat("\n============================================================\n")
cat(" RINGKASAN AKHIR J2_06\n")
cat("============================================================\n")
print(comparison_tbl)
cat(sprintf("\n%d/%d perbandingan signifikan (p<0.10).\n", n_sig, nrow(comparison_tbl)))
cat("\nOK J2_06_Paired_Bootstrap_Compare_Rho.R selesai.\n")
