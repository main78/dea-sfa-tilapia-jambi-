# ============================================================
# J2_04_Divergence_Predictors_RQ4.R
# RQ4: WHAT SYSTEMATICALLY PREDICTS THE MAGNITUDE OF DEA-SFA DIVERGENCE?
# ============================================================
# Journal-2 pipeline (rebuilt, clean) - Ma'in (NIM P3C123010), FEB UNJA
# Step 4 of 4 (final). Run AFTER J2_02.
#
# Regresses |TE_DEA - TE_SFA(Model B)| on AWQI, cluster sample size, FCR,
# SR, and farming system (WAD). A cluster-level skewness covariate was
# tested in a preliminary pilot and found to be collinear with AWQI
# (r = 0.816); it is dropped from the main specification and retained
# only as part of the robustness comparison below, so that no unstable
# intermediate result is reported as if it were final.
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/word/tex: J2_04_RQ4_Divergence_Predictors
#   plots/J2_04_*.png
#   rds/J2_04_Results.rds
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

set.seed(2026)

cat("\n============================================================\n")
cat(" J2_04: RQ4 - PREDICTORS OF DEA-SFA DIVERGENCE\n")
cat("============================================================\n")

# ============================================================
# STEP 1: LOAD DATA + BUILD DIVERGENCE & PREDICTORS
# ============================================================
cat("\n--- STEP 1: Load Data & Build Variables ---\n")

j2_02 <- readRDS(file.path(DIR_RDS, "J2_02_Results.rds"))
data_sfa <- j2_02$data_sfa
stopifnot(all(c("TE_DEA", "TE_SFA_modelB", "klaster", "WAD", "FCR", "SR") %in% names(data_sfa)))

data_sfa$Divergence     <- data_sfa$TE_DEA - data_sfa$TE_SFA_modelB
data_sfa$Abs_Divergence <- abs(data_sfa$Divergence)
cat(sprintf("  Divergence (TE_DEA - TE_SFA Model B): mean=%.4f, sd=%.4f\n",
            mean(data_sfa$Divergence, na.rm = TRUE), sd(data_sfa$Divergence, na.rm = TRUE)))

# AWQI per cluster x farming-system combination (Tabel 5.86, dissertation - final values)
awqi_map <- data.frame(klaster = c(1, 1, 2, 2, 3, 3), WAD = c(0, 1, 0, 1, 0, 1),
                        AWQI = c(60.85, 74.15, 45.50, 65.73, 90.67, 77.31))
data_sfa <- merge(data_sfa, awqi_map, by = c("klaster", "WAD"), all.x = TRUE)

n_klaster_map <- data.frame(klaster = c(1, 2, 3), Cluster_n = c(193, 51, 136))
data_sfa <- merge(data_sfa, n_klaster_map, by = "klaster", all.x = TRUE)

skew_fun <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); mean((x - m)^3, na.rm = TRUE) / s^3 }
skew_tbl <- data_sfa %>% group_by(klaster) %>% summarise(Skewness = skew_fun(ln_Y), .groups = "drop")
data_sfa <- merge(data_sfa, skew_tbl, by = "klaster", all.x = TRUE)

data_sfa$FCR_z <- as.numeric(scale(data_sfa$FCR))
data_sfa$SR_z  <- as.numeric(scale(data_sfa$SR))

# ============================================================
# STEP 2: COLLINEARITY CHECK
# ============================================================
cat("\n--- STEP 2: Collinearity Check Among Predictors ---\n")

pred_vars <- c("AWQI", "Cluster_n", "Skewness", "FCR_z", "SR_z", "WAD")
cor_mat <- cor(data_sfa[, pred_vars], use = "complete.obs")
cat("  Correlation matrix:\n"); print(round(cor_mat, 3))
high_pairs <- which(abs(cor_mat) > 0.8 & abs(cor_mat) < 1, arr.ind = TRUE)
if (nrow(high_pairs) > 0) {
  cat("\n  WARNING: AWQI and Skewness are highly correlated (|r| > 0.8).\n")
  cat("  Skewness is therefore EXCLUDED from the main specification (Step 3)\n")
  cat("  and retained only in the robustness check (Step 4).\n")
}

# ============================================================
# STEP 3: MAIN SPECIFICATION (Skewness excluded)
# ============================================================
cat("\n--- STEP 3: Main OLS Specification (Skewness Excluded) ---\n")

mod_main <- lm(Abs_Divergence ~ AWQI + Cluster_n + FCR_z + SR_z + WAD, data = data_sfa)
smry_main <- summary(mod_main)
print(smry_main)

extract_coef <- function(mod) {
  s <- summary(mod)
  df <- as.data.frame(s$coefficients)
  df$Variable <- rownames(df)
  names(df) <- c("Estimate", "SE", "t_value", "p_value", "Variable")
  df <- df[df$Variable != "(Intercept)", c("Variable", "Estimate", "SE", "p_value")]
  df$Significance <- sapply(df$p_value, sig_label)
  rownames(df) <- NULL
  df
}
coef_main <- extract_coef(mod_main)
cat("\n  Main specification coefficients:\n"); print(coef_main)
cat(sprintf("\n  R-squared = %.4f | Adj. R-squared = %.4f | n = %d\n",
            smry_main$r.squared, smry_main$adj.r.squared, nrow(data_sfa)))

# ============================================================
# STEP 4: ROBUSTNESS CHECK (Skewness included)
# ============================================================
cat("\n--- STEP 4: Robustness Check (Skewness Included) ---\n")

mod_robust <- lm(Abs_Divergence ~ AWQI + Cluster_n + Skewness + FCR_z + SR_z + WAD, data = data_sfa)
smry_robust <- summary(mod_robust)
coef_robust <- extract_coef(mod_robust)
cat("  Robustness-check coefficients (with Skewness):\n"); print(coef_robust)

robust_compare <- merge(
  coef_main[, c("Variable", "Estimate", "p_value")],
  coef_robust[, c("Variable", "Estimate", "p_value")],
  by = "Variable", suffixes = c("_Main", "_Robustness"), all = TRUE
)
robust_compare$Pct_Change <- round(100 * (robust_compare$Estimate_Robustness - robust_compare$Estimate_Main) /
                                      abs(robust_compare$Estimate_Main), 1)
robust_compare$Stable <- ifelse(is.na(robust_compare$Pct_Change), TRUE, abs(robust_compare$Pct_Change) < 20)
cat("\n  Main vs. robustness-check coefficient comparison:\n"); print(robust_compare)

model_fit_tbl <- data.frame(
  Model = c("Main (Skewness excluded)", "Robustness check (Skewness included)"),
  R_squared = c(smry_main$r.squared, smry_robust$r.squared),
  Adj_R_squared = c(smry_main$adj.r.squared, smry_robust$adj.r.squared)
)
cat("\n  Model fit comparison:\n"); print(model_fit_tbl)

# ============================================================
# STEP 5: PLOTS
# ============================================================
cat("\n--- STEP 5: Plots ---\n")

data_sfa$Cluster <- factor(data_sfa$klaster, 1:3,
                             c("K1 (0-100 m a.s.l.)", "K2 (>100-500 m a.s.l.)", "K3 (>500 m a.s.l.)"))
COL_K <- c("K1 (0-100 m a.s.l.)" = "#2E75B6", "K2 (>100-500 m a.s.l.)" = "#E69138",
           "K3 (>500 m a.s.l.)" = "#1A6B3C")

p1 <- ggplot(data_sfa, aes(x = Cluster, y = Abs_Divergence, fill = Cluster)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8, alpha = 0.85) +
  scale_fill_manual(values = COL_K, guide = "none") +
  labs(x = "Elevation cluster (sample size: K1 = 193, K2 = 51, K3 = 136)",
       y = "Absolute divergence |TE(DEA) - TE(SFA)|",
       caption = "Distribution across elevation clusters; n = 380.") +
  theme_j2_en()
save_png_j2(p1, "J2_04_Boxplot_Divergence_by_Cluster", w = 7, h = 5.2)

p2 <- ggplot(data_sfa, aes(x = AWQI, y = Abs_Divergence, color = Cluster)) +
  geom_point(alpha = 0.55, size = 1.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "#000000",
              fill = "grey75", linewidth = 0.6, alpha = 0.25) +
  scale_color_manual(values = COL_K) +
  labs(x = "Aquaculture Water Quality Index (AWQI)", y = "Absolute divergence |TE(DEA) - TE(SFA)|",
       caption = "OLS trend line with 95% confidence band; n = 380.") +
  theme_j2_en()
save_png_j2(p2, "J2_04_Scatter_Divergence_vs_AWQI", w = 7.5, h = 5.5)

coef_plot_df <- coef_main
coef_plot_df$CI_lo <- coef_plot_df$Estimate - 1.96 * coef_plot_df$SE
coef_plot_df$CI_hi <- coef_plot_df$Estimate + 1.96 * coef_plot_df$SE
p3 <- ggplot(coef_plot_df, aes(x = reorder(Variable, Estimate), y = Estimate)) +
  geom_hline(yintercept = 0, color = "#000000", linewidth = 0.4) +
  geom_pointrange(aes(ymin = CI_lo, ymax = CI_hi), color = "#2E75B6", linewidth = 0.6, size = 0.45) +
  coord_flip() +
  labs(x = NULL, y = "Coefficient (effect on absolute divergence)",
       caption = "OLS regression; 95% confidence intervals; n = 380.") +
  theme_j2_en()
save_png_j2(p3, "J2_04_Coefplot_Main_Specification", w = 7.5, h = 5)

# ============================================================
# STEP 6: EXPORT
# ============================================================
all_tables <- list(
  Main_Coefficients = list(df = coef_main,
    caption = "Table 8. OLS determinants of absolute DEA-SFA divergence (main specification)."),
  Robustness_Comparison = list(df = robust_compare,
    caption = "Table 9. Coefficient stability: main specification vs. robustness check including cluster skewness."),
  Model_Fit_Comparison = list(df = model_fit_tbl,
    caption = "Table 10. Model fit comparison, main specification vs. robustness check.")
)
finalize_outputs(all_tables, "J2_04_RQ4_Divergence_Predictors",
                  "RQ4 - Predictors of DEA-SFA Divergence Magnitude")

# ============================================================
# STEP 7: SAVE RDS
# ============================================================
saveRDS(list(data_sfa = data_sfa, mod_main = mod_main, mod_robust = mod_robust,
             coef_main = coef_main, robust_compare = robust_compare, model_fit_tbl = model_fit_tbl),
        file.path(DIR_RDS, "J2_04_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_04_Results.rds")))

cat("\nOK J2_04_Divergence_Predictors_RQ4.R complete.\n")
cat("\n============================================================\n")
cat(" JOURNAL-2 PIPELINE (REBUILT) COMPLETE: J2_01 -> J2_02 -> J2_03 -> J2_04\n")
cat(" All outputs are under: Journal_2_rev/\n")
cat("============================================================\n")
