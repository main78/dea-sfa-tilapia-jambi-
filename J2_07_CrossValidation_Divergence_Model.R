# ============================================================
# J2_07_CrossValidation_Divergence_Model.R
# OUT-OF-SAMPLE VALIDATION OF THE RQ4 DIVERGENCE-PREDICTOR REGRESSION
# (addresses reviewer's "no cross-validation / holdout sample" critique)
# ============================================================
# Journal-2 pipeline (rebuilt) - Ma'in (NIM P3C123010), FEB UNJA
# Run AFTER J2_04 (reads J2_04_Results.rds; refits the SAME OLS specification
# repeatedly on resampled training/test splits - does not change the model).
#
# MOTIVATION:
#   Section 4.5 reports R-squared = 0.4545 for the main divergence-predictor
#   model, computed in-sample (training set = testing set = all 380 farms).
#   A reviewer could reasonably ask whether this R-squared reflects genuine
#   out-of-sample predictive power or merely in-sample overfitting.
#
# METHOD:
#   Repeated k-fold cross-validation (k = 10 folds, R = 50 repeats = 500
#   total fold-fits). In each repeat, the 380 farms are randomly partitioned
#   into 10 folds; for each fold, the model is fit on the other 9 folds
#   (training set) and R-squared is computed on the held-out fold (test set)
#   by correlating predicted and observed absolute divergence. The same
#   procedure is run for (a) the main specification (AWQI, ClusterN, FCR_z,
#   SR_z, WAD) and (b) the robustness specification (same + Skewness), to
#   check whether the stability ranking found in-sample (Table 9) also holds
#   out-of-sample.
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/word/tex: J2_07_RQ4_CrossValidation_OutOfSample
#   plots/J2_07_*.png
#   rds/J2_07_Results.rds
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

set.seed(2026)

cat("\n============================================================\n")
cat(" J2_07: K-FOLD CROSS-VALIDATION OF THE RQ4 DIVERGENCE MODEL\n")
cat("============================================================\n")

# ============================================================
# STEP 1: LOAD DATA (same data_sfa used in J2_04 - no refitting of SFA/DEA)
# ============================================================
cat("\n--- STEP 1: Load Data from J2_04 ---\n")

j2_04 <- readRDS(file.path(DIR_RDS, "J2_04_Results.rds"))
data_div <- j2_04$data_sfa  # must contain Abs_Divergence, AWQI, Cluster_n, FCR_z, SR_z, WAD, Skewness

required_vars <- c("Abs_Divergence", "AWQI", "Cluster_n", "FCR_z", "SR_z", "WAD", "Skewness")
stopifnot(all(required_vars %in% names(data_div)))
data_div <- data_div[complete.cases(data_div[, required_vars]), ]
n <- nrow(data_div)
cat(sprintf("  OK. n=%d complete cases.\n", n))

specs <- list(
  Main       = Abs_Divergence ~ AWQI + Cluster_n + FCR_z + SR_z + WAD,
  Robustness = Abs_Divergence ~ AWQI + Cluster_n + FCR_z + SR_z + WAD + Skewness
)

# ============================================================
# STEP 2: IN-SAMPLE R-SQUARED (cross-check against manuscript Table 8/10)
# ============================================================
cat("\n--- STEP 2: In-Sample R-squared (cross-check vs. Table 8/10) ---\n")

insample_r2 <- sapply(specs, function(f) {
  m <- lm(f, data = data_div)
  summary(m)$r.squared
})
cat("  In-sample R-squared per specification:\n")
print(round(insample_r2, 4))
cat("\n  (Compare against manuscript Table 8/10: Main=0.4545, Robustness=0.4559.)\n")

# ============================================================
# STEP 3: REPEATED 10-FOLD CROSS-VALIDATION (R = 50 repeats)
# ============================================================
cat("\n--- STEP 3: Repeated 10-Fold Cross-Validation (k=10, R=50) ---\n")

k <- 10
R <- 50
cv_r2 <- data.frame(Repeat = integer(), Fold = integer(), Spec = character(), R2_OOS = numeric())

t_start <- Sys.time()
for (rep_i in seq_len(R)) {
  fold_id <- sample(rep(1:k, length.out = n))

  for (fold in seq_len(k)) {
    train <- data_div[fold_id != fold, ]
    test  <- data_div[fold_id == fold, ]

    for (spec_name in names(specs)) {
      m_train <- lm(specs[[spec_name]], data = train)
      pred <- predict(m_train, newdata = test)
      obs  <- test$Abs_Divergence

      # Out-of-sample R-squared: 1 - SS_res/SS_tot, using the TRAINING mean
      # as the naive benchmark (standard out-of-sample R2 definition)
      ss_res <- sum((obs - pred)^2)
      ss_tot <- sum((obs - mean(train$Abs_Divergence))^2)
      r2_oos <- 1 - ss_res / ss_tot

      cv_r2 <- rbind(cv_r2, data.frame(Repeat = rep_i, Fold = fold, Spec = spec_name, R2_OOS = r2_oos))
    }
  }
  if (rep_i %% 10 == 0) {
    elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 2)
    cat(sprintf("  ... repeat %d/%d done | %.2f min elapsed\n", rep_i, R, elapsed))
  }
}
cat(sprintf("\n  Cross-validation complete. %d total fold-fits per specification.\n", k * R))

# ============================================================
# STEP 4: SUMMARY - IN-SAMPLE VS. OUT-OF-SAMPLE R-SQUARED
# ============================================================
cat("\n--- STEP 4: Summary ---\n")

cv_summary <- do.call(rbind, lapply(names(specs), function(spec_name) {
  vals <- cv_r2$R2_OOS[cv_r2$Spec == spec_name]
  data.frame(
    Specification     = spec_name,
    InSample_R2       = round(insample_r2[spec_name], 4),
    Mean_OOS_R2       = round(mean(vals), 4),
    Median_OOS_R2     = round(median(vals), 4),
    SD_OOS_R2         = round(sd(vals), 4),
    Pct_Negative_OOS  = round(100 * mean(vals < 0), 1),
    CI_Lower_95       = round(quantile(vals, 0.025), 4),
    CI_Upper_95       = round(quantile(vals, 0.975), 4)
  )
}))
rownames(cv_summary) <- NULL

cat("\n  In-sample vs. out-of-sample R-squared (500 fold-fits per specification):\n")
print(cv_summary)

shrinkage <- cv_summary$InSample_R2 - cv_summary$Mean_OOS_R2
cat(sprintf("\n  R-squared shrinkage (in-sample minus mean out-of-sample): Main=%.4f, Robustness=%.4f\n",
            shrinkage[1], shrinkage[2]))
cat("  A small, stable shrinkage indicates the in-sample R-squared reflects genuine\n")
cat("  predictive structure rather than overfitting; a large or highly variable\n")
cat("  shrinkage (or frequent negative out-of-sample R-squared) would indicate overfitting.\n")

# ============================================================
# STEP 5: PLOT
# ============================================================
cat("\n--- STEP 5: Plot ---\n")

p1 <- ggplot(cv_r2, aes(x = Spec, y = R2_OOS, fill = Spec)) +
  geom_boxplot(width = 0.45, alpha = 0.8, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  geom_point(data = data.frame(Spec = names(specs), R2_OOS = insample_r2),
             aes(x = Spec, y = R2_OOS), color = "#C0504D", size = 3, shape = 18) +
  scale_fill_manual(values = c(Main = "#2E75B6", Robustness = "#70AD47")) +
  labs(x = NULL, y = "Out-of-sample R-squared (10-fold CV, 50 repeats)",
       caption = wrap_caption("Boxplots show 500 out-of-sample fold R-squared values per specification; red diamond marks the in-sample R-squared for comparison.")) +
  theme_j2_en() + theme(legend.position = "none")
save_png_j2(p1, "J2_07_OutOfSample_R2_Boxplot", w = 7, h = 5.5)

# ============================================================
# STEP 6: EXPORT
# ============================================================
all_tables <- list(
  CV_Summary = list(
    df = cv_summary,
    caption = "Table 13. In-sample versus repeated 10-fold cross-validated out-of-sample R-squared for the RQ4 divergence-predictor model."
  )
)
finalize_outputs(all_tables, "J2_07_RQ4_CrossValidation_OutOfSample",
                  "Out-of-Sample Cross-Validation of the Divergence-Predictor Model")

# ============================================================
# STEP 7: SAVE RDS
# ============================================================
saveRDS(list(cv_r2 = cv_r2, cv_summary = cv_summary, insample_r2 = insample_r2),
        file.path(DIR_RDS, "J2_07_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_07_Results.rds")))

cat("\n============================================================\n")
cat(" RINGKASAN AKHIR J2_07\n")
cat("============================================================\n")
print(cv_summary)
cat("\nOK J2_07_CrossValidation_Divergence_Model.R selesai.\n")
