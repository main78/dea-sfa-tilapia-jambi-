# ============================================================
# J2_02_Model_Specification_RQ2.R
# RQ2: WHICH SFA SPECIFICATION IS NUMERICALLY STABLE AND CONVERGENT?
# ============================================================
# Journal-2 pipeline (rebuilt, clean) - Ma'in (NIM P3C123010), FEB UNJA
# Step 2 of 4. Run AFTER J2_01.
#
# Compares two ways of letting SFA accommodate cross-cluster technology
# heterogeneity without re-introducing the small-K2-subsample fragility
# documented in J2_01:
#   Model A - pooled frontier with FULL slope interactions (D_K2, D_K3 x
#             every input). TE is computed MANUALLY via the closed-form
#             Battese & Coelli (1988) estimator, bypassing frontier's
#             efficiencies() function (which is unreliable whenever the
#             MLE covariance matrix is not positive semidefinite).
#   Model B - pooled frontier with cluster INTERCEPT shifters only
#             (D_K2, D_K3), common input slopes.
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/word/tex: J2_02_RQ2_Model_Specification
#   plots/J2_02_*.png
#   rds/J2_02_Results.rds   <- read by J2_03 and J2_04
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

if (!requireNamespace("frontier", quietly = TRUE)) {
  install.packages("frontier", repos = "https://cran.r-project.org")
}
suppressPackageStartupMessages(library(frontier))

set.seed(2026)

cat("\n============================================================\n")
cat(" J2_02: RQ2 - SFA MODEL SPECIFICATION COMPARISON\n")
cat("============================================================\n")

# ============================================================
# STEP 1: LOAD J2_01 RESULTS
# ============================================================
cat("\n--- STEP 1: Load J2_01 Results ---\n")

j2_01 <- readRDS(file.path(DIR_RDS, "J2_01_Results.rds"))
data_sfa <- j2_01$data_sfa
fits_all <- j2_01$fits_all
cat(sprintf("  OK. n=%d.\n", nrow(data_sfa)))

data_sfa$D_K2 <- as.integer(data_sfa$klaster == 2)
data_sfa$D_K3 <- as.integer(data_sfa$klaster == 3)

# ============================================================
# STEP 2: FIT MODEL A (FULL SLOPE INTERACTIONS) - MANUAL BC88 TE
# ============================================================
cat("\n--- STEP 2: Fit Model A (full slope interactions) ---\n")

formula_A <- as.formula(
  ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6 + D_K2 + D_K3 +
    D_K2:ln_X1 + D_K2:ln_X2 + D_K2:ln_X3 + D_K2:ln_X6 +
    D_K3:ln_X1 + D_K3:ln_X2 + D_K3:ln_X3 + D_K3:ln_X6
)

cov_warning_A <- FALSE
fit_A <- withCallingHandlers(
  tryCatch(sfa(formula_A, data = data_sfa, ineffDecrease = TRUE),
           error = function(e) { cat(sprintf("  FAILED: %s\n", e$message)); NULL }),
  warning = function(w) {
    cat(sprintf("  [WARNING] %s\n", w$message))
    if (grepl("positive semidefinite", w$message)) cov_warning_A <<- TRUE
    invokeRestart("muffleWarning")
  }
)

model_A <- list(label = "Model A: full slope interactions", n = nrow(data_sfa),
                  converged = FALSE, gamma = NA_real_, loglik = NA_real_,
                  boundary_issue = FALSE, cov_warning = cov_warning_A, coef = NULL, TE = NULL)

# Battese & Coelli (1988) closed-form TE estimator - bypasses efficiencies()
compute_TE_BC88 <- function(e, sigmaSq, gamma) {
  sigma_u2    <- gamma * sigmaSq
  sigma_v2    <- (1 - gamma) * sigmaSq
  sigma_star2 <- sigma_u2 * sigma_v2 / sigmaSq
  sigma_star  <- sqrt(sigma_star2)
  mu_star     <- -e * sigma_u2 / sigmaSq
  z           <- mu_star / sigma_star
  exp(-mu_star + 0.5 * sigma_star2) * pnorm(z - sigma_star) / pnorm(z)
}

if (!is.null(fit_A)) {
  smry_A <- tryCatch(summary(fit_A), error = function(e) NULL)
  if (!is.null(smry_A)) {
    mle_A <- as.data.frame(smry_A$mleParam); mle_A$Parameter <- rownames(mle_A)
    gamma_A <- mle_A[grepl("^gamma$", mle_A$Parameter, ignore.case = TRUE), "Estimate"]
    gamma_A <- if (length(gamma_A) == 1) gamma_A else NA_real_
    sigmaSq_A <- mle_A[grepl("sigmaSq", mle_A$Parameter, ignore.case = TRUE), "Estimate"]
    sigmaSq_A <- if (length(sigmaSq_A) == 1) sigmaSq_A else NA_real_

    model_A$converged       <- TRUE
    model_A$gamma            <- gamma_A
    model_A$loglik           <- tryCatch(as.numeric(logLik(fit_A)), error = function(e) NA_real_)
    model_A$boundary_issue    <- !is.na(gamma_A) && (gamma_A < 0.05 || gamma_A > 0.95)
    model_A$coef             <- mle_A

    X_mat <- model.matrix(formula_A, data = data_sfa)
    beta_named <- setNames(mle_A$Estimate, mle_A$Parameter)
    beta_match <- beta_named[colnames(X_mat)]
    fitted_vals <- as.numeric(X_mat %*% beta_match)
    e_resid <- data_sfa$ln_Y - fitted_vals
    model_A$TE <- compute_TE_BC88(e_resid, sigmaSq_A, gamma_A)

    cat(sprintf("  OK. gamma=%.4f | logLik=%.2f | boundary=%s | cov_warning=%s\n",
                gamma_A, model_A$loglik, ifelse(model_A$boundary_issue, "YES", "no"),
                ifelse(model_A$cov_warning, "YES", "no")))
    cat(sprintf("  TE (manual BC88): mean=%.4f, sd=%.4f, n_valid=%d/%d\n",
                mean(model_A$TE, na.rm = TRUE), sd(model_A$TE, na.rm = TRUE),
                sum(!is.na(model_A$TE)), length(model_A$TE)))
  }
}

# ============================================================
# STEP 3: FIT MODEL B (INTERCEPT SHIFTERS ONLY)
# ============================================================
cat("\n--- STEP 3: Fit Model B (intercept shifters only) ---\n")

formula_B <- as.formula(ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6 + D_K2 + D_K3)

cov_warning_B <- FALSE
fit_B <- withCallingHandlers(
  tryCatch(sfa(formula_B, data = data_sfa, ineffDecrease = TRUE),
           error = function(e) { cat(sprintf("  FAILED: %s\n", e$message)); NULL }),
  warning = function(w) {
    cat(sprintf("  [WARNING] %s\n", w$message))
    if (grepl("positive semidefinite", w$message)) cov_warning_B <<- TRUE
    invokeRestart("muffleWarning")
  }
)

model_B <- list(label = "Model B: intercept shifters only", n = nrow(data_sfa),
                  converged = FALSE, gamma = NA_real_, loglik = NA_real_,
                  boundary_issue = FALSE, cov_warning = cov_warning_B, coef = NULL, TE = NULL)

if (!is.null(fit_B)) {
  smry_B <- tryCatch(summary(fit_B), error = function(e) NULL)
  if (!is.null(smry_B)) {
    mle_B <- as.data.frame(smry_B$mleParam); mle_B$Parameter <- rownames(mle_B)
    gamma_B <- mle_B[grepl("^gamma$", mle_B$Parameter, ignore.case = TRUE), "Estimate"]
    gamma_B <- if (length(gamma_B) == 1) gamma_B else NA_real_
    model_B$converged       <- TRUE
    model_B$gamma            <- gamma_B
    model_B$loglik           <- tryCatch(as.numeric(logLik(fit_B)), error = function(e) NA_real_)
    model_B$boundary_issue    <- !is.na(gamma_B) && (gamma_B < 0.05 || gamma_B > 0.95)
    model_B$coef             <- mle_B
    model_B$TE                <- tryCatch(as.numeric(efficiencies(fit_B)), error = function(e) NULL)

    if (is.null(model_B$TE)) {  # fallback, same manual estimator
      sigmaSq_B <- mle_B[grepl("sigmaSq", mle_B$Parameter, ignore.case = TRUE), "Estimate"]
      X_matB <- model.matrix(formula_B, data = data_sfa)
      beta_B <- setNames(mle_B$Estimate, mle_B$Parameter)
      fitted_B <- as.numeric(X_matB %*% beta_B[colnames(X_matB)])
      model_B$TE <- compute_TE_BC88(data_sfa$ln_Y - fitted_B, sigmaSq_B, gamma_B)
    }
    cat(sprintf("  OK. gamma=%.4f | logLik=%.2f | boundary=%s | cov_warning=%s\n",
                gamma_B, model_B$loglik, ifelse(model_B$boundary_issue, "YES", "no"),
                ifelse(model_B$cov_warning, "YES", "no")))
  }
}

data_sfa$TE_SFA_modelA <- model_A$TE
data_sfa$TE_SFA_modelB <- model_B$TE

# ============================================================
# STEP 4: COMPARE ALL SPECIFICATIONS AGAINST DEA
# ============================================================
cat("\n--- STEP 4: Compare All Specifications vs. DEA ---\n")

compare_rho <- function(x, y, label) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 5) return(data.frame(Specification = label, n = sum(ok), Rho = NA_real_, p_value = NA_real_))
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  data.frame(Specification = label, n = sum(ok), Rho = round(unname(ct$estimate), 4),
             p_value = round(ct$p.value, 6))
}

rho_spec_tbl <- rbind(
  compare_rho(data_sfa$TE_SFA_cluster, data_sfa$TE_DEA, "Per-cluster SFA (J2_01)"),
  compare_rho(data_sfa$TE_SFA_pooled,  data_sfa$TE_DEA, "Pooled SFA, no shifters (J2_01)"),
  compare_rho(data_sfa$TE_SFA_modelB,  data_sfa$TE_DEA, "Model B: intercept shifters only"),
  compare_rho(data_sfa$TE_SFA_modelA,  data_sfa$TE_DEA, "Model A: full slope interactions")
)
cat("\n  Spearman rho vs. DEA, all specifications:\n"); print(rho_spec_tbl)

spec_diag_tbl <- data.frame(
  Specification = c("Model B: intercept shifters only", "Model A: full slope interactions"),
  n = c(model_B$n, model_A$n),
  Gamma = round(c(model_B$gamma, model_A$gamma), 4),
  Boundary_Issue = c(model_B$boundary_issue, model_A$boundary_issue),
  Covariance_Warning = c(model_B$cov_warning, model_A$cov_warning)
)
cat("\n  Numerical diagnostics, Model A vs. Model B:\n"); print(spec_diag_tbl)

best_spec <- rho_spec_tbl$Specification[which.max(rho_spec_tbl$Rho)]
cat(sprintf("\n  Highest-converging specification: %s\n", best_spec))

# ============================================================
# STEP 5: PLOTS
# ============================================================
cat("\n--- STEP 5: Plots ---\n")

p1 <- ggplot(rho_spec_tbl, aes(x = reorder(Specification, Rho), y = Rho)) +
  geom_col(fill = "#2E75B6", width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", Rho)), hjust = -0.15, size = 3.4, color = "#000000") +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "Spearman rho vs. DEA meta-frontier TE",
       caption = "Spearman rank correlation across four candidate SFA specifications; n = 380.") +
  theme_j2_en()
save_png_j2(p1, "J2_02_Barplot_Rho_AllSpecifications", w = 8.5, h = 5)

p2 <- ggplot(data_sfa, aes(TE_SFA_modelB, TE_DEA, color = Cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  scale_color_manual(values = c("K1 (0-100 m a.s.l.)" = "#2E75B6", "K2 (>100-500 m a.s.l.)" = "#E69138",
                                 "K3 (>500 m a.s.l.)" = "#1A6B3C")) +
  labs(x = "Technical efficiency, SFA (Model B)", y = "Technical efficiency, DEA (meta-frontier)",
       caption = "Model B: pooled frontier with cluster-specific intercept shifters; dashed line is the 45-degree reference.") +
  theme_j2_en()
save_png_j2(p2, "J2_02_Scatter_ModelB_vs_DEA", w = 7.5, h = 6)

# ============================================================
# STEP 6: EXPORT
# ============================================================
all_tables <- list(
  Rho_All_Specifications = list(df = rho_spec_tbl,
    caption = "Table 4. Spearman rank correlation with DEA meta-frontier TE, across SFA specifications."),
  Model_A_vs_B_Diagnostics = list(df = spec_diag_tbl,
    caption = "Table 5. Numerical diagnostics: Model A (full slope interactions) vs. Model B (intercept shifters only).")
)
finalize_outputs(all_tables, "J2_02_RQ2_Model_Specification",
                  "RQ2 - SFA Model Specification Comparison")

# ============================================================
# STEP 7: SAVE RDS FOR J2_03 AND J2_04
# ============================================================
saveRDS(list(data_sfa = data_sfa, model_A = model_A, model_B = model_B,
             rho_spec_tbl = rho_spec_tbl, spec_diag_tbl = spec_diag_tbl),
        file.path(DIR_RDS, "J2_02_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_02_Results.rds")))

cat("\nOK J2_02_Model_Specification_RQ2.R complete.\n")
