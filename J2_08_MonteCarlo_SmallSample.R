# ============================================================
# J2_08_MonteCarlo_SmallSample.R
# MONTE CARLO SIMULATION (OPTION B - single-factor, sample size only)
# Does the small-sample fragility documented for K2 (n=51) reflect a
# general property of SFA estimation, or a Jambi-specific quirk?
# ============================================================
# Journal-2 pipeline (rebuilt) - Ma'in (NIM P3C123010), FEB UNJA
# STANDALONE: does NOT read any real survey data. The data-generating
# process (DGP) below is entirely synthetic, with known true parameters,
# so that estimated technical efficiency can be compared directly against
# TRUE technical efficiency (exp(-u_i)) rather than against DEA. This
# isolates the SFA-side estimation mechanism cleanly, without any
# confound from DEA's own small-sample behaviour.
#
# DESIGN:
#   For n in {50, 100, 200} - chosen to bracket the real clusters
#   (K2 n=51, K3 n=136, K1 n=193) - simulate S = 1,000 independent
#   samples from a Cobb-Douglas frontier with TRUE gamma = 0.80 (close to
#   the empirically estimated gamma for K1/K3 in the real data). For each
#   simulated sample, fit a half-normal Cobb-Douglas SFA by maximum
#   likelihood and record:
#     (a) whether the model converged at all,
#     (b) whether the estimated gamma sits at the boundary (>0.99),
#     (c) the Spearman rank correlation between estimated TE (via the
#         Battese-Coelli, 1988, closed-form point estimator) and the
#         TRUE simulated TE = exp(-u_i).
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/word/tex: J2_08_MonteCarlo_SmallSample_OptionB
#   plots/J2_08_*.png
#   rds/J2_08_Results.rds
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")
library(frontier)

set.seed(2026)

cat("\n============================================================\n")
cat(" J2_08: MONTE CARLO SIMULATION - SAMPLE SIZE AND SFA FRAGILITY\n")
cat("============================================================\n")

# ============================================================
# STEP 1: DGP PARAMETERS (fixed, known truth)
# ============================================================
cat("\n--- STEP 1: Data-Generating Process ---\n")

true_beta   <- c(beta0 = 1.00, beta1 = 0.20, beta2 = 0.25, beta3 = 0.30, beta4 = 0.15)
sigma_v     <- 0.15                                  # noise SD
true_gamma  <- 0.80                                  # target variance ratio
sigma_u     <- sqrt(true_gamma * sigma_v^2 / (1 - true_gamma))  # implied inefficiency SD

cat(sprintf("  True beta:    %s\n", paste(names(true_beta), round(true_beta, 3), sep = "=", collapse = ", ")))
cat(sprintf("  sigma_v = %.4f | sigma_u = %.4f | implied true gamma = %.4f\n",
            sigma_v, sigma_u, sigma_u^2 / (sigma_u^2 + sigma_v^2)))

# Battese & Coelli (1988) closed-form point estimator for individual TE,
# implemented directly (self-contained; does not assume a specific
# helper signature from J2_00_Helpers.R)
compute_TE_BC88_local <- function(resid, sigma_u, sigma_v) {
  sigma2   <- sigma_u^2 + sigma_v^2
  sigma_st <- sqrt(sigma_u^2 * sigma_v^2 / sigma2)
  mu_st    <- -resid * sigma_u^2 / sigma2
  TE <- exp(-mu_st + 0.5 * sigma_st^2) *
    pnorm(mu_st / sigma_st - sigma_st) / pnorm(mu_st / sigma_st)
  pmin(pmax(TE, 1e-6), 1)
}

# ============================================================
# STEP 2: SIMULATION LOOP
# ============================================================
cat("\n--- STEP 2: Simulation (S = 1,000 replications per n) ---\n")

n_levels <- c(50, 100, 200)
S <- 1000

results <- data.frame()
t_start <- Sys.time()

for (n in n_levels) {
  cat(sprintf("\n  n = %d (mimics %s in the real data):\n", n,
              ifelse(n == 50, "K2", ifelse(n == 100, "between K3 and K1", "K1/K3"))))

  for (s in seq_len(S)) {
    # --- simulate inputs (mean-centered logs, mimicking real pipeline convention) ---
    X <- matrix(rnorm(n * 4, mean = 0, sd = 0.5), nrow = n, ncol = 4)
    colnames(X) <- paste0("lnX", 1:4)

    u_true <- abs(rnorm(n, mean = 0, sd = sigma_u))   # half-normal inefficiency
    v_true <- rnorm(n, mean = 0, sd = sigma_v)        # two-sided noise
    TE_true <- exp(-u_true)

    y <- true_beta["beta0"] + X %*% true_beta[2:5] + v_true - u_true
    sim_df <- data.frame(y = as.numeric(y), X)

    fit <- tryCatch(
      sfa(y ~ lnX1 + lnX2 + lnX3 + lnX4, data = sim_df, ineffDecrease = TRUE),
      error = function(e) NULL
    )

    if (is.null(fit) || is.null(fit$mleParam)) {
      results <- rbind(results, data.frame(n = n, rep = s, converged = FALSE,
                                            gamma_hat = NA, boundary = NA, rho_true = NA))
      next
    }

    converged <- isTRUE(fit$convergence == 0) || isTRUE(fit$code %in% c(1, 2))
    gamma_hat <- tryCatch(unname(fit$mleParam["gamma"]), error = function(e) NA_real_)

    resid_hat <- tryCatch(residuals(fit, asInData = TRUE), error = function(e) NULL)
    if (is.null(resid_hat) || is.na(gamma_hat)) {
      results <- rbind(results, data.frame(n = n, rep = s, converged = converged,
                                            gamma_hat = gamma_hat, boundary = NA, rho_true = NA))
      next
    }

    sigma2_hat  <- tryCatch(unname(fit$mleParam["sigmaSq"]), error = function(e) NA_real_)
    sigma_u_hat <- sqrt(gamma_hat * sigma2_hat)
    sigma_v_hat <- sqrt((1 - gamma_hat) * sigma2_hat)

    TE_hat <- compute_TE_BC88_local(resid_hat, sigma_u_hat, sigma_v_hat)
    ok <- complete.cases(TE_hat, TE_true)
    rho_true <- if (sum(ok) >= 10) cor(TE_hat[ok], TE_true[ok], method = "spearman") else NA_real_

    results <- rbind(results, data.frame(
      n = n, rep = s, converged = converged,
      gamma_hat = gamma_hat, boundary = gamma_hat > 0.99, rho_true = rho_true
    ))
  }
  elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 2)
  cat(sprintf("    done | %.2f min elapsed (cumulative)\n", elapsed))
}

cat(sprintf("\n  Simulation complete. %d total fits attempted (%d per n level).\n",
            nrow(results), S))

# ============================================================
# STEP 3: SUMMARY TABLE
# ============================================================
cat("\n--- STEP 3: Summary by Sample Size ---\n")

summary_tbl <- do.call(rbind, lapply(n_levels, function(n_val) {
  sub <- results[results$n == n_val, ]
  conv_rate <- mean(sub$converged, na.rm = TRUE)
  sub_conv  <- sub[sub$converged %in% TRUE, ]
  data.frame(
    n                     = n_val,
    Convergence_Rate      = round(100 * conv_rate, 1),
    Pct_Boundary_Gamma    = round(100 * mean(sub_conv$boundary, na.rm = TRUE), 1),
    Mean_Gamma_Hat        = round(mean(sub_conv$gamma_hat, na.rm = TRUE), 4),
    Mean_Rho_vs_TrueTE    = round(mean(sub_conv$rho_true, na.rm = TRUE), 4),
    Median_Rho_vs_TrueTE  = round(median(sub_conv$rho_true, na.rm = TRUE), 4),
    SD_Rho_vs_TrueTE      = round(sd(sub_conv$rho_true, na.rm = TRUE), 4)
  )
}))
cat("\n  Monte Carlo summary (S = 1,000 replications per n; true gamma = 0.80):\n")
print(summary_tbl)

cat("\n  Interpretation: if Pct_Boundary_Gamma falls and Mean_Rho_vs_TrueTE rises\n")
cat("  monotonically as n increases from 50 to 200, this confirms that the\n")
cat("  fragility documented for the real K2 cluster (n=51) is a GENERAL property\n")
cat("  of small-sample SFA estimation, not a quirk specific to the Jambi tilapia data.\n")

# ============================================================
# STEP 4: PLOTS
# ============================================================
cat("\n--- STEP 4: Plots ---\n")

results$n_label <- factor(results$n, levels = n_levels,
                           labels = paste0("n=", n_levels, c(" (~K2)", "", " (~K1/K3)")))

p1 <- ggplot(results[results$converged %in% TRUE, ], aes(x = n_label, y = rho_true, fill = n_label)) +
  geom_boxplot(width = 0.45, alpha = 0.8, outlier.alpha = 0.25) +
  scale_fill_manual(values = c("#C0504D", "#ED9B40", "#2E75B6")) +
  labs(x = NULL, y = "Spearman rho: estimated TE vs. true simulated TE",
       caption = wrap_caption("S = 1,000 Monte Carlo replications per sample size; true gamma = 0.80; half-normal Cobb-Douglas SFA.")) +
  theme_j2_en() + theme(legend.position = "none")
save_png_j2(p1, "J2_08_MonteCarlo_Rho_by_SampleSize", w = 7, h = 5.5)

boundary_df <- summary_tbl[, c("n", "Pct_Boundary_Gamma")]
p2 <- ggplot(boundary_df, aes(x = factor(n), y = Pct_Boundary_Gamma)) +
  geom_col(fill = "#2E75B6", width = 0.55) +
  geom_text(aes(label = paste0(Pct_Boundary_Gamma, "%")), vjust = -0.5, size = 4) +
  labs(x = "Sample size (n)", y = "Replications with boundary gamma (%)",
       caption = wrap_caption("S = 1,000 Monte Carlo replications per sample size; boundary defined as gamma_hat > 0.99.")) +
  theme_j2_en()
save_png_j2(p2, "J2_08_MonteCarlo_PctBoundary_by_SampleSize", w = 6.5, h = 5)

# Third plot: outright non-convergence rate by n - emphasized in the manuscript text as a
# more severe failure mode than boundary-degenerate convergence alone, and therefore given
# its own figure rather than left as a table-only statistic.
conv_df <- summary_tbl[, c("n", "Convergence_Rate")]
conv_df$NonConvergence <- 100 - conv_df$Convergence_Rate
p3 <- ggplot(conv_df, aes(x = factor(n), y = NonConvergence)) +
  geom_col(fill = "#C0504D", width = 0.55) +
  geom_text(aes(label = paste0(round(NonConvergence, 1), "%")), vjust = -0.5, size = 4) +
  labs(x = "Sample size (n)", y = "Non-convergence rate (%)",
       caption = wrap_caption("S = 1,000 Monte Carlo replications per sample size; non-convergence = model fitting failed entirely.")) +
  theme_j2_en()
save_png_j2(p3, "J2_08_MonteCarlo_NonConvergence_by_SampleSize", w = 6.5, h = 5)

# ============================================================
# STEP 5: EXPORT
# ============================================================
all_tables <- list(
  MonteCarlo_Summary = list(
    df = summary_tbl,
    caption = "Table 14. Monte Carlo simulation: SFA estimation behaviour as a function of sample size (S = 1,000 replications; true gamma = 0.80)."
  )
)
finalize_outputs(all_tables, "J2_08_MonteCarlo_SmallSample_OptionB",
                  "Monte Carlo Simulation: Sample Size and SFA Estimation Fragility")

# ============================================================
# STEP 6: SAVE RDS
# ============================================================
saveRDS(list(results = results, summary_tbl = summary_tbl,
             true_beta = true_beta, sigma_u = sigma_u, sigma_v = sigma_v, true_gamma = true_gamma),
        file.path(DIR_RDS, "J2_08_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_08_Results.rds")))

cat("\n============================================================\n")
cat(" RINGKASAN AKHIR J2_08\n")
cat("============================================================\n")
print(summary_tbl)
cat("\nOK J2_08_MonteCarlo_SmallSample.R selesai.\n")
