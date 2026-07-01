# ============================================================
# J2_03_BC95_Inefficiency_Effects_RQ3.R
# RQ3: ARE INEFFICIENCY DETERMINANTS CONSISTENT ACROSS DEA AND SFA?
# ============================================================
# Journal-2 pipeline (rebuilt, clean) - Ma'in (NIM P3C123010), FEB UNJA
# Step 3 of 4. Run AFTER J2_02.
#
# Fits a Battese & Coelli (1995) inefficiency-effects model:
#   ln_Y ~ ln_X1+ln_X2+ln_X3+ln_X6  |  Density_z+FCR_z+SR_z+Laborint_z
# using the simplest pooled frontier (parsimony favoured per RQ2; adding
# cluster dummies here was verified in pilot diagnostics to leave the
# inefficiency-effects coefficients materially unchanged).
#
# NOTE ON INFERENCE: analytic (Hessian-based) standard errors for this
# class of model are well known to be unreliable in small/unbalanced
# panels and were confirmed unreliable here in pilot diagnostics (a
# near-uniform spread of standard errors close to 1.0 across ALL
# parameters, including otherwise precisely identified input
# elasticities). Inference is therefore based on B=300 non-parametric
# bootstrap resampling throughout, not on the package's analytic SEs.
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/word/tex: J2_03_RQ3_BC95_Inefficiency_Effects
#   plots/J2_03_*.png
#   rds/J2_03_Results.rds   <- read by downstream manuscript assembly (optional)
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

if (!requireNamespace("frontier", quietly = TRUE)) {
  install.packages("frontier", repos = "https://cran.r-project.org")
}
suppressPackageStartupMessages(library(frontier))
suppressPackageStartupMessages(library(Formula))
suppressPackageStartupMessages(library(tidyr))

set.seed(2026)

# Print warnings immediately as they occur, rather than deferring them to
# a batch summary at the end of the script. With B=300 bootstrap fits, many
# of which legitimately warn about gamma near the parameter boundary, the
# deferred end-of-script warning printer can itself crash (an internal R/
# package printing bug, unrelated to this script's logic); printing warnings
# immediately avoids that code path entirely.
old_warn_opt <- options(warn = 1)

cat("\n============================================================\n")
cat(" J2_03: RQ3 - BATTESE & COELLI (1995) INEFFICIENCY-EFFECTS MODEL\n")
cat("============================================================\n")

# ============================================================
# STEP 1: LOAD DATA + DETERMINANTS
# ============================================================
cat("\n--- STEP 1: Load Data ---\n")

j2_02 <- readRDS(file.path(DIR_RDS, "J2_02_Results.rds"))
data_sfa <- j2_02$data_sfa

stopifnot(all(c("Density", "FCR", "SR", "Laborint") %in% names(data_sfa)))

data_sfa$Density_z  <- as.numeric(scale(data_sfa$Density))
data_sfa$FCR_z      <- as.numeric(scale(data_sfa$FCR))
data_sfa$SR_z       <- as.numeric(scale(data_sfa$SR))
data_sfa$Laborint_z <- as.numeric(scale(data_sfa$Laborint))

cat(sprintf("  OK. n=%d.\n", nrow(data_sfa)))

formula_bc95 <- Formula::as.Formula(
  ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6 | Density_z + FCR_z + SR_z + Laborint_z
)

get_est <- function(mle_df, varname) {
  row <- mle_df[grepl(paste0(varname, "$"), mle_df$Parameter), ]
  if (nrow(row) != 1) return(NA_real_)
  row$Estimate
}

# ============================================================
# STEP 2: POINT ESTIMATE (FULL SAMPLE)
# ============================================================
cat("\n--- STEP 2: Point Estimate (Full Sample) ---\n")

fit_point <- suppressWarnings(tryCatch(
  sfa(formula_bc95, data = data_sfa, ineffDecrease = TRUE, maxit = 5000),
  error = function(e) NULL
))
if (is.null(fit_point)) stop("Point-estimate BC95 model failed to fit.")

mle_point <- as.data.frame(suppressWarnings(summary(fit_point))$mleParam); mle_point$Parameter <- rownames(mle_point)

point_est <- c(
  Density  = get_est(mle_point, "Z_Density_z"),
  FCR      = get_est(mle_point, "Z_FCR_z"),
  SR       = get_est(mle_point, "Z_SR_z"),
  Laborint = get_est(mle_point, "Z_Laborint_z")
)
point_gamma <- get_est(mle_point, "^gamma$")
cat(sprintf("  Point estimates (effect on u): Density=%.4f, FCR=%.4f, SR=%.4f, Laborint=%.4f, gamma=%.4f\n",
            point_est["Density"], point_est["FCR"], point_est["SR"], point_est["Laborint"], point_gamma))

# ============================================================
# STEP 3: BOOTSTRAP RESAMPLING (B=300)
# ============================================================
cat("\n--- STEP 3: Bootstrap Resampling (B=300) ---\n")

B <- 300
n <- nrow(data_sfa)
boot_results <- vector("list", B)
t_start <- Sys.time()

for (b in seq_len(B)) {
  d_b <- data_sfa[sample.int(n, n, replace = TRUE), ]
  mod_b <- suppressWarnings(tryCatch(sfa(formula_bc95, data = d_b, ineffDecrease = TRUE, maxit = 3000),
                                       error = function(e) NULL))
  if (is.null(mod_b)) { boot_results[[b]] <- NULL; next }
  smry_b <- tryCatch(suppressWarnings(summary(mod_b)), error = function(e) NULL)
  if (is.null(smry_b)) { boot_results[[b]] <- NULL; next }
  mle_b <- as.data.frame(smry_b$mleParam); mle_b$Parameter <- rownames(mle_b)
  gamma_b <- get_est(mle_b, "^gamma$")
  if (is.na(gamma_b)) { boot_results[[b]] <- NULL; next }
  boot_results[[b]] <- c(Density = get_est(mle_b, "Z_Density_z"), FCR = get_est(mle_b, "Z_FCR_z"),
                          SR = get_est(mle_b, "Z_SR_z"), Laborint = get_est(mle_b, "Z_Laborint_z"),
                          gamma = gamma_b)
  if (b %% 50 == 0) {
    elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
    cat(sprintf("  ... replication %d/%d done | %.1f min elapsed\n", b, B, elapsed))
  }
}

boot_valid <- boot_results[!vapply(boot_results, is.null, logical(1))]
conv_rate <- length(boot_valid) / B
cat(sprintf("\n  Bootstrap convergence rate: %d/%d (%.1f%%)\n", length(boot_valid), B, conv_rate * 100))
if (length(boot_valid) < 30) stop("Too few valid bootstrap replications (<30).")

boot_df <- do.call(rbind, lapply(boot_valid, function(x) as.data.frame(t(x))))
gamma_pct_boundary <- mean(boot_df$gamma > 0.95 | boot_df$gamma < 0.05, na.rm = TRUE) * 100
cat(sprintf("  Pct. bootstrap replications with gamma at boundary: %.1f%%\n", gamma_pct_boundary))

# ============================================================
# STEP 4: BOOTSTRAP INFERENCE SUMMARY
# ============================================================
cat("\n--- STEP 4: Bootstrap Inference Summary ---\n")

bootstrap_p <- function(boot_vec, point) {
  p1 <- if (point >= 0) mean(boot_vec < 0, na.rm = TRUE) else mean(boot_vec > 0, na.rm = TRUE)
  min(1, 2 * p1)
}

build_boot_row <- function(varname) {
  vec <- boot_df[[varname]]
  pe  <- point_est[varname]
  data.frame(Variable = varname, Point_Estimate = round(pe, 4),
             Bootstrap_SE = round(sd(vec, na.rm = TRUE), 4),
             CI_Lower_95 = round(quantile(vec, 0.025, na.rm = TRUE), 4),
             CI_Upper_95 = round(quantile(vec, 0.975, na.rm = TRUE), 4),
             Bootstrap_p_value = round(bootstrap_p(vec, pe), 4))
}
bc95_boot_tbl <- rbind(build_boot_row("Density"), build_boot_row("FCR"),
                         build_boot_row("SR"), build_boot_row("Laborint"))
bc95_boot_tbl$Significance <- sapply(bc95_boot_tbl$Bootstrap_p_value, sig_label)
rownames(bc95_boot_tbl) <- NULL
cat("\n  Bootstrap inference (effect on u):\n"); print(bc95_boot_tbl)

# ============================================================
# STEP 5: COMPARISON WITH DEA TWO-STAGE TRUNCATED REGRESSION (M1)
# ============================================================
cat("\n--- STEP 5: Comparison with DEA M1 (Tabel 5.30, dissertation) ---\n")

dea_m1 <- data.frame(
  Variable = c("Density", "FCR", "SR", "Laborint"),
  Beta_DEA_TE = c(0.0011, -0.7020, 0.0019, 0.0017),
  SE_DEA_TE   = c(0.0004, 0.0313, 0.0004, 0.0002),
  p_DEA_TE    = c(0.002, 0.000, 0.000, 0.000)
)
dea_m1$Significance_DEA <- sapply(dea_m1$p_DEA_TE, sig_label)
dea_m1$Sign_DEA <- ifelse(dea_m1$Beta_DEA_TE > 0, "+", "-")

sfa_final <- bc95_boot_tbl[, c("Variable", "Point_Estimate", "Bootstrap_SE", "Bootstrap_p_value", "Significance")]
names(sfa_final)[names(sfa_final) == "Significance"] <- "Significance_SFA"
sfa_final$Implied_Effect_on_TE <- -sfa_final$Point_Estimate
sfa_final$Sign_SFA_on_TE <- ifelse(sfa_final$Implied_Effect_on_TE > 0, "+", "-")

comparison_tbl <- merge(dea_m1, sfa_final, by = "Variable")
comparison_tbl$Direction_Concordance <- ifelse(comparison_tbl$Sign_DEA == comparison_tbl$Sign_SFA_on_TE,
                                                "Concordant", "Discordant")
comparison_tbl$Significance_Concordance <- ifelse(
  (comparison_tbl$Significance_DEA != "ns") == (comparison_tbl$Significance_SFA != "ns"),
  "Concordant", "Discordant"
)
comparison_tbl <- comparison_tbl[, c("Variable", "Beta_DEA_TE", "SE_DEA_TE", "Significance_DEA", "Sign_DEA",
                                       "Point_Estimate", "Bootstrap_SE", "Significance_SFA",
                                       "Implied_Effect_on_TE", "Sign_SFA_on_TE",
                                       "Direction_Concordance", "Significance_Concordance")]
names(comparison_tbl) <- c("Variable", "DEA_Beta_on_TE", "DEA_SE", "DEA_Significance", "DEA_Sign",
                             "SFA_Delta_on_u", "SFA_Bootstrap_SE", "SFA_Significance",
                             "SFA_Implied_Effect_on_TE", "SFA_Sign_on_TE",
                             "Direction_Concordance", "Significance_Concordance")
cat("\n  DEA (two-stage) vs. SFA-BC95 (bootstrap), comparison table:\n"); print(comparison_tbl)

n_concordant <- sum(comparison_tbl$Direction_Concordance == "Concordant")
fcr_row <- comparison_tbl[comparison_tbl$Variable == "FCR", ]
cat(sprintf("\n  Direction concordance: %d/4 variables. FCR paradox replicated: %s (SFA significance: %s)\n",
            n_concordant, fcr_row$Direction_Concordance, fcr_row$SFA_Significance))

# ============================================================
# STEP 6: PLOTS
# ============================================================
cat("\n--- STEP 6: Plots ---\n")

boot_long <- boot_df[, c("Density", "FCR", "SR", "Laborint")] %>%
  tidyr::pivot_longer(everything(), names_to = "Variable", values_to = "Delta")
point_df <- data.frame(Variable = names(point_est), Point = as.numeric(point_est))

p1 <- ggplot(boot_long, aes(x = Delta)) +
  geom_histogram(fill = "#2E75B6", alpha = 0.75, bins = 40, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  geom_vline(data = point_df, aes(xintercept = Point), color = "#C0504D", linewidth = 0.8) +
  facet_wrap(~Variable, scales = "free", nrow = 1) +
  labs(x = "Bootstrap estimate of delta (effect on u)", y = "Frequency",
       caption = wrap_caption("B = 300 non-parametric bootstrap replications; vertical red line marks the full-sample point estimate.")) +
  theme_j2_en()
save_png_j2(p1, "J2_03_Bootstrap_Delta_Histograms", w = 12, h = 4.2)

p2 <- ggplot(comparison_tbl, aes(x = Variable)) +
  geom_hline(yintercept = 0, color = "#000000", linewidth = 0.4) +
  geom_point(aes(y = DEA_Beta_on_TE, color = "DEA (two-stage)"), size = 2.8,
             position = position_nudge(x = -0.15)) +
  geom_point(aes(y = SFA_Implied_Effect_on_TE, color = "SFA-BC95 (bootstrap)"), size = 2.8,
             position = position_nudge(x = 0.15)) +
  scale_color_manual(values = c("DEA (two-stage)" = "#2E75B6", "SFA-BC95 (bootstrap)" = "#1A6B3C")) +
  labs(x = NULL, y = "Coefficient / implied effect on technical efficiency", color = NULL,
       caption = wrap_caption("Comparison of direction, not magnitude (the two models scale a different dependent variable).")) +
  theme_j2_en() + theme(legend.position = "bottom")
save_png_j2(p2, "J2_03_Direction_Comparison_DEA_vs_SFA", w = 8, h = 5.5)

# ============================================================
# STEP 7: EXPORT
# ============================================================
all_tables <- list(
  Bootstrap_Inference = list(df = bc95_boot_tbl,
    caption = "Table 6. Bootstrap inference for Battese and Coelli (1995) inefficiency-effects coefficients (B = 300)."),
  DEA_vs_SFA_Comparison = list(df = comparison_tbl,
    caption = "Table 7. Direction and significance concordance between DEA two-stage and SFA-BC95 inefficiency determinants.")
)
finalize_outputs(all_tables, "J2_03_RQ3_BC95_Inefficiency_Effects",
                  "RQ3 - Consistency of Inefficiency Determinants: DEA vs. SFA-BC95")

# ============================================================
# STEP 8: SAVE RDS
# ============================================================
saveRDS(list(data_sfa = data_sfa, boot_df = boot_df, bc95_boot_tbl = bc95_boot_tbl,
             comparison_tbl = comparison_tbl, conv_rate = conv_rate,
             gamma_pct_boundary = gamma_pct_boundary),
        file.path(DIR_RDS, "J2_03_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_03_Results.rds")))

cat("\nOK J2_03_BC95_Inefficiency_Effects_RQ3.R complete.\n")

# Clean up complex model objects (frontier-class S4 objects with nested
# formula/call structures) before returning control to the console. This
# is believed to prevent RStudio's Environment-pane object preview from
# triggering an internal print-formatting bug (Error in is.name(callee) ...)
# unrelated to this script's own logic - all outputs above are unaffected
# either way, since they are saved before this point.
suppressWarnings(rm(fit_point, mod_b, smry_b))

invisible(options(old_warn_opt))
