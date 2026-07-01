# ============================================================
# J2_01_DEA_SFA_Convergence_RQ1.R
# RQ1: DOES SFA CONVERGE WITH THE DEA BOOTSTRAP META-FRONTIER?
# ============================================================
# Journal-2 pipeline (rebuilt, clean) - Ma'in (NIM P3C123010), FEB UNJA
# Step 1 of 4. Run this FIRST.
#
# Fits SFA (Cobb-Douglas, half-normal) three ways - per cluster (K1, K2,
# K3), pooled, and compares each against the DEA bootstrap meta-frontier
# TE (TE_meta_bc). Runs a B=500 bootstrap to test whether per-cluster
# estimates (especially K2, n=51) are numerically stable.
#
# OUTPUT (all under Journal_2_rev/):
#   xlsx/J2_01_RQ1_Convergence.xlsx
#   tex/J2_01_RQ1_Convergence.tex
#   word/J2_01_RQ1_Convergence.docx
#   plots/J2_01_*.png
#   rds/J2_01_Results.rds   <- read by J2_02
# ============================================================

source("00_setup.R")
source("J2_00_Helpers.R")

if (!requireNamespace("frontier", quietly = TRUE)) {
  install.packages("frontier", repos = "https://cran.r-project.org")
}
suppressPackageStartupMessages(library(frontier))

set.seed(2026)

cat("\n============================================================\n")
cat(" J2_01: RQ1 - DEA-SFA CONVERGENCE (PER-CLUSTER, POOLED, BOOTSTRAP)\n")
cat("============================================================\n")

# ============================================================
# STEP 1: LOAD DATA
# ============================================================
cat("\n--- STEP 1: Load Data ---\n")

data_clean <- as.data.frame(readRDS("output/rds/data_clean.rds"), stringsAsFactors = FALSE)
d_te <- as.data.frame(readRDS("output/rds/02_TE_results.rds"), stringsAsFactors = FALSE)

data_full <- merge(data_clean, d_te[, c("DMU", "TE_meta_bc", "TE_group_bc")], by = "DMU", all.x = TRUE)

data_full$ln_Y  <- log(pmax(data_full$Y_OUT,   0.001))
data_full$ln_X1 <- log(pmax(data_full$X1_SIZE, 0.001))
data_full$ln_X2 <- log(pmax(data_full$X2_SEED, 0.001))
data_full$ln_X3 <- log(pmax(data_full$X3_FEED, 0.001))
data_full$ln_X6 <- log(pmax(data_full$X6_LABR, 0.001))

cols_needed <- c("DMU", "klaster", "WAD", "Density", "FCR", "SR", "Laborint",
                  "ln_Y", "ln_X1", "ln_X2", "ln_X3", "ln_X6", "TE_meta_bc", "TE_group_bc")
data_sfa <- data_full[complete.cases(data_full[, cols_needed]), ]
data_sfa$Cluster <- factor(data_sfa$klaster, 1:3,
                             c("K1 (0-100 m a.s.l.)", "K2 (>100-500 m a.s.l.)", "K3 (>500 m a.s.l.)"))
data_sfa$TE_DEA <- data_sfa$TE_meta_bc

n_total <- nrow(data_sfa)
n_k1 <- sum(data_sfa$klaster == 1); n_k2 <- sum(data_sfa$klaster == 2); n_k3 <- sum(data_sfa$klaster == 3)
cat(sprintf("  n total=%d | K1=%d | K2=%d | K3=%d\n", n_total, n_k1, n_k2, n_k3))

# ============================================================
# STEP 2: FIT SFA - PER CLUSTER + POOLED
# ============================================================
cat("\n--- STEP 2: Fit SFA (half-normal Cobb-Douglas) ---\n")

fit_sfa_safe <- function(d, label) {
  cat(sprintf("  Fitting: %s (n=%d)...\n", label, nrow(d)))
  out <- list(label = label, n = nrow(d), converged = FALSE,
              gamma = NA_real_, loglik = NA_real_, boundary_issue = FALSE,
              coef = NULL, TE = NULL)
  mod <- tryCatch(sfa(ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6, data = d, ineffDecrease = TRUE),
                   error = function(e) { cat(sprintf("    FAILED: %s\n", e$message)); NULL })
  if (is.null(mod)) return(out)
  smry <- tryCatch(summary(mod), error = function(e) NULL)
  if (is.null(smry)) return(out)
  mle <- as.data.frame(smry$mleParam); mle$Parameter <- rownames(mle)
  gamma_val <- mle[grepl("^gamma$", mle$Parameter, ignore.case = TRUE), "Estimate"]
  gamma_val <- if (length(gamma_val) == 1) gamma_val else NA_real_
  sigmaSq_val <- mle[grepl("sigmaSq", mle$Parameter, ignore.case = TRUE), "Estimate"]
  sigmaSq_val <- if (length(sigmaSq_val) == 1) sigmaSq_val else NA_real_
  out$converged <- TRUE
  out$gamma <- gamma_val
  out$loglik <- tryCatch(as.numeric(logLik(mod)), error = function(e) NA_real_)
  out$boundary_issue <- !is.na(gamma_val) && (gamma_val < 0.05 || gamma_val > 0.95)
  out$coef <- mle

  # Primary: package efficiencies(). NEVER fail silently - print the
  # actual error/mismatch so the cause is always visible in the log.
  te_pkg <- tryCatch(as.numeric(efficiencies(mod)), error = function(e) {
    cat(sprintf("    NOTE: efficiencies() raised an error (%s).\n", e$message)); NULL
  })
  if (is.null(te_pkg) || length(te_pkg) != nrow(d)) {
    if (!is.null(te_pkg)) {
      cat(sprintf("    NOTE: efficiencies() returned length %d, expected %d.\n", length(te_pkg), nrow(d)))
    }
    cat("    Falling back to manual Battese & Coelli (1988) TE estimator.\n")
    X_mat <- model.matrix(ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6, data = d)
    beta_named <- setNames(mle$Estimate, mle$Parameter)
    beta_match <- beta_named[colnames(X_mat)]
    fitted_vals <- as.numeric(X_mat %*% beta_match)
    e_resid <- d$ln_Y - fitted_vals
    out$TE <- compute_TE_BC88(e_resid, sigmaSq_val, gamma_val)
  } else {
    out$TE <- te_pkg
  }

  cat(sprintf("    OK. gamma=%.4f | logLik=%.2f | boundary=%s | TE n=%d/%d\n",
              gamma_val, out$loglik, ifelse(out$boundary_issue, "YES", "no"),
              sum(!is.na(out$TE)), nrow(d)))
  out
}

fit_k1     <- fit_sfa_safe(data_sfa[data_sfa$klaster == 1, ], "K1")
fit_k2     <- fit_sfa_safe(data_sfa[data_sfa$klaster == 2, ], "K2")
fit_k3     <- fit_sfa_safe(data_sfa[data_sfa$klaster == 3, ], "K3")
fit_pooled <- fit_sfa_safe(data_sfa, "Pooled")
fits_all <- list(K1 = fit_k1, K2 = fit_k2, K3 = fit_k3, Pooled = fit_pooled)

convergence_tbl <- do.call(rbind, lapply(names(fits_all), function(nm) {
  f <- fits_all[[nm]]
  data.frame(Specification = nm, n = f$n, Converged = f$converged,
             Gamma = round(f$gamma, 4), LogLik = round(f$loglik, 2),
             Boundary_Issue = f$boundary_issue)
}))
cat("\n  Convergence summary:\n"); print(convergence_tbl)

# ============================================================
# STEP 3: BOOTSTRAP STABILITY (B=500, per cluster)
# ============================================================
cat("\n--- STEP 3: Bootstrap Stability (B=500 per cluster) ---\n")

run_bootstrap_stability <- function(d, label, B = 500) {
  cat(sprintf("  Bootstrapping %s (n=%d, B=%d)...\n", label, nrow(d), B))
  n <- nrow(d); res <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    mod_b <- tryCatch(suppressWarnings(sfa(ln_Y ~ ln_X1 + ln_X2 + ln_X3 + ln_X6,
                                            data = d[idx, ], ineffDecrease = TRUE)),
                       error = function(e) NULL)
    if (is.null(mod_b)) { res[[b]] <- NULL; next }
    smry_b <- tryCatch(summary(mod_b), error = function(e) NULL)
    if (is.null(smry_b)) { res[[b]] <- NULL; next }
    mle_b <- as.data.frame(smry_b$mleParam); mle_b$Parameter <- rownames(mle_b)
    gamma_b <- mle_b[grepl("^gamma$", mle_b$Parameter, ignore.case = TRUE), "Estimate"]
    if (length(gamma_b) != 1 || is.na(gamma_b)) { res[[b]] <- NULL; next }
    res[[b]] <- gamma_b
  }
  gamma_vec <- unlist(res[!vapply(res, is.null, logical(1))])
  conv_rate <- length(gamma_vec) / B
  pct_boundary <- mean(gamma_vec > 0.95 | gamma_vec < 0.05, na.rm = TRUE) * 100
  cv_gamma <- sd(gamma_vec, na.rm = TRUE) / mean(gamma_vec, na.rm = TRUE)
  cat(sprintf("    Convergence rate=%.1f%% | mean gamma=%.4f | pct at boundary=%.1f%% | CV=%.4f\n",
              conv_rate * 100, mean(gamma_vec, na.rm = TRUE), pct_boundary, abs(cv_gamma)))
  list(label = label, gamma_vec = gamma_vec, conv_rate = conv_rate,
       pct_boundary = pct_boundary, cv_gamma = abs(cv_gamma))
}

boot_k1 <- run_bootstrap_stability(data_sfa[data_sfa$klaster == 1, ], "K1")
boot_k2 <- run_bootstrap_stability(data_sfa[data_sfa$klaster == 2, ], "K2")
boot_k3 <- run_bootstrap_stability(data_sfa[data_sfa$klaster == 3, ], "K3")
boots_all <- list(K1 = boot_k1, K2 = boot_k2, K3 = boot_k3)

stability_tbl <- do.call(rbind, lapply(names(boots_all), function(nm) {
  b <- boots_all[[nm]]
  data.frame(Cluster = nm, n = c(n_k1, n_k2, n_k3)[match(nm, c("K1", "K2", "K3"))],
             Bootstrap_Convergence_Rate_Pct = round(b$conv_rate * 100, 1),
             Mean_Gamma = round(mean(b$gamma_vec, na.rm = TRUE), 4),
             Pct_At_Boundary = round(b$pct_boundary, 1),
             CV_Gamma = round(b$cv_gamma, 4))
}))
cat("\n  Bootstrap stability summary:\n"); print(stability_tbl)

# ============================================================
# STEP 4: RANK CORRELATION (SPEARMAN) vs DEA
# ============================================================
cat("\n--- STEP 4: Spearman Rank Correlation vs DEA ---\n")

data_sfa$TE_SFA_cluster <- NA_real_
check_and_assign <- function(data_sfa, fit_obj, klaster_id) {
  if (!fit_obj$converged) {
    cat(sprintf("  NOTE: K%d model did not converge - TE left as NA.\n", klaster_id))
    return(data_sfa)
  }
  n_expected <- sum(data_sfa$klaster == klaster_id)
  if (is.null(fit_obj$TE) || length(fit_obj$TE) != n_expected) {
    stop(sprintf(
      "K%d: TE vector length (%s) does not match number of K%d observations (%d). Check fit_sfa_safe().",
      klaster_id, ifelse(is.null(fit_obj$TE), "NULL", length(fit_obj$TE)), klaster_id, n_expected))
  }
  data_sfa$TE_SFA_cluster[data_sfa$klaster == klaster_id] <- fit_obj$TE
  data_sfa
}
data_sfa <- check_and_assign(data_sfa, fit_k1, 1)
data_sfa <- check_and_assign(data_sfa, fit_k2, 2)
data_sfa <- check_and_assign(data_sfa, fit_k3, 3)
data_sfa$TE_SFA_pooled <- if (fit_pooled$converged) fit_pooled$TE else NA_real_

compare_rho <- function(x, y, label) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 5) return(data.frame(Comparison = label, n = sum(ok), Rho = NA_real_, p_value = NA_real_))
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  data.frame(Comparison = label, n = sum(ok), Rho = round(unname(ct$estimate), 4),
             p_value = round(ct$p.value, 6))
}

rho_tbl <- rbind(
  compare_rho(data_sfa$TE_SFA_cluster, data_sfa$TE_DEA, "Overall: SFA per-cluster vs DEA"),
  compare_rho(data_sfa$TE_SFA_pooled,  data_sfa$TE_DEA, "Overall: SFA pooled vs DEA")
)
for (k in c(1, 2, 3)) {
  d_k <- data_sfa[data_sfa$klaster == k, ]
  rho_tbl <- rbind(rho_tbl,
    compare_rho(d_k$TE_SFA_cluster, d_k$TE_DEA, sprintf("K%d: SFA per-cluster vs DEA", k)),
    compare_rho(d_k$TE_SFA_pooled,  d_k$TE_DEA, sprintf("K%d: SFA pooled vs DEA", k)))
}
cat("\n  Spearman rank correlation:\n"); print(rho_tbl)

# ============================================================
# STEP 5: PLOTS (English, no titles, black footnote)
# ============================================================
cat("\n--- STEP 5: Plots ---\n")

COL_K <- c("K1 (0-100 m a.s.l.)" = "#2E75B6", "K2 (>100-500 m a.s.l.)" = "#E69138",
           "K3 (>500 m a.s.l.)" = "#1A6B3C")

gamma_boot_df <- do.call(rbind, lapply(names(boots_all), function(nm) {
  data.frame(Cluster = nm, Gamma = boots_all[[nm]]$gamma_vec)
}))

p1 <- ggplot(gamma_boot_df, aes(x = Cluster, y = Gamma, fill = Cluster)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.6) +
  geom_hline(yintercept = c(0.05, 0.95), linetype = "dashed", color = "#000000", linewidth = 0.4) +
  scale_fill_manual(values = c("K1" = "#2E75B6", "K2" = "#E69138", "K3" = "#1A6B3C"), guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Elevation cluster", y = "Bootstrap estimate of gamma",
       caption = wrap_caption("Half-normal SFA, B = 500 within-cluster bootstrap replications; dashed lines mark the boundary of the parameter space (0.05, 0.95)")) +
  theme_j2_en()
save_png_j2(p1, "J2_01_Bootstrap_Gamma_Distribution", w = 7, h = 5.2)

p2 <- ggplot(data_sfa, aes(TE_SFA_cluster, TE_DEA, color = Cluster)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  facet_wrap(~Cluster, nrow = 1) +
  scale_color_manual(values = COL_K, guide = "none") +
  labs(x = "Technical efficiency, SFA (cluster-specific)", y = "Technical efficiency, DEA (meta-frontier)",
       caption = wrap_caption("Half-normal SFA per elevation cluster vs. bootstrap bias-corrected DEA meta-frontier; dashed line is the 45-degree reference.")) +
  theme_j2_en()
save_png_j2(p2, "J2_01_Scatter_ClusterSFA_vs_DEA", w = 11, h = 4.5)

rho_plot_df <- rho_tbl[grepl("^K[123]:", rho_tbl$Comparison), ]
rho_plot_df$Cluster <- substr(rho_plot_df$Comparison, 1, 2)
rho_plot_df$Specification <- ifelse(grepl("per-cluster", rho_plot_df$Comparison), "Per-cluster SFA", "Pooled SFA")

p3 <- ggplot(rho_plot_df, aes(x = Cluster, y = Rho, fill = Specification)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Rho)), position = position_dodge(width = 0.6),
            vjust = -0.4, size = 3.2, color = "#000000") +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "#000000", linewidth = 0.4) +
  scale_fill_manual(values = c("Per-cluster SFA" = "#1A6B3C", "Pooled SFA" = "#8C8C8C")) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(x = "Elevation cluster", y = "Spearman rho vs. DEA meta-frontier TE",
       caption = wrap_caption("Spearman rank correlation; dashed line marks the conventional convergence threshold (rho = 0.70).")) +
  theme_j2_en()
save_png_j2(p3, "J2_01_Barplot_Rho_PerCluster_vs_Pooled", w = 7.5, h = 5.2)

# ============================================================
# STEP 6: EXPORT (xlsx + tex + docx via shared helper)
# ============================================================
all_tables <- list(
  Convergence_Summary = list(df = convergence_tbl,
    caption = "Table 1. Convergence diagnostics of half-normal SFA (per cluster and pooled)."),
  Bootstrap_Stability = list(df = stability_tbl,
    caption = "Table 2. Bootstrap stability of the variance-ratio parameter (B = 500)."),
  Rank_Correlation     = list(df = rho_tbl,
    caption = "Table 3. Spearman rank correlation between SFA and DEA meta-frontier technical efficiency.")
)
finalize_outputs(all_tables, "J2_01_RQ1_Convergence",
                  "RQ1 - DEA-SFA Convergence: Per-Cluster, Pooled, and Bootstrap Stability")

# ============================================================
# STEP 7: SAVE RDS FOR J2_02
# ============================================================
saveRDS(list(data_sfa = data_sfa, fits_all = fits_all, boots_all = boots_all,
             convergence_tbl = convergence_tbl, stability_tbl = stability_tbl, rho_tbl = rho_tbl),
        file.path(DIR_RDS, "J2_01_Results.rds"))
cat(sprintf("\n  OK  %s\n", file.path(DIR_RDS, "J2_01_Results.rds")))

cat("\nOK J2_01_DEA_SFA_Convergence_RQ1.R complete.\n")
