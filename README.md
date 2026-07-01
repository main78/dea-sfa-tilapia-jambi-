# dea-sfa-tilapia-jambi-
# When and Why Do DEA and SFA Diverge? Evidence from Spatially Heterogeneous Tilapia Aquaculture

This repository contains the R scripts used to produce all tables and figures reported in the manuscript *"When and Why Do DEA and SFA Diverge? Evidence from Spatially Heterogeneous Tilapia Aquaculture"*, submitted to the *Journal of Productivity Analysis*.

## What this repository contains

- **Scripts only.** Farm-level survey data are not included, since they contain potentially identifying information about individual smallholder respondents (see the manuscript's Data Availability statement for how to request de-identified data).
- All scripts are written in R and depend on the [`frontier`](https://CRAN.R-project.org/package=frontier) package for Stochastic Frontier Analysis estimation, plus standard `tidyverse`/`ggplot2` packages for data handling and visualization.

## Script-to-section mapping

| Script | Research Question | Manuscript Section(s) | Output |
|---|---|---|---|
| `J2_00_Helpers.R` | — (shared utilities) | — | Theme, plotting, table-export, and caption-wrapping helper functions sourced by every other script |
| `J2_01_DEA_SFA_Convergence_RQ1.R` | RQ1 | Sections 4.1–4.2 | Convergence diagnostics, bootstrap stability, Spearman rho (Tables 1–3, Figs. 1–3) |
| `J2_02_Model_Specification_RQ2.R` | RQ2 | Section 4.3 | Model A / Model B comparison (Tables 4–5, Figs. 4–5) |
| `J2_03_BC95_Inefficiency_Effects_RQ3.R` | RQ3 | Section 4.4 | Battese & Coelli (1995) bootstrap inference, DEA vs. SFA-BC95 comparison (Tables 6–7, Figs. 6–7) |
| `J2_04_Divergence_Predictors_RQ4.R` | RQ4 | Section 4.5 | OLS divergence-predictor model, robustness check (Tables 8–10, Figs. 8–10) |
| `J2_05_Check_Divergence_By_Cluster.R` | RQ4 (supplementary) | Section 4.5 | Mean absolute divergence by elevation cluster (Table 11) |
| `J2_06_Paired_Bootstrap_Compare_Rho.R` | RQ1/RQ2 (robustness) | Section 4.3 | Formal paired-bootstrap test of pairwise rho differences (Table 12) |
| `J2_07_CrossValidation_Divergence_Model.R` | RQ4 (robustness) | Section 4.5 | Repeated 10-fold cross-validation of the divergence model (Table 13) |
| `J2_08_MonteCarlo_SmallSample.R` | RQ1 (robustness) | Sections 3.10, 4.6 | Monte Carlo simulation confirming sample-size fragility is a general estimator property (Table 14, Figs. 11–13) |

## Reproducing the analysis

Scripts are designed to be run in numerical order (`J2_01` through `J2_08`); each later script reads the `.rds` output of an earlier one rather than re-fitting models from scratch. `J2_08` is fully self-contained and synthetic (no survey data required).

```r
# Example: from the repository root, with working directory set accordingly
source("J2_00_Helpers.R")
source("J2_01_DEA_SFA_Convergence_RQ1.R")
source("J2_02_Model_Specification_RQ2.R")
# ... etc.
```

## Software requirements

- R (>= 4.2)
- R packages: `frontier`, `ggplot2`, `dplyr`, `openxlsx`, `officer`, `stargazer` or equivalent (see individual script headers for exact dependencies)

## Citation

If you use this code, please cite the manuscript:

> [Author(s)]. (forthcoming). When and Why Do DEA and SFA Diverge? Evidence from Spatially Heterogeneous Tilapia Aquaculture. *Journal of Productivity Analysis*.

## License

This repository is released under the MIT License (see `LICENSE`). The underlying survey data are not covered by this license and remain subject to the data-sharing terms described in the manuscript.

## Contact

Questions about this code should be directed to the corresponding author (see manuscript title page).
