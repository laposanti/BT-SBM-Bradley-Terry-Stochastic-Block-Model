# Bradley-Terry meets Stochastic Block Models

This repository contains code and data for the paper:  
**"Bradley-Terry meets Stochastic Block Models: Clustering Players from Pairwise Comparisons"**  
by Lapo Santi, Nial Friel [University College Dublin].

---

## 🔍 Overview

We propose a Bayesian Bradley-Terry model with block structure over players, combining discrete clustering with pairwise ranking data.  
This repo includes:

- an MCMC sampler implementing the full posterior inference,
- simulation studies,
- empirical application to ATP tennis data,
- and code to reproduce all figures and tables.

---

## 📁 Folder Structure

| Folder        | Description |
|---------------|-------------|
| `main/`       | MCMC sampler and core model implementation |
| `analysis/`   | Scripts to load posterior samples, produce plots, and replicate results in the paper |
| `simulations/`| Scripts to generate synthetic data and assess recovery of block structures |
| `images/`     | Precomputed plots for the paper (can be regenerated) |
| `data/`       | Raw input data (ranking histories, tournament results, etc.) |

---

## ▶️ How to Run

1. **Clone the repo:**

    ```bash
    git clone https://github.com/[your-username]/bradley-terry-sbm.git
    cd bradley-terry-sbm
    ```

2. **Install dependencies in R:**

    ```R
    install.packages(c("coda", "ggplot2", "tidygraph", "ggside", "loo"))
    ```

3. **Run the full pipeline:**

    This is a three-step process:

    **a) Load the core MCMC engine**

    ```R
    source("main/mcmc_bt_sbm.R")
    ```

    This defines the MCMC function `run_bt_sbm_mcmc()`, but **does not execute anything** when sourced.

    **b) Run the analysis on real data**

    ```R
    source("analysis/run_analysis.R")
    ```

    This driver script loads the ATP data, sets priors, runs the MCMC sampler, and saves the output to `results/mcmc_atp_results.rds`.

    **c) Generate figures and summaries**

    ```R
    source("analysis/postprocessing.R")
    ```

    This reads the saved output and produces all plots and tables used in the paper.

---

## 📊 Reproducing the Figures

The following table maps each figure in the paper to its corresponding output file and code section.

| Description                                                  | Script / Object                             | Output file |
|--------------------------------------------------------------|---------------------------------------------|-------------|
| Posterior adjacency matrix (block-ordered)                   | `postprocessing.R` / `geom_adjacency_fixed` | [`images/geom_adjacency_fixed.png`](./images/geom_adjacency_fixed.png) |
| Assignment-probabilities heatmap (players × clusters)        | `postprocessing.R` / `ass_prob_plot`        | [`images/ass_prob_plot.png`](./images/ass_prob_plot.png) |
| Player skill (λ) uncertainty — median + 90% HPD (log10)      | `postprocessing.R` / `plot_lambda`          | [`images/lambda_uncertainty.png`](./images/lambda_uncertainty.png) |
| P(Top block) by season — jittered points                     | `postprocessing.R` / `p_top_across_time`    | [`images/Ptop_across_time1.png`](./images/Ptop_across_time1.png) |
| Shannon entropy across seasons (mean with 90% band)          | `postprocessing.R` / `entropy_plot`         | [`images/entropy_plot1.png`](./images/entropy_plot1.png) |
| Nº of players in top block by season (bar chart)             | `postprocessing.R` / `num_block_plot`       | [`images/num_block_plot1.png`](./images/num_block_plot1.png) |

> All outputs are saved to the `images/` folder unless otherwise noted.  
> You can customize the output location by modifying the save paths in `postprocessing.R`.

---
