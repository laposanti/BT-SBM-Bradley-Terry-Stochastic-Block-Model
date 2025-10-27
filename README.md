# **"Bradley–Terry meets Stochastic Block Models: Clustering Players from Pairwise Comparisons"**  

Lapo Santi, Nial Friel — University College Dublin

This repo reproduces the results of the paper and relies on the **[`BTSBM`]([https://laposanti](https://github.com/laposanti/BTSBM))** R package for model implementation and MCMC.
 

---

## 🔍 What’s inside (lean & focused)

- **Scripts at repo root**

  - `RUN_MCMC.R` — run the code across the 22 seasons of male ATP tennis tournaments. 
  
  - `multiple_seasons_analysis.R` — post_processes all seasons, saving tables and plots
  
  - `single_season_analysis.R` — post_processes just the selected season, saving tables and plots
  
  - `Model comparison` - compares the BT model with the BT-SBM proposed model.


- **`Other folder**:

  - `simulation/` — it contains the code to run the simulation study (Appendix D)
  
  - `results/` — it collects the raw outputs of the scripts.

  - `images/` — it collects the figures used in the paper (PNG format).
  
  - `tables/` — tables reported in the paper (LaTeX format).

> If `results/`, `images/` or `tables/` are missing, the scripts will create them.


---


### 1) Install & load

```r
# install once (adjust to your GitHub origin if needed)
# install.packages("devtools")
devtools::install_github("laposanti/BTSBM")

library(BTSBM)
```

## ▶️ Fit the model to the data


Place yourself in the repo root and run:

```r
source("RUN_MCMC.R")
```

What it does:

- Iterates over all ATP seasons provided by BTSBM::ATP_2000_2022.

- Prints the season under analysis (with quick stats).

- Saves one file:

`results/augmented_multiple_seasonsGN2.rds`

This will take approximately 35 minutes.

## ▶️ Single-Season Analysis

You can reproduce figures/tables for one season in isolation (e.g. useful for paper insets or diagnostics).

### 📊 Outputs — Single-Season Analysis 

| Description | Script / Object | Preview | Output file |
|---|---|----|---|
| Posterior adjacency matrix - Figure 3 | `postprocessing.R` / `geom_adjacency_fixed` | <a href="./images/adjacency_reordered.png"><img src="./images/adjacency_reordered.png" width="160" alt="Block-ordered adjacency"></a> | [`images/adjacency_reordered.png`](./images/adjacency_reordered.png) |
| Assignment-probabilities heatmap - Figure 4 | `postprocessing.R` / `ass_prob_plot` | <a href="./images/assignment_uncertainty.png"><img src="./images/assignment_uncertainty.png" width="160" alt="Assignment probabilities heatmap"></a> | [`images/assignment_uncertainty.png`](./images/assignment_uncertainty.png) |
| Player skill (λ) uncertainty — Figure 5 | `postprocessing.R` / `plot_lambda` | <a href="./images/lambda_uncertainty.png"><img src="./images/lambda_uncertainty.png" width="160" alt="Lambda uncertainty"></a> | [`images/lambda_uncertainty.png`](./images/lambda_uncertainty.png) |

---

## ▶️ Model Comparison Analysis

```r
source("Model Comparison.R")  # saves the results as csv and reproduces the plot in Fig. 7
```

### 📊 Outputs — Model Comparison Analysis

| Description | Script / Object | Preview | Output file |
|---|---|----|---|
| Model comparison plot - Figure 6 | `Model Comparison.R` / `DELPD_plot.png` | <a href="./images/DELPD_plot.png"><img src="./images/DELPD_plot.png" width="160" alt="Delta ELPD between BT and BT-SBM"></a> | [`images/DELPD_plot.png`](./images/DELPD_plot.png) |
| Model comparison table | `Model Comparison.R` / `model_comparison.csv` | <a href="./tables/model_comparison.csv"><img src="./tables/model_comparison.csv" width="160" alt="Delta ELPD between BT and BT-SBM"></a> | [`tables/model_comparison.csv`](./tables/model_comparison.csv) |

### ▶️ Multiple-Seasons Analysis

- Make sure you have the current directory correctly set in the project root

- Run the `Multiple_seasons_analysis.R` file. 

Outputs (preview below) are written to `images/` and `tables`



### 📊 Outputs — Multiple Seasons 

The following table maps each figure in the paper to its generating code and output file, with a live thumbnail preview.

| Description  | Preview | Output file |
|---|----|---|
| Table with the posterior K probabilities — Table 4   | <a href="./tables/post_numb_block_across_years_table.tex"><img src="./tables/post_numb_block_across_years_table.tex" width="160" alt="Posterior of K across seasons"></a> | [`./tables/post_numb_block_across_years_table.tex`](./tables/post_numb_block_across_years_table.tex) |
| Nº of players in top block by season — Figure 7  | <a href="./images/num_block_plot.png"><img src="./images/num_block_plot.png" width="160" alt="# players in top block"></a> | [`images/num_block_plot.png`](./images/num_block_plot.png) |
| P(Top block) by season — Figure 8 |  <a href="./images/Ptop_across_time.png"><img src="./images/Ptop_across_time.png" width="160" alt="P(top block) by season"></a> | [`images/Ptop_across_time.png`](./images/Ptop_across_time.png) |
| Shannon entropy across seasons — Figure 9 |  | <a href="./images/entropy_plot.png"><img src="./images/entropy_plot.png" width="160" alt="Entropy across seasons"></a> | [`images/entropy_plot.png`](./images/entropy_plot.png) |

> All plots are saved to the `images/` folder, while the table in the `./tables` folder.
> You can customize the output location by modifying the save paths in `postprocessing.R`.
> From table 4, if you take the 2017/2018 row you also get table 2.

---

### ▶️ Appendix- Prior sentivity check.

- Make sure you have the current directory correctly set in the project root

- Run the `simulation.R` file. 

Outputs (preview below) are written to `simulation/`.



### 📊 Output — Prior sensitivity

The following table maps each figure in the paper to its generating code and output file, with a live thumbnail preview.

| Description  | Preview | Output file |
|---|----|---|
| Prior sensitivity table — Table 5  | <a href="./images/num_block_plot.png"><img src="./tables/prior_sensitivity.csv" width="160" alt="#table"></a> | [`./tables/prior_sensitivity.csv`](./tables/prior_sensitivity.csv) |
| Prior sensitivity plot (b = exp(psi(a))) — Figure 10  | <a href="./tables/post_numb_block_across_years_table.tex"><img src="./tables/post_numb_block_across_years_table.tex" width="160" alt="Posterior of K across seasons"></a> | [`./tables/post_numb_block_across_years_table.tex`](./tables/post_numb_block_across_years_table.tex) |
| Prior sensitivity plot (b = 1)  | <a href="./tables/post_numb_block_across_years_table.tex"><img src="./tables/post_numb_block_across_years_table.tex" width="160" alt="Posterior of K across seasons"></a> | [`./tables/post_numb_block_across_years_table.tex`](./tables/post_numb_block_across_years_table.tex) |

> All plots are saved to the `images/` folder, while the table in the `./tables` folder.
> You can customize the output location by modifying the save paths in `postprocessing.R`.

---


### ▶️ Appendix- simulation study.

- Make sure you have the current directory correctly set in the project root

- Run the `simulation.R` file. 

Outputs (preview below) are written to `simulation/`.



### 📊 Output — Simulation study

The following table maps each figure in the paper to its generating code and output file, with a live thumbnail preview.

| Description  | Preview | Output file |
|---|----|---|
| ARI plot — Figure 12  | <a href="./tables/post_numb_block_across_years_table.tex"><img src="./tables/post_numb_block_across_years_table.tex" width="160" alt="Posterior of K across seasons"></a> | [`./tables/post_numb_block_across_years_table.tex`](./tables/post_numb_block_across_years_table.tex) |
| Contingency table — Table 6  | <a href="./images/num_block_plot.png"><img src="./images/num_block_plot.png" width="160" alt="# players in top block"></a> | [`images/num_block_plot.png`](./images/num_block_plot.png) |

> All plots are saved to the `images/` folder, while the table in the `./tables` folder.
> You can customize the output location by modifying the save paths in `postprocessing.R`.

---


