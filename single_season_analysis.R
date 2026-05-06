library(mcclust)
library(mcclust.ext)
library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)
library(ggrepel)
library(kableExtra)
library(reshape2)
library(stringr)
library(ggside)
library(BTSBM)

# Used for combining plots and plot margins
library(cowplot)
library(grid)
library(knitr)
library(tibble)








# -------------------------------
# Single seasons Analysis
# -------------------------------

# ============================================================
# 0. SETUP
# ============================================================

# Set working directory (change as needed)
current_wd <- "/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model/"
setwd(current_wd)

if (!dir.exists("./images"))  dir.create("./images",  recursive = TRUE)
if (!dir.exists("./results")) dir.create("./results", recursive = TRUE)

source("plotting_functions_temp.R")

# --- Load data and MCMC output ---
data <- readRDS("./data/ATP_2000_2026_SN_extended.rds")
res  <- readRDS("./raw_output_ext/MCMC_raw_output_ext.rds")

# Season of interest
season <- "2017"
w_ij   <- data[[season]]$Y_ij
pl_df  <- data[[season]]$players_df

player_names <- pl_df$player_label
if (is.null(player_names) || length(player_names) != nrow(w_ij)) {
  # fallback if the dataset uses a different name column
  player_names <- pl_df$player_slug
}

# Run relabelling: sorts block labels by decreasing lambda, computes credible balls
# (slow, ~1 min)
inf_i <- relabel_by_lambda(res[[season]]$x_samples, res[[season]]$lambda_samples)





# ============================================================
# 1. Exploratory adjacency (Fig 1)
# ============================================================

p_expl_adjacency <- exploratory_adjacency(
  w_ij,
  players_df = pl_df,
  players_id_col = "player_label",
  players_rank_col = "last_rank"
) + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))

ggsave(
  filename = "./images/exploratory_reorderdered_bw.pdf",
  plot = p_expl_adjacency,
  height = 8,
  width = 10
)


# ============================================================
# 2. Reordered adjacency heatmaps
#    - point estimates (minVI + Binder)
#    - K=4 representative partition
#    - credible ball bounds (upper/lower/horizontal)
# ============================================================

plot_margin0 <- ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))

# 2.1 Point estimate partitions
p_adj_minVI <- plot_block_adjacency(
  w_ij = w_ij,
  x_hat = inf_i$minVI_partition,
  players_df = pl_df,
  bw_preview = FALSE
) + plot_margin0

ggsave(
  filename = "./images/reordered_heatmap_point_estimatebw.pdf",
  plot = p_adj_minVI,
  height = 8,
  width = 11
)


# 2.2 Alternative point estimate with K = 4 (SALSO)
if (requireNamespace("salso", quietly = TRUE)) {
  x_hat_K4 <- salso::salso(
    inf_i$x_samples_relabel,
    loss = salso::VI(a = .5),
    maxNClusters = 4
  )

  p_adj_K4 <- plot_block_adjacency(
    w_ij = w_ij,
    x_hat = x_hat_K4,
    players_df = pl_df,
    bw_preview = FALSE
  ) + plot_margin0

  ggsave(
    filename = "./images/reordered_heatmap_point_estimateK4bw.pdf",
    plot = p_adj_K4,
    height = 8,
    width = 11
  )
}

# 2.3 Credible ball bounds
p_cb_upper <- plot_block_adjacency(
  fit = inf_i,
  w_ij = w_ij,
  players_df = pl_df,
  x_hat = inf_i$credible_ball_upper_partition
) + plot_margin0

p_cb_lower <- plot_block_adjacency(
  fit = inf_i,
  w_ij = w_ij,
  players_df = pl_df,
  x_hat = inf_i$credible_ball_lower_partition
) + plot_margin0

p_cb_horiz <- plot_block_adjacency(
  fit = inf_i,
  w_ij = w_ij,
  players_df = pl_df,
  x_hat = inf_i$credible_ball_horiz_partition
) + plot_margin0

ggsave(filename = "./images/reordered_heatmap_v_ubbw.pdf", plot = p_cb_upper, height = 8, width = 11)
ggsave(filename = "./images/reordered_heatmap_v_lbbw.pdf", plot = p_cb_lower, height = 8, width = 11)
ggsave(filename = "./images/reordered_heatmap_horizbw.pdf", plot = p_cb_horiz, height = 8, width = 11)

# Optional: a 3-panel combined figure for the bounds (axes/legend stripped)
strip_xy_legend <- ggplot2::theme(
  legend.position = "none",
  axis.title.y   = ggplot2::element_blank(),
  axis.title.x   = ggplot2::element_blank(),
  axis.text.y    = ggplot2::element_blank(),
  axis.text.x    = ggplot2::element_blank(),
  axis.ticks.y   = ggplot2::element_blank(),
  axis.ticks.x   = ggplot2::element_blank(),
  strip.text.y   = ggplot2::element_blank(),
  strip.text.x   = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank()
)

add_bottom_title <- function(p, label, fontsize = 11, pad_bottom = 0.08) {
  cowplot::ggdraw() +
    cowplot::draw_plot(p, x = 0, y = pad_bottom, width = 1, height = 1 - pad_bottom) +
    cowplot::draw_label(
      label,
      x = 0.5,
      y = pad_bottom * 0.5,
      vjust = 1,
      hjust = 0.5,
      fontface = "bold",
      size = fontsize
    )
}

p_ub <- add_bottom_title(p_cb_upper + strip_xy_legend, "Vertical upper bound")
p_lb <- add_bottom_title(p_cb_lower + strip_xy_legend, "Vertical lower bound")
p_hz <- add_bottom_title(p_cb_horiz + strip_xy_legend, "Horizontal bound")

combined_bounds <- cowplot::plot_grid(p_ub, p_lb, p_hz, nrow = 1, align = "hv", axis = "tb")
ggsave("./images/combined_bounds_bw.pdf", combined_bounds, width = 12, height = 4)


# ============================================================
# 3. Assignment probabilities uncertainty (Fig 5)
# ============================================================

p_assignment <- plot_assignment_probabilities(
  inf_i,
  w_ij,
  players_df = pl_df,
  max_n_clust = 4
)

ggsave(filename = "./images/plot_ass.pdf", plot = p_assignment, height = 8, width = 11)


# ============================================================
# 4. Lambda uncertainty plots (conditional/unconditional)
# ============================================================

p_lambda_cond <- plot_lambda_uncertainty(
  inf_i,
  w_ij = w_ij,
  max_n_clust = 4,
  prob = 0.9,
  conditional = TRUE
)

p_lambda_uncond <- plot_lambda_uncertainty(
  inf_i,
  w_ij = w_ij,
  prob = 0.9,
  conditional = FALSE
)

ggsave(filename = "./images/conditional_lambda_plot.pdf", plot = p_lambda_cond, height = 9, width = 9)
ggsave(filename = "./images/unconditional_lambda_plot.pdf", plot = p_lambda_uncond, height = 9, width = 9)


# ============================================================
# 5. Posterior ranks via midranks of lambda draws
# ============================================================

rank_res <- compute_expected_wins_rank_posterior(
  inf_i$lambda_samples_relabel,
  player_names = player_names
)

p_rank_intervals <- plot_rank_intervals(rank_res)
ggsave(filename = "./images/plot_rank.pdf", plot = p_rank_intervals, height = 9, width = 9)

if (exists("plot_rank_distribution_heatmap")) {
  p_rank_heat <- plot_rank_distribution_heatmap(rank_res$rank_prob, rank_res$summary, max_rank = 40)
  ggsave(filename = "./images/rank_distribution_heatmap_midrank.pdf", plot = p_rank_heat, height = 10, width = 10)
}



