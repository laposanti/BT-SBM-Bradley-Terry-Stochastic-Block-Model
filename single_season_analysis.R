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

# -------------------------------
# Single seasons Analysis
# -------------------------------
#change this
current_wd = "./Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model/"
setwd(current_wd)

# ensure images dir exists
if (!dir.exists("./images")) dir.create("./images", recursive = TRUE)

# Setup
first_year <- 1999
yr <- 18  # season index (e.g. 2000/01 + 18 = 2018/2019)
model <- 'GN'
season_label <- paste0(first_year + yr, "/", first_year + yr + 1)
print(season_label)

# Load the output of gibbs_bt_sbm()
fit <- readRDS(paste0("./raw outputs/augmented_multiple_seasonsGN2.rds"))[[yr]]
x_samples <- fit$x_samples
lambda_samples <- fit$lambda_samples

# External data for this season
tennis_years = BTSBM::ATP_2000_2022
w_ij <- tennis_years[[yr]]$Y_ij #pairwise success matrix
pl_df <- tennis_years[[yr]]$players_df #info about players

#Relabeling posterior draws
inf_i <- relabel_by_lambda(x_samples, lambda_samples = lambda_samples)

BTSBM::pretty_table_K_distribution(inf_i)

count_K = function(x) length(unique(x))
coda::HPDinterval(mcmc(apply(x_samples,1,count_K)))

#Reordered heatmap
reordered_heatmap <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$minVI_partition)
reordered_heatmap

#Uncertainty over the assignment
ass_prob_plot <- plot_assignment_probabilities(inf_i,w_ij = w_ij,max_n_clust = 4)
ass_prob_plot

#Posterior Lambdas, and their uncertainty
plot_lambda <- plot_lambda_uncertainty(inf_i,w_ij = w_ij,max_n_clust = 4)
plot_lambda

ggsave(filename = "./images/reordered_heatmap.png",reordered_heatmap, height = 8, width = 10)
ggsave(filename = "./images/ass_prob_plot.png",ass_prob_plot, height = 8, width = 10)
ggsave(filename = "./images/lambda_uncertainty2.png",plot_lambda, height = 6, width = 7)


#-------------------------------
# Saving also the credible balls

#vertical upper bound (coarsest partition)
reordered_heatmap_vert_ub <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$credible_ball_upper_partition)+theme( plot.margin = unit(c(0, 0, 0, 0), "pt"))  # remove outer padding
#vertical lower bound (finest partition)
reordered_heatmap_vert_lb <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$credible_ball_lower_partition)+theme( plot.margin = unit(c(0, 0, 0, 0), "pt"))
#Horizontal bound (complexity unaware)
reordered_heatmap_horiz <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$credible_ball_horiz_partition)+theme( plot.margin = unit(c(0, 0, 0, 0), "pt"))


ggsave(filename = "./images/reordered_heatmap_v_ub.png",reordered_heatmap_vert_ub, height = 8, width = 10)
ggsave(filename = "./images/reordered_heatmap_v_lb.png",reordered_heatmap_vert_lb, height = 8, width = 10)
ggsave(filename = "./images/reordered_heatmap_horiz.png",reordered_heatmap_horiz, height = 8, width = 10)
