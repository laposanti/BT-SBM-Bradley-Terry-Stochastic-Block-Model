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
current_wd<- "/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model/"
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
fit <- readRDS(paste0("results/augmented_multiple_seasonsGN2.rds"))[[yr]]
x_samples <- fit$x_samples
lambda_samples <- fit$lambda_samples

# External data for this season
tennis_years = BTSBM::ATP_2000_2022
w_ij <- tennis_years[[yr]]$Y_ij #pairwise success matrix
pl_df <- tennis_years[[yr]]$players_df #info about players

#Relabeling posterior draws
inf_i <- relabel_by_lambda(res_i$x_samples, lambda_samples = res_i$lambda_samples)

#Producing the plots

#Reordered heatmap
geom_adjacency_fixed<- plot_block_adjacency(inf_i,w_ij = Y_ij)
geom_adjacency_fixed

#Uncertainty over the assignment
ass_prob_plot<- BTSBM::plot_assignment_probabilities(inf_i,w_ij = Y_ij,max_n_clust = 4)
ass_prob_plot

#Posterior Lambdas, and their uncertainty
plot_lambda<- BTSBM::plot_lambda_uncertainty(inf_i,w_ij = Y_ij,max_n_clust = 4)
plot_lambda

ggsave(filename = "./images/ass_prob_plot2.png",ass_prob_plot, height = 8, width = 10)
ggsave(filename = "./images/geom_adjacency_fixed2.png",geom_adjacency_fixed, height = 8, width = 10)
ggsave(filename = "./images/lambda_uncertainty2.png",plot_lambda, height = 6, width = 7)


