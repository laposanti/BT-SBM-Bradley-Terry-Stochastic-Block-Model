# --------------------------------
# Run the code on each Tennis year
# --------------------------------

# Set your project root 
#setwd('/.../current folder')

# Packages
library(BTSBM)          # core package
library(mcclust)        # comp.psm, vi.dist
library(mcclust.ext)    # minVI, minbinder
library(mclust)         # adjustedRandIndex
library(foreach)        # parallel loops             "
library(LaplacesDemon)  # WAIC (optional)
library(filelock)
library(MASS)           # mvrnorm

# Data
tennis_years <- BTSBM::ATP_2000_2022

# Derive season labels
years <- names(tennis_years)
if (is.null(years) || !all(nzchar(years))) {
  years <- as.character(2000 + seq_along(tennis_years) - 1)
}

# Ensure results directory exists
results_dir <- file.path(getwd(), "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  message(sprintf("Created results directory at: %s", results_dir))
}

# Storage
res_list <- vector("list", length(tennis_years))
names(res_list) <- years

# Reproducibility
set.seed(1234)

# Run per season
for (i in seq_along(tennis_years)) {
  yr <- years[i]
  Y_ij <- tennis_years[[i]]$Y_ij
  N_ij <- tennis_years[[i]]$N_ij
  
  # Progress print with a couple of quick stats
  message(sprintf("▶ Season %s \n",yr))
  
  res <- BTSBM::gibbs_bt_sbm(
    w_ij = Y_ij,
    n_ij = N_ij,          # K x K matrices for pairwise data
    a = 0.1,
    b = 1,                # Gamma(a, b) hyperparams for each block rate
    prior = "GN",         # one of c("DP", "PY", "DM", "GN")
    alpha_PY = 1,         # DP/PY param (ignored for GN)
    sigma_PY = 0,         # PY discount (ignored for GN)
    # beta_DM, H_DM removed (not used under GN)
    gamma_GN = 0.8,       # Gnedin parameter
    n_iter = 500,
    burnin = 500,
    init_x = NULL,        # optional init for cluster labels
    store_z = FALSE,      # store all latent Z?
    verbose = TRUE
  )
  
  res_list[[i]] <- res
}

# Save computations
out_file <- file.path(results_dir, "augmented_multiple_seasonsGN1.rds")
saveRDS(res_list, out_file)
message(sprintf("Saved results to: %s", out_file))
