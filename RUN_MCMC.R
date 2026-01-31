# --------------------------------
# Run the code on each Tennis year
# --------------------------------

# Set your project root 
#setwd('/.../current folder')
setwd("/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model/")
# Packages
library(BTSBM)          # core package
library(mcclust)        # comp.psm, vi.dist
library(mcclust.ext)    # minVI, minbinder
library(mclust)         # adjustedRandIndex
library(filelock)
library(MASS)           # mvrnorm

# Data
tennis_years <- readRDS("./data/ATP_2000_2025_SN_extended.rds")

# Derive season labels
years <- names(tennis_years)
if (is.null(years) || !all(nzchar(years))) {
  years <- as.character(2000 + seq_along(tennis_years) - 1)
}

# Ensure results directory exists
results_dir <- file.path(getwd(), "raw_output_ext")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  message(sprintf("Created results directory at: %s", results_dir))
}

# --- helpers -------------------------------------------------
fmt_dur <- function(sec) {
  if (sec < 90) return(sprintf("%.1f s", sec))
  if (sec < 5400) return(sprintf("%.1f min", sec/60))
  sprintf("%.2f h", sec/3600)
}


n <- 105
gamma <- 0.8
meanK <- BTSBM::gnedin_K_mean(n,gamma)
varK  <- BTSBM::gnedin_K_var(n,gamma)

c(mean = meanK, variance = varK)
# mean ≈ 2.364270161, variance ≈ 45.95131612
 

# --- storage -------------------------------------------------
res_list <- vector("list", length(tennis_years))
names(res_list) <- years
season_times <- numeric(length(tennis_years))

# --- reproducibility ----------------------------------------
set.seed(1234)

for (i in seq_along(tennis_years)) {
  yr  <- years[i]
  w_ij <- tennis_years[[i]]$Y_ij
  
  message(sprintf("▶ Season %s", yr))
  
  res <- gibbs_bt_sbm(
    w_ij     = w_ij,
    a        = 2,
    prior    = "GN",
    gamma_GN = 0.8,
    T_iter   = 30000,
    T_burn   = 5000,
    init_x   = NULL,
    store_z  = FALSE,
    verbose  = TRUE
  )

  res_list[[i]] <- res
}

# Save computations
out_file <- file.path(results_dir, "MCMC_raw_output_ext.rds")

saveRDS(res_list, out_file)
message(sprintf("Saved results to: %s", out_file))
