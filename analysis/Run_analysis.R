# --------------------------------
# Run the code on each Tennis year
# --------------------------------

setwd('/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM/')
#assuming one is running from the project root `bradley-terry-sbm/`

#importing the functions
source("main/mcmc_bt_sbm.R")
source("main/utils.R")

#importing the data
tennis_years = readRDS("./data/2000_2022_data.rds")

# ------------------------------------
# Prior hyperparam setting

H_dm <- 20 
beta_dm <- 1/105
sigma_dm <- 0  

probs_gnedin <- HGnedin(105, 1:105, gamma = .9)
round(sum(1:105*probs_gnedin))
round(expected_cl_py(105, sigma = sigma_dm, theta = beta_dm*H_dm, H = H_dm))

# ------------------------------------
# Running code for each year

#storer of the results
res_list= list()

for(yr in 1:length(tennis_years)){
  Y_ij = tennis_years[[yr]]$Y_ij
  N_ij = tennis_years[[yr]]$N_ij
  
  res <- gibbs_bt_urn_log(
    w_ij=Y_ij, n_ij=N_ij,     # K x K matrices for pairwise data
    a=0.1, 
    b=1,                      # Gamma(a,b) hyperparams for each block rate
    prior = "GN",             # which prior among c("DP","PY","DM","GN")
    alpha_PY = 1,             # parameter for DP/PY
    sigma_PY = 0,             # discount for PY
    beta_DM = beta_dm,        # parameter for DM
    H_DM = H_dm,              # max # clusters for DM
    gamma_GN = .99,           # parameter for Gnedin
    n_iter = 50000,
    burnin = 20000,
    init_x = NULL,            # optional initialization for cluster labels
    store_z = FALSE,          # whether to store all latent Z
    verbose = T
  ) 
  
  res_list[[yr]]<-res
}

#store computations
saveRDS(res_list, "./results/augmented_multiple_seasonsGN.rds")


