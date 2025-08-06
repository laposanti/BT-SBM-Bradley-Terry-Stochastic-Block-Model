# ------------------------------------------------
# 1) Load packages
# ------------------------------------------------
library(mcclust)       # for comp.psm, vi.dist
library(mcclust.ext)   # for minVI and minbinder
library(mclust)        # for adjustedRandIndex
library(doParallel)    # for parallel backend
library(foreach)       # for parallel loops
library(LaplacesDemon) # for WAIC (if not installed, install.packages("LaplacesDemon"))
library(filelock) 
library(MASS)  # needed for mvrnorm
setwd('/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM/')


data_real = readRDS('/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM/data/2000_2022_data.rds')
n_rep = length(data_real)
## ----------------------------------------------------------
## 1.  Helper functions (define ONCE)
## ----------------------------------------------------------


sample_Nij <- function(n_players, counts, density) {
  pairs <- utils::combn(n_players, 2)
  N_out <- matrix(0L, n_players, n_players)
  for (k in seq_len(ncol(pairs))) {
    n_ij <- rpois(1,5)
    N_out[i, j] <- N_out[j, i] <- n_ij
  }
  N_out
}


lambda_to_theta <- function(lambda) {
  outer(lambda, lambda, function(a, b) a / (a + b))
}

## ----------------------------------------------------------
## 2.  Generate datasets for K = 3..10
## ----------------------------------------------------------
for (K in 3:10) {
  K_data <- vector("list", n_rep)
  
  for (s in seq_len(n_rep)) {
    set.seed(123 + s)
    
    n       <- 105
    ## --- simulate match topology ----------
    N_sim <- sample_Nij(n)
    
    ## --- latent block structure ----------
    x_star      <- sample.int(K, n, replace = TRUE)   # random labels
    
    lambda_star <- seq(0.1,2,length.out=K)
    theta_star  <- t(lambda_to_theta(lambda_star))
    
    ## --- generate win counts -------------
    w_ij <- matrix(0L, n, n)
    idx  <- which(upper.tri(N_sim) & N_sim > 0L, arr.ind = TRUE)
    for (k in seq_len(nrow(idx))) {
      i <- idx[k, 1]; j <- idx[k, 2]
      nij <- N_sim[i, j]
      pij <- theta_star[x_star[i], x_star[j]]
      wij <- rbinom(1, nij, pij)
      w_ij[i, j] <- wij
      w_ij[j, i] <- nij - wij
    }
    
    ## --- stash ---------------------------
    K_data[[s]] <- list(
      Y_ij       = w_ij,
      N_ij       = N_sim,
      lambda_star= lambda_star,
      theta_star = theta_star,
      z_star     = x_star,
      K          = K,
      uid        = sprintf("s%03d-K%d", s, K)
    )
  }
  
  saveRDS(K_data,
          file = sprintf("./data/simulated_data_K%d.RDS", K))
}

# ------------------------------------------------
# 2) Setup parallel cluster
# ------------------------------------------------
n_cores <- max(1, parallel::detectCores() - 1)
cl      <- makePSOCKcluster(n_cores)
registerDoParallel(cl)

# ------------------------------------------------
# 3) Define your data generator as before
#----------------------------------------
generate_data <- function(i, B) {
  data_path <- paste0("./data/simulated_data_K", B, ".RDS")
  data_list <- readRDS(data_path) 
  data_i    <- data_list[[i]]
  return(list(
    Y_ij   = data_i$Y_ij,
    N_ij   = data_i$N_ij,
    z_star = data_i$z_star
  ))
}
# --- Setup paths

out_dir  <- normalizePath("results", mustWork = TRUE)
out_file <- file.path(out_dir, "augmented_simulation_results_real.csv")
lockfile <- paste0(out_file, ".lock")
if (file.exists(out_file)) file.remove(out_file)

# --- Cluster
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)


# Export any objects/functions defined in the global environment *before* foreach.
# If your source file defines them, we source once per worker instead.
all_combos <- expand.grid(i = 1:n_rep, B = 8:10)

results_df <- foreach(myCombo = iter(all_combos, by = "row"),
                      .combine = rbind,
                      .packages = c("mcclust","mcclust.ext","mclust",
                                    "LaplacesDemon","fossil","coda",
                                    "readr","filelock")) %dopar% {
                                      
                                      ## One-time per worker initialization
                                      if (!exists(".INIT_DONE", envir = .GlobalEnv)) {
                                        source("./main/mcmc_bt_sbm.R", local = FALSE)
                                        assign(".INIT_DONE", TRUE, envir = .GlobalEnv)
                                      }
                                      
                                      i     <- myCombo$i
                                      B     <- myCombo$B
                                      prior <- "GN"
                                      
                                      ## Wrap the *entire* work so we always write *something*
                                      out_row <- tryCatch({
                                        
                                        # ---- 1. Data
                                        sim_data <- generate_data(i = i, B = B)
                                        Y_ij <- sim_data$Y_ij
                                        N_ij <- sim_data$N_ij
                                        z_true <- sim_data$z_star
                                        
                                        
                                        # ---- 2. MCMC
                                        res <- gibbs_bt_sbm(
                                          w_ij     = Y_ij,
                                          n_ij     = N_ij,
                                          a        = 0.1, b = 1,
                                          prior    = prior,
                                          gamma_GN = 0.8,
                                          burnin   = 5000, n_iter = 10000,
                                          init_x   = NULL, store_z = FALSE, verbose = T
                                        )
                                        
                                        # ---- 3. Summaries
                                        x_samples         <- res$x_samples
                                        lambda_samples    <- res$lambda_samples
                                        inf.help          <- inference_helper(x_samples, lambda_samples)
                                        x_samples_relabel <- inf.help$x_samples_relabel
                                        lambda_relabeled  <- inf.help$lambda_samples_relabel
                                        
                                        vector_cl <- apply(x_samples_relabel, 1, function(x) length(unique(x)))
                                        mean_cl   <- round(median(vector_cl),0)
                                        HPD_int   <- coda::HPDinterval(coda::mcmc(vector_cl))
                                        
                                        minVI_cl        <- inf.help$minVI_partition
                                        binder_cl       <- inf.help$partition_binder
                                        vi_distance     <- mcclust::vi.dist(minVI_cl, z_true)
                                        ari_value       <- fossil::adj.rand.index(minVI_cl, z_true)
                                        binder_distance <- mcclust::vi.dist(binder_cl, z_true)
                                        binder_ari      <- fossil::adj.rand.index(binder_cl, z_true)
                                        
                                        # ---- WAIC (unchanged, could vectorize)
                                        n_samples <- nrow(x_samples_relabel)
                                        n_edges   <- sum(N_ij > 0 & upper.tri(N_ij))
                                        ll_matrix <- matrix(0, nrow = n_samples, ncol = n_edges)
                                        for (t in seq_len(n_samples)) {
                                          idx <- 1
                                          for (col_j in 1:105) for (row_i in 1:col_j) if (N_ij[row_i, col_j] > 0) {
                                            lam_i <- lambda_relabeled[t, row_i]; lam_j <- lambda_relabeled[t, col_j]
                                            p_ij  <- lam_i / (lam_i + lam_j)
                                            ll_matrix[t, idx] <- Y_ij[row_i, col_j] * log(p_ij) +
                                              (N_ij[row_i, col_j] - Y_ij[row_i, col_j]) * log(1 - p_ij)
                                            idx <- idx + 1
                                          }
                                        }
                                        waic_val <- LaplacesDemon::WAIC(ll_matrix)$WAIC
                                        ess_lambda <- coda::effectiveSize(lambda_relabeled)
                                        mean_ess   <- mean(ess_lambda)
                                        
                                        data.frame(
                                          dataset_id  = i,
                                          B           = B,
                                          prior       = prior,
                                          VI_K        = length(unique(minVI_cl)),
                                          mean_num_cl = mean_cl,
                                          lower_ci_cl = HPD_int[1, "lower"],
                                          upper_ci_cl = HPD_int[1, "upper"],
                                          ARI_minVI   = ari_value,
                                          VI_minVI    = vi_distance,
                                          ARI_binder  = binder_ari,
                                          VI_binder   = binder_distance,
                                          waic        = waic_val,
                                          mean_ess    = mean_ess,
                                          error       = NA_character_
                                        )
                                        
                                      }, error = function(e) {
                                        # Return an error row so you can see what failed later
                                        data.frame(
                                          dataset_id  = i,
                                          B           = B,
                                          prior       = prior,
                                          VI_K        = NA,
                                          mean_num_cl = NA,
                                          lower_ci_cl = NA,
                                          upper_ci_cl = NA,
                                          ARI_minVI   = NA,
                                          VI_minVI    = NA,
                                          ARI_binder  = NA,
                                          VI_binder   = NA,
                                          waic        = NA,
                                          mean_ess    = NA,
                                          error       = as.character(e$message)
                                        )
                                      })
                                      
                                      # ---- Thread-safe write
                                      # Lock a *separate* file so we never block ourselves by truncation races
                                      lock_obj <- filelock::lock(lockfile)
                                      on.exit(filelock::unlock(lock_obj), add = TRUE)
                                      
                                      write_header <- !file.exists(out_file) || file.size(out_file) == 0
                                      readr::write_csv(out_row, out_file, append = TRUE, col_names = write_header)
                                      
                                      out_row
                                    }

stopCluster(cl)
write.csv(results_df, file = "./results/augmented_simulation_results_real.csv", row.names = FALSE)


table(results_df$B,results_df$VI_K)



table(results_df$B,results_df$VI_K)

## ----------------------------------------------------------
## 2.  Helper: compute ESS & R‑hat for a list of chains
## ----------------------------------------------------------
diagnostics_lambda <- function(lambda_chain_list) {
  if (length(lambda_chain_list) < 2)
    stop("Need at least two chains for Gelman–Rubin diagnostics")
  
  ## —‑‑ find the shortest surviving chain ‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑
  m_vec   <- vapply(lambda_chain_list, nrow, integer(1))
  m_min   <- min(m_vec)                 # common length
  
  ## —‑‑ coerce to mcmc objects with identical mcpar ‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑
  mcmc_list <- mcmc.list(lapply(lambda_chain_list, function(mat) {
    mcmc(mat[seq_len(m_min), , drop = FALSE], start = 1, thin = 1)
  }))
  
  ess_vec  <- effectiveSize(mcmc_list)               # ESS per λ_j
  gel <- gelman.diag(mcmc_list[,1:10],
                     autoburnin   = FALSE,
                     multivariate = F)
  
  
  rhat_mean         <- mean(gel$psrf[,'Point est.'])
  ## after computing ess_vec and rhat_mat
  n_par <- length(ess_vec)
  param_names <- colnames(lambda_chain_list[[1]])
  if (is.null(param_names))
    param_names <- paste0("lambda_", seq_len(n_par))  # fallback names
  
  data.frame(
    param = 1,
    ESS   = mean(ess_vec, na.rm = T),
    R_mean  = rhat_mean
  )
}

run_one_chain <- function(seed, Y_ij, N_ij, ...) {
  set.seed(seed)
  gibbs_bt_sbm(
    w_ij   = Y_ij,
    n_ij   = N_ij,
    a      = 0.1, b = 1,
    prior  = "GN",
    gamma_GN = 0.8,
    burnin = 2000, n_iter = 5000,
    init_x = NULL, store_z = FALSE, verbose = T,
    ...
  )
}
generate_data <- function(i, B) {
  data_path <- paste0("./data/simulated_data_K", B, ".RDS")
  data_list <- readRDS(data_path) 
  data_i    <- data_list[[i]]
  return(list(
    Y_ij   = data_i$Y_ij,
    N_ij   = data_i$N_ij,
    z_star = data_i$z_star
  ))
}


# --- Setup paths

out_dir  <- normalizePath("results", mustWork = TRUE)
out_file <- file.path(out_dir, "augmented_ simulation_diagnostics.csv")
lockfile <- paste0(out_file, ".lock")
if (file.exists(out_file)) file.remove(out_file)




results_dummy <- foreach(B = 3:10,
                         .packages = c("mcclust", "mcclust.ext", "mclust",
                                       "LaplacesDemon", "fossil", "coda",
                                       "abind", "filelock")) %do% {
                                         if (!exists(".INIT_DONE", envir = .GlobalEnv)) {
                                           source("./main/mcmc_bt_sbm.R", local = FALSE)
                                           assign(".INIT_DONE", TRUE, envir = .GlobalEnv)
                                         }
                                         i   <- 1
                                         dat <- generate_data(i, B)
                                         Y_ij <- dat$Y_ij ; N_ij <- dat$N_ij ; z_true <- dat$z_star
                                         n_chains=3
                                         chain_seeds <- 1234 + seq_len(n_chains) + 1000 * i + 10 * B
                                         chain_out   <- lapply(chain_seeds, run_one_chain, Y_ij = Y_ij, N_ij = N_ij)
                                         
                                         lambda_all <- do.call(rbind, lapply(chain_out, `[[`, "lambda_samples"))
                                         x_all      <- do.call(rbind, lapply(chain_out, `[[`, "x_samples"))
                                         chain_id   <- rep(seq_along(chain_out), times = vapply(chain_out, function(ch) nrow(ch$x), integer(1)))
                                         n_players  <- ncol(x_all)
                                         
                                         relabel_one <- function(x_row, lam_row) {
                                           occ     <- sort(unique(x_row))
                                           lam_occ <- lam_row[occ]
                                           ord     <- order(lam_occ, decreasing = TRUE)
                                           map     <- setNames(seq_along(ord), occ[ord])
                                           x_new   <- map[as.character(x_row)]
                                           lam_new <- lam_occ[ord][x_new]
                                           list(x = x_new, lambda = lam_new)
                                         }
                                         
                                         relab_list <- Map(relabel_one,
                                                           split(x_all, row(x_all)),
                                                           split(lambda_all, row(lambda_all)))
                                         
                                         x_relab      <- t(vapply(relab_list, `[[`, integer(n_players), 1L))
                                         lambda_relab <- t(vapply(relab_list, `[[`, numeric(n_players),  2L))
                                         
                                         K_hat <- length(unique(minVI(mcclust::comp.psm(x_relab))$cl))
                                         
                                         keep_by_chain <- tapply(seq_len(nrow(x_relab)), chain_id, function(idx) {
                                           idx[vapply(idx, function(r) length(unique(x_relab[r, ])) == K_hat, logical(1L))]
                                         })
                                         
                                         lambda_keep_list <- mapply(function(idx, mat = lambda_relab)
                                           mat[idx, , drop = FALSE], keep_by_chain, SIMPLIFY = FALSE)
                                         lambda_keep_list <- Filter(function(mat) nrow(mat) > 0, lambda_keep_list)
                                         
                                         if (length(lambda_keep_list) >= 2) {
                                           diag_tbl <- diagnostics_lambda(lambda_chain_list = lambda_keep_list)
                                           mean_ess <-  diag_tbl$ESS
                                           mean_rhat <- diag_tbl$R_mean
                                         } else {
                                           mean_ess <- NA_real_
                                           max_rhat <- NA_real_
                                           warning("Only one chain left after filtering – R-hat not defined")
                                         }
                                         
                                         out_row <- data.frame(
                                           dataset_id  = i,
                                           B           = B,
                                           mean_ess    = mean_ess,
                                           mean_rhat <- diag_tbl$R_mean
                                         )
                                         
                                         out_file <- file.path("results", "augmented_simulation_diagnostics.csv")
                                         lockfile <- paste0(out_file, ".lock")
                                         
                                         # Write output row safely
                                         lock_obj <- filelock::lock(lockfile)
                                         on.exit(filelock::unlock(lock_obj), add = TRUE)
                                         
                                         write_header <- !file.exists(out_file) || file.size(out_file) == 0
                                         readr::write_csv(out_row, out_file, append = TRUE, col_names = write_header)
                                         
                                         out_row
                                       }


time_df = data.frame(K=3:10, time=numeric(8))

for(K in 3:10){
  
  tm = system.time({
    res <- gibbs_bt_sbm(
      w_ij     = Y_ij,
      n_ij     = N_ij,
      a        = 0.1, 
      b        = 1,
      prior    = 'GN',
      gamma_GN = 0.8,
      burnin   = 500, 
      n_iter   = 3000,
      init_x   = NULL, 
      store_z  = FALSE, 
      verbose  = TRUE
    )
  })
  time_df[which(time_df$K==K),"time"]<- tm[3]
}

library(dplyr)
library(kableExtra)
time_df$time = round(time_df$time*10/60,2)
time_df|>
  kableExtra::kable()

results_df = read.csv('./results/augmented_simulation_results.csv')

results_df|>
  group_by(B)|>
  summarise(mean_VI = round(mean(VI_minVI),2),
            mean_ARI = round(mean(ARI_minVI),2))|>
  kable(format = 'latex',booktabs=T)



