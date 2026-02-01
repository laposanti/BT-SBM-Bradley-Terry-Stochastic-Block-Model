# ================================================================
# Bradley–Terry SBM — Simulation + MCMC Harness (cleaned)
# ================================================================

# ------------------------------------------------
# 0) Packages
# ------------------------------------------------
suppressPackageStartupMessages({
  library(mcclust)        # comp.psm, vi.dist
  library(mcclust.ext)    # minVI, minbinder
  library(mclust)         # adjustedRandIndex
  library(doParallel)     # parallel backend
  library(foreach)        # parallel loops
  library(LaplacesDemon)  # WAIC
  library(filelock)
  library(MASS)
  library(fossil)         # adj.rand.index
  library(coda)           # ESS, Gelman diag
  library(readr)          # write_csv
  library(abind)
  library(dplyr)
  library(kableExtra)
  library(tidyr)
  library(ggplot2)
})

proj_dir <- normalizePath('/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model', mustWork = TRUE)

# Helper: absolute path builders
data_file_K <- function(B) file.path(proj_dir, "data", sprintf("simulated_data_K%d.RDS", B))
out_file    <- file.path(proj_dir, "results", "augmented_simulation_results1.csv")
lockfile    <- paste0(out_file, ".lock")

# (Re)create results file
if (file.exists(out_file)) file.remove(out_file)


# ------------------------------------------------
# 1) Load real data meta (for n_rep) — safe if file exists
# ------------------------------------------------
data_path_real <- './data/2000_2022_data.rds'
if (file.exists(data_path_real)) {
  data_real <- readRDS(data_path_real)
  n_rep <- length(data_real)
} 

# Ensure results dir exists
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("data"))    dir.create("data", recursive = TRUE)

N_ij1 = data_real[[1]]$N_ij
N_ij2 = data_real[[2]]$N_ij
N_ij3 = data_real[[3]]$N_ij

# ---- Gnedin prior: sequential sampler for labels (unconditional) ----
# gamma \in (0,1)
sample_gnedin_labels <- function(n, gamma = 0.5, check = TRUE) {
  stopifnot(n >= 1, is.finite(gamma), gamma > 0, gamma < 1)
  
  x <- integer(n)
  H <- 0L
  sizes <- integer(0)
  
  for (i in seq_len(n)) {
    if (H == 0L) { 
      H <- 1L; sizes <- 1L; x[i] <- 1L
      next 
    }
    
    n_curr <- i - 1L  # should equal sum(sizes)
    
    if (check) {
      stopifnot(length(sizes) == H)
      stopifnot(sum(sizes) == n_curr)
      stopifnot(all(sizes >= 1L))
    }
    
    existing_weights <- (sizes + 1) * (n_curr - H + gamma)  # (n_j+1)(n-k+γ)
    new_weight <- H * (H - gamma)                           # k(k-γ)
    w <- c(existing_weights, new_weight)
    
    if (check) {
      total <- n_curr * (n_curr + gamma)
      if (!all(is.finite(w)) || any(w < 0)) stop("Non-finite or negative weights.")
      if (abs(sum(w) - total) > 1e-8 * max(1, total)) {
        stop("Weights do not sum to n(n+gamma), something is inconsistent.")
      }
    }
    
    a <- sample.int(H + 1L, size = 1L, prob = w)
    if (a <= H) {
      sizes[a] <- sizes[a] + 1L
      x[i] <- a
    } else {
      H <- H + 1L
      sizes <- c(sizes, 1L)
      x[i] <- H
    }
  }
  
  list(x = x, K = H, sizes = sizes)
}

x_samples = matrix(NA, 1000,105)
K_samples = numeric(1000)
for(i in 1:1000){
  x_cur_i = sample_gnedin_labels(105,0.6)$x
  K = length(unique(x_cur_i))
  
  x_samples[i,1:105] = x_cur_i
  K_samples[i] = K
}
draw_gnedin_given_K <- function(n, K_target, gamma = 0.6, max_tries = 1e6) {
  for (t in seq_len(max_tries)) {
    z <- sample_gnedin_labels(n, gamma, check = FALSE)$x
    if (length(unique(z)) == K_target) return(z)
  }
  stop("Failed to sample GN partition with target K. Increase max_tries or adjust gamma.")
}



# Controlled separation: log lambda_k = delta*(K-k)
make_lambda_ladder <- function(K_true, delta = 0.5) {
  log_lambda <- delta * (K_true - seq_len(K_true))
  lambda <- exp(log_lambda)
  lambda
}

assign_lambda_by_cluster_size <- function(z, delta = 0.5) {
  
  
  lambda_by_label
}


simulate_W_from_blocks <- function(N, z_true, lambda_true, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(N)
  W <- matrix(0L, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      nij <- N[i, j]
      if (nij > 0) {
        p_ij <- lambda_true[z_true[i]] / (lambda_true[z_true[i]] + lambda_true[z_true[j]])
        w_ij <- rbinom(1, size = nij, prob = p_ij)
        W[i, j] <- w_ij
        W[j, i] <- nij - w_ij
      }
    }
  }
  diag(W) <- 0L
  W
}

make_lambda_from_padj <- function(K, p_adj = 0.75, base = 1) {
  stopifnot(K >= 2, p_adj > 0.5, p_adj < 1)
  delta <- qlogis(p_adj)              # delta = logit(p_adj)
  lambda <- base * exp(delta * (0:(K-1)))
  lambda
}


lambda_true
bt_p <- function(lambda) {
  outer(lambda, lambda, function(a, b) a / (a + b))
}

summarise_tier_probs <- function(lambda) {
  K <- length(lambda)
  P <- bt_p(lambda)
  adj <- sapply(2:K, function(k) P[k, k-1])  # prob tier k beats tier k-1
  list(
    lambda = lambda,
    p_adj = adj,
    p_strongest_vs_weakest = P[K, 1]
  )
}

summarise_tier_probs(make_lambda_from_padj(5, 0.70))
summarise_tier_probs(make_lambda_from_padj(5, 0.80))
# ~ 1, 3, 9, 27, 81   (ratio 3 each step)






# Simulate data for K = 3,5,7
# Using N_ij1, N_ij2, N_ij3 from real
# To do so, first generate 1000 samples from GN, and randomly pick one with K=3,5,7
# Importantly, smaller blocks need to have the largest lambda, to be better identified
# Sample Lambda given K, x
# generate the wins
# fit the model BTSBM::fit_bit_sbm

n_players <- nrow(N_ij1)
z_true <- draw_gnedin_given_K(n_players, K_target = 5, gamma = 0.6)


lambda_true <- assign_lambda_by_cluster_size(z_true, delta = 0.7)
W <- simulate_W_from_blocks(N_ij1, z_true, lambda_true)

# --- conditional GN draw (via rejection), returning sizes ---
draw_gnedin_sizes_given_K <- function(n, K_target, gamma = 0.6, max_tries = 1e6) {
  for (t in seq_len(max_tries)) {
    z <- sample_gnedin_labels(n, gamma, check = FALSE)$x
    k <- length(unique(z))
    if (k == K_target) {
      sizes <- tabulate(z, nbins = k)
      return(sizes)
    }
  }
  stop("Failed to sample GN partition with target K. Increase max_tries or adjust gamma.")
}

# --- build z by assigning most active players to smallest clusters ---
# activity: numeric length n, e.g. rowSums(N)
# This enforces: argmax(activity) is in a minimum-size cluster.
make_z_activity_to_smallest <- function(sizes, activity, seed = NULL, shuffle_within_cluster = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(activity)
  stopifnot(sum(sizes) == n)
  
  K <- length(sizes)
  
  # clusters ordered by size ascending (smallest first)
  ord_sizes <- order(sizes, decreasing = FALSE)
  sizes_sorted <- sizes[ord_sizes]
  
  # players ordered by activity descending (most active first)
  ord_players <- order(activity, decreasing = TRUE)
  
  z_sorted <- integer(n)
  idx <- 1L
  
  for (c in seq_len(K)) {
    m <- sizes_sorted[c]
    players_c <- ord_players[idx:(idx + m - 1L)]
    
    if (shuffle_within_cluster && m > 1L) players_c <- sample(players_c)
    
    z_sorted[players_c] <- c
    idx <- idx + m
  }
  
  # Now map cluster labels back to the original "unsorted sizes" label space if you want,
  # but you probably *don't*: keeping labels ordered by size is useful (ties, lambda ordering).
  z_sorted
}

# --- assert your constraint explicitly ---
check_top_player_in_smallest_cluster <- function(z, activity) {
  sizes <- tabulate(z, nbins = length(unique(z)))
  min_size <- min(sizes)
  top_player <- which.max(activity)
  cl_top <- z[top_player]
  if (sizes[cl_top] != min_size) stop("Constraint violated: top-activity player not in a smallest cluster.")
  TRUE
}


simulate_dataset_design_based <- function(N, K_target, gamma = 0.8, delta = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.matrix(N), nrow(N) == ncol(N))
  
  n <- nrow(N)
  activity <- rowSums(N)  # total matches per player
  
  # 1) GN sizes (conditional on K)
  sizes <- draw_gnedin_sizes_given_K(n, K_target = K_target, gamma = gamma)
  
  # 2) Assign most active players to smallest clusters
  z_true <- make_z_activity_to_smallest(sizes, activity, seed = seed, shuffle_within_cluster = TRUE)
  check_top_player_in_smallest_cluster(z_true, activity)
  
  # 3) Lambda: since z labels are already ordered by cluster size (smallest = 1),
  # you can directly do: smallest cluster gets largest lambda (ladder)
  
  sizes <- tabulate(z_true, nbins = K_target)
  
  # order clusters by size ascending: smallest gets rank 1
  ord <- order(sizes, decreasing = FALSE)
  
  lambda_rank <- rev(make_lambda_from_padj(K = K_target, p_adj = 0.9))
  lambda_by_label <- numeric(K_target)
  lambda_by_label[ord] <- lambda_rank
  
  # 4) Outcomes
  W <- simulate_W_from_blocks(N, z_true, lambda_by_label)
  
  list(N = N, W = W, z_true = z_true, lambda_true = lambda_by_label, activity = activity, sizes = sizes)
}



sim_data = simulate_dataset_design_based(N_ij1, K_target = 3, gamma = 0.71, delta = 0.7)
table(sim_data$z_true)

res = gibbs_bt_sbm(w_ij = sim_data$W, T_iter = 5000, T_burn = 1000, a =2, gamma=0.71, verbose = T, prior = 'GN')

# --- Summaries ---
x_samples      <- res$x_samples
lambda_samples <- res$lambda_samples

inf.help          <- relabel_by_lambda(x_samples, lambda_samples)
x_samples_relabel <- inf.help$x_samples_relabel
lambda_relabeled  <- inf.help$lambda_samples_relabel

vector_cl <- apply(x_samples_relabel, 1, function(x) length(unique(x)))
mean_cl   <- as.integer(round(median(vector_cl), 0))
HPD_int   <- coda::HPDinterval(coda::mcmc(vector_cl))

minVI_cl        <- inf.help$minVI_partition
binder_cl       <- inf.help$partition_binder
vi_distance     <- mcclust::vi.dist(minVI_cl, sim_data$z_true)
ari_value       <- fossil::adj.rand.index(minVI_cl, sim_data$z_true)
binder_distance <- mcclust::vi.dist(binder_cl, sim_data$z_true)
binder_ari      <- fossil::adj.rand.index(binder_cl, sim_data$z_true)

vi_distance
ari_value


## -------------------------------
## Helpers
## -------------------------------

wins_to_Pi <- function(W) {
  stopifnot(is.matrix(W), nrow(W) == ncol(W))
  idx <- which(W > 0, arr.ind = TRUE)
  n <- sum(W[idx])
  Pi <- matrix(NA_integer_, nrow = n, ncol = 2)
  pos <- 1L
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]
    j <- idx[k, 2]
    m <- W[i, j]
    Pi[pos:(pos + m - 1L), ] <- cbind(rep.int(i, m), rep.int(j, m))
    pos <- pos + m
  }
  Pi
}

rank_from_strength <- function(strengths, player_names) {
  tibble::tibble(
    player = player_names,
    strength = strengths,
    rank = rank(-strengths, ties.method = "average")
  )
}

## Force geometric mean 1: mean(log(v)) = 0
norm_geo1 <- function(v) {
  v <- as.numeric(v)
  v / exp(mean(log(v)))
}

## Nice comparison metrics
compare_lambda <- function(est, truth) {
  est <- as.numeric(est)
  truth <- as.numeric(truth)
  
  tibble::tibble(
    rmse = sqrt(mean((est - truth)^2)),
    mae  = mean(abs(est - truth)),
    spearman = suppressWarnings(cor(est, truth, method = "spearman")),
    kendall  = suppressWarnings(cor(est, truth, method = "kendall")),
    pearson  = suppressWarnings(cor(est, truth, method = "pearson"))
  )
}

## -------------------------------
## 1) Simulate x (cluster labels)
## -------------------------------
set.seed(1)

cluster_sizes <- c(3, 5, 7)
K_true <- length(cluster_sizes)
n <- sum(cluster_sizes)

x_true <- rep(seq_len(K_true), times = cluster_sizes)
x_true <- sample(x_true)  # shuffle node order
table(x_true)

player_names <- paste0("p", seq_len(n))

## -------------------------------
## 2) Simulate well-separated lambda by cluster
##    and DOUBLE CHECK the scale: mean(log lambda)=0
## -------------------------------
## Choose log-worth levels (cluster-wise) with decent separation,
## and then enforce sum(log lambda)=0 exactly.

## Start with something intentionally separated.
## We solve 3*a1 + 5*a2 + 7*a3 = 0 by choosing a2 = 0, a1 = -2, a3 = 6/7.
log_levels <- c(-2, 0, 6/7)  # cluster 1 is weak, cluster 3 strong
stopifnot(abs(sum(cluster_sizes * log_levels)) < 1e-10)

log_lambda_true <- log_levels[x_true]
lambda_true <- exp(log_lambda_true)

## Sanity checks:
cat("Sanity check: mean(log lambda_true) =", mean(log(lambda_true)), "\n")
cat("Sanity check: geometric mean =", exp(mean(log(lambda_true))), "\n")
cat("Cluster-wise lambda_true means:\n")
print(tapply(lambda_true, x_true, mean))

## -------------------------------
## 3) Simulate pairwise win-count matrix w_ij
## -------------------------------
simulate_wij_bt <- function(lambda, mu_games = 12, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(lambda)
  W <- matrix(0L, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      n_games <- rpois(1, mu_games)
      if (n_games == 0) next
      
      p_i <- lambda[i] / (lambda[i] + lambda[j])
      w_i <- rbinom(1, size = n_games, prob = p_i)
      w_j <- n_games - w_i
      
      W[i, j] <- w_i
      W[j, i] <- w_j
    }
  }
  diag(W) <- 0L
  W
}

w_ij <- simulate_wij_bt(lambda_true, mu_games = 12, seed = 2)

## Optional quick signal check: do strong cluster nodes tend to win?
cat("Total wins by node (first 6):\n")
print(head(rowSums(w_ij), 6))

Pi <- wins_to_Pi(w_ij)

## -------------------------------
## 4) Fit models + time them
## -------------------------------

T_iter <- 3000
T_burn <- 600

## 4.1 Plain BT (BTSBM)
t_bt <- system.time({
  fit_bt <- BTSBM::gibbs_bt_simple(
    w_ij = w_ij,
    T_iter = T_iter,
    T_burn = T_burn,
    verbose = FALSE
  )
})

## 4.2 Rank-clustered (rankclust: Erosheva/Pearce)
## First run BTL to initialise nu0, then RCBTL
t_rcbtl <- system.time({
  res_btl <- mcmc_BTL(
    Pi = Pi,
    J = n,
    a_gamma = 1,
    b_gamma = 1,
    num_iters = T_iter,
    burn_prop = T_burn / T_iter,
    chains = 1,
    groupwise = TRUE,
    seed = 1
  )
  
  nu0 <- colMeans(res_btl[, grep("^omega", names(res_btl))])
  
  res_rcbtl <- mcmc_RCBTL(
    Pi = Pi,
    J = n,
    a_gamma = 1,
    b_gamma = 1,
    lambda = 2,   # prior intensity for rank-clusters in their model (not your worth)
    nu0 = nu0,
    num_iters = T_iter,
    nu_reps = 2,
    burn_prop = T_burn / T_iter,
    thin = 1,
    chains = 1,
    groupwise = TRUE,
    seed = 1
  )
})

## 4.3 BTSBM (BT + latent clusters)
t_btsbm <- system.time({
  fit_btsbm <- BTSBM::gibbs_bt_sbm(
    w_ij = as.matrix(w_ij),
    T_iter = T_iter,
    T_burn = T_burn,
    a = 2,
    gamma = 0.8,
    prior = "GN",
    verbose = FALSE
  )
})

timings <- tibble::tibble(
  model = c("BT (BTSBM::gibbs_bt_simple)", "RCBTL (rankclust)", "BT-SBM (BTSBM::gibbs_bt_sbm)"),
  user_time = c(t_bt["user.self"], t_rcbtl["user.self"], t_btsbm["user.self"]),
  sys_time  = c(t_bt["sys.self"],  t_rcbtl["sys.self"],  t_btsbm["sys.self"]),
  elapsed   = c(t_bt["elapsed"],   t_rcbtl["elapsed"],   t_btsbm["elapsed"])
)

print(timings)

## -------------------------------
## 5) Extract lambdas and NORMALISE SCALE consistently
##    Target: mean(log(lambda)) = 0  (geometric mean 1)
## -------------------------------

## BT
bt_lambda_hat <- colMeans(fit_bt$lambda_samples)
bt_lambda_hat <- norm_geo1(bt_lambda_hat)

## RCBTL (rankclust calls them omega)
rcbtl_lambda_hat <- colMeans(res_rcbtl[, grep("^omega", names(res_rcbtl))])
rcbtl_lambda_hat <- norm_geo1(rcbtl_lambda_hat)

## BTSBM
## NOTE: if you relabel for identifiability of clusters, do it before summarising lambdas.
## Your earlier code used relabel_by_lambda; keep it if it’s needed for x, but lambda itself
## is per-node, so we mostly just scale-normalise.
x_samples      <- fit_btsbm$x_samples
lambda_samples <- fit_btsbm$lambda_samples

inf.help          <- relabel_by_lambda(x_samples, lambda_samples)
x_samples_relabel <- inf.help$x_samples_relabel
lambda_relabeled  <- inf.help$lambda_samples_relabel




btsbm_lambda_hat <- colMeans(lambda_relabeled)
btsbm_lambda_hat <- norm_geo1(btsbm_lambda_hat)

## Truth (already mean log = 0 by construction)
lambda_true_scaled <- norm_geo1(lambda_true)

## Double-check all scales:
scale_checks <- tibble::tibble(
  model = c("truth", "bt", "rcbtl", "btsbm"),
  mean_log = c(mean(log(lambda_true_scaled)),
               mean(log(bt_lambda_hat)),
               mean(log(rcbtl_lambda_hat)),
               mean(log(btsbm_lambda_hat))),
  geom_mean = c(exp(mean(log(lambda_true_scaled))),
                exp(mean(log(bt_lambda_hat))),
                exp(mean(log(rcbtl_lambda_hat))),
                exp(mean(log(btsbm_lambda_hat))))
)
print(scale_checks)

## -------------------------------
## 6) Compare lambda estimates
## -------------------------------
perf <- bind_rows(
  compare_lambda(bt_lambda_hat,    lambda_true_scaled) %>% mutate(model = "BT"),
  compare_lambda(rcbtl_lambda_hat, lambda_true_scaled) %>% mutate(model = "RCBTL"),
  compare_lambda(btsbm_lambda_hat, lambda_true_scaled) %>% mutate(model = "BT-SBM")
) %>% dplyr::select(model, everything())

print(perf)

## A compact per-node comparison table (optional)
lambda_table <- tibble::tibble(
  player = player_names,
  x_true = x_true,
  lambda_true = lambda_true_scaled,
  bt = bt_lambda_hat,
  rcbtl = rcbtl_lambda_hat,
  btsbm = btsbm_lambda_hat
) %>%
  arrange(desc(lambda_true))

print(lambda_table)

## If you want rank comparisons too:
bt_ranks    <- rank_from_strength(bt_lambda_hat, player_names)    %>% rename(bt_strength = strength, bt_rank = rank)
rcbtl_ranks <- rank_from_strength(rcbtl_lambda_hat, player_names) %>% rename(rcbtl_strength = strength, rcbtl_rank = rank)
btsbm_ranks <- rank_from_strength(btsbm_lambda_hat, player_names) %>% rename(btsbm_strength = strength, btsbm_rank = rank)

rank_comp <- bt_ranks %>%
  inner_join(rcbtl_ranks, by = "player") %>%
  inner_join(btsbm_ranks, by = "player") %>%
  left_join(tibble::tibble(player = player_names, x_true = x_true, lambda_true = lambda_true_scaled), by = "player") %>%
  arrange(bt_rank)

print()
View(rank_comp)





library(knitr)
library(kableExtra)
library(dplyr)

perf_tex <- perf %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename(
    Model = model,
    RMSE = rmse,
    MAE = mae,
    Spearman = spearman,
    Kendall = kendall,
    Pearson = pearson
  )

tab1 <- kable(
  perf_tex,
  format = "latex",
  booktabs = TRUE,
  caption = "Comparison of posterior mean strength estimates across models. All strengths are normalised to geometric mean 1.",
  align = "lccccc"
) %>%
  kable_styling(latex_options = c("hold_position"))

cat(tab1, file = "results/performance_table.tex")
lambda_tex <- lambda_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  rename(
    Player = player,
    Cluster = x_true,
    `True $\\lambda$` = lambda_true,
    BT = bt,
    RCBTL = rcbtl,
    `BT--SBM` = btsbm
  )

tab2 <- kable(
  lambda_tex,
  format = "latex",
  booktabs = TRUE,
  caption = "True and estimated player strengths ($\\lambda$) under each model. Values scaled to geometric mean 1.",
  align = "lcccccc"
) %>%
  kable_styling(latex_options = c("scale_down"))

cat(tab2, file = "results/lambda_table.tex")
rank_tex <- rank_comp %>%
  select(player, bt_rank, rcbtl_rank, btsbm_rank, x_true) %>%
  rename(
    Player = player,
    `BT rank` = bt_rank,
    `RCBTL rank` = rcbtl_rank,
    `BT--SBM rank` = btsbm_rank,
    Cluster = x_true
  )

tab3 <- kable(
  rank_tex,
  format = "latex",
  booktabs = TRUE,
  caption = "Player ranks induced by posterior mean strengths.",
  align = "lcccc"
)

cat(tab3, file = "results/rank_table.tex")
library(dplyr)
library(knitr)
library(kableExtra)

timings_tex <- timings %>%
  mutate(
    elapsed_min = elapsed / 60,
    user_min    = user_time / 60,
    sys_min     = sys_time / 60
  ) %>%
  transmute(
    Model = model,
    `User (min)`    = round(user_min, 2),
    `System (min)`  = round(sys_min, 2),
    `Elapsed (min)` = round(elapsed_min, 2)
  )

tab_time <- kable(
  timings_tex,
  format = "latex",
  booktabs = TRUE,
  caption = "Computation time for each algorithm. Times are reported in minutes.",
  align = "lccc"
) %>%
  kable_styling(latex_options = c("hold_position"))

cat(tab_time, file = "results/timings_table.tex")







