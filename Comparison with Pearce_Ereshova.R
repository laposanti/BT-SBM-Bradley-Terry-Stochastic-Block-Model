## ------------------------------------------------------------
## Simulate a dataset using the design-based DGP from new_sim_study.R
## (fixed N_ij design matrix + GN sizes conditional on K + p_adj separation),
## then fit:
##  - BT (BTSBM::gibbs_bt_simple)
##  - Rank-clustered BT (rankclust: mcmc_RCBTL)
##  - BT-SBM (BTSBM::gibbs_bt_sbm)
## Then time + compare lambda estimates on the SAME SCALE:
##   mean(log(lambda)) = 0  (geometric mean = 1)
## ------------------------------------------------------------

library(rankclust)
library(BTSBM)
library(coda)
library(parallel)
library(mcclust)
library(mcclust.ext)
library(dplyr)
library(tidyr)


relabel_by_lambda = function (x_samples, lambda_samples) 
{
  stopifnot(is.matrix(x_samples))
  S <- nrow(x_samples)
  N <- ncol(x_samples)
  is_list_format <- is.list(lambda_samples)
  get_lambda_vec <- function(iter) {
    if (is_list_format) {
      v <- lambda_samples[[iter]]
      if (!is.numeric(v)) 
        stop("lambda_samples[[iter]] must be numeric.")
      v
    }
    else {
      lambda_samples[iter, ]
    }
  }
  x_relabeled <- matrix(NA_integer_, S, N)
  lambda_per_item <- matrix(NA_real_, S, N)
  cluster_lambda_ordered <- vector("list", S)
  n_clusters_each_iter <- integer(S)
  top_block_count_per_iter <- integer(S)
  for (iter in seq_len(S)) {
    xi <- as.integer(x_samples[iter, ])
    occ_raw <- sort(unique(xi))
    xi_seq <- match(xi, occ_raw)
    K <- max(xi_seq)
    lam_vec_full <- get_lambda_vec(iter)
    lam_occ <- rep(NA_real_, length(occ_raw))
    ok_idx <- occ_raw <= length(lam_vec_full)
    lam_occ[ok_idx] <- lam_vec_full[occ_raw[ok_idx]]
    ord <- order(lam_occ, decreasing = TRUE, na.last = TRUE)
    occ_ord <- occ_raw[ord]
    lam_ord <- lam_occ[ord]
    if (anyNA(lam_ord)) 
      lam_ord[is.na(lam_ord)] <- .Machine$double.xmin
    raw_to_ord_id <- integer(max(occ_ord))
    raw_to_ord_id[occ_ord] <- seq_len(K)
    xi_new <- raw_to_ord_id[occ_raw[xi_seq]]
    x_relabeled[iter, ] <- xi_new
    lambda_per_item[iter, ] <- lam_ord[xi_new]
    cluster_lambda_ordered[[iter]] <- lam_ord
    n_clusters_each_iter[iter] <- K
    top_block_count_per_iter[iter] <- sum(xi_new == 1L)
  }
  modal_K <- as.integer(names(which.max(table(n_clusters_each_iter))))
  psm <- mcclust::comp.psm(x_samples)
  partition_binder <- mcclust.ext::minbinder.ext(psm, cls.draw = x_samples, 
                                                 method = "all")$cl[1, ]
  partition_minVI <- mcclust.ext::minVI(psm, cls.draw = x_samples, 
                                        method = "all")$cl[1, ]
  x_ball <- mcclust.ext::credibleball(c.star = partition_minVI, 
                                      cls.draw = x_samples, c.dist = "VI")
  relabel_partition_by_item_mean_lambda <- function(z, lambda_item_mean) {
    stopifnot(length(z) == length(lambda_item_mean))
    z <- as.integer(z)
    labs <- sort(unique(z))
    cl_means <- vapply(labs, function(k) mean(lambda_item_mean[z == 
                                                                 k], na.rm = TRUE), numeric(1))
    ord <- order(cl_means, decreasing = TRUE)
    new_ids <- seq_along(labs)
    names(new_ids) <- labs[ord]
    z_new <- new_ids[as.character(z)]
    as.integer(z_new)
  }
  lambda_item_mean <- colMeans(lambda_per_item, na.rm = TRUE)
  partition_minVI = relabel_partition_by_item_mean_lambda(partition_minVI, 
                                                          lambda_item_mean)
  partition_binder = relabel_partition_by_item_mean_lambda(partition_binder, 
                                                           lambda_item_mean)
  get_part <- function(obj, name1, name2) {
    if (!is.null(obj[[name1]])) 
      obj[[name1]]
    else obj[[name2]]
  }
  c_lower_raw <- get_part(x_ball, "c.lower", "c.lowervert")
  c_upper_raw <- get_part(x_ball, "c.upper", "c.uppervert")
  c_horiz_raw <- x_ball$c.horiz
  
  pick_row <- function(obj, centre) {
    if (is.vector(obj) && length(obj) == N) return(as.integer(obj))
    if (is.matrix(obj) && ncol(obj) == N) {
      d <- apply(obj, 1, function(z) mcclust::vi.dist(as.integer(z), as.integer(centre)))
      return(as.integer(obj[which.min(d), ]))
    }
    stop("Unexpected credibleball partition format.")
  }
  
  c_lower_vec <- pick_row(c_lower_raw, partition_minVI)
  c_upper_vec <- pick_row(c_upper_raw, partition_minVI)
  c_horiz_vec <- pick_row(c_horiz_raw, partition_minVI)
  
  c_lower_rl <- relabel_partition_by_item_mean_lambda(c_lower_vec, lambda_item_mean)
  c_upper_rl <- relabel_partition_by_item_mean_lambda(c_upper_vec, lambda_item_mean)
  c_horiz_rl <- relabel_partition_by_item_mean_lambda(c_horiz_vec, lambda_item_mean)
  
  K_VI_upper <- length(unique(c_upper_rl))
  K_VI_lower <- length(unique(c_lower_rl))
  K_VI_horiz <- length(unique(c_lower_rl))
  Kmax <- N
  assignment_probs <- matrix(0, nrow = N, ncol = Kmax)
  for (k in seq_len(Kmax)) {
    assignment_probs[, k] <- colMeans(x_relabeled == k, na.rm = TRUE)
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(Kmax))
  rownames(assignment_probs) <- paste0("Item_", seq_len(N))
  assignment_probs_df <- as.data.frame(assignment_probs)
  bc_tab <- table(n_clusters_each_iter)
  block_count_df <- data.frame(num_blocks = as.integer(names(bc_tab)), 
                               count = as.vector(bc_tab), prob = as.vector(bc_tab)/sum(bc_tab))
  list(x_samples_relabel = x_relabeled, lambda_samples_relabel = lambda_per_item, 
       cluster_lambda_ordered = cluster_lambda_ordered, co_clustering = psm, 
       minVI_partition = partition_minVI, partition_binder = partition_binder, 
       n_clusters_each_iter = n_clusters_each_iter, block_count_distribution = block_count_df, 
       item_cluster_assignment_probs = assignment_probs_df, 
       avg_top_block_count = mean(top_block_count_per_iter), 
       top_block_count_per_iter = top_block_count_per_iter, 
       credible_ball_lower_partition = c_lower_rl, credible_ball_upper_partition = c_upper_rl, 
       credible_ball_horiz_partition = c_horiz_rl, K_VI_lower = K_VI_lower, 
       K_VI_upper = K_VI_upper, K_VI_horiz = K_VI_horiz)
}


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

## -------------------------------
## Design-based DGP (from new_sim_study.R)
## -------------------------------
sample_gnedin_labels <- function(n, gamma = 0.5, check = TRUE) {
  stopifnot(n >= 1, is.finite(gamma), gamma > 0, gamma < 1)

  x <- integer(n)
  H <- 0L
  sizes <- integer(0)

  for (i in seq_len(n)) {
    if (H == 0L) {
      H <- 1L
      sizes <- 1L
      x[i] <- 1L
      next
    }

    n_curr <- i - 1L
    if (check) {
      stopifnot(length(sizes) == H, sum(sizes) == n_curr, all(sizes >= 1L))
    }

    existing_weights <- (sizes + 1) * (n_curr - H + gamma)
    new_weight <- H * (H - gamma)
    w <- c(existing_weights, new_weight)

    if (check) {
      total <- n_curr * (n_curr + gamma)
      if (!all(is.finite(w)) || any(w < 0)) stop("Non-finite or negative weights.")
      if (abs(sum(w) - total) > 1e-8 * max(1, total)) stop("Weight sum check failed.")
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

draw_gnedin_sizes_given_K <- function(n, K_target, gamma = 0.6, max_tries = 1e6) {
  for (t in seq_len(max_tries)) {
    z <- sample_gnedin_labels(n, gamma, check = FALSE)$x
    if (length(unique(z)) == K_target) {
      return(tabulate(z, nbins = K_target))
    }
  }
  stop("Failed to draw GN sizes given K. Increase max_tries or adjust gamma.")
}

make_z_activity_to_smallest <- function(sizes, activity, seed = NULL, shuffle_within_cluster = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(sum(sizes) == length(activity))

  K <- length(sizes)
  ord_sizes <- order(sizes, decreasing = FALSE)
  sizes_sorted <- sizes[ord_sizes]
  ord_players <- order(activity, decreasing = TRUE)

  z <- integer(length(activity))
  idx <- 1L
  for (c in seq_len(K)) {
    m <- sizes_sorted[c]
    players_c <- ord_players[idx:(idx + m - 1L)]
    if (shuffle_within_cluster && m > 1L) players_c <- sample(players_c)
    z[players_c] <- c
    idx <- idx + m
  }
  z
}

check_top_player_in_smallest_cluster <- function(z, activity) {
  sizes <- tabulate(z, nbins = length(unique(z)))
  min_size <- min(sizes)
  top_player <- which.max(activity)
  if (sizes[z[top_player]] != min_size) stop("Constraint violated: top-activity player not in a smallest cluster.")
  TRUE
}

make_lambda_from_padj <- function(K, p_adj = 0.85, base = 1) {
  stopifnot(K >= 2, p_adj > 0.5, p_adj < 1)
  delta <- qlogis(p_adj)
  base * exp(delta * (0:(K - 1)))
}

simulate_W_from_blocks <- function(N, z_true, lambda_by_label, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(N)
  W <- matrix(0L, n, n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      nij <- N[i, j]
      if (nij > 0) {
        li <- lambda_by_label[z_true[i]]
        lj <- lambda_by_label[z_true[j]]
        p_ij <- li / (li + lj)
        w_ij <- rbinom(1, size = nij, prob = p_ij)
        W[i, j] <- w_ij
        W[j, i] <- nij - w_ij
      }
    }
  }
  diag(W) <- 0L
  W
}

simulate_dataset_design_based <- function(N, K_target, gamma_true = 0.71, p_adj = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(is.matrix(N), nrow(N) == ncol(N))

  n <- nrow(N)
  activity <- rowSums(N)

  sizes <- draw_gnedin_sizes_given_K(n, K_target = K_target, gamma = gamma_true)
  z_true <- make_z_activity_to_smallest(sizes, activity, seed = seed, shuffle_within_cluster = TRUE)
  check_top_player_in_smallest_cluster(z_true, activity)

  # IMPORTANT: label 1 is smallest cluster -> give it the BIGGEST lambda
  lambda_by_label <- rev(make_lambda_from_padj(K_target, p_adj = p_adj))
  W <- simulate_W_from_blocks(N, z_true, lambda_by_label, seed = seed)
  lambda_true_players <- norm_geo1(lambda_by_label[z_true])

  list(
    N = N, W = W,
    z_true = z_true,
    lambda_by_label = lambda_by_label,
    lambda_true_players = lambda_true_players,
    sizes = tabulate(z_true, nbins = K_target),
    activity = activity,
    gamma_true = gamma_true,
    p_adj = p_adj
  )
}

## Nice comparison metrics
compare_lambda <- function(est, truth) {
  est <- as.numeric(est)
  truth <- as.numeric(truth)
  stopifnot(length(est) == length(truth))

  eps <- .Machine$double.eps
  est_pos <- pmax(est, eps)
  truth_pos <- pmax(truth, eps)
  rel_err <- (est_pos - truth_pos) / truth_pos
  log_ratio <- log(est_pos / truth_pos)
  
  tibble::tibble(
    rmse = sqrt(mean((est - truth)^2)),
    mae  = mean(abs(est - truth)),
    rmse_rel = sqrt(mean(rel_err^2)),
    mae_rel = mean(abs(rel_err)),
    mean_abs_log_ratio = mean(abs(log_ratio)),
    rms_log_ratio = sqrt(mean(log_ratio^2)),
    spearman = suppressWarnings(cor(est, truth, method = "spearman")),
    kendall  = suppressWarnings(cor(est, truth, method = "kendall")),
    pearson  = suppressWarnings(cor(est, truth, method = "pearson"))
  )
}

## -------------------------------
## 1) Simulate data (DGP from new_sim_study.R)
## -------------------------------
seed <- as.integer(Sys.getenv("SEED", "123"))      # same default seed base as new_sim_study.R
gamma_true <- as.numeric(Sys.getenv("GAMMA_TRUE", "0.71"))
p_adj <- as.numeric(Sys.getenv("P_ADJ", "0.85"))

N_list <- readRDS("./N_list_for_simulation.rds")
if (!is.list(N_list) || length(N_list) < 1L) stop("N_list_for_simulation.rds must contain a non-empty list of N designs")

design_name <- Sys.getenv("SIM_DESIGN", names(N_list)[1])
if (!design_name %in% names(N_list)) design_name <- names(N_list)[1]

K_true <- as.integer(Sys.getenv("SIM_K_TRUE", "3"))
K_true <- max(2L, K_true)

N <- as.matrix(N_list[[design_name]])
if (!is.matrix(N) || nrow(N) != ncol(N)) stop("Selected design is not a square matrix: ", design_name)

sim <- simulate_dataset_design_based(
  N,
  K_target = K_true,
  gamma_true = gamma_true,
  p_adj = p_adj,
  seed = seed
)

w_ij <- sim$W
x_true <- sim$z_true
lambda_true <- sim$lambda_true_players
n <- nrow(w_ij)

player_names <- paste0("p", seq_len(n))

cat("Using DGP from new_sim_study.R\n")
cat("Design:", design_name, "  n=", n, "  K_true=", K_true, "  seed=", seed, "\n", sep = "")
cat("True cluster sizes:\n")
print(sim$sizes)
cat("Sanity check: mean(log lambda_true) =", mean(log(lambda_true)), "\n")
cat("Sanity check: geometric mean =", exp(mean(log(lambda_true))), "\n")
cat("Cluster-wise lambda_true means:\n")
print(tapply(lambda_true, x_true, mean))

cat("Total wins by node (first 6):\n")
print(head(rowSums(w_ij), 6))

Pi <- wins_to_Pi(w_ij)

## -------------------------------
## 4) Fit models + time them
## -------------------------------

T_iter <- as.integer(Sys.getenv("T_ITER", "3000"))
T_burn <- as.integer(Sys.getenv("T_BURN", "600"))
if (!is.finite(T_iter) || T_iter < 10L) stop("T_ITER must be a positive integer")
if (!is.finite(T_burn) || T_burn < 0L || T_burn >= T_iter) stop("T_BURN must be in [0, T_ITER)")

## --- Pearce & Erosheva (arXiv:2406.19563v2) recommendation ---
## RJ partition mixing can be slow => run multiple independent chains and
## check convergence across chains (Gelman-Rubin), focusing on label-invariant
## summaries (e.g., K) and worth parameters rather than raw labels.
##
## The paper does not mandate parallel execution of chains; here we optionally
## parallelize independent chains across cores for wall-clock speed.
n_chains <- as.integer(Sys.getenv("N_CHAINS", "4"))
if (!is.finite(n_chains) || n_chains < 1L) stop("N_CHAINS must be a positive integer")

slurm_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
cores_avail <- if (!is.na(slurm_cores) && slurm_cores >= 1L) slurm_cores else suppressWarnings(parallel::detectCores())
if (is.na(cores_avail) || cores_avail < 1L) cores_avail <- 1L

n_cores <- min(n_chains, max(1L, as.integer(cores_avail) - 1L))
use_parallel <- isTRUE(as.logical(Sys.getenv("USE_PARALLEL", if (n_cores > 1L) "TRUE" else "FALSE")))
chain_seeds <- 1000L + seq_len(n_chains)

run_independent_chains <- function(n_chains, FUN, seeds, n_cores = 1L, use_parallel = TRUE) {
  stopifnot(length(seeds) >= n_chains)
  idx <- seq_len(n_chains)
  if (isTRUE(use_parallel) && n_cores > 1L) {
    parallel::mclapply(
      idx,
      function(ch) FUN(seed = as.integer(seeds[ch]), chain_id = ch),
      mc.cores = n_cores
    )
  } else {
    lapply(idx, function(ch) FUN(seed = as.integer(seeds[ch]), chain_id = ch))
  }
}

norm_geo1_rows <- function(M, eps = .Machine$double.eps) {
  M <- as.matrix(M)
  M <- pmax(M, eps)
  M / exp(rowMeans(log(M)))
}

as_mcmc_safe <- function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  coda::as.mcmc(M)
}

gelman_report <- function(mcmc_list, label) {
  if (length(mcmc_list) < 2L) {
    message(label, ": need >=2 chains for Gelman-Rubin")
    return(invisible(NULL))
  }
  ml <- coda::mcmc.list(mcmc_list)
  out <- coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)
  cat("\n--- Gelman-Rubin (", label, ") ---\n", sep = "")
  print(out)
  invisible(out)
}

## 4.1 Plain BT (BTSBM) -- multiple independent chains
t_bt <- system.time({
  fit_bt_list <- run_independent_chains(
    n_chains = n_chains,
    seeds = chain_seeds,
    n_cores = n_cores,
    use_parallel = use_parallel,
    FUN = function(seed, chain_id) {
      set.seed(seed)
      BTSBM::gibbs_bt_simple(
        w_ij = w_ij,
        T_iter = T_iter,
        T_burn = T_burn,
        verbose = FALSE
      )
    }
  )
})

## 4.2 Rank-clustered (rankclust: Pearce/Erosheva)
## Run multiple independent chains and diagnose convergence across chains.
## Note: rankclust's internal sampler initializes g(0),nu(0) at random in code;
## the paper recommends singleton initialization g(0)={1,...,J}, but the current
## rankclust API does not expose g(0) directly.
t_rcbtl <- system.time({
  res_rcbtl_list <- run_independent_chains(
    n_chains = n_chains,
    seeds = chain_seeds + 10000L,
    n_cores = n_cores,
    use_parallel = use_parallel,
    FUN = function(seed, chain_id) {
      set.seed(seed)
      res_btl <- mcmc_BTL(
        Pi = Pi,
        J = n,
        a_gamma = 1,
        b_gamma = 1,
        num_iters = T_iter,
        burn_prop = T_burn / T_iter,
        chains = 1,
        groupwise = TRUE,
        seed = seed
      )

      nu0 <- colMeans(res_btl[, grep("^omega", names(res_btl)), drop = FALSE])

      res <- mcmc_RCBTL(
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
        seed = seed,
        normalize_omega = TRUE
      )

      res$chain <- factor(chain_id)
      res
    }
  )
})

## 4.3 BT-SBM (BTSBM::gibbs_bt_sbm) -- multiple independent chains
## Use singleton initialization (all items separate) as recommended in Pearce/Erosheva
## for RJ partition samplers (easier to merge than to recover missed splits).
t_btsbm <- system.time({
  fit_btsbm_list <- run_independent_chains(
    n_chains = n_chains,
    seeds = chain_seeds + 20000L,
    n_cores = n_cores,
    use_parallel = use_parallel,
    FUN = function(seed, chain_id) {
      set.seed(seed)
      BTSBM::gibbs_bt_sbm(
        w_ij = as.matrix(w_ij),
        T_iter = T_iter,
        T_burn = T_burn,
        init_x = seq_len(n),
        a = 2,
        gamma_GN = 0.8,
        prior = "GN",
        verbose = FALSE
      )
    }
  )
})

timings <- tibble::tibble(
  model = c(
    paste0("BT (BTSBM::gibbs_bt_simple) x", n_chains, " chains"),
    paste0("RCBTL (rankclust) x", n_chains, " chains"),
    paste0("BT-SBM (BTSBM::gibbs_bt_sbm) x", n_chains, " chains")
  ),
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
bt_draws_all <- do.call(rbind, lapply(fit_bt_list, function(f) f$lambda_samples))
bt_lambda_hat <- colMeans(bt_draws_all)
bt_lambda_hat <- norm_geo1(bt_lambda_hat)

## RCBTL (rankclust calls them omega)
res_rcbtl <- dplyr::bind_rows(res_rcbtl_list)
rcbtl_omega_cols <- grep("^omega", names(res_rcbtl), value = TRUE)
rcbtl_lambda_hat <- colMeans(as.matrix(res_rcbtl[, rcbtl_omega_cols, drop = FALSE]))
rcbtl_lambda_hat <- norm_geo1(rcbtl_lambda_hat)

## BTSBM
## NOTE: if you relabel for identifiability of clusters, do it before summarising lambdas.
## Your earlier code used relabel_by_lambda; keep it if it’s needed for x, but lambda itself
## is per-node, so we mostly just scale-normalise.
btsbm_chain_summaries <- lapply(fit_btsbm_list, function(fit) {
  relabel_by_lambda(fit$x_samples, fit$lambda_samples)
})
lambda_relabeled_all <- do.call(rbind, lapply(btsbm_chain_summaries, function(inf) inf$lambda_samples_relabel))
btsbm_lambda_hat <- colMeans(lambda_relabeled_all)
btsbm_lambda_hat <- norm_geo1(btsbm_lambda_hat)

## -------------------------------
## 5b) Convergence diagnostics across chains (Gelman-Rubin)
## -------------------------------

## BT: worth parameters (log-scale, per-iter geo-mean normalized)
bt_mcmc_list <- lapply(fit_bt_list, function(f) {
  M <- log(norm_geo1_rows(f$lambda_samples))
  as_mcmc_safe(M)
})
gelman_report(bt_mcmc_list, "BT log(lambda) (geo-mean normalized)")

## RCBTL: K and omega (omega already normalized in rankclust if normalize_omega=TRUE)
rcbtl_mcmc_omega <- lapply(res_rcbtl_list, function(df) {
  omega_cols <- grep("^omega", names(df), value = TRUE)
  M <- as.matrix(df[, omega_cols, drop = FALSE])
  M <- log(norm_geo1_rows(M))
  as_mcmc_safe(M)
})
gelman_report(rcbtl_mcmc_omega, "RCBTL log(omega) (geo-mean normalized)")

rcbtl_mcmc_K <- lapply(res_rcbtl_list, function(df) {
  as_mcmc_safe(matrix(as.numeric(df$K), ncol = 1, dimnames = list(NULL, "K")))
})
gelman_report(rcbtl_mcmc_K, "RCBTL K")

## BT-SBM: worth parameters in player space via relabel_by_lambda(), and K
btsbm_mcmc_lambda <- lapply(btsbm_chain_summaries, function(inf) {
  M <- log(norm_geo1_rows(inf$lambda_samples_relabel))
  as_mcmc_safe(M)
})
gelman_report(btsbm_mcmc_lambda, "BT-SBM log(lambda_player) (geo-mean normalized)")

btsbm_mcmc_K <- lapply(btsbm_chain_summaries, function(inf) {
  as_mcmc_safe(matrix(as.numeric(inf$n_clusters_each_iter), ncol = 1, dimnames = list(NULL, "K")))
})
gelman_report(btsbm_mcmc_K, "BT-SBM K")

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
) %>% select(model, everything())

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

print(rank_comp)
