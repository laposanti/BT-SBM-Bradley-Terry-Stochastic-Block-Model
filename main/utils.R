
# A helper function for name formatting
clean_players_names = function(name) {
  name_parts <- unlist(strsplit(name, " "))
  if (length(name_parts) == 1) return(name) 
  first_initial <- substr(name_parts[1], 1, 1)
  surname <- name_parts[length(name_parts)]
  paste(surname, first_initial, ".", sep = " ")
}


library(mcclust)
library(mcclust.ext)

inference_helper <- function(x_samples, lambda_samples){
  n_iter <- nrow(x_samples)
  K <- ncol(x_samples)
  
  new_x_samples <- matrix(NA, n_iter, K)
  new_lambda_samples <- matrix(NA, n_iter, K)
  
  for (iter in seq_len(n_iter)) {
    occupied <- sort(unique(x_samples[iter, ]))
    current_lambda <- lambda_samples[iter, occupied]
    ord <- order(current_lambda, decreasing = TRUE)
    sorted_lambda <- current_lambda[ord]
    
    label_map <- rep(NA_integer_, K)
    for (r in seq_along(ord)) {
      old_label <- occupied[ord[r]]
      label_map[old_label] <- r
    }
    
    new_x_samples[iter, ] <- label_map[x_samples[iter, ]]
    for (i in seq_len(K)) {
      new_label_i <- new_x_samples[iter, i]
      if (!is.na(new_label_i) && new_label_i >= 1 && new_label_i <= length(sorted_lambda)) {
        new_lambda_samples[iter, i] <- sorted_lambda[new_label_i]
      } else {
        new_lambda_samples[iter, i] <- NA
      }
    }
  }
  
  co_clust <- comp.psm(new_x_samples)
  partition_binder <- minbinder.ext(psm = co_clust)$cl
  partition_minVI <- minVI(psm = co_clust)$cl
  partition_expected <- apply(new_x_samples, 2, median)
  
  lambda_means <- colMeans(new_lambda_samples)
  lambda_medians <- apply(new_lambda_samples, 2, median)
  
  n_clust_vec <- apply(x_samples, 1, function(z) length(unique(z)))
  avg_n_clust <- mean(n_clust_vec)
  
  assignment_probs <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (k in seq_len(K)) {
      assignment_probs[i, k] <- mean(new_x_samples[, i] == k)
    }
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(K))
  rownames(assignment_probs) <- paste0("Player_", seq_len(K))
  assignment_probs_df <- as.data.frame(assignment_probs)
  
  count_cl <- function(z) length(unique(z[!is.na(z)]))
  vector_cl <- apply(new_x_samples, 1, count_cl)
  block_count_freq <- table(vector_cl)
  block_count_df <- data.frame(
    num_blocks = as.numeric(names(block_count_freq)),
    count = as.vector(block_count_freq),
    prob = as.vector(block_count_freq) / sum(block_count_freq)
  )
  
  top_block_count_per_iter <- rowSums(new_x_samples == 1, na.rm = TRUE)
  avg_top_block_count <- mean(top_block_count_per_iter)
  
  return(list(
    partition_binder = partition_binder,
    partition_expected = partition_expected,
    x_samples_relabel = new_x_samples,
    lambda_samples_relabel = new_lambda_samples,
    co_clustering = co_clust,
    minVI_partition = partition_minVI,
    n_clusters_each_iter = n_clust_vec,
    avg_n_clusters = avg_n_clust,
    lambda_means = lambda_means,
    lambda_medians = lambda_medians,
    player_block_assignment_probs = assignment_probs_df,
    block_count_distribution = block_count_df,
    top_block_count_per_iter = top_block_count_per_iter,
    avg_top_block_count = avg_top_block_count
  ))
}


HGnedin <- function(V, h, gamma=0.5){
  exp(lchoose(V, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(V+ gamma - h) - lgamma(V +gamma))
}


urn_DM<- function(v_minus, beta_DM, H_DM){
  # Dirichlet-Multinomial with maximum H_DM clusters
  # Existing cluster weight = n_c + beta_DM
  # New cluster weight = beta_DM*(H_DM - H)*(H_DM>H)
  # so if we haven't used all H_DM clusters yet, there's a positive new cluster weight; 
  # otherwise it's zero
  
  H <- length(v_minus)
  # from your snippet: c(v_minus+beta_DM, beta_DM*(H_DM-H)*(H_DM>H))
  # (H_DM>H) is a logical that acts like 1 if true, 0 if false
  return(c(v_minus + beta_DM, beta_DM * (H_DM - H) * (H_DM > H)))
}



urn_DP<- function(v_minus, alpha_PY){
  # v_minus is a vector of cluster sizes (excluding the node being moved)
  # length(v_minus) = number of currently active clusters
  # Return:
  #  c( (existing cluster weights), (new cluster weight) )
  
  # Dirichlet Process: existing weight = n_c, new weight = alpha
  # But in the snippet, they do "c(v_minus, alpha_PY)" i.e. n_c, alpha
  return(c(v_minus, alpha_PY))
}

urn_GN <- function(v_minus, gamma) {
  ## corrected off‑by‑one: n = sum(v) + 1 when node i joins
  H  <- length(v_minus)
  n_ <- sum(v_minus) + 1               # +1 for the focal node
  c((v_minus + 1) * (n_ - H + gamma),
    H^2 - H * gamma)
}

urn_PY<- function(v_minus, alpha_PY, sigma_PY){
  # Pitman-Yor( alpha, sigma ), with 0 < sigma < 1
  # Existing cluster weight = n_c - sigma
  # New cluster weight = alpha + sigma * (#clusters)
  # The snippet's code does "c(v_minus - sigma_PY, alpha_PY + length(v_minus)*sigma_PY)"
  
  H <- length(v_minus)
  return(c(v_minus - sigma_PY, alpha_PY + H*sigma_PY))
}


########################################################################
# Post-processing
########################################################################



inference_helper <- function(x_samples, lambda_samples){
  n_iter <- nrow(x_samples)
  K <- ncol(x_samples)    # number of players
  stopifnot(nrow(lambda_samples) == n_iter, ncol(lambda_samples) == K)
  
  # =========================================================
  # 1) RELABEL EACH ITERATION BY DESCENDING LAMBDA
  # =========================================================
  new_x_samples <- matrix(NA, n_iter, K)
  new_lambda_samples <- matrix(NA, n_iter, K)
  
  for (iter in seq_len(n_iter)) {
    #finding occupied clusters, in order, like 1 2 3,...
    occupied <- sort(unique(x_samples[iter, ])) 
    #retrieving the lambdas associated to occupied clusters,
    current_lambda <- lambda_samples[iter, occupied]
    
    ord <- order(current_lambda, decreasing = TRUE) 
    #new lambdas
    sorted_lambda <- current_lambda[ord]
    
    label_map <- rep(NA_integer_, K)
    for (r in seq_along(sorted_lambda)) {
      old_label <- occupied[ord[r]]
      label_map[old_label] <- r
    }
    
    new_x_samples[iter, ] <- label_map[x_samples[iter, ]]
    for (i in seq_len(K)) {
      new_label_i <- new_x_samples[iter, i]  
      new_lambda_samples[iter, i] <- sorted_lambda[new_label_i]
    }
    
  }
  # =========================================================
  # 2) BUILD A CO-CLUSTERING MATRIX & MINVI PARTITION
  # =========================================================
  co_clust = mcclust::comp.psm(new_x_samples)
  
  # If you're using a function like comp.psm, you can overwrite co_clust:
  # (Example usage -- you might want to focus on a subset of iterations, etc.)
  # co_clust <- comp.psm(new_x_samples[2000:4000, ])
  # Adjust to your actual iteration range as you see fit:
  # co_clust = comp.psm(new_x_samples[2000:4000,])
  binder = minbinder.ext(psm = co_clust)
  partition_binder = binder$cl
  vi_results <- minVI(psm = co_clust)
  partition_minVI <- vi_results$cl
  partition_expected <- apply(new_x_samples, 2, median)
  # =========================================================
  # 3) EXTRACT POINT ESTIMATES OF LAMBDA
  # =========================================================
  lambda_means <- colMeans(new_lambda_samples)
  lambda_medians <- apply(new_lambda_samples, 2, median)
  
  # The average number of occupied clusters per iteration:
  n_clust_vec <- integer(n_iter)
  for (iter in seq_len(n_iter)) {
    n_clust_vec[iter] <- length(unique(new_x_samples[iter, ]))
  }
  avg_n_clust <- mean(n_clust_vec)
  
  # =========================================================
  # 4) COMPUTE ASSIGNMENT PROBABILITIES
  #    i.e., P(player i is in cluster k), for k = 1..K
  # =========================================================
  assignment_probs <- matrix(0, nrow = K, ncol = K)
  for (player in seq_len(K)) {
    for (k in seq_len(K)) {
      assignment_probs[player, k] <- mean(new_x_samples[, player] == k)
    }
  }
  # Turn into a data frame for convenience
  assignment_probs_df <- as.data.frame(assignment_probs)
  colnames(assignment_probs_df) <- paste0("Cluster_", seq_len(K))
  rownames(assignment_probs_df) <- paste0("Player_", seq_len(K))
  
  # =========================================================
  # 5) COMPUTE POSTERIOR PROBABILITY FOR EACH NUMBER OF BLOCKS
  #    (based on new_x_samples or old x_samples; shown with new_x_samples here)
  # =========================================================
  # We'll define a helper that counts the number of distinct clusters in one iteration row:
  count_cl <- function(row_labels) length(unique(row_labels[!is.na(row_labels)]))
  
  vector_cl <- apply(new_x_samples, 1, count_cl)
  # Frequency table for each possible # blocks
  block_count_freq <- table(vector_cl)
  # Convert to data frame with probabilities
  block_count_df <- data.frame(
    num_blocks = as.numeric(names(block_count_freq)),
    count      = as.vector(block_count_freq),
    prob       = as.vector(block_count_freq) / sum(block_count_freq)
  )
  
  # =========================================================
  # 6) NUMBER OF PLAYERS ASSIGNED TO THE TOP BLOCK
  #    (Block #1 has the highest lambda by construction)
  # =========================================================
  top_block_count_per_iter <- rowSums(new_x_samples == 1, na.rm = TRUE)
  avg_top_block_count <- mean(top_block_count_per_iter)
  
  # Finally, return everything (including the new outputs).
  return(list(
    partition_binder              = partition_binder,
    partition_expected            = partition_expected,
    x_samples_relabel             = new_x_samples,
    lambda_samples_relabel        = new_lambda_samples,
    co_clustering                 = co_clust,
    minVI_partition               = partition_minVI,
    n_clusters_each_iter          = n_clust_vec,
    avg_n_clusters                = avg_n_clust,
    lambda_means                  = lambda_means,
    lambda_medians                = lambda_medians,
    player_block_assignment_probs = assignment_probs_df,  # (1) a posteriori P for each player-block
    block_count_distribution      = block_count_df,        # (2) prob for each # blocks
    top_block_count_per_iter      = top_block_count_per_iter, 
    avg_top_block_count           = avg_top_block_count    # (3) number of players in top block
  ))
}




# Helper to compute the expected number of clusters under different gibbs-type p.

expected_cl_py <- function(n, sigma, theta, H){
  n <- as.integer(n)
  stopifnot(sigma >= 0, sigma < 1, theta > - sigma, n > 0, H > 1)
  if(H == Inf){
    if(sigma==0) {
      out <- theta * sum(1/(theta - 1 + 1:n))
    } else {
      out <- 1/sigma*exp(lgamma(theta + sigma + n) - lgamma(theta + sigma) - lgamma(theta + n) + lgamma(theta + 1)) - theta/sigma
    }
  } else if(H < Inf){
    if(sigma==0) {
      index <- 0:(n-1)
      out <- H - H*exp(sum(log(index + theta*(1 - 1/H)) - log(theta+ index)))
    }
  }
  return(out)
}

gibbs_bt_simple <- function(w_ij, n_ij,
                            a = 0.01, b = 1,
                            n_iter = 5000, burnin = 1000,
                            verbose = TRUE) {
  K <- nrow(w_ij)
  stopifnot(K == ncol(w_ij),
            K == nrow(n_ij), K == ncol(n_ij))
  
  # pre‑computation
  w_i <- rowSums(w_ij)
  
  # initial values
  lambda <- rgamma(K, a, b)
  Z      <- matrix(0, K, K)
  
  n_keep <- n_iter - burnin
  lambda_store <- matrix(NA_real_, n_keep, K)
  keep_idx <- 0L
  
  for (iter in seq_len(n_iter)) {
    # ----- latent Z_{ij} -----
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        nij <- n_ij[i, j]
        if (nij > 0) {
          rate  <- lambda[i] + lambda[j]
          Z_ij  <- rgamma(1, nij, rate)
          Z[i,j] <- Z[j,i] <- Z_ij
        } else {
          Z[i,j] <- Z[j,i] <- 0
        }
      }
    }
    
    # ----- update λ_i -----
    for (i in seq_len(K)) {
      shape <- a + w_i[i]
      rate  <- b + sum(Z[i, ])
      lambda[i] <- rgamma(1, shape, rate)
    }
    
    # store
    if (iter > burnin) {
      keep_idx <- keep_idx + 1L
      lambda_store[keep_idx, ] <- lambda
    }
    
    if (verbose && iter %% 1000 == 0)
      cat("simple BT – iter", iter, "\n")
  }
  
  invisible(list(lambda_samples = lambda_store))
}

# ============================================================
# 2. Log‑likelihood matrix for LOO — *simple* model
#    Each unordered pair (i<j) is one observation with
#    y = wins_i_vs_j, n = total matches.
# ============================================================
make_bt_simple_loo <- function(w_ij, n_ij, lambda_samples) {
  K   <- nrow(w_ij)
  S   <- nrow(lambda_samples)
  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D   <- nrow(idx)                    # #data points
  ll  <- matrix(NA_real_, S, D)
  
  for (s in seq_len(S)) {
    lam <- lambda_samples[s, ]
    for (d in seq_len(D)) {
      i <- idx[d,1]; j <- idx[d,2]
      n <- n_ij[i,j]; w <- w_ij[i,j]
      p <- lam[i] / (lam[i] + lam[j])
      ll[s, d] <- dbinom(w, n, p, log = TRUE)
    }
  }
  
  structure(list(ll = ll, obs_idx = idx), class = "bt_ll_matrix")
}

# ============================================================
# 3. Log‑likelihood matrix for LOO — *clustered* model
#    Needs *both* λ_k and assignments x_i.
#    x_samples must be relabelled consistently with λ_samples.
# ============================================================
make_bt_cluster_loo <- function(w_ij, n_ij, lambda_samples, x_samples) {
  stopifnot(nrow(lambda_samples) == nrow(x_samples))
  K   <- ncol(x_samples)
  S   <- nrow(x_samples)
  idx <- which(upper.tri(n_ij) & n_ij > 0, arr.ind = TRUE)
  D   <- nrow(idx)
  ll  <- matrix(NA_real_, S, D)
  
  for (s in seq_len(S)) {
    lam <- lambda_samples[s, ]        # length K (after relabeling): λ_1, …, λ_K
    xs  <- x_samples[s, ]             # assignments 1..K or NA for empty
    for (d in seq_len(D)) {
      i <- idx[d,1]; j <- idx[d,2]
      n <- n_ij[i,j]; w <- w_ij[i,j]
      
      # cluster labels for this draw
      ki <- xs[i]; kj <- xs[j]
      # skip if either is NA (should not happen if relabeled)
      if (is.na(ki) || is.na(kj)) {
        ll[s, d] <- NA_real_
        next
      }
      p <- lam[ki] / (lam[ki] + lam[kj])
      ll[s, d] <- dbinom(w, n, p, log = TRUE)
    }
  }
  
  structure(list(ll = ll, obs_idx = idx), class = "bt_ll_matrix")
}

# -----------------------------------------------------------------
# 4. Convenience wrapper – returns loo objects for both models
# -----------------------------------------------------------------
compare_bt_models_loo <- function(w_ij, n_ij,
                                  simple_draws,
                                  cluster_draws) {
  if (!requireNamespace("loo", quietly = TRUE))
    stop("Package 'loo' not installed.")
  
  ll_simple  <- simple_draws$ll
  loo_simple <- loo::loo(ll_simple)
  
  ll_cluster  <- cluster_draws$ll
  loo_cluster <- loo::loo(ll_cluster)
  
  comp <- loo::compare_models(loo_simple, loo_cluster)
  list(simple = loo_simple,
       cluster = loo_cluster,
       comparison = comp)
}
# -------------------------------------------------------------------
#  make_bt_loo_cluster()
#     • x_draws       – M × K matrix of cluster labels (each row = one draw)
#     • lambda_draws  – M × K matrix of cluster-specific rates
#     • w_ij, n_ij    – K × K win counts and match counts as in the sampler
# Returns: M × P matrix where P is the number of distinct pairs (i<j)
# -------------------------------------------------------------------
make_bt_loo_cluster <- function(x_draws,
                                lambda_draws,
                                w_ij,
                                n_ij) {
  stopifnot(dim(x_draws) == dim(lambda_draws),
            nrow(w_ij)  == ncol(w_ij),
            nrow(w_ij)  == nrow(n_ij),
            ncol(w_ij)  == ncol(n_ij))
  
  M <- nrow(x_draws)               # posterior draws
  K <- ncol(x_draws)               # number of players
  
  # Index of the upper-triangular (i < j) pairs – these are our “data points”
  pair_idx <- which(upper.tri(n_ij), arr.ind = TRUE)
  P <- nrow(pair_idx)
  
  log_lik <- matrix(0, nrow = M, ncol = P)
  
  for (m in seq_len(M)) {
    # draw-specific cluster labels and rates
    x_m   <- x_draws[m, ]
    lam_m <- lambda_draws[m, ]
    
    # pre-compute player-level λ to avoid double look-ups inside the loop
    lambda_player <- lam_m[x_m]    # length K
    
    for (p in seq_len(P)) {
      i <- pair_idx[p, 1]
      j <- pair_idx[p, 2]
      nij <- n_ij[i, j]
      if (nij == 0) {
        log_lik[m, p] <- 0          # no information in that cell
      } else {
        wij <- w_ij[i, j]
        pij <- lambda_player[i] / (lambda_player[i] + lambda_player[j])
        log_lik[m, p] <- dbinom(wij, size = nij, prob = pij, log = TRUE)
      }
    }
  }
  
  # Nice column names for loo output: "(1,2)", "(1,3)", …
  colnames(log_lik) <- apply(pair_idx, 1, paste, collapse = ",")
  
  log_lik
}

shannon_entropy <- function(p) {
  probs <- p / sum(p)
  probs <- probs[probs > 0]  # remove zero probabilities to avoid log(0)
  -sum(probs * log(probs))
}


