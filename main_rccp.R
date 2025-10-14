# ----------------------- Main Gibbs sampler (fixed) ----------------------------------
# Key changes vs your original:
# - Separate N (items) from L (label capacity). Never change N; L can grow.
# - Robust label creation: grow lambda_curr only; assignments x_curr stay length N.
# - All tabulations over labels use nbins = L.
# - Store lambdas in a ragged list (no column mismatch) + keep L_hist and K_occ.
# - Optional: learn a via slice sampling on log(a). If learn_b=TRUE, we skip
#   lambda re-centering; if learn_b=FALSE and center_lambda=TRUE, we center.
# - Return traces for a and b; return sufficient stats per iteration if requested.
library(Rcpp) 
Rcpp::sourceCpp(code='
#include <Rcpp.h>

// We use fully qualified types; attribute must sit directly above the function.
// [[Rcpp::export]]
Rcpp::List updateZ_and_rowsums(const Rcpp::IntegerMatrix& n_ij,
                               const Rcpp::IntegerVector& x,
                               const Rcpp::NumericVector& lambda) {
  const int K = n_ij.nrow();
  Rcpp::NumericMatrix Z(K, K);
  Rcpp::NumericVector rowSumsZ(K);

  // lambda per node via its cluster label (1..K in R)
  Rcpp::NumericVector lam_i(K);
  for (int i = 0; i < K; ++i) {
    int lab = x[i];                  // 1-based label
    if (lab <= 0 || lab > lambda.size() || Rcpp::NumericVector::is_na(lambda[lab - 1])) {
      lam_i[i] = NA_REAL;
    } else {
      lam_i[i] = lambda[lab - 1];
    }
  }

  for (int i = 0; i < K; ++i) {
    for (int j = i + 1; j < K; ++j) {
      int nij = n_ij(i, j);
      if (nij > 0) {
        double rate = lam_i[i] + lam_i[j];
        if (!R_finite(rate) || rate <= 0.0) {
          Z(i, j) = 0.0; Z(j, i) = 0.0;
        } else {
          // R::rgamma(shape, scale). We have shape=nij, *rate*=rate => scale=1/rate
          double draw = R::rgamma((double)nij, 1.0 / rate);
          Z(i, j) = draw; Z(j, i) = draw;
          rowSumsZ[i] += draw; rowSumsZ[j] += draw;
        }
      } else {
        Z(i, j) = 0.0; Z(j, i) = 0.0;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("Z") = Z,
                            Rcpp::Named("rowSumsZ") = rowSumsZ);
}
')
# --- small utilities ---------------------------------------------------------------
.safe_log <- function(x) ifelse(x > 0, log(x), -Inf)

# Univariate slice sampler on the log-scale (for log(a))
.slice_on_log <- function(
    logpost,           # function(loga) -> log posterior up to additive const
    loga0,             # current value
    w = 1.0,           # initial bracket width
    m = 20L,           # max stepping-out steps
    lower = -Inf,      # bounds on loga
    upper =  Inf
) {
  f0 <- logpost(loga0)
  if (!is.finite(f0)) stop("Initial loga has -Inf logpost.")
  # Draw vertical level
  y <- f0 - rexp(1)
  # Initial interval
  L <- loga0 - runif(1, 0, w)
  R <- L + w
  # Step out
  J <- floor(runif(1) * m); K <- (m - 1L) - J
  while (J > 0L && L > lower && logpost(L) > y) { L <- L - w; J <- J - 1L }
  while (K > 0L && R < upper && logpost(R) > y) { R <- R + w; K <- K - 1L }
  # Shrinkage
  for (iter in 1:100) {
    loga1 <- runif(1, max(L, lower), min(R, upper))
    f1 <- logpost(loga1)
    if (is.finite(f1) && f1 >= y) return(loga1)
    if (loga1 < loga0) L <- loga1 else R <- loga1
  }
  warning("Slice sampler reached max shrink steps; returning last draw.")
  return(loga1)
}
.safe_log <- function(x) ifelse(x > 0, log(x), -Inf)
.logsumexp <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x - m))) }

# Compute integrated log-likelihood for a NEW cluster (Gamma-Poisson integral)
.new_cluster_integral_log <- function(wi, Zi, a_now, b_now) {
  (a_now * log(b_now)) + lgamma(a_now + wi) - lgamma(a_now) - (a_now + wi) * log(b_now + Zi)
}
.sample_from_logweights <- function(logw) {
  # returns an index in 1..length(logw); never passes NA probs to sample()
  if (all(!is.finite(logw))) {
    # fall back to uniform over all positions
    return(sample.int(length(logw), 1L))
  }
  lse <- .logsumexp(logw)
  if (!is.finite(lse)) {
    return(sample.int(length(logw), 1L))
  }
  p <- exp(logw - lse)
  # clean numerical dust
  p[!is.finite(p)] <- 0
  s <- sum(p)
  if (s <= 0) {
    return(sample.int(length(p), 1L))
  } else {
    p <- p / s
    return(sample.int(length(p), 1L, prob = p))
  }
}
# --- anchored collapsed score for NEW cluster (drops the a-only constant) ---
.new_cluster_integral_log_anchored <- function(wi, Zi, a_now, b_now) {
  # equals: logGamma(a+wi) - (a+wi)*log(b+Zi)
  lgamma(a_now + wi) - (a_now + wi) * log(b_now + Zi)
}

# ----------------------- Main Gibbs sampler ----------------------------------
gibbs_bt_urn_rcpp <- function(
    w_ij,
    a = 1, b = 1,
    prior = c("DP", "PY", "DM", "GN"),
    alpha_PY = NA_real_,
    sigma_PY = NA_real_,
    beta_DM  = NA_real_,
    H_DM     = NA_integer_,
    gamma_GN = NA_real_,
    # hyperpriors / learning switches
    learn_b = F, e_b = 1, f_b = 1,     # b ~ Gamma(e_b, f_b)
    learn_a = F, a_shape = 1, a_rate = 1,  # a ~ Gamma(a_shape, a_rate), via slice on log(a)
    n_iter = 2000, burnin = 1000,
    init_x = NULL,
    store_z = FALSE,
    store_suff = FALSE,     # store sufficient stats (S1,S2) per iter (list)
    verbose = TRUE
) {
  stopifnot(is.matrix(w_ij))
  N <- nrow(w_ij)
  stopifnot(N == ncol(w_ij))
  if (any(w_ij < 0)) stop("Negative entries not allowed.")
  n_ij <- w_ij + t(w_ij)
  if (any(w_ij > n_ij)) stop("w_ij must be <= n_ij elementwise.")
  if (burnin >= n_iter) stop("burnin < n_iter required.")
  if (a <= 0 || b <= 0) stop("Gamma(a,b) needs a>0, b>0.")
  
  prior <- match.arg(prior)
  if (prior == "DP" && !is.finite(alpha_PY)) stop("DP requires alpha_PY.")
  if (prior == "PY" && (!is.finite(alpha_PY) || !is.finite(sigma_PY) || sigma_PY <= 0 || sigma_PY >= 1))
    stop("PY requires alpha_PY and sigma_PY in (0,1).")
  if (prior == "DM" && (!is.finite(beta_DM) || !is.finite(H_DM) || H_DM < 1))
    stop("DM requires beta_DM>0 and H_DM>=1.")
  if (prior == "GN" && !is.finite(gamma_GN)) stop("GN requires gamma_GN.")
  
  urn_fun <- switch(
    prior,
    DP = function(v) urn_DP(v, alpha_PY),
    PY = function(v) urn_PY(v, alpha_PY, sigma_PY),
    DM = function(v) urn_DM(v, beta_DM, H_DM),
    GN = function(v) urn_GN(v, gamma_GN)
  )
  
  # precompute row wins
  w_i <- rowSums(w_ij)
  
  # --- init labels and rates ------------------------------------------------------
  L <- N                            # label capacity (can grow)
  if (is.null(init_x)) x_curr <- sample.int(L, N, replace = TRUE) else {
    stopifnot(length(init_x) == N)
    x_curr <- as.integer(init_x)
    L <- max(L, max(x_curr))
  }
  a_curr <- a; 
  b_eff <- if (learn_b) {
    b
  } else{
    as.numeric(exp(digamma(a_curr)))   # b := exp(psi(a))
  }
  lambda_curr <- rep(NA_real_, L)
  csize0 <- tabulate(x_curr, nbins = L)
  occ0 <- which(csize0 > 0L)
  lambda_curr[occ0] <- rgamma(length(occ0), shape = a_curr, rate = b_eff)
  
  # latent Z
  Z_curr <- matrix(0, N, N)
  if (store_z) z_store <- array(0, dim = c(n_iter - burnin, N, N)) else z_store <- NULL
  
  # storage
  S <- n_iter - burnin
  x_samples <- matrix(NA_integer_, S, N)
  lambda_list <- vector("list", S)         # ragged storage for lambdas
  L_hist <- integer(S)
  K_occ <- integer(S)
  a_trace <- if (learn_a) numeric(S) else NULL
  b_trace <- if (learn_b) numeric(S) else NULL
  suff_list <- if (store_suff) vector("list", S) else NULL
  
  # helpers
  logpost_loga <- function(loga, lam_occ) {
    a_now <- exp(loga)
    if (!is.finite(a_now) || a_now <= 0) return(-Inf)
    Kocc <- length(lam_occ)
    val <- sum((a_now - 1) * log(lam_occ)) - Kocc * lgamma(a_now) + (a_shape - 1) * loga - a_rate * a_now
    as.numeric(val)
  }
  
  save_i <- 0L
  for (iter in seq_len(n_iter)) {
    # ---- 1) Update Z and row sums in C++ (expects lambda length >= max(x))
    outZ <- updateZ_and_rowsums(n_ij = structure(n_ij, .Dim = dim(n_ij)),
                                x = as.integer(x_curr),
                                lambda = lambda_curr)
    Z_curr <- outZ$Z
    Z_row_sum <- outZ$rowSumsZ
    
    # ---- 2) Single-site updates for x_i
    for (i in seq_len(N)) {
      csize_minus <- tabulate(x_curr[-i], nbins = L)
      occupied <- which(csize_minus > 0L)
      H <- length(occupied)
      v_minus <- csize_minus[occupied]
      
      prior_vec <- if (H > 0L) urn_fun(v_minus) else c(0, urn_fun(integer(0))[1L])
      new_weight <- prior_vec[H + 1L]
      existing_weights <- if (H > 0L) prior_vec[seq_len(H)] else numeric(0)
      
      Zi <- Z_row_sum[i]
      wi <- w_i[i]
      
      logp <- rep(-Inf, H + 1L)
      
      if (H > 0L) {
        lam_h <- lambda_curr[occupied]
        llh <- wi * .safe_log(lam_h) - lam_h * Zi
        lprior <- .safe_log(existing_weights)
        logp[seq_len(H)] <- lprior + llh
      }
      # new (anchored):
      if (new_weight > 0) {
        logp[H + 1L] <- .safe_log(new_weight) +
          .new_cluster_integral_log_anchored(wi, Zi, a_curr, b_eff)
      }
      
      
      choice <- .sample_from_logweights(logp)
      
      if (choice <= H) {
        x_curr[i] <- occupied[choice]
      } else {
        # create / reuse empty label id
        csize_all <- tabulate(x_curr[-i], nbins = L)
        empties <- which(csize_all == 0L)
        if (length(empties) == 0L) {
          L <- L + 1L
          lambda_curr <- c(lambda_curr, NA_real_)
          new_label <- L
        } else new_label <- empties[1L]
        x_curr[i] <- new_label
        
        shape_k <- a_curr + wi
        rate_k  <- b_eff + Zi
        lambda_curr[new_label] <- rgamma(1L, shape = shape_k, rate = rate_k)
      }
    }
    
    # ---- 3) Update lambda_k for occupied clusters (conjugate)
    csize <- tabulate(x_curr, nbins = L)
    occ <- which(csize > 0L)
    if (length(occ)) {
      for (k in occ) {
        members <- which(x_curr == k)
        shape_k <- a_curr + sum(w_i[members])
        rate_k  <- b_eff + sum(Z_row_sum[members])
        lambda_curr[k] <- rgamma(1L, shape = shape_k, rate = rate_k)
      }
    }
    lambda_curr[csize == 0L] <- NA_real_
    
    # ---- 4) Update b (conjugate) if requested
    if (learn_b) {
      Kocc <- length(occ)
      sumlam <- sum(lambda_curr[occ])
      b_eff <- rgamma(1, shape = e_b + a_curr * Kocc, rate = f_b + sumlam)
    }
    
  
    # ---- 6) Optionally learn a via slice sampling on log(a)
    occupied <- which(csize > 0L)
    if (length(occupied) > 0L) {
      g_mean <- exp(mean(log(lambda_curr[occupied])))
      lambda_curr[occupied] <- lambda_curr[occupied] / g_mean
    }
    # ---- 7) Store
    if (iter > burnin) {
      save_i <- save_i + 1L
      x_samples[save_i, ] <- x_curr
      lambda_list[[save_i]] <- lambda_curr
      L_hist[save_i] <- L
      K_occ[save_i] <- length(occ)
      if (!is.null(a_trace)) a_trace[save_i] <- a_curr
      if (!is.null(b_trace)) b_trace[save_i] <- b_eff
      if (store_z) z_store[save_i, , ] <- Z_curr
      if (store_suff) {
        # sufficient stats for each occupied label this iter
        S1 <- S2 <- rep(NA_real_, L)
        for (k in occ) {
          members <- which(x_curr == k)
          S1[k] <- sum(w_i[members])
          S2[k] <- sum(Z_row_sum[members])
        }
        suff_list[[save_i]] <- list(S1 = S1, S2 = S2)
      }
    }
    
    if (verbose && iter %% 1000L == 0L) {
      cat("iter", iter, "occupied =", length(unique(x_curr)), "a=", round(a_curr,3), "b=", round(b_eff,3), "L=", L, "\n")
    }
  }
  
  list(
    x_samples = x_samples,
    lambda_samples = lambda_list,  # ragged storage
    L_hist = L_hist,
    K_occ = K_occ,
    z_samples = if (store_z) z_store else NULL,
    suff_stats = if (store_suff) suff_list else NULL,
    hyper = list(a_final = a_curr, b_final = b_eff, a_trace = a_trace, b_trace = b_trace)
  )
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

urn_GN<- function(v_minus, gamma_GN){
  # Gnedin's process
  # snippet code: c( (v_minus+1)*(sum(v_minus)-H+gamma_GN), H^2 - H*gamma_GN )
  
  H <- length(v_minus)
  s_vminus <- sum(v_minus)
  existing_weights <- (v_minus + 1) * (s_vminus - H + gamma_GN)
  new_weight <- H^2 - H * gamma_GN
  
  return(c(existing_weights, new_weight))
}

urn_PY<- function(v_minus, alpha_PY, sigma_PY){
  # Pitman-Yor( alpha, sigma ), with 0 < sigma < 1
  # Existing cluster weight = n_c - sigma
  # New cluster weight = alpha + sigma * (#clusters)
  # The snippet's code does "c(v_minus - sigma_PY, alpha_PY + length(v_minus)*sigma_PY)"
  
  H <- length(v_minus)
  return(c(v_minus - sigma_PY, alpha_PY + H*sigma_PY))
}

# ----------------------- Helper: relabel by lambda (handles ragged lambda list) ---
# Usage:
# out <- relabel_by_lambda(x_samples, lambda_samples)
# where:
# - x_samples: S x N integer matrix of item -> label ids per iteration
# - lambda_samples: either (a) list length S, each a numeric vector of length L_t
# with NA for empty labels; or (b) S x L matrix (back-compat)
# Returns a list with relabeled x, per-item cluster lambdas, partitions, and summaries.

relabel_by_lambda <- function(x_samples, lambda_samples) {
  stopifnot(is.matrix(x_samples))
  S <- nrow(x_samples); N <- ncol(x_samples)
  
  # Support both new ragged-list and old matrix formats
  is_list_format <- is.list(lambda_samples)
  if (!is_list_format && !is.matrix(lambda_samples)) {
    stop("lambda_samples must be a list (ragged) or a matrix (legacy).")
  }
  if (!is_list_format && nrow(lambda_samples) != S) {
    stop("lambda_samples (matrix) must have same rows as x_samples.")
  }
  
  new_x_samples <- matrix(NA_integer_, S, N)
  new_lambda_samples <- matrix(NA_real_, S, N) # per-item cluster λ after relabeling
  
  for (iter in seq_len(S)) {
    xi <- x_samples[iter, ]
    occ_labels <- sort(unique(xi[!is.na(xi)]))
    H <- length(occ_labels)
    
    # Get lambda vector for *this* iteration aligned to occ_labels
    if (is_list_format) {
      lam_vec_full <- lambda_samples[[iter]]
      if (is.null(lam_vec_full)) lam_vec_full <- numeric(0)
    } else {
      lam_vec_full <- lambda_samples[iter, ]
    }
    # Pull the lambdas for occupied labels; allow NAs
    lam_occ <- rep(NA_real_, H)
    if (length(lam_vec_full)) {
      # Some labels may exceed length(lam_vec_full) if capacity grew; guard
      max_idx <- min(length(lam_vec_full), max(occ_labels))
      idx_ok <- occ_labels <= max_idx
      lam_occ[idx_ok] <- lam_vec_full[occ_labels[idx_ok]]
    }
    
    # Order occupied labels by decreasing lambda (NAs go last)
    ord <- order(lam_occ, decreasing = TRUE, na.last = TRUE)
    sorted_lambda <- lam_occ[ord]
    new_ids <- seq_len(H)
    
    # Build label map old_label -> new sequential label 1..H (according to ord)
    label_map <- rep(NA_integer_, max(c(occ_labels, 1L)))
    label_map[occ_labels[ord]] <- new_ids
    
    # Apply map to items
    xi_mapped <- ifelse(!is.na(xi) & xi <= length(label_map), label_map[xi], NA_integer_)
    new_x_samples[iter, ] <- xi_mapped
    
    # For each item, attach its cluster lambda after relabeling
    # (i.e., item takes sorted_lambda[new_label])
    li <- rep(NA_real_, N)
    ok <- !is.na(xi_mapped) & xi_mapped >= 1L & xi_mapped <= length(sorted_lambda)
    li[ok] <- sorted_lambda[xi_mapped[ok]]
    new_lambda_samples[iter, ] <- li
  }
  
  # Posterior similarity matrix and partitions
  co_clust <- mcclust::comp.psm(new_x_samples)
  partition_binder <- mcclust.ext::minbinder.ext(psm = co_clust)$cl
  partition_minVI <- mcclust.ext::minVI(psm = co_clust)$cl
  partition_expected <- apply(new_x_samples, 2, stats::median, na.rm = TRUE)
  
  # Lambda summaries (per-item assigned λ after relabeling)
  lambda_means <- colMeans(new_lambda_samples, na.rm = TRUE)
  lambda_medians <- apply(new_lambda_samples, 2, stats::median, na.rm = TRUE)
  
  # Cluster count diagnostics
  n_clust_vec <- apply(new_x_samples, 1, function(z) length(unique(z[!is.na(z)])))
  avg_n_clust <- mean(n_clust_vec)
  
  # Assignment probabilities: P(item i -> cluster k) with k up to N
  Kmax <- N
  assignment_probs <- matrix(0, nrow = N, ncol = Kmax)
  for (i in seq_len(N)) {
    for (k in seq_len(Kmax)) {
      assignment_probs[i, k] <- mean(new_x_samples[, i] == k, na.rm = TRUE)
    }
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(Kmax))
  rownames(assignment_probs) <- paste0("Item_", seq_len(N))
  assignment_probs_df <- as.data.frame(assignment_probs)
  
  # Distribution of number of blocks per iteration
  block_count_freq <- table(n_clust_vec)
  block_count_df <- data.frame(
    num_blocks = as.numeric(names(block_count_freq)),
    count = as.vector(block_count_freq),
    prob = as.vector(block_count_freq) / sum(block_count_freq)
  )
  
  # Size of top cluster per iteration (label 1 after relabeling)
  top_block_count_per_iter <- rowSums(new_x_samples == 1L, na.rm = TRUE)
  avg_top_block_count <- mean(top_block_count_per_iter)
  
  list(
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
    item_cluster_assignment_probs = assignment_probs_df,
    block_count_distribution = block_count_df,
    top_block_count_per_iter = top_block_count_per_iter,
    avg_top_block_count = avg_top_block_count
  )
}
