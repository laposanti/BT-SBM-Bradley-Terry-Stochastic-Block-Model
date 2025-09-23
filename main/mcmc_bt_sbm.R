
########################################################################
# MAIN MCMC FUNCTION
########################################################################
compress_by_lambda_desc <- function(x, lambda) {
  L <- length(lambda)
  csize <- tabulate(x, nbins = L)
  occupied <- which(csize > 0L)
  H <- length(occupied)
  if (H == 0L) stop("No occupied clusters; x is empty?")
  
  lam_occ <- lambda[occupied]
  ord <- order(lam_occ, decreasing = TRUE)     # sort by lambda desc
  new_order <- occupied[ord]                   # old labels in new rank order
  
  # old label -> new label map: top lambda -> 1, next -> 2, ...
  map <- integer(L); map[new_order] <- seq_len(H)
  
  x_new <- map[x]
  
  lambda_new <- rep(NA_real_, L)
  lambda_new[seq_len(H)] <- lambda[new_order]
  
  list(x = x_new, lambda = lambda_new)
}

gibbs_bt_sbm <- function(w_ij, n_ij,
                         a = 0.1, b = 1,
                         prior = c("DP","PY","DM","GN"),
                         alpha_PY = 1, sigma_PY = 0,
                         beta_DM  = 1, H_DM = 20,
                         gamma_GN = 0.5,
                         n_iter = 4000, burnin = 2000,
                         init_x  = NULL, store_z = FALSE,
                         verbose = TRUE) {
  
  stopifnot(all(dim(w_ij) == dim(n_ij)),
            nrow(w_ij) == ncol(w_ij))
  N <- nrow(w_ij)
  stopifnot(N == ncol(w_ij))
  L <- N
  prior <- match.arg(prior)
  
  ## fast helpers ----------------------------------------------------
  w_i        <- rowSums(w_ij)                 # total wins per player
  sumZi_row  <- function(mat) rowSums(mat)    # later alias
  
  ## starting values -------------------------------------------------
  if (is.null(init_x)) x <- sample.int(L, N, TRUE) else x <- init_x
  lambda <- rgamma(L, a, b)
  Z      <- matrix(0, N, N)
  
  ## storage ---------------------------------------------------------
  keep  <- n_iter - burnin
  xsamp <- matrix(NA_integer_, keep, N)
  lsamp <- matrix(NA_real_,    keep, N)
  zsamp <- if (store_z) array(0, c(keep, N, N)) else NULL
  keep_id <- 0L
  
  ## choose urn function --------------------------------------------
  urn_fun <- switch(prior,
                    DP = function(v) urn_DP(v, alpha_PY),
                    PY = function(v) urn_PY(v, alpha_PY, sigma_PY),
                    DM = function(v) urn_DM(v, beta_DM, H_DM),
                    GN = function(v) urn_GN(v, gamma_GN))
  
  ## helper for the collapsed marginal (new cluster) ----------------
  log_marg_new <- function(wi, Zi) {
    (a * log(b)) + lgamma(a + wi) -
      lgamma(a)   - (a + wi) * log(b + Zi)
  }
  
  ## Gibbs loop ------------------------------------------------------
  for (it in seq_len(n_iter)) {
    
    ## --- 1. latent Z ----------------------------------------------
    for (i in 1:(N-1)) {
      li <- lambda[x[i]]
      for (j in (i+1):N) {
        nij <- n_ij[i, j]
        if (nij > 0) {
          lj <- lambda[x[j]]
          Zij <- rgamma(1, nij, li + lj)
          Z[i, j] <- Z[j, i] <- Zij
        } else Z[i, j] <- Z[j, i] <- 0
      }
    }
    
    ## --- 2. single‑site update for each player --------------------
    for (i in seq_len(N)) {
      csize    <- tabulate(x[-i], nbins = L)
      occupied <- which(csize > 0L)
      H        <- length(occupied)
      wm <- w_i[i];  Zi <- sum(Z[i, ])
      
      prior_w <- urn_fun(csize[occupied])
      new_w   <- prior_w[H + 1]
      old_w   <- prior_w[1:H]
      
      logp <- numeric(H + 1)
      
      ## likelihood * prior for existing clusters
      for (h in seq_len(H)) {
        k  <- occupied[h]
        lk <- lambda[k]
        logp[h] <- if (old_w[h] == 0) -Inf else
          log(old_w[h]) + wm * log(lk) - lk * Zi
      }
      
      ## new cluster option
      logp[H + 1] <- if (new_w == 0) -Inf else
        log(new_w) + log_marg_new(wm, Zi)
      
      pr <- exp(logp - max(logp))
      pr <- pr / sum(pr)
      dest <- sample.int(H + 1, 1, prob = pr)
      
      # inside the single-site update for item i:
      csize    <- tabulate(x[-i], nbins = L)   # counts AFTER removing i
      occupied <- which(csize > 0L)
      H        <- length(occupied)
      if (dest <= H) {
        ## move to existing cluster
        x[i] <- occupied[dest]
      } else {
        ## create / reuse empty label
        empty <- which(csize == 0L)
        stopifnot(length(empty) >= 1L)     # should always hold with nbins = L = N
        newlab <- empty[1L]
        x[i]   <- newlab
        lambda[newlab] <- rgamma(1, a + wm, b + Zi)
      }
    }
    

    ## --- 3. update lambda for each occupied block ----------------------
    csize <- tabulate(x, nbins = L)
    Zrs   <- rowSums(Z)
    for (k in which(csize > 0L)) {
      mem      <- which(x == k)
      shape_k  <- a + sum(w_i[mem])
      rate_k   <- b + sum(Zrs[mem])
      lambda[k] <- rgamma(1, shape_k, rate_k)
    }
    
    cmp <- compress_by_lambda_desc(x, lambda)
    x      <- cmp$x
    lambda <- cmp$lambda
    
    occ <- sort(unique(x))
    H   <- length(occ)
    b   <- rgamma(1, shape = 1 + H * a,
                  rate  = 1 + sum(lambda[occ]))
    ## --- 5. store -------------------------------------------------
    if (it > burnin) {
      keep_id <- keep_id + 1L
      xsamp[keep_id, ] <- x
      
      # mark empty labels as NA in storage only (no effect on the chain)
      csize <- tabulate(x, nbins = L)
      lam_store <- lambda
      lam_store[csize == 0L] <- NA_real_
      lsamp[keep_id, ] <- lam_store
      
      if (store_z) zsamp[keep_id, , ] <- Z
    }
    
    if (verbose && it %% 500 == 0)
      cat("iter", it, "  #blocks:",
          length(unique(x)), "\n")
  }
  
  list(x_samples = xsamp,
       lambda_samples = lsamp,
       z_samples = zsamp)
}

