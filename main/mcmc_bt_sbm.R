
########################################################################
# MAIN MCMC FUNCTION
########################################################################

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
  K <- nrow(w_ij)
  prior <- match.arg(prior)
  
  ## fast helpers ----------------------------------------------------
  w_i        <- rowSums(w_ij)                 # total wins per player
  sumZi_row  <- function(mat) rowSums(mat)    # later alias
  
  ## starting values -------------------------------------------------
  if (is.null(init_x))
    x <- sample.int(K, K, TRUE) else x <- init_x
  lambda <- rgamma(K, a, b)
  Z      <- matrix(0, K, K)
  
  ## storage ---------------------------------------------------------
  keep  <- n_iter - burnin
  xsamp <- matrix(NA_integer_, keep, K)
  lsamp <- matrix(NA_real_,    keep, K)
  zsamp <- if (store_z) array(0, c(keep, K, K)) else NULL
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
    for (i in 1:(K-1)) {
      li <- lambda[x[i]]
      for (j in (i+1):K) {
        nij <- n_ij[i, j]
        if (nij > 0) {
          lj <- lambda[x[j]]
          Zij <- rgamma(1, nij, li + lj)
          Z[i, j] <- Z[j, i] <- Zij
        } else Z[i, j] <- Z[j, i] <- 0
      }
    }
    
    ## --- 2. single‑site update for each player --------------------
    for (i in seq_len(K)) {
      
      csize <- tabulate(x[-i], nbins = K)
      occupied <- which(csize > 0L)
      H  <- length(occupied)
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
      
      if (dest <= H) {
        ## move to existing cluster
        x[i] <- occupied[dest]
      } else {
        ## create / reuse empty label
        empty <- which(csize == 0L)
        if (length(empty) == 0L) {
          K <- K + 1L
          lambda <- c(lambda, NA_real_)
          empty  <- K
          Z      <- rbind(Z, 0);  Z <- cbind(Z, 0)
          n_ij   <- rbind(n_ij, 0); n_ij <- cbind(n_ij, 0)
        }
        newlab     <- empty[1]
        x[i]       <- newlab
        lambda[newlab] <- rgamma(1, a + wm, b + Zi)
      }
    }
    
    ## --- 3. update lambda for each occupied block ----------------------
    csize <- tabulate(x, nbins = K)
    Zrs   <- rowSums(Z)
    for (k in which(csize > 0L)) {
      mem      <- which(x == k)
      shape_k  <- a + sum(w_i[mem])
      rate_k   <- b + sum(Zrs[mem])
      lambda[k] <- rgamma(1, shape_k, rate_k)
    }
    
    ## --- 4. identifiability – geometric‑mean = 1 ------------------
    occ <- which(csize > 0L)
    if (length(occ) > 0L) {
      g <- exp(mean(log(lambda[occ])))
      lambda[occ] <- lambda[occ] / g
    }
    lambda[csize == 0L] <- NA_real_
    
    ## --- 5. store -------------------------------------------------
    if (it > burnin) {
      keep_id                <- keep_id + 1L
      xsamp[keep_id, ]       <- x
      lsamp[keep_id, ]       <- lambda
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

