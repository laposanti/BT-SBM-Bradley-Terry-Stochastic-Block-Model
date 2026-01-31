## Beta (Github) version
# install.packages("devtools") # uncomment if you haven't installed 'devtools' before
devtools::install_github("pearce790/rankclust",force = TRUE)


library(rankclust)
wins_to_Pi <- function(W) {
  stopifnot(is.matrix(W), nrow(W) == ncol(W))
  idx <- which(W > 0, arr.ind = TRUE)
  n   <- sum(W[idx])
  Pi  <- matrix(NA_integer_, nrow = n, ncol = 2)
  pos <- 1L
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]; j <- idx[k, 2]
    m <- W[i, j]
    Pi[pos:(pos + m - 1L), ] <- cbind(rep.int(i, m), rep.int(j, m))
    pos <- pos + m
  }
  Pi
}

Pi <- wins_to_Pi(w_ij)

library(rankclust)

# (optional but recommended) fit plain BTL first to initialise nu0
res_btl <- mcmc_BTL(
  Pi = Pi, J = nrow(w_ij),
  a_gamma = 1, b_gamma = 1,
  num_iters = 3000, burn_prop = 0.2,
  chains = 1, groupwise = TRUE,
  seed = 1
)
# Note: The above code runs MCMC for the Ranked Choice BTL model using the 'rankclust' package.

nu0 <- colMeans(res_btl[, grep("^omega", names(res_btl))])

res_rcbtl <- mcmc_RCBTL(
  Pi = Pi, J = nrow(w_ij),
  a_gamma = 1, b_gamma = 1,
  lambda = 2, nu0 = nu0,
  num_iters = 3000, nu_reps = 2,
  burn_prop = 0.2, thin = 1,
  chains = 1, groupwise = TRUE,
  seed = 1
)

G <- as.matrix(res_rcbtl[, grep("^G", names(res_rcbtl))])
J <- ncol(G)

P_same <- matrix(0, J, J)
for (s in seq_len(nrow(G))) {
  g <- G[s, ]
  P_same <- P_same + outer(g, g, "==")
}
P_same <- P_same / nrow(G)


BTSBM::gibbs_bt_sbm(A = as.matrix(w_ij),
                      num_iters = 3000,
                      burn_in = 600,
                      alpha = 2,
                      beta = 2,
                      chains = 1,
                      groupwise = T)