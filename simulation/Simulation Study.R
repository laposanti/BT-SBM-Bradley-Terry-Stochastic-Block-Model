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


setwd(proj_dir)
source("./main_rccp.R")
# ------------------------------------------------
# 1) Load real data meta (for n_rep) — safe if file exists
# ------------------------------------------------
data_path_real <- './data/2000_2022_data.rds'
if (file.exists(data_path_real)) {
  data_real <- readRDS(data_path_real)
  n_rep <- length(data_real)
} else {
  # Fallback: choose a small number of replicates if the real file is missing
  n_rep <- 5L
}

# Ensure results dir exists
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
if (!dir.exists("data"))    dir.create("data", recursive = TRUE)

# ------------------------------------------------
# 2) Helpers
# ------------------------------------------------

# 2.1 Sparse-ish match topology
sample_Nij <- function(n_players, mean_matches = 10, p_edge = 0.05) {
  N <- matrix(0L, n_players, n_players)
  pairs <- utils::combn(n_players, 2)
  m     <- ncol(pairs)
  keep  <- which(runif(m) < p_edge)
  for (k in keep) {
    i <- pairs[1, k]; j <- pairs[2, k]
    n_ij <- rpois(1, mean_matches)
    if (n_ij > 0L) { N[i, j] <- N[j, i] <- n_ij }
  }
  diag(N) <- 0L
  N
}
# ---- Gnedin prior: sequential sampler for labels (unconditional) ----
# gamma \in (0,1)
sample_gnedin_labels <- function(n, gamma = 0.5) {
  stopifnot(n >= 1, is.finite(gamma), gamma > 0, gamma < 1)
  x <- integer(n)
  H <- 0L
  sizes <- integer(0)

  for (i in seq_len(n)) {
    n_curr <- sum(sizes)  # i-1
    if (H == 0L) { H <- 1L; sizes <- 1L; x[i] <- 1L; next }

    existing_weights <- (sizes + 1) * (n_curr - H + gamma)
    new_weight <- H*H - gamma*H
    if (new_weight < 0) new_weight <- 0

    w <- c(existing_weights, new_weight)
    if (!any(is.finite(w)) || sum(w) <= 0) w <- c(rep(1, H), 0)

    a <- sample.int(H + 1L, size = 1L, prob = w)
    if (a <= H) { sizes[a] <- sizes[a] + 1L; x[i] <- a }
    else        { H <- H + 1L; sizes <- c(sizes, 1L); x[i] <- H }
  }
  list(x = x, K = H, sizes = sizes)
}

# ---- Rejection sampler: draw until K(x) == K_target (exact conditioning) ----
sample_gnedin_labels_given_K <- function(n, K_target, gamma = 0.5,
                                         max_tries = 2000, verbose = FALSE) {
  best <- NULL
  for (t in seq_len(max_tries)) {
    out <- sample_gnedin_labels(n, gamma)
    if (out$K == K_target) return(out)
    # keep best (closest) attempt just to help debugging if it fails
    if (is.null(best) || abs(out$K - K_target) < abs(best$K - K_target)) best <- out
    if (verbose && t %% 100 == 0) message(sprintf("[Gnedin K|K*] try=%d got K=%d", t, out$K))
  }
  stop(sprintf("Failed to draw K(x) = %d after %d tries (closest was K=%d). ",
               K_target, max_tries, best$K),
       call. = FALSE)
}

lambda_to_theta <- function(lambda) {
  outer(lambda, lambda, function(a, b) a / (a + b))
}

# your existing match-topology helper is fine
sample_Nij <- function(n_players, mean_matches = 5, p_edge = 0.02) {
  N <- matrix(0L, n_players, n_players)
  pairs <- utils::combn(n_players, 2)
  m     <- ncol(pairs)
  keep  <- which(runif(m) < p_edge)
  for (k in keep) {
    i <- pairs[1, k]; j <- pairs[2, k]
    n_ij <- rpois(1, mean_matches)
    if (n_ij > 0L) { N[i, j] <- N[j, i] <- n_ij }
  }
  diag(N) <- 0L
  N
}

# --------- SIMULATE for K = 3..10 using Gnedin-conditioned labels ----------
for (K in 3:10) {
  K_data <- vector("list", n_rep)
  for (s in seq_len(n_rep)) {


    set.seed(123 + 1000*K + s)

    n    <- 150
    gamma_gn <- 0.1    # tweak if hitting K is too rare
    N_sim <- sample_Nij(n, mean_matches = 5, p_edge = 0.5)

    # ---- latent labels from Gnedin PRIOR conditioned on K(x)=K ----
    gn <- sample_gnedin_labels_given_K(n, K_target = K, gamma = gamma_gn,
                                       max_tries = 5000, verbose = FALSE)
    x_star <- rep(1:K,length.out=n)    # already 1..K in construction
    stopifnot(length(unique(x_star)) == K)

    # ---- block strengths and win-prob matrix ----
    make_lambda_geometric <- function(K, base=1, ratio=2.3) base * ratio^(0:(K-1))
    lambda_star <- make_lambda_geometric(K, base=0.08, ratio=2.3)

    theta_star  <- lambda_to_theta(lambda_star)

    # ---- generate outcome counts ----
    w_ij <- matrix(0L, n, n)
    idx  <- which(upper.tri(N_sim) & N_sim > 0L, arr.ind = TRUE)
    for (krow in seq_len(nrow(idx))) {
      i <- idx[krow, 1]; j <- idx[krow, 2]
      nij <- N_sim[i, j]
      pij <- theta_star[x_star[i], x_star[j]]
      wij <- rbinom(1, nij, pij)
      w_ij[i, j] <- wij
      w_ij[j, i] <- nij - wij
    }

    K_data[[s]] <- list(
      Y_ij        = w_ij,
      N_ij        = N_sim,
      lambda_star = lambda_star,
      theta_star  = theta_star,
      z_star      = x_star,
      K           = K,
      uid         = sprintf("s%03d-K%d", s, K)
    )
  }
  saveRDS(K_data, file = sprintf("./data/simulated_data_K%d.RDS", K))
}


# ------------------------------------------------
# 4) Parallel MCMC over (i, B)
# ------------------------------------------------


# Absolute project root (avoid relative path issues on workers)
proj_dir <- normalizePath('/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model', mustWork = TRUE)

# Helper: absolute path builders
data_file_K <- function(B) file.path(proj_dir, "data", sprintf("simulated_data_K%d.RDS", B))
out_file    <- file.path(proj_dir, "results", "augmented_simulation_results1.csv")
lockfile    <- paste0(out_file, ".lock")

# (Re)create results file
if (file.exists(out_file)) file.remove(out_file)

# Generate grid
all_combos <- expand.grid(B = 3:10, i = 1:n_rep)

# Define here so we can export it
generate_data <- function(i, B) {
  data_list <- readRDS(data_file_K(B))
  data_i    <- data_list[[i]]
  list(Y_ij = data_i$Y_ij, N_ij = data_i$N_ij, z_star = data_i$z_star)
}



# 3) Run foreach WITHOUT per-iteration init flags
results_df <- foreach(myCombo = iter(all_combos, by = "row"),
                      .combine  = rbind,
                      .packages = c("mcclust","mcclust.ext","mclust",
                                    "LaplacesDemon","fossil","coda",
                                    "readr","filelock",'Rcpp')) %do% {

                                      i <- myCombo$i
                                      B <- myCombo$B
                                      prior <- "GN"

                                      tryCatch({
                                        # --- Data ---
                                        sim_data <- generate_data(i = i, B = B)
                                        Y_ij <- sim_data$Y_ij
                                        N_ij <- sim_data$N_ij
                                        z_true <- sim_data$z_star

                                        # --- MCMC ---
                                        res <- gibbs_bt_urn_rcpp(
                                          w_ij     = Y_ij,
                                          a        = 2,
                                          b        = 1,
                                          prior    = "GN",
                                          gamma_GN = 0.8,
                                          burnin   = 5000,
                                          n_iter   = 20000,
                                          init_x   = NULL,
                                          store_z  = TRUE,
                                          verbose  = TRUE
                                        )

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
                                        vi_distance     <- mcclust::vi.dist(minVI_cl, z_true)
                                        ari_value       <- fossil::adj.rand.index(minVI_cl, z_true)
                                        binder_distance <- mcclust::vi.dist(binder_cl, z_true)
                                        binder_ari      <- fossil::adj.rand.index(binder_cl, z_true)

                                        # --- WAIC ---
                                        n_samples <- nrow(x_samples_relabel)
                                        edge_idx  <- which(upper.tri(N_ij) & N_ij > 0L, arr.ind = TRUE)
                                        n_edges   <- nrow(edge_idx)
                                        ll_matrix <- matrix(0, nrow = n_samples, ncol = n_edges)

                                        for (t in seq_len(n_samples)) {
                                          lam_t <- lambda_relabeled[t, ]
                                          for (e in seq_len(n_edges)) {
                                            i1 <- edge_idx[e, 1]; j1 <- edge_idx[e, 2]
                                            p_ij <- lam_t[i1] / (lam_t[i1] + lam_t[j1])
                                            y_ij <- Y_ij[i1, j1]; n_ij <- N_ij[i1, j1]
                                            ll_matrix[t, e] <- y_ij * log(p_ij) + (n_ij - y_ij) * log1p(-p_ij)
                                          }
                                        }
                                        waic_val   <- LaplacesDemon::WAIC(ll_matrix)$WAIC
                                        ess_lambda <- coda::effectiveSize(coda::mcmc(lambda_relabeled))
                                        mean_ess   <- mean(ess_lambda, na.rm = TRUE)

                                        row <- data.frame(
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

                                        # Thread-safe append (keep if you want immediate on-disk results)
                                        lk <- filelock::lock(lockfile)
                                        on.exit(filelock::unlock(lk), add = TRUE)
                                        write_header <- !file.exists(out_file) || file.size(out_file) == 0
                                        readr::write_csv(row, out_file, append = TRUE, col_names = write_header)

                                        row
                                      }, error = function(e) {
                                        row <- data.frame(
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
                                        lk <- filelock::lock(lockfile)
                                        on.exit(filelock::unlock(lk), add = TRUE)
                                        write_header <- !file.exists(out_file) || file.size(out_file) == 0
                                        readr::write_csv(row, out_file, append = TRUE, col_names = write_header)
                                        row
                                      })
                                    }

parallel::stopCluster(cl)

# Also write the combined object (in case you later want to re-aggregate)
write.csv(results_df, file = file.path(proj_dir, "results", "augmented_simulation_results1.csv"), row.names = FALSE)
print(table(results_df$B, results_df$VI_K))


### Plotting the results -- Adjusted Rand Index (ARI) based on MINVI partition point estimate
ARI_plot = results_df %>%
  ggplot(aes(x = factor(B), y = ARI_minVI)) +
  geom_boxplot(aes(fill = after_stat(..middle..))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_gradient(limits = c(0,1),
    low = "white",
    high = "#2E8B57",
    name = "Median ARI"
  ) +
  theme_minimal()+
  labs(y = "Median ARI", x= "K")

ggsave(filename = "./images/ARI_plot.png",ARI_plot , height = 8, width = 10)


# Contingency table of true vs estitmated K
tab_hard <- results_df %>%
  mutate(K_hat = round(mean_num_cl)) %>%
  count(B, K_hat, name = "n") %>%                   # contingency counts
  complete(B, K_hat, fill = list(n = 0)) %>%
  arrange(B, K_hat) %>%
  pivot_wider(names_from = K_hat, values_from = n, values_fill = 0)



latex_tab = kable(tab_hard, format = 'latex', booktabs = TRUE)
writeLines(latex_tab, "./results/table_K_sim_study.csv")

res_ARIVI = results_df%>%
  group_by(B)%>%
  summarise(mean_ARI = mean(ARI_minVI),
            mean_VI = mean(VI_minVI))
latex_tab_ARIVI = kable(res_ARIVI, format = 'latex', booktabs = TRUE)
writeLines(latex_tab_ARIVI, "./results/table_latex_tab_ARIVI_sim_study.csv")



# ------------------------------------------------
# 5) Diagnostics (4 chains, relabel, Gelman–Rubin + Geweke)
# ------------------------------------------------

diagnostics_lambda <- function(lambda_chain_list) {
  stopifnot(length(lambda_chain_list) >= 2)
  m_min <- min(vapply(lambda_chain_list, nrow, integer(1)))
  to_mcmc <- function(mat) coda::mcmc(mat[seq_len(m_min), , drop = FALSE], start = 1, thin = 1)
  ml <- coda::mcmc.list(lapply(lambda_chain_list, to_mcmc))

  # Limit to first 50 dimensions for stability (as in your earlier code)
  p <- min(ncol(lambda_chain_list[[1]]), 50L)
  ml_sub <- coda::mcmc.list(lapply(ml, function(ch) ch[, seq_len(p)]))

  gel <- coda::gelman.diag(ml_sub, autoburnin = FALSE, multivariate = FALSE)
  ess_vec <- coda::effectiveSize(ml_sub)

  # Geweke: compute per chain, then average |z| across params and chains
  geweke_per_chain <- lapply(ml_sub, function(ch) {
    gz <- coda::geweke.diag(ch, frac1 = 0.1, frac2 = 0.5)
    # gz$z is a named vector of z-scores for each parameter
    mean(abs(gz$z), na.rm = TRUE)
  })

  data.frame(
    ESS_mean = mean(ess_vec, na.rm = TRUE),
    Rhat_mean = mean(gel$psrf[, "Point est."], na.rm = TRUE),
    Geweke_mean_absz = mean(unlist(geweke_per_chain), na.rm = TRUE)
  )
}

diag_out_file <- file.path("results", "augmented_simulation_diagnostics.csv")
if (file.exists(diag_out_file)) file.remove(diag_out_file)

# One dataset per K (to mirror your earlier intent); 4 chains, relabel, GR + Geweke
results_diag <- foreach(B = 3:10,
                        .combine = rbind,
                        .packages = c("mcclust","mcclust.ext","mclust",
                                      "LaplacesDemon","fossil","coda",
                                      "abind","filelock","readr")) %do% {
                                        i <- 1
                                        dat <- generate_data(i, B)
                                        Y_ij <- dat$Y_ij;

                                        n_chains <- 4L
                                        chain_seeds <- 4321 + seq_len(n_chains) + 1000L * i + 10L * B

                                        chain_out = list()
                                        for(chain in 1:n_chains){
                                          res <- gibbs_bt_urn_rcpp(
                                            w_ij     = Y_ij,
                                            a        = 4,
                                            b        = 1,
                                            prior    = "GN",
                                            gamma_GN = 0.8,
                                            burnin   = 5000,
                                            n_iter   = 20000,
                                            init_x   = NULL,
                                            store_z  = TRUE,
                                            verbose  = TRUE
                                          )
                                          chain_out[[chain]]<- res
                                        }

                                        # ---------- DIAGNOSTICS: use cluster-level λ, NA-free ----------
                                        # (inside your foreach over B)
                                        x_list   <- lapply(chain_out, `[[`, "x_samples")
                                        lam_list <- lapply(chain_out, `[[`, "lambda_samples")  # <- LIST per chain

                                        # rows per chain and pooling
                                        n_by_chain <- vapply(x_list, nrow, integer(1))
                                        x_pool     <- do.call(rbind, x_list)
                                        # lambda is a list-of-lists: flatten to a single list of per-iter lambda vectors
                                        lam_pool   <- do.call(c, lam_list)

                                        relab        <- relabel_by_lambda(x_pool, lam_pool)
                                        x_rel_pool   <- relab$x_samples_relabel
                                        lam_cl_pool  <- relab$lambda_cluster  # list length S_total; each numeric length H_iter

                                        # split back by chain
                                        idx_split <- split(seq_len(nrow(x_rel_pool)),
                                                           rep(seq_along(n_by_chain), times = n_by_chain))
                                        x_rel_list   <- lapply(idx_split, function(idx) x_rel_pool[idx, , drop = FALSE])
                                        lam_cl_list  <- lapply(idx_split, function(idx) lam_cl_pool[idx])

                                        # K_hat from pooled relabeled x
                                        K_hat <- length(unique(mcclust.ext::minVI(mcclust::comp.psm(x_rel_pool))$cl))

                                        # keep only iterations with exactly K_hat clusters
                                        keep_idx_list <- lapply(x_rel_list, function(xrel) {
                                          which(apply(xrel, 1, function(r) length(unique(r)) == K_hat))
                                        })

                                        # build per-chain matrices (iters x K_hat) from ordered cluster λ vectors
                                        lam_chain_cl <- mapply(function(keep_idx, vec_list) {
                                          if (length(keep_idx) == 0) return(NULL)
                                          rows <- lapply(keep_idx, function(j) {
                                            v <- vec_list[[j]]
                                            if (!is.numeric(v)) return(NULL)
                                            if (length(v) != K_hat) return(NULL)   # skip if not exactly K_hat
                                            if (!all(is.finite(v))) return(NULL)
                                            v
                                          })
                                          rows <- Filter(Negate(is.null), rows)
                                          if (length(rows) < 2L) return(NULL)
                                          M <- do.call(rbind, rows)
                                          storage.mode(M) <- "double"
                                          M
                                        }, keep_idx_list, lam_cl_list, SIMPLIFY = FALSE)

                                        lam_chain_cl <- Filter(Negate(is.null), lam_chain_cl)

                                        if (length(lam_chain_cl) >= 2) {
                                          diag_tbl   <- diagnostics_lambda(lam_chain_cl)  # your helper: list -> mcmc.list -> GR/ESS/Geweke
                                          ess_mean   <- diag_tbl$ESS_mean
                                          rhat_mean  <- diag_tbl$Rhat_mean
                                          gz_meanabs <- diag_tbl$Geweke_mean_absz
                                        } else {
                                          ess_mean <- rhat_mean <- gz_meanabs <- NA_real_
                                          warning(sprintf("K=%d: fewer than 2 usable chains for diagnostics after filtering.", B))
                                        }

                                        out_row <- data.frame(
                                          dataset_id = i,
                                          B          = B,
                                          ESS_mean   = ess_mean,
                                          Rhat_mean  = rhat_mean,
                                          Geweke_mean_absz = gz_meanabs
                                        )

                                        lk <- filelock::lock(paste0(diag_out_file, ".lock"))
                                        on.exit(filelock::unlock(lk), add = TRUE)
                                        write_header <- !file.exists(diag_out_file) || file.size(diag_out_file) == 0
                                        readr::write_csv(out_row, diag_out_file, append = TRUE, col_names = write_header)

                                        out_row
                                      }

# ------------------------------------------------
# 6) Timing vs dataset size (n = 100, 500, 1000, 1000)
# ------------------------------------------------

# helper to simulate a BT dataset of size n (independent of earlier K-sim pipeline)
simulate_bt_data <- function(n, mean_matches = 5, p_edge = 0.2, K) {
  # reuse your sample_Nij helper already defined above
  N_ij <- sample_Nij(n_players = n, mean_matches = mean_matches, p_edge = p_edge)
  # player-specific lambda, then generate outcomes
  x_star <- rep(1:K,length.out=n)    # already 1..K in construction
  stopifnot(length(unique(x_star)) == K)

  # ---- block strengths and win-prob matrix ----
  make_lambda_geometric <- function(K, base=1, ratio=2.3) base * ratio^(0:(K-1))
  lambda_star <- make_lambda_geometric(K, base=0.08, ratio=2.3)

  theta_star  <- lambda_to_theta(lambda_star)
  # ---- generate outcome counts ----
  w_ij <- matrix(0L, n, n)
  idx  <- which(upper.tri(N_ij) & N_ij > 0L, arr.ind = TRUE)
  for (krow in seq_len(nrow(idx))) {
    i <- idx[krow, 1]; j <- idx[krow, 2]
    nij <- N_ij[i, j]
    pij <- theta_star[x_star[i], x_star[j]]
    wij <- rbinom(1, nij, pij)
    w_ij[i, j] <- wij
    w_ij[j, i] <- nij - wij
  }
  return(w_ij)
}

sizes_n <- c(100L, 500L, 1000L, 5000L)
p_edge <- c(0.5,0.10,0.05,0.01)
time_df <- data.frame(n = sizes_n, p_edge=p_edge, E = integer(length(sizes_n)), time_sec = NA_real_)

for (ii in seq_along(sizes_n)) {
  n <- sizes_n[ii]
  w_ij <- simulate_bt_data(n, mean_matches = 5, p_edge = p_edge[ii],K = 5)
  time_df$E[ii] = sum(w_ij)
  tm <- system.time({
    invisible(
      gibbs_bt_urn_rcpp(
        w_ij     = w_ij,
        a        = 4, b = 1,
        prior    = "GN", gamma_GN = 0.8,
        burnin   = 1000,
        n_iter   = 10000,
        init_x   = NULL,
        store_z  = FALSE,
        verbose  = T
      )
    )
  })
  time_df$time_sec[ii] <- unname(tm[["elapsed"]])
  print(round(time_df$time_sec/60, 2))
}

time_df <- time_df %>% mutate(time_min = round(time_sec/60, 2))
print(kable(time_df, align = "r",format = 'latex',booktabs = T))


