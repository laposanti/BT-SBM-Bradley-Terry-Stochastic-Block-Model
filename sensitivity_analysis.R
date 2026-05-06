# ============================================================
# BT–SBM sensitivity over a in {1,2,3,4} on ONE dataset
#   - one simulated dataset (you can swap in your own W, z_true)
#   - fit ONLY BT–SBM (BTSBM::gibbs_bt_sbm), GN prior
#   - compute VI to truth (per-draw + minVI partition)
#   - store posterior K draws
#   - boxplots for VI and K by a
# ============================================================

suppressPackageStartupMessages({
  library(BTSBM)        # provides gibbs_bt_sbm
  library(mcclust)      # comp.psm, vi.dist
  library(mcclust.ext)  # minVI (and minbinder.ext if you want)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

# ----------------------------
# Helpers
# ----------------------------
compress_labels <- function(z) {
  z <- as.integer(z)
  u <- sort(unique(z))
  match(z, u)
}

get_minVI_partition <- function(psm, cls_draw) {
  # mcclust.ext::minVI has slightly different return formats across versions
  obj <- mcclust.ext::minVI(psm, cls.draw = cls_draw, method = "all")
  if (!is.null(obj$cl)) {
    return(as.integer(obj$cl[1, ]))
  }
  if (!is.null(obj$cstar)) {
    return(as.integer(obj$cstar))
  }
  if (!is.null(obj$c.star)) {
    return(as.integer(obj$c.star))
  }
  stop("Could not extract minVI partition from mcclust.ext::minVI output.")
}

mode_int <- function(x) {
  tx <- table(x)
  as.integer(names(tx)[which.max(tx)])
}

# ----------------------------
# One-dataset DGP (simple, transparent)
#   Uses a match design N_ij (symmetric), simulates W_ij
# ----------------------------
simulate_one_dataset <- function(N, K_true = 5, p_adj = 0.85, seed = 1) {
  set.seed(seed)
  stopifnot(is.matrix(N), nrow(N) == ncol(N))
  n <- nrow(N)
  
  # Ground-truth partition: roughly equal sizes, permuted
  z_true <- rep(seq_len(K_true), length.out = n)
  z_true <- sample(z_true, size = n, replace = FALSE)
  z_true <- as.integer(z_true)
  
  # Block strengths (tiered): p_adj controls adjacent win prob in expectation
  # lambda_k proportional to exp(delta*(k-1)), delta = logit(p_adj)
  delta <- qlogis(p_adj)
  lambda_by_block <- exp(delta * (0:(K_true - 1L)))
  lambda_by_block <- rev(lambda_by_block) # block 1 strongest, purely for readability
  
  W <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      nij <- N[i, j]
      if (nij > 0) {
        li <- lambda_by_block[z_true[i]]
        lj <- lambda_by_block[z_true[j]]
        p_ij <- li / (li + lj)
        w_ij <- rbinom(1L, size = nij, prob = p_ij)
        W[i, j] <- w_ij
        W[j, i] <- nij - w_ij
      }
    }
  }
  diag(W) <- 0L
  
  list(N = N, W = W, z_true = z_true, K_true = K_true,
       lambda_by_block = lambda_by_block)
}

# ----------------------------
# Main sensitivity runner
# ----------------------------
run_sensitivity_a <- function(
    N,
    a_grid = 1:4,
    gamma_fit = 0.71,
    T_iter = 10000,
    T_burn = 2000,
    thin = 1,
    K_true = 5,
    p_adj = 0.85,
    seed_data = 1,
    seed_mcmc = 100,
    out_dir = "results/sensitivity_a"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1) Make ONE dataset (swap this out if you want to use real data)
  dat <- simulate_one_dataset(N, K_true = K_true, p_adj = p_adj, seed = seed_data)
  W <- dat$W
  z_true <- dat$z_true
  
  # Save dataset for reproducibility
  saveRDS(dat, file.path(out_dir, "dataset_used.rds"))
  
  all_draws <- list()
  summaries <- list()
  k_posteriors <- list()
  
  for (a_val in a_grid) {
    cat("Fitting BT–SBM with a =", a_val, "...\n")
    set.seed(seed_mcmc + as.integer(a_val))
    
    t0 <- proc.time()[["elapsed"]]
    fit <- BTSBM::gibbs_bt_sbm(
      w_ij   = W,
      T_iter = T_iter,
      T_burn = T_burn,
      a      = a_val,
      gamma  = gamma_fit,
      prior  = "GN",
      verbose = FALSE
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    
    # x_samples: iterations x n_items
    x_raw <- fit$x_samples
    stopifnot(is.matrix(x_raw), ncol(x_raw) == length(z_true))
    
    # Optionally thin (after burn already applied inside BTSBM sampler)
    keep_idx <- seq(1, nrow(x_raw), by = thin)
    x_raw <- x_raw[keep_idx, , drop = FALSE]
    
    # Compress labels per draw (removes gaps like {1,3,7})
    x_draw <- t(apply(x_raw, 1, compress_labels))
    
    # Posterior K draws
    K_draw <- apply(x_draw, 1, function(z) length(unique(z)))
    
    # VI per draw against truth
    VI_draw <- apply(x_draw, 1, function(z) mcclust::vi.dist(as.integer(z), z_true))
    
    # minVI posterior point estimate
    psm <- mcclust::comp.psm(x_draw)
    z_minVI <- get_minVI_partition(psm, cls_draw = x_draw)
    vi_minVI <- mcclust::vi.dist(z_minVI, z_true)
    
    # Summaries you’ll likely want in a table
    sum_row <- tibble(
      a = a_val,
      gamma = gamma_fit,
      T_iter = T_iter,
      T_burn = T_burn,
      thin = thin,
      elapsed_sec = elapsed,
      K_true = dat$K_true,
      K_mode = mode_int(K_draw),
      K_median = as.integer(round(median(K_draw))),
      pr_K_true = mean(K_draw == dat$K_true),
      vi_minVI = vi_minVI,
      vi_mean = mean(VI_draw),
      vi_median = median(VI_draw)
    )
    
    # Store draw-level output for boxplots
    draws_df <- tibble(
      a = a_val,
      iter = seq_along(K_draw),
      K = K_draw,
      VI = VI_draw
    )
    
    # Store posterior K distribution (as a tidy table)
    k_tab <- as.data.frame(table(K_draw), stringsAsFactors = FALSE) |>
      transmute(a = a_val,
                K = as.integer(K_draw),
                count = Freq,
                prob = Freq / sum(Freq))
    
    all_draws[[as.character(a_val)]] <- draws_df
    summaries[[as.character(a_val)]] <- sum_row
    k_posteriors[[as.character(a_val)]] <- k_tab
    
    # Save per-a objects too (handy when something crashes later)
    saveRDS(list(fit = fit, x_draw = x_draw, draws = draws_df, summary = sum_row, K_tab = k_tab),
            file.path(out_dir, paste0("fit_a", a_val, ".rds")))
  }
  
  draws_all <- bind_rows(all_draws)
  summary_all <- bind_rows(summaries)
  k_post_all <- bind_rows(k_posteriors)
  
  # Write tables
  write_csv(summary_all, file.path(out_dir, "summary_by_a.csv"))
  write_csv(draws_all, file.path(out_dir, "draws_VI_K_by_a.csv"))
  write_csv(k_post_all, file.path(out_dir, "posterior_K_distribution_by_a.csv"))
  
  # ----------------------------
  # Plots
  # ----------------------------
  p_vi <- ggplot(draws_all, aes(x = factor(a), y = VI)) +
    geom_boxplot(outlier.alpha = 0.25) +
    theme_bw(base_size = 12) +
    labs(x = "a", y = "VI distance to truth",
         title = "BT–SBM sensitivity: VI per MCMC draw")
  
  p_k <- ggplot(draws_all, aes(x = factor(a), y = K)) +
    geom_boxplot(outlier.alpha = 0.25) +
    theme_bw(base_size = 12) +
    labs(x = "a", y = "Posterior K draws",
         title = "BT–SBM sensitivity: posterior K per MCMC draw")
  
  p_kbar <- ggplot(k_post_all, aes(x = factor(K), y = prob)) +
    geom_col() +
    facet_wrap(~ a, nrow = 1) +
    theme_bw(base_size = 12) +
    labs(x = "K", y = "Posterior probability",
         title = "Posterior distribution of K by a")
  
  ggsave(file.path(out_dir, "boxplot_VI_by_a.png"), p_vi, width = 7.5, height = 4.5, dpi = 300)
  ggsave(file.path(out_dir, "boxplot_K_by_a.png"),  p_k,  width = 7.5, height = 4.5, dpi = 300)
  ggsave(file.path(out_dir, "bar_posteriorK_by_a.png"), p_kbar, width = 10.5, height = 3.8, dpi = 300)
  
  invisible(list(summary = summary_all, draws = draws_all, posteriorK = k_post_all,
                 plots = list(VI = p_vi, K = p_k, posteriorK = p_kbar)))
}

# ============================================================
# Example run
#   - uses your N_list file if available, otherwise stop with a clear message
# ============================================================
if (file.exists("./N_list_for_simulation.rds")) {
  N_list <- readRDS("./N_list_for_simulation.rds")
  design_name <- names(N_list)[1]
  N <- N_list[[design_name]]
  
  res <- run_sensitivity_a(
    N = N,
    a_grid = 1:4,
    gamma_fit = 0.8,
    T_iter = 10000,
    T_burn = 2000,
    thin = 2,
    K_true = 5,
    p_adj = 0.85,
    seed_data = 123,
    seed_mcmc = 999,
    out_dir = file.path("results", "sensitivity_a")
  )
  
  print(res$summary)
} else {
  stop("Cannot find ./N_list_for_simulation.rds. Provide N (match design) or adapt the data-loading section.")
}


res$summary%>%
  select(a, K_true, K_mode, pr_K_true, vi_mean) %>%
  kable(booktabs = TRUE, digits = 3,
        format = 'latex') 



