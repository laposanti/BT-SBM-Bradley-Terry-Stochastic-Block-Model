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
if (!requireNamespace("fossil", quietly = TRUE)) {
  stop("Package 'fossil' is required for ARI. Install it via install.packages('fossil').")
}
library(dplyr)
library(tidyr)


## -------------------------------
## Experiment configuration
## -------------------------------

seed0 <- as.integer(Sys.getenv("SEED", "123"))      # base seed
gamma_true <- as.numeric(Sys.getenv("GAMMA_TRUE", "0.71"))
p_adj <- as.numeric(Sys.getenv("P_ADJ", "0.85"))

## User request: default to 10,000 iterations.
T_iter <- as.integer(Sys.getenv("T_ITER", "15000"))
T_burn <- as.integer(Sys.getenv("T_BURN", "5000"))
if (!is.finite(T_iter) || T_iter < 10L) stop("T_ITER must be a positive integer")
if (!is.finite(T_burn) || T_burn < 0L || T_burn >= T_iter) stop("T_BURN must be in [0, T_ITER)")

## RCBTL prior hyperparameters (rankclust)
a_gamma_rcbtl <- as.numeric(Sys.getenv("A_GAMMA_RCBTL", "5"))
b_gamma_rcbtl <- as.numeric(Sys.getenv("B_GAMMA_RCBTL", "3"))
if (!is.finite(a_gamma_rcbtl) || a_gamma_rcbtl <= 0) stop("A_GAMMA_RCBTL must be > 0")
if (!is.finite(b_gamma_rcbtl) || b_gamma_rcbtl <= 0) stop("B_GAMMA_RCBTL must be > 0")

## Optional: save raw fits (can be large)
save_raw_fits <- isTRUE(as.logical(Sys.getenv("SAVE_RAW_FITS", "TRUE")))
raw_fits_dir <- NA_character_

## User request: do a 5-run test first.
n_runs <- as.integer(Sys.getenv("N_RUNS", "1"))
if (!is.finite(n_runs) || n_runs < 1L) stop("N_RUNS must be a positive integer")

out_dir <- Sys.getenv("OUT_DIR", "./results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (save_raw_fits) {
  raw_fits_dir <- file.path(out_dir, "raw_fits")
  dir.create(raw_fits_dir, recursive = TRUE, showWarnings = FALSE)
}

N_list <- readRDS("./N_list_for_simulation.rds")
if (!is.list(N_list) || length(N_list) < 1L) stop("N_list_for_simulation.rds must contain a non-empty list of N designs")

design_name <- Sys.getenv("SIM_DESIGN", names(N_list)[1])
if (!design_name %in% names(N_list)) design_name <- names(N_list)[1]

K_true <- as.integer(Sys.getenv("SIM_K_TRUE", "3"))
K_true <- max(2L, K_true)

N <- as.matrix(N_list[[design_name]])
if (!is.matrix(N) || nrow(N) != ncol(N)) stop("Selected design is not a square matrix: ", design_name)

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

gelman_summary_df <- function(mcmc_list, run_id, seed_run, model, quantity) {
  if (length(mcmc_list) < 2L) {
    return(tibble::tibble(
      run = run_id,
      seed = seed_run,
      model = model,
      quantity = quantity,
      n_chains = length(mcmc_list),
      n_params = NA_integer_,
      mpsrf = NA_real_,
      max_psrf_point = NA_real_,
      max_psrf_upper = NA_real_,
      median_psrf_upper = NA_real_
    ))
  }
  
  ## Align chains to a common iteration length (min over chains)
  mats <- lapply(mcmc_list, function(x) {
    M <- as.matrix(x)
    storage.mode(M) <- "double"
    M
  })
  n_iters <- vapply(mats, nrow, integer(1))
  min_iters <- min(n_iters)
  mats <- lapply(mats, function(M) M[seq_len(min_iters), , drop = FALSE])
  mcmc_list_aligned <- lapply(mats, coda::as.mcmc)
  
  ml <- coda::mcmc.list(mcmc_list_aligned)
  gd <- coda::gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)
  psrf <- gd$psrf
  point <- psrf[, 1]
  upper <- psrf[, 2]
  
  mpsrf_val <- NA_real_
  if (nrow(psrf) > 1L) {
    mpsrf_val <- tryCatch(
      as.numeric(coda::gelman.diag(ml, autoburnin = FALSE, multivariate = TRUE)$mpsrf),
      error = function(e) NA_real_
    )
  }
  
  tibble::tibble(
    run = run_id,
    seed = seed_run,
    model = model,
    quantity = quantity,
    n_chains = length(mcmc_list),
    n_params = nrow(psrf),
    mpsrf = mpsrf_val,
    max_psrf_point = suppressWarnings(max(point, na.rm = TRUE)),
    max_psrf_upper = suppressWarnings(max(upper, na.rm = TRUE)),
    median_psrf_upper = suppressWarnings(stats::median(upper, na.rm = TRUE))
  )
}

## -------------------------------
## Run one simulation + fit + summarise
## -------------------------------

run_one <- function(run_id, seed_run) {
  sim <- simulate_dataset_design_based(
    N,
    K_target = K_true,
    gamma_true = gamma_true,
    p_adj = p_adj,
    seed = seed_run
  )
  
  w_ij <- sim$W
  x_true <- sim$z_true
  lambda_true <- sim$lambda_true_players
  n <- nrow(w_ij)
  Pi <- wins_to_Pi(w_ij)
  
  cat("\n==============================\n")
  cat("Run ", run_id, "/", n_runs, "  seed=", seed_run, "  design=", design_name, "  n=", n, "  K_true=", K_true, "\n", sep = "")
  
  ## 4.1 Plain BT (BTSBM)
  t_bt <- system.time({
    fit_bt_list <- run_independent_chains(
      n_chains = n_chains,
      seeds = chain_seeds + 100000L * run_id,
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
  
  ## 4.2 Rank-clustered (rankclust)
  t_rcbtl <- system.time({
    res_rcbtl_list <- run_independent_chains(
      n_chains = n_chains,
      seeds = chain_seeds + 100000L * run_id + 10000L,
      n_cores = n_cores,
      use_parallel = use_parallel,
      FUN = function(seed, chain_id) {
        set.seed(seed)
        res_btl <- mcmc_BTL(
          Pi = Pi,
          J = n,
          a_gamma = a_gamma_rcbtl,
          b_gamma = b_gamma_rcbtl,
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
          a_gamma = a_gamma_rcbtl,
          b_gamma = b_gamma_rcbtl,
          lambda = 2,
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
        
        ## rankclust sometimes returns a matrix; `cluster_ggplots()` assumes
        ## a data.frame-like object with `chain` and `iteration` columns.
        res_df <- as.data.frame(res)
        res_df$chain <- factor(chain_id)
        res_df
      }
    )
  })
  
  ## 4.3 BT-SBM (BTSBM)
  t_btsbm <- system.time({
    fit_btsbm_list <- run_independent_chains(
      n_chains = n_chains,
      seeds = chain_seeds + 100000L * run_id + 20000L,
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
  
  if (save_raw_fits) {
    stamp <- paste0(
      design_name,
      "_K", K_true,
      "_iters", T_iter,
      "_burn", T_burn,
      "_run", run_id,
      "_seed", seed_run
    )
    saveRDS(
      list(
        model = "RCBTL",
        design = design_name,
        K_true = K_true,
        run = run_id,
        seed = seed_run,
        T_iter = T_iter,
        T_burn = T_burn,
        n_chains = n_chains,
        a_gamma = a_gamma_rcbtl,
        b_gamma = b_gamma_rcbtl,
        fits = res_rcbtl_list
      ),
      file = file.path(raw_fits_dir, paste0("rawfit_RCBTL_", stamp, ".rds"))
    )
    saveRDS(
      list(
        model = "BT-SBM",
        design = design_name,
        K_true = K_true,
        run = run_id,
        seed = seed_run,
        T_iter = T_iter,
        T_burn = T_burn,
        n_chains = n_chains,
        fits = fit_btsbm_list
      ),
      file = file.path(raw_fits_dir, paste0("rawfit_BT-SBM_", stamp, ".rds"))
    )
  }
  
  ## -------------------------------
  ## Summaries: lambda on common scale
  ## -------------------------------
  
  bt_draws_all <- do.call(rbind, lapply(fit_bt_list, function(f) f$lambda_samples))
  bt_lambda_hat <- norm_geo1(colMeans(bt_draws_all))
  
  res_rcbtl <- dplyr::bind_rows(res_rcbtl_list)
  rcbtl_omega_cols <- grep("^omega", names(res_rcbtl), value = TRUE)
  rcbtl_lambda_hat <- norm_geo1(colMeans(as.matrix(res_rcbtl[, rcbtl_omega_cols, drop = FALSE])))
  
  btsbm_chain_summaries <- lapply(fit_btsbm_list, function(fit) {
    relabel_by_lambda(fit$x_samples, fit$lambda_samples)
  })
  lambda_relabeled_all <- do.call(rbind, lapply(btsbm_chain_summaries, function(inf) inf$lambda_samples_relabel))
  btsbm_lambda_hat <- norm_geo1(colMeans(lambda_relabeled_all))
  
  lambda_true_scaled <- norm_geo1(lambda_true)
  
  bt_lambda_metrics <- compare_lambda(bt_lambda_hat, lambda_true_scaled)
  rcbtl_lambda_metrics <- compare_lambda(rcbtl_lambda_hat, lambda_true_scaled)
  btsbm_lambda_metrics <- compare_lambda(btsbm_lambda_hat, lambda_true_scaled)
  
  ## -------------------------------
  ## Gelman-Rubin diagnostics across chains
  ## (summarized to keep CSV small)
  ## -------------------------------
  
  bt_mcmc_list <- lapply(fit_bt_list, function(f) {
    M <- log(norm_geo1_rows(f$lambda_samples))
    as_mcmc_safe(M)
  })
  
  rcbtl_mcmc_omega <- lapply(res_rcbtl_list, function(df) {
    omega_cols <- grep("^omega", names(df), value = TRUE)
    M <- as.matrix(df[, omega_cols, drop = FALSE])
    M <- log(norm_geo1_rows(M))
    as_mcmc_safe(M)
  })
  
  rcbtl_mcmc_K <- lapply(res_rcbtl_list, function(df) {
    as_mcmc_safe(matrix(as.numeric(df$K), ncol = 1, dimnames = list(NULL, "K")))
  })
  
  btsbm_mcmc_lambda <- lapply(btsbm_chain_summaries, function(inf) {
    M <- log(norm_geo1_rows(inf$lambda_samples_relabel))
    as_mcmc_safe(M)
  })
  
  btsbm_mcmc_K <- lapply(btsbm_chain_summaries, function(inf) {
    as_mcmc_safe(matrix(as.numeric(inf$n_clusters_each_iter), ncol = 1, dimnames = list(NULL, "K")))
  })
  
  gelman_out <- dplyr::bind_rows(
    gelman_summary_df(bt_mcmc_list, run_id, seed_run, model = "BT", quantity = "log_lambda"),
    gelman_summary_df(rcbtl_mcmc_omega, run_id, seed_run, model = "RCBTL", quantity = "log_omega"),
    gelman_summary_df(rcbtl_mcmc_K, run_id, seed_run, model = "RCBTL", quantity = "K"),
    gelman_summary_df(btsbm_mcmc_lambda, run_id, seed_run, model = "BT-SBM", quantity = "log_lambda"),
    gelman_summary_df(btsbm_mcmc_K, run_id, seed_run, model = "BT-SBM", quantity = "K")
  )
  
  ## -------------------------------
  ## Summaries: partitions via PSM -> minVI, then ARI/VI vs truth
  ## -------------------------------
  
  ## RCBTL partition draws
  rcbtl_G_cols <- grep("^G[0-9]+$", names(res_rcbtl), value = TRUE)
  if (length(rcbtl_G_cols) < 1L) {
    stop("RCBTL output has no G1..GJ columns; cannot compute partitions.")
  }
  rcbtl_cls_draw <- as.matrix(res_rcbtl[, rcbtl_G_cols, drop = FALSE])
  storage.mode(rcbtl_cls_draw) <- "integer"
  rcbtl_hat <- relabel_consecutive(mcclust.ext::minVI(rcbtl_cls_draw))
  rcbtl_part_metrics <- partition_metrics(rcbtl_hat, x_true)
  rcbtl_K_draws <- if ("K" %in% names(res_rcbtl)) as.numeric(res_rcbtl$K) else apply(rcbtl_cls_draw, 1, function(z) length(unique(z)))
  rcbtl_K_sum <- summarise_K_draws(rcbtl_K_draws)
  
  ## BTSBM partition draws
  btsbm_cls_draw <- do.call(rbind, lapply(fit_btsbm_list, function(f) f$x_samples))
  storage.mode(btsbm_cls_draw) <- "integer"
  btsbm_hat <- relabel_consecutive(mcclust.ext::minVI(btsbm_cls_draw))
  btsbm_part_metrics <- partition_metrics(btsbm_hat, x_true)
  btsbm_K_draws <- apply(btsbm_cls_draw, 1, function(z) length(unique(z)))
  btsbm_K_sum <- summarise_K_draws(btsbm_K_draws)
  
  ## -------------------------------
  ## Output row(s) for CSV
  ## -------------------------------
  
  timing_df <- tibble::tibble(
    model = c("BT", "RCBTL", "BT-SBM"),
    elapsed = c(t_bt["elapsed"], t_rcbtl["elapsed"], t_btsbm["elapsed"])
  )
  
  base_cols <- tibble::tibble(
    run = run_id,
    seed = seed_run,
    design = design_name,
    n = n,
    K_true = K_true,
    T_iter = T_iter,
    T_burn = T_burn,
    n_chains = n_chains
  )
  
  out <- dplyr::bind_rows(
    dplyr::bind_cols(
      base_cols,
      tibble::tibble(model = "BT"),
      bt_lambda_metrics,
      tibble::tibble(ari = NA_real_, vi = NA_real_, K_mean = NA_real_, K_ci_low = NA_real_, K_ci_high = NA_real_),
      tibble::tibble(elapsed = timing_df$elapsed[timing_df$model == "BT"])
    ),
    dplyr::bind_cols(
      base_cols,
      tibble::tibble(model = "RCBTL"),
      rcbtl_lambda_metrics,
      rcbtl_part_metrics,
      rcbtl_K_sum,
      tibble::tibble(elapsed = timing_df$elapsed[timing_df$model == "RCBTL"])
    ),
    dplyr::bind_cols(
      base_cols,
      tibble::tibble(model = "BT-SBM"),
      btsbm_lambda_metrics,
      btsbm_part_metrics,
      btsbm_K_sum,
      tibble::tibble(elapsed = timing_df$elapsed[timing_df$model == "BT-SBM"])
    )
  )
  
  partitions_out <- dplyr::bind_rows(
    partition_to_long_df(
      partition = x_true,
      run_id = run_id,
      seed_run = seed_run,
      model = "truth",
      design_name = design_name,
      K_true = K_true,
      T_iter = T_iter,
      T_burn = T_burn
    ),
    partition_to_long_df(
      partition = rcbtl_hat,
      run_id = run_id,
      seed_run = seed_run,
      model = "RCBTL",
      design_name = design_name,
      K_true = K_true,
      T_iter = T_iter,
      T_burn = T_burn
    ),
    partition_to_long_df(
      partition = btsbm_hat,
      run_id = run_id,
      seed_run = seed_run,
      model = "BT-SBM",
      design_name = design_name,
      K_true = K_true,
      T_iter = T_iter,
      T_burn = T_burn
    )
  )
  
  list(summary = out, partitions = partitions_out, gelman = gelman_out)
}

## -------------------------------
## Main: run n_runs experiments and write CSV
## -------------------------------

all_runs <- lapply(seq_len(n_runs), function(r) {
  run_one(run_id = r, seed_run = seed0 + 1000L * (r - 1L))
})

all_summaries <- dplyr::bind_rows(lapply(all_runs, `[[`, "summary"))
all_partitions <- dplyr::bind_rows(lapply(all_runs, `[[`, "partitions"))
all_gelman <- dplyr::bind_rows(lapply(all_runs, `[[`, "gelman"))

out_csv <- file.path(
  out_dir,
  paste0(
    "comparison_pearce_ereshova_",
    design_name,
    "_K",
    K_true,
    "_iters",
    T_iter,
    "_runs",
    n_runs,
    ".csv"
  )
)

utils::write.csv(all_summaries, out_csv, row.names = FALSE)
cat("\nWrote summary CSV to: ", out_csv, "\n", sep = "")

out_csv_partitions <- sub("\\.csv$", "_partitions.csv", out_csv)
utils::write.csv(all_partitions, out_csv_partitions, row.names = FALSE)
cat("Wrote partitions CSV to: ", out_csv_partitions, "\n", sep = "")

out_csv_gelman <- sub("\\.csv$", "_gelman.csv", out_csv)
utils::write.csv(all_gelman, out_csv_gelman, row.names = FALSE)
cat("Wrote Gelman-Rubin CSV to: ", out_csv_gelman, "\n", sep = "")

## Keep a compact printout in the console
print(dplyr::arrange(all_summaries, run, model))
