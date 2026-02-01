# ================================================================
# BT–SBM Simulation + Multi-Model Comparison (SERVER-FRIENDLY)
#   DGP: design-based N_ij + GN sizes|K + activity->smallest cluster
#   Fits: BT, RCBTL (Pearce/Erosheva), BT-SBM (+ optional bt_sbm)
#   Output: chunk CSVs + merged single big CSV
# ================================================================

# -----------------------------
# 0) Packages
# -----------------------------
suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(mcclust)
  library(mcclust.ext)
  library(fossil)
  library(coda)
})

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# -----------------------------
# 1) Paths (server friendly)
# -----------------------------
proj_dir <- normalizePath(Sys.getenv("PROJECT_DIR", getwd()), mustWork = FALSE)
data_dir <- file.path(proj_dir, "data")
res_dir  <- file.path(proj_dir, "results")
chunks_dir <- file.path(res_dir, "chunks")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)

data_path_real <- file.path(data_dir, "2000_2022_data.rds")
stopifnot(file.exists(data_path_real))
data_real <- readRDS(data_path_real)

# Use 3 observed designs (edit freely)
N_list <- list(
  design1 = data_real[[1]]$N_ij,
  design2 = data_real[[2]]$N_ij,
  design3 = data_real[[3]]$N_ij
)

# single merged output
out_file <- file.path(res_dir, "simulation_comparison_all_models.csv")

# overwrite policy
if (isTRUE(as.logical(Sys.getenv("OVERWRITE_RESULTS", "TRUE")))) {
  if (file.exists(out_file)) file.remove(out_file)
  # remove old chunks
  old_chunks <- list.files(chunks_dir, full.names = TRUE, pattern = "\\.csv$")
  if (length(old_chunks) > 0) file.remove(old_chunks)
}

# -----------------------------
# 2) Parallel setup
# -----------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
n_cores <- max(1L, n_cores)

cl <- makeCluster(n_cores)
registerDoParallel(cl)

on.exit({
  try(stopCluster(cl), silent = TRUE)
}, add = TRUE)

# Load packages on each worker (avoid library() inside foreach)
clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(mcclust)
    library(mcclust.ext)
    library(fossil)
    library(coda)
  })
  NULL
})

# If these packages exist, load them on workers too
if (has_pkg("BTSBM")) clusterEvalQ(cl, { suppressPackageStartupMessages(library(BTSBM)); NULL })
if (has_pkg("rankclust")) clusterEvalQ(cl, { suppressPackageStartupMessages(library(rankclust)); NULL })

# -----------------------------
# 3) Utilities
# -----------------------------
time_fit <- function(expr) {
  err <- NA_character_
  res <- NULL
  tt <- system.time({
    res <- tryCatch(expr, error = function(e) { err <<- conditionMessage(e); NULL })
  })
  list(res = res, time = tt, error = err)
}

norm_geo1 <- function(v) {
  v <- as.numeric(v)
  v / exp(mean(log(v)))
}

compare_lambda <- function(est, truth) {
  est <- as.numeric(est)
  truth <- as.numeric(truth)
  tibble(
    rmse = sqrt(mean((est - truth)^2)),
    mae  = mean(abs(est - truth)),
    spearman = suppressWarnings(cor(est, truth, method = "spearman")),
    kendall  = suppressWarnings(cor(est, truth, method = "kendall")),
    pearson  = suppressWarnings(cor(est, truth, method = "pearson"))
  )
}

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

mode_int <- function(x) {
  tx <- table(x)
  as.integer(names(tx)[which.max(tx)])
}

# -----------------------------
# 4) GN sampler + conditional sizes
# -----------------------------
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

# Assign most active players to smallest clusters (label 1 = smallest)
make_z_activity_to_smallest <- function(sizes, activity, seed = NULL, shuffle_within_cluster = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(sum(sizes) == length(activity))
  
  K <- length(sizes)
  ord_sizes <- order(sizes, decreasing = FALSE)  # smallest first
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

# -----------------------------
# 5) DGP (strengths + outcomes)
# -----------------------------
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

# -----------------------------
# 6) Relabel BT-SBM draws (pragmatic)
# -----------------------------
relabel_by_lambda_draws <- function(x_samples, lambda_samples) {
  stopifnot(is.matrix(x_samples), is.matrix(lambda_samples))
  stopifnot(nrow(x_samples) == nrow(lambda_samples), ncol(x_samples) == ncol(lambda_samples))
  
  S <- nrow(x_samples)
  x_out <- x_samples
  
  for (s in seq_len(S)) {
    x <- x_samples[s, ]
    lam <- lambda_samples[s, ]
    labs <- sort(unique(x))
    cl_mean <- sapply(labs, function(k) mean(lam[x == k]))
    ord <- order(cl_mean, decreasing = TRUE)
    map <- setNames(seq_along(labs), labs[ord])
    x_out[s, ] <- as.integer(map[as.character(x)])
  }
  
  list(x_samples_relabel = x_out, lambda_samples_relabel = lambda_samples)
}

summarise_btsbm_clustering <- function(fit_btsbm, z_true, K_true) {
  x_samples <- fit_btsbm$x_samples
  lambda_samples <- fit_btsbm$lambda_samples
  stopifnot(is.matrix(x_samples), is.matrix(lambda_samples))
  
  rel <- relabel_by_lambda_draws(x_samples, lambda_samples)
  x_rel <- rel$x_samples_relabel
  
  K_draws <- apply(x_rel, 1, function(x) length(unique(x)))
  
  psm <- mcclust::comp.psm(x_rel)
  minVI <- mcclust.ext::minVI(psm)$cl
  binder <- mcclust.ext::minbinder(psm)$cl
  
  hpd <- coda::HPDinterval(coda::mcmc(K_draws))
  
  tibble(
    K_median = as.integer(round(median(K_draws), 0)),
    K_mode = mode_int(K_draws),
    pr_K_true = mean(K_draws == K_true),
    HPD_low = as.integer(hpd[, 1]),
    HPD_high = as.integer(hpd[, 2]),
    ari_minVI = fossil::adj.rand.index(minVI, z_true),
    vi_minVI  = mcclust::vi.dist(minVI, z_true),
    ari_binder = fossil::adj.rand.index(binder, z_true),
    vi_binder  = mcclust::vi.dist(binder, z_true)
  )
}

# -----------------------------
# 7) Fit wrappers
# -----------------------------
get_fun <- function(pkg, fname) {
  if (has_pkg(pkg) && exists(fname, where = asNamespace(pkg), mode = "function")) {
    return(get(fname, envir = asNamespace(pkg)))
  }
  if (exists(fname, mode = "function")) return(get(fname, mode = "function"))
  NULL
}

fit_bt <- function(W, T_iter, T_burn, verbose = FALSE) {
  fun <- get_fun("BTSBM", "gibbs_bt_simple")
  if (is.null(fun)) stop("Cannot find BTSBM::gibbs_bt_simple (or gibbs_bt_simple).")
  fun(w_ij = as.matrix(W), T_iter = T_iter, T_burn = T_burn, verbose = verbose)
}

fit_btsbm <- function(W, T_iter, T_burn, a, gamma, prior = "GN", verbose = FALSE) {
  fun <- get_fun("BTSBM", "gibbs_bt_sbm")
  if (is.null(fun)) stop("Cannot find BTSBM::gibbs_bt_sbm (or gibbs_bt_sbm).")
  fun(w_ij = as.matrix(W), T_iter = T_iter, T_burn = T_burn, a = a, gamma = gamma, prior = prior, verbose = verbose)
}

fit_optional_bt_sbm <- function(W, ...) {
  cand <- c("fit_bt_sbm", "bt_sbm")
  for (nm in cand) {
    fun <- get_fun("BTSBM", nm)
    if (!is.null(fun)) return(list(name = nm, fit = fun(w_ij = as.matrix(W), ...)))
  }
  NULL
}

maybe_source_rankclust_code <- function() {
  # If mcmc_BTL/mcmc_RCBTL aren't in scope, source user file.
  if (!exists("mcmc_BTL", mode = "function") || !exists("mcmc_RCBTL", mode = "function")) {
    src <- Sys.getenv("RANKCLUST_MCMC_PATH", "")
    if (nzchar(src) && file.exists(src)) source(src)
  }
}

fit_rankclust_rcbtl <- function(W, T_iter, T_burn, seed = 1) {
  maybe_source_rankclust_code()
  if (!exists("mcmc_BTL", mode = "function") || !exists("mcmc_RCBTL", mode = "function")) {
    stop("mcmc_BTL/mcmc_RCBTL not found. Set env var RANKCLUST_MCMC_PATH to a file that defines them.")
  }
  
  Pi <- wins_to_Pi(as.matrix(W))
  n <- nrow(W)
  burn_prop <- T_burn / T_iter
  
  res_btl <- mcmc_BTL(
    Pi = Pi, J = n,
    a_gamma = 1, b_gamma = 1,
    num_iters = T_iter,
    burn_prop = burn_prop,
    chains = 1,
    groupwise = TRUE,
    seed = seed
  )
  
  omega_cols <- grep("^omega", names(res_btl))
  if (length(omega_cols) == 0) stop("mcmc_BTL output missing omega columns.")
  
  nu0 <- colMeans(res_btl[, omega_cols, drop = FALSE])
  
  res_rcbtl <- mcmc_RCBTL(
    Pi = Pi, J = n,
    a_gamma = 1, b_gamma = 1,
    lambda = 2,
    nu0 = nu0,
    num_iters = T_iter,
    nu_reps = 2,
    burn_prop = burn_prop,
    thin = 1,
    chains = 1,
    groupwise = TRUE,
    seed = seed
  )
  
  list(res_btl = res_btl, res_rcbtl = res_rcbtl)
}

# -----------------------------
# 8) Extract lambdas
# -----------------------------
extract_bt_lambda_hat <- function(fit_bt) {
  if (is.null(fit_bt$lambda_samples)) stop("BT fit missing lambda_samples.")
  norm_geo1(colMeans(fit_bt$lambda_samples))
}

extract_btsbm_lambda_hat <- function(fit_btsbm) {
  if (is.null(fit_btsbm$lambda_samples)) stop("BT-SBM fit missing lambda_samples.")
  norm_geo1(colMeans(fit_btsbm$lambda_samples))
}

extract_rankclust_lambda_hat <- function(fit_rankclust) {
  res_rcbtl <- fit_rankclust$res_rcbtl
  omega_cols <- grep("^omega", names(res_rcbtl))
  if (length(omega_cols) == 0) stop("RCBTL output missing omega columns.")
  norm_geo1(colMeans(res_rcbtl[, omega_cols, drop = FALSE]))
}

# -----------------------------
# 9) One replicate runner
# -----------------------------
run_one <- function(design_name, N, rep_id,
                    K_true, gamma_true, p_adj,
                    T_iter, T_burn,
                    a_fit, gamma_fit,
                    seed_base = 1,
                    verbose = FALSE) {
  
  seed <- seed_base + 100000 * rep_id + 1000 * K_true + match(design_name, names(N_list)) * 10
  sim <- simulate_dataset_design_based(N, K_target = K_true, gamma_true = gamma_true, p_adj = p_adj, seed = seed)
  
  W <- sim$W
  z_true <- sim$z_true
  lambda_true <- sim$lambda_true_players
  n <- nrow(N)
  
  base_row <- tibble(
    design = design_name,
    rep = rep_id,
    n_players = n,
    K_true = K_true,
    gamma_true = gamma_true,
    p_adj = p_adj,
    seed = seed
  )
  
  rows <- list()
  
  # ---- BT
  bt_fit <- time_fit(fit_bt(W, T_iter = T_iter, T_burn = T_burn, verbose = verbose))
  if (is.na(bt_fit$error)) {
    bt_hat <- extract_bt_lambda_hat(bt_fit$res)
    perf <- compare_lambda(bt_hat, lambda_true)
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = "BT",
             time_user = unname(bt_fit$time["user.self"]),
             time_sys = unname(bt_fit$time["sys.self"]),
             time_elapsed = unname(bt_fit$time["elapsed"]),
             fit_error = NA_character_),
      perf
    )
  } else {
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = "BT",
             time_user = unname(bt_fit$time["user.self"]),
             time_sys = unname(bt_fit$time["sys.self"]),
             time_elapsed = unname(bt_fit$time["elapsed"]),
             fit_error = bt_fit$error),
      tibble(rmse = NA_real_, mae = NA_real_, spearman = NA_real_, kendall = NA_real_, pearson = NA_real_)
    )
  }
  
  # ---- RCBTL
  rc_fit <- time_fit(fit_rankclust_rcbtl(W, T_iter = T_iter, T_burn = T_burn, seed = seed))
  if (is.na(rc_fit$error)) {
    rc_hat <- extract_rankclust_lambda_hat(rc_fit$res)
    perf <- compare_lambda(rc_hat, lambda_true)
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = "RCBTL",
             time_user = unname(rc_fit$time["user.self"]),
             time_sys = unname(rc_fit$time["sys.self"]),
             time_elapsed = unname(rc_fit$time["elapsed"]),
             fit_error = NA_character_),
      perf
    )
  } else {
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = "RCBTL",
             time_user = unname(rc_fit$time["user.self"]),
             time_sys = unname(rc_fit$time["sys.self"]),
             time_elapsed = unname(rc_fit$time["elapsed"]),
             fit_error = rc_fit$error),
      tibble(rmse = NA_real_, mae = NA_real_, spearman = NA_real_, kendall = NA_real_, pearson = NA_real_)
    )
  }
  
  # ---- BT-SBM
  btsbm_fit <- time_fit(fit_btsbm(W, T_iter = T_iter, T_burn = T_burn, a = a_fit, gamma = gamma_fit, prior = "GN", verbose = verbose))
  if (is.na(btsbm_fit$error)) {
    btsbm_hat <- extract_btsbm_lambda_hat(btsbm_fit$res)
    perf <- compare_lambda(btsbm_hat, lambda_true)
    cl_summ <- summarise_btsbm_clustering(btsbm_fit$res, z_true = z_true, K_true = K_true)
    
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = "BT-SBM",
             time_user = unname(btsbm_fit$time["user.self"]),
             time_sys = unname(btsbm_fit$time["sys.self"]),
             time_elapsed = unname(btsbm_fit$time["elapsed"]),
             fit_error = NA_character_),
      perf,
      cl_summ
    )
  } else {
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = "BT-SBM",
             time_user = unname(btsbm_fit$time["user.self"]),
             time_sys = unname(btsbm_fit$time["sys.self"]),
             time_elapsed = unname(btsbm_fit$time["elapsed"]),
             fit_error = btsbm_fit$error),
      tibble(rmse = NA_real_, mae = NA_real_, spearman = NA_real_, kendall = NA_real_, pearson = NA_real_,
             K_median = NA_integer_, K_mode = NA_integer_, pr_K_true = NA_real_,
             HPD_low = NA_integer_, HPD_high = NA_integer_,
             ari_minVI = NA_real_, vi_minVI = NA_real_, ari_binder = NA_real_, vi_binder = NA_real_)
    )
  }
  
  # ---- Optional baseline if present
  opt <- fit_optional_bt_sbm(W, T_iter = T_iter, T_burn = T_burn)
  if (!is.null(opt)) {
    # only timing (and lambda recovery if lambda_samples exists)
    opt_name <- paste0("OPTIONAL_", opt$name)
    # many optional fits won't return timing; wrap with time_fit if needed
    rmse <- mae <- spearman <- kendall <- pearson <- NA_real_
    if (!is.null(opt$fit$lambda_samples)) {
      opt_hat <- norm_geo1(colMeans(opt$fit$lambda_samples))
      perf <- compare_lambda(opt_hat, lambda_true)
      rmse <- perf$rmse; mae <- perf$mae; spearman <- perf$spearman; kendall <- perf$kendall; pearson <- perf$pearson
    }
    rows[[length(rows) + 1]] <- bind_cols(
      base_row,
      tibble(model = opt_name,
             time_user = NA_real_, time_sys = NA_real_, time_elapsed = NA_real_,
             fit_error = NA_character_,
             rmse = rmse, mae = mae, spearman = spearman, kendall = kendall, pearson = pearson)
    )
  }
  
  bind_rows(rows)
}

# -----------------------------
# 10) Merge chunks to big CSV
# -----------------------------
merge_chunks_to_csv <- function(chunk_files, out_file) {
  chunk_files <- chunk_files[file.exists(chunk_files)]
  if (length(chunk_files) == 0) stop("No chunk files found to merge.")
  
  first <- TRUE
  for (f in chunk_files) {
    df <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
    write.table(df, file = out_file, sep = ",",
                row.names = FALSE,
                col.names = first,
                append = !first)
    first <- FALSE
  }
  invisible(out_file)
}

# -----------------------------
# 11) Simulation grid + run
# -----------------------------
K_grid <- c(3, 5, 7)
R_reps <- as.integer(Sys.getenv("N_REPS", "10"))

gamma_true <- 0.71
p_adj <- 0.85

T_iter <- as.integer(Sys.getenv("T_ITER", "5000"))
T_burn <- as.integer(Sys.getenv("T_BURN", "1000"))

a_fit <- 2
gamma_fit <- 0.71

seed_base <- 123
verbose <- FALSE

grid <- expand.grid(
  design = names(N_list),
  rep = seq_len(R_reps),
  K_true = K_grid,
  stringsAsFactors = FALSE
)

# Export needed objects/functions to workers (important for PSOCK)
clusterExport(
  cl,
  varlist = c(
    "N_list", "simulate_dataset_design_based", "run_one",
    "sample_gnedin_labels", "draw_gnedin_sizes_given_K", "make_z_activity_to_smallest",
    "check_top_player_in_smallest_cluster", "make_lambda_from_padj", "simulate_W_from_blocks",
    "time_fit", "norm_geo1", "compare_lambda", "wins_to_Pi", "mode_int",
    "relabel_by_lambda_draws", "summarise_btsbm_clustering",
    "has_pkg", "get_fun", "fit_bt", "fit_btsbm", "fit_optional_bt_sbm",
    "maybe_source_rankclust_code", "fit_rankclust_rcbtl",
    "extract_bt_lambda_hat", "extract_btsbm_lambda_hat", "extract_rankclust_lambda_hat"
  ),
  envir = environment()
)

chunk_files <- foreach(ii = seq_len(nrow(grid)),
                       .packages = c("readr", "dplyr", "tidyr", "mcclust", "mcclust.ext", "fossil", "coda"),
                       .errorhandling = "pass") %dopar% {
                         
                         design_name <- grid$design[ii]
                         rep_id <- grid$rep[ii]
                         K_true <- grid$K_true[ii]
                         N <- N_list[[design_name]]
                         
                         out <- tryCatch(
                           run_one(
                             design_name = design_name,
                             N = N,
                             rep_id = rep_id,
                             K_true = K_true,
                             gamma_true = gamma_true,
                             p_adj = p_adj,
                             T_iter = T_iter,
                             T_burn = T_burn,
                             a_fit = a_fit,
                             gamma_fit = gamma_fit,
                             seed_base = seed_base,
                             verbose = verbose
                           ),
                           error = function(e) {
                             # still write an informative row
                             tibble(
                               design = design_name, rep = rep_id, n_players = nrow(N),
                               K_true = K_true, gamma_true = gamma_true, p_adj = p_adj,
                               seed = NA_integer_, model = "WORKER_ERROR",
                               time_user = NA_real_, time_sys = NA_real_, time_elapsed = NA_real_,
                               fit_error = conditionMessage(e),
                               rmse = NA_real_, mae = NA_real_, spearman = NA_real_, kendall = NA_real_, pearson = NA_real_
                             )
                           }
                         )
                         
                         # write chunk file unique per task
                         chunk_file <- file.path(
                           chunks_dir,
                           sprintf("chunk_%s_K%d_rep%03d_pid%d.csv", design_name, K_true, rep_id, Sys.getpid())
                         )
                         readr::write_csv(out, chunk_file)
                         chunk_file
                       }

# Merge on master
chunk_files <- unlist(chunk_files)
merge_chunks_to_csv(chunk_files, out_file)

message("Done. Merged results written to: ", out_file)

# Optionally clean chunks after merge
if (isTRUE(as.logical(Sys.getenv("CLEAN_CHUNKS", "FALSE")))) {
  file.remove(chunk_files[file.exists(chunk_files)])
}

# -----------------------------
# 12) (COMMENTED) LaTeX printing
# -----------------------------
# library(knitr)
# library(kableExtra)
# ... later ...


# 
# 
# library(knitr)
# library(kableExtra)
# library(dplyr)
# 
# perf_tex <- perf %>%
#   mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
#   rename(
#     Model = model,
#     RMSE = rmse,
#     MAE = mae,
#     Spearman = spearman,
#     Kendall = kendall,
#     Pearson = pearson
#   )
# 
# tab1 <- kable(
#   perf_tex,
#   format = "latex",
#   booktabs = TRUE,
#   caption = "Comparison of posterior mean strength estimates across models. All strengths are normalised to geometric mean 1.",
#   align = "lccccc"
# ) %>%
#   kable_styling(latex_options = c("hold_position"))
# 
# cat(tab1, file = "results/performance_table.tex")
# lambda_tex <- lambda_table %>%
#   mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
#   rename(
#     Player = player,
#     Cluster = x_true,
#     `True $\\lambda$` = lambda_true,
#     BT = bt,
#     RCBTL = rcbtl,
#     `BT--SBM` = btsbm
#   )
# 
# tab2 <- kable(
#   lambda_tex,
#   format = "latex",
#   booktabs = TRUE,
#   caption = "True and estimated player strengths ($\\lambda$) under each model. Values scaled to geometric mean 1.",
#   align = "lcccccc"
# ) %>%
#   kable_styling(latex_options = c("scale_down"))
# 
# cat(tab2, file = "results/lambda_table.tex")
# rank_tex <- rank_comp %>%
#   select(player, bt_rank, rcbtl_rank, btsbm_rank, x_true) %>%
#   rename(
#     Player = player,
#     `BT rank` = bt_rank,
#     `RCBTL rank` = rcbtl_rank,
#     `BT--SBM rank` = btsbm_rank,
#     Cluster = x_true
#   )
# 
# tab3 <- kable(
#   rank_tex,
#   format = "latex",
#   booktabs = TRUE,
#   caption = "Player ranks induced by posterior mean strengths.",
#   align = "lcccc"
# )
# 
# cat(tab3, file = "results/rank_table.tex")
# library(dplyr)
# library(knitr)
# library(kableExtra)
# 
# timings_tex <- timings %>%
#   mutate(
#     elapsed_min = elapsed / 60,
#     user_min    = user_time / 60,
#     sys_min     = sys_time / 60
#   ) %>%
#   transmute(
#     Model = model,
#     `User (min)`    = round(user_min, 2),
#     `System (min)`  = round(sys_min, 2),
#     `Elapsed (min)` = round(elapsed_min, 2)
#   )
# 
# tab_time <- kable(
#   timings_tex,
#   format = "latex",
#   booktabs = TRUE,
#   caption = "Computation time for each algorithm. Times are reported in minutes.",
#   align = "lccc"
# ) %>%
#   kable_styling(latex_options = c("hold_position"))
# 
# cat(tab_time, file = "results/timings_table.tex")
# 
# 
# 
# 
# 
# 
# 
