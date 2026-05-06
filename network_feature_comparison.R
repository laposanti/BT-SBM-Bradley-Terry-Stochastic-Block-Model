#!/usr/bin/env Rscript

# ==============================================================
# Compare network features: simulated W vs real W
#
# For each of the 3 observed designs in N_list_for_simulation.rds,
# we:
#  1) find the matching real-data season in data/kallog_2000_2026_top105_data.rds
#     by exact equality of N_ij
#  2) compute simple network features for:
#       - real wins matrix W_real (= Y_ij)
#       - simulated wins matrix W_sim generated using the same N
#
# Outputs:
#  - results/network_features_sim_vs_real.csv
# ==============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# -----------------------------
# Config
# -----------------------------
proj_dir <- normalizePath(Sys.getenv("PROJECT_DIR", getwd()), mustWork = FALSE)
res_dir <- file.path(proj_dir, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# N-list: historically singular, but allow plural filename too.
# You can also override with N_LIST_FILE (path can be relative to PROJECT_DIR).
n_list_override <- Sys.getenv("N_LIST_FILE", "")
n_list_candidates <- c(
  if (nzchar(n_list_override)) file.path(proj_dir, n_list_override) else character(0),
  file.path(proj_dir, "N_list_for_simulations.rds"),
  file.path(proj_dir, "N_list_for_simulation.rds")
)
n_list_file <- n_list_candidates[file.exists(n_list_candidates)][1]
real_data_file <- file.path(proj_dir, "data", "ATP_2000_2026_SN_extended.rds")

# Simulation DGP params (do not affect N; only how wins are split)
K_true <- as.integer(Sys.getenv("K_TRUE", "3"))
gamma_true <- as.numeric(Sys.getenv("GAMMA_TRUE", "0.71"))
p_adj <- as.numeric(Sys.getenv("P_ADJ", "0.85"))
seed_base <- as.integer(Sys.getenv("SEED_BASE", "123"))

# Optional manual mapping if auto-match fails.
# Example: DESIGN_YEAR_MAP="design1:2018,design2:2019,design3:2022"
manual_map <- Sys.getenv("DESIGN_YEAR_MAP", "")

# -----------------------------
# Minimal simulation helpers (mirrors RUN_SIMULATION_STUDY.R)
# -----------------------------

norm_geo1 <- function(v) {
  v <- as.numeric(v)
  v <- pmax(v, .Machine$double.xmin)
  v / exp(mean(log(v)))
}

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
# Feature extraction
# -----------------------------

network_features <- function(N, W) {
  stopifnot(is.matrix(N), is.matrix(W), nrow(N) == ncol(N), all(dim(N) == dim(W)))
  n <- nrow(N)

  ut <- upper.tri(N, diag = FALSE)
  n_pairs <- n * (n - 1) / 2

  edge_present <- ut & (N > 0)
  n_edges <- sum(edge_present)
  density <- n_edges / n_pairs
  sparsity <- 1 - density

  total_matches <- sum(N[ut])
  total_wins <- sum(W) # should equal total_matches if W is consistent with N

  # sanity: each match contributes exactly one win across both directions
  if (!isTRUE(all.equal(as.numeric(total_wins), as.numeric(total_matches), tolerance = 1e-8))) {
    warning(
      "Inconsistent W vs N: sum(W) != sum(N[upper.tri]). ",
      "Check that W[i,j] + W[j,i] == N[i,j] for all i<j."
    )
  }

  # win-structure metrics (these distinguish real vs simulated W)
  nij <- N[ut]
  wij <- W[ut]
  p_ij <- ifelse(nij > 0, wij / nij, NA_real_)
  pair_margin <- abs(p_ij - 0.5)
  one_sided <- (nij > 0) & (wij == 0 | wij == nij)

  # binary entropy (in nats) for each played pair; 0 at p in {0,1}, max at 0.5
  entropy_bin <- function(p) {
    p <- pmin(pmax(p, .Machine$double.eps), 1 - .Machine$double.eps)
    -(p * log(p) + (1 - p) * log(1 - p))
  }

  # player win rates
  matches_i <- rowSums(N)
  wins_i <- rowSums(W)
  win_rate_i <- ifelse(matches_i > 0, wins_i / matches_i, NA_real_)

  gini <- function(x) {
    x <- x[is.finite(x) & !is.na(x)]
    if (length(x) == 0) return(NA_real_)
    if (all(x == 0)) return(0)
    x <- sort(as.numeric(x))
    n <- length(x)
    (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
  }

  directed_win_edges <- sum(W > 0) - sum(diag(W) > 0)
  players_with_matches <- sum(matches_i > 0)
  players_with_wins <- sum(wins_i > 0)

  tibble(
    n_players = n,
    n_pairs = n_pairs,
    n_edges = n_edges,
    density = density,
    sparsity = sparsity,
    total_matches = total_matches,
    total_wins = total_wins,
    directed_win_edges = directed_win_edges,
    one_sided_pairs = sum(one_sided),
    one_sided_fraction = ifelse(n_edges > 0, sum(one_sided) / n_edges, NA_real_),
    mean_pair_margin = mean(pair_margin, na.rm = TRUE),
    sd_pair_margin = sd(pair_margin, na.rm = TRUE),
    mean_pair_entropy = mean(entropy_bin(p_ij), na.rm = TRUE),
    sd_win_rate = sd(win_rate_i, na.rm = TRUE),
    gini_wins = gini(wins_i),
    players_with_matches = players_with_matches,
    players_with_wins = players_with_wins
  )
}

# -----------------------------
# Find real-data year for each design
# -----------------------------

parse_manual_map <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) return(list())
  parts <- strsplit(x, ",")[[1]]
  out <- list()
  for (p in parts) {
    kv <- strsplit(trimws(p), ":")[[1]]
    if (length(kv) != 2) next
    out[[trimws(kv[1])]] <- trimws(kv[2])
  }
  out
}

find_matching_year_by_N <- function(N_target, real_list) {
  N_target <- as.matrix(N_target)
  for (yr in names(real_list)) {
    obj <- real_list[[yr]]
    if (!is.list(obj) || is.null(obj$N_ij)) next
    Nij <- as.matrix(obj$N_ij)
    if (!all(dim(Nij) == dim(N_target))) next

    # exact match (after integer coercion)
    if (isTRUE(all.equal(as.integer(Nij), as.integer(N_target), check.attributes = FALSE))) {
      return(yr)
    }
  }
  NA_character_
}

# -----------------------------
# Main
# -----------------------------

if (is.na(n_list_file) || !nzchar(n_list_file) || !file.exists(n_list_file)) {
  stop(
    "Missing N-list file. Looked for one of: ",
    paste(basename(n_list_candidates), collapse = ", "),
    "\nYou can set N_LIST_FILE to override (relative to PROJECT_DIR or absolute)."
  )
}
if (!file.exists(real_data_file)) stop("Missing: ", real_data_file)

N_list <- readr::read_rds(n_list_file)
real <- readRDS(real_data_file)

design_names <- names(N_list)
if (is.null(design_names) || any(!nzchar(design_names))) {
  design_names <- paste0("design", seq_along(N_list))
  names(N_list) <- design_names
}

manual <- parse_manual_map(manual_map)

rows <- list()

for (design in names(N_list)) {
  N <- as.matrix(N_list[[design]])

  yr <- manual[[design]]
  if (is.null(yr) || !nzchar(yr)) {
    yr <- find_matching_year_by_N(N, real)
  }

  if (is.na(yr) || !nzchar(yr) || is.null(real[[yr]])) {
    stop(
      "Could not find matching real-data year for ", design, ".\n",
      "Provide mapping via DESIGN_YEAR_MAP, e.g. DESIGN_YEAR_MAP=\"design1:2018,design2:2019,design3:2022\""
    )
  }

  W_real <- as.matrix(real[[yr]]$Y_ij)
  N_real <- as.matrix(real[[yr]]$N_ij)

  if (!isTRUE(all.equal(as.integer(N_real), as.integer(N), check.attributes = FALSE))) {
    stop("Matched year ", yr, " for ", design, " but N_ij differs from N_list.")
  }

  seed <- seed_base + match(design, names(N_list)) * 10
  sim <- simulate_dataset_design_based(N, K_target = K_true, gamma_true = gamma_true, p_adj = p_adj, seed = seed)
  W_sim <- sim$W

  f_real <- network_features(N_real, W_real) %>%
    mutate(design = design, source = "real", year = yr)
  f_sim <- network_features(N, W_sim) %>%
    mutate(design = design, source = "simulated", year = yr,
           K_true = K_true, gamma_true = gamma_true, p_adj = p_adj, seed = seed)

  rows[[length(rows) + 1]] <- bind_rows(f_real, f_sim)
}

out <- bind_rows(rows) %>%
  relocate(design, year, source)

out_file <- file.path(res_dir, "network_features_sim_vs_real.csv")
write_csv(out, out_file)

cat("Wrote: ", out_file, "\n", sep = "")
print(out, n = Inf, width = Inf)



# an Adjacency matrix that instead of having Item 1, it has I.1 on columns
BTSBM::plot_block_adjacency(w_ij = sim$W, x_hat = sim$z_true, clean_fun = function(x) gsub("Item ", "Pl ", x))
