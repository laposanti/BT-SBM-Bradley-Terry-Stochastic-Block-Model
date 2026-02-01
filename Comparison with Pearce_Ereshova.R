## Beta (Github) version
# install.packages("devtools") # uncomment if you haven't installed 'devtools' before
# devtools::install_github("pearce790/rankclust", force = TRUE)

library(rankclust)
library(BTSBM)
library(dplyr)
library(tidyr)
library(readr)

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

co_cluster_prob <- function(x_samples) {
  n <- ncol(x_samples)
  p_same <- matrix(0, n, n)
  for (s in seq_len(nrow(x_samples))) {
    g <- x_samples[s, ]
    p_same <- p_same + outer(g, g, "==")
  }
  p_same / nrow(x_samples)
}

rank_from_strength <- function(strengths, player_names) {
  tibble::tibble(
    player = player_names,
    strength = strengths,
    rank = rank(-strengths, ties.method = "average")
  )
}

# -------------------------------
# Load 2017 season data (same source as single_season_analysis.R)
# -------------------------------
first_year <- 1999
season_start_year <- 2017
season_index <- season_start_year - first_year
season_label <- paste0(season_start_year, "/", season_start_year + 1)

if (!dir.exists("./results")) dir.create("./results", recursive = TRUE)

# External data for this season
# (path matches single_season_analysis.R)
tennis_years <- readRDS("./data/kallog_2000_2026_top105_data.rds")

w_ij <- tennis_years[[season_index]]$Y_ij
players_df <- tennis_years[[season_index]]$players_df
player_names <- players_df$player_slug
if (is.null(player_names) || length(player_names) != nrow(w_ij)) {
  player_names <- rownames(w_ij)
}

Pi <- wins_to_Pi(w_ij)

# -------------------------------
# Fit models
# -------------------------------
# 1) Plain BT
fit_bt <- BTSBM::gibbs_bt_simple(
  w_ij = w_ij,
  T_iter = 3000,
  T_burn = 600,
  verbose = TRUE
)

# 2) Ereshova rank-clustered (rankclust)
res_btl <- mcmc_BTL(
  Pi = Pi,
  J = nrow(w_ij),
  a_gamma = 1,
  b_gamma = 1,
  num_iters = 3000,
  burn_prop = 0.2,
  chains = 1,
  groupwise = TRUE,
  seed = 1
)

nu0 <- colMeans(res_btl[, grep("^omega", names(res_btl))])

res_rcbtl <- mcmc_RCBTL(
  Pi = Pi,
  J = nrow(w_ij),
  a_gamma = 1,
  b_gamma = 1,
  lambda = 2,
  nu0 = nu0,
  num_iters = 3000,
  nu_reps = 2,
  burn_prop = 0.2,
  thin = 1,
  chains = 1,
  groupwise = TRUE,
  seed = 1
)

# 3) BTSBM
fit_btsbm <- BTSBM::gibbs_bt_sbm(
  A = as.matrix(w_ij),
  num_iters = 3000,
  burn_in = 600,
  alpha = 2,
  beta = 2,
  chains = 1,
  groupwise = TRUE
)

# -------------------------------
# Save raw results
# -------------------------------
raw_results <- list(
  season = season_label,
  season_index = season_index,
  w_ij = w_ij,
  players_df = players_df,
  bt = fit_bt,
  rcbtl = res_rcbtl,
  btsbm = fit_btsbm
)

saveRDS(raw_results, file = "./results/pearce_ereshova_2017_raw_results.rds")

# -------------------------------
# Post-hoc analysis
# -------------------------------
# Strength summaries
bt_strength <- colMeans(fit_bt$lambda_samples)
rcbtl_strength <- colMeans(res_rcbtl[, grep("^omega", names(res_rcbtl))])
btsbm_strength <- colMeans(fit_btsbm$lambda_samples)

bt_ranks <- rank_from_strength(bt_strength, player_names)
rcbtl_ranks <- rank_from_strength(rcbtl_strength, player_names)
btsbm_ranks <- rank_from_strength(btsbm_strength, player_names)

rank_comparison <- bt_ranks %>%
  rename(bt_strength = strength, bt_rank = rank) %>%
  inner_join(rcbtl_ranks %>% rename(rcbtl_strength = strength, rcbtl_rank = rank),
             by = "player") %>%
  inner_join(btsbm_ranks %>% rename(btsbm_strength = strength, btsbm_rank = rank),
             by = "player")

rank_correlations <- tibble::tibble(
  metric = c("spearman", "kendall"),
  bt_vs_rcbtl = c(
    cor(rank_comparison$bt_strength, rank_comparison$rcbtl_strength, method = "spearman"),
    cor(rank_comparison$bt_strength, rank_comparison$rcbtl_strength, method = "kendall")
  ),
  bt_vs_btsbm = c(
    cor(rank_comparison$bt_strength, rank_comparison$btsbm_strength, method = "spearman"),
    cor(rank_comparison$bt_strength, rank_comparison$btsbm_strength, method = "kendall")
  ),
  rcbtl_vs_btsbm = c(
    cor(rank_comparison$rcbtl_strength, rank_comparison$btsbm_strength, method = "spearman"),
    cor(rank_comparison$rcbtl_strength, rank_comparison$btsbm_strength, method = "kendall")
  )
)

# Top-10 overlap
get_top_players <- function(rank_df, n = 10) {
  rank_df %>% arrange(rank) %>% slice_head(n = n) %>% pull(player)
}

bt_top <- get_top_players(bt_ranks)
rcbtl_top <- get_top_players(rcbtl_ranks)
btsbm_top <- get_top_players(btsbm_ranks)

_top_overlap <- function(x, y) length(intersect(x, y))

rank_top10_overlap <- tibble::tibble(
  bt_rcbtl = _top_overlap(bt_top, rcbtl_top),
  bt_btsbm = _top_overlap(bt_top, btsbm_top),
  rcbtl_btsbm = _top_overlap(rcbtl_top, btsbm_top)
)

# Co-clustering comparison (rank-clustered vs BTSBM)
rcbtl_groups <- as.matrix(res_rcbtl[, grep("^G", names(res_rcbtl))])
rcbtl_co_cluster <- co_cluster_prob(rcbtl_groups)

btsbm_co_cluster <- co_cluster_prob(fit_btsbm$x_samples)

co_cluster_similarity <- tibble::tibble(
  frobenius_distance = sqrt(sum((rcbtl_co_cluster - btsbm_co_cluster)^2)),
  pearson_correlation = cor(as.vector(rcbtl_co_cluster), as.vector(btsbm_co_cluster))
)

posthoc_results <- list(
  rank_comparison = rank_comparison,
  rank_correlations = rank_correlations,
  rank_top10_overlap = rank_top10_overlap,
  co_cluster_similarity = co_cluster_similarity
)

saveRDS(posthoc_results, file = "./results/pearce_ereshova_2017_posthoc.rds")

readr::write_csv(rank_comparison, "./results/pearce_ereshova_2017_rankings.csv")
readr::write_csv(rank_correlations, "./results/pearce_ereshova_2017_rank_correlations.csv")
readr::write_csv(rank_top10_overlap, "./results/pearce_ereshova_2017_top10_overlap.csv")
readr::write_csv(co_cluster_similarity, "./results/pearce_ereshova_2017_cluster_similarity.csv")

message("Saved raw and post-hoc results for season ", season_label)
