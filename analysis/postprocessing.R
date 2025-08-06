# -------------------------------
# Post-processing and Summary Plots for BT-SBM
# -------------------------------
library(mcclust)
library(mcclust.ext)
library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)
library(ggrepel)
library(kableExtra)
library(reshape2)
library(stringr)
library(ggside)
# Load helper functions
source("main/utils.R")
source("main/mcmc_bt_sbm.R")

# Load full MCMC results across seasons
res_list <- readRDS("results/augmented_multiple_seasonsGN.rds")
first_year <- 1999

# Initialize containers
top_block_counts_across_years <- data.frame(season = character(), avg_top_block_cnt = numeric())
prob_assignment_across_years  <- data.frame()
post_numb_block_across_years  <- data.frame(season = character(), num_blocks = integer(), count = integer(), prob = numeric())
avg_strength_each_player      <- data.frame()
entropy_per_season <- data.frame(season = character(),
                                 quantile005 = numeric(), 
                                 quantile095 = numeric(), 
                                 mean_entropy = numeric())

# -------------------------------
# Main Loop Over Seasons
# -------------------------------
for (yr in seq_along(res_list)) {
  season_label <- paste0(first_year + yr, "/", first_year + yr + 1)
  res_i <- res_list[[yr]]
  x_samples <- res_i$x_samples
  lambda_samples <- res_i$lambda_samples
  
  inf_i <- inference_helper(x_samples, lambda_samples)
  x_relabeled <- inf_i$x_samples_relabel
  lambdas_reordered <- inf_i$lambda_samples_relabel
  T_iter <- nrow(lambdas_reordered)
  n_players <- ncol(lambdas_reordered)
  
  # Player-level skill computation and entropy per iteration
  pl_lambda <- matrix(0, nrow = T_iter, ncol = n_players)
  entropy_container <- numeric(T_iter)
  
  for (i in 1:T_iter) {
    x_cur <- x_relabeled[i, ]
    p1 <- mean(x_cur == 1)
    p2 <- 1 - p1
    entropy_container[i] <- -sum(c(p1, p2) * log(c(p1, p2) + 1e-10))
    
    lambda_cur <- lambdas_reordered[i, ]
    unique_lambdas <- unique(lambda_cur)
    new_unique_lambdas <- unique_lambdas / (prod(unique_lambdas)^(1 / length(unique_lambdas)))
    pl_lambda[i, ] <- new_unique_lambdas[x_cur]
  }
  
  HPD_entropy <- HPDinterval(as.mcmc(entropy_container))
  entropy_per_season <- rbind(entropy_per_season, data.frame(
    season = season_label,
    quantile005 = HPD_entropy[1],
    quantile095 = HPD_entropy[2],
    mean_entropy = mean(entropy_container)
  ))
  
  avg_strength_each_player <- rbind(avg_strength_each_player, data.frame(
    season = rep(season_label, n_players),
    mean_str = apply(pl_lambda, 2, median, na.rm = TRUE),
    lower_quantile = apply(pl_lambda, 2, quantile, probs = 0.025),
    upper_quantile = apply(pl_lambda, 2, quantile, probs = 0.975)
  ))
  
  top_block_counts_across_years <- rbind(top_block_counts_across_years, data.frame(
    season = season_label,
    avg_top_block_cnt = inf_i$avg_top_block_count
  ))
  
  block_assignment_i <- inf_i$player_block_assignment_probs
  block_assignment_i$season <- season_label
  block_assignment_i$pl_name    <- rownames(tennis_years[[yr]]$Y_ij)
  block_assignment_i$eos_ranking <- tennis_years[[yr]]$players_df$last_rank
  
  # Clean underscores
  block_assignment_i$pl_name <- gsub(pattern = "_", 
                                     replacement = " ", 
                                     x = block_assignment_i$pl_name)
  
  prob_assignment_across_years <- rbind(prob_assignment_across_years, block_assignment_i)
  
  post_numb_blocks_i <- inf_i$block_count_distribution
  post_numb_blocks_i$season <- season_label
  post_numb_block_across_years <- rbind(post_numb_block_across_years, post_numb_blocks_i)
  
  message("Processed season: ", season_label)
}

# Save outputs
# write.csv(avg_strength_each_player, "results/avg_strength_each_player.csv", row.names = FALSE)
# write.csv(prob_assignment_across_years, "results/prob_assignment_across_years.csv", row.names = FALSE)
# write.csv(entropy_per_season, "results/entropy_per_season.csv", row.names = FALSE)
# write.csv(post_numb_block_across_years, "results/post_num_blocks.csv", row.names = FALSE)

# Generate plots and tables 
# ------------------------------------------------------------------------
# Probability on the number of clusters for each season
# Table in kable that has rows = season, columns = #blocks, 
# and cells = probability of that #blocks
# ------------------------------------------------------------------------

post_numb_block_across_years_wide <- post_numb_block_across_years %>%
  filter(num_blocks < 7) %>%
  select(-count)%>%
  pivot_wider(
    names_from  = num_blocks,
    values_from = prob
  )

# Identify the most likely number of blocks per season
num_blocks_season <- post_numb_block_across_years %>%
  group_by(season) %>%
  summarise(num_blocks = num_blocks[ which.max(prob) ]) %>%
  mutate(num_blocks = factor(num_blocks, ordered = TRUE))

post_numb_block_across_years_table = cbind(post_numb_block_across_years_wide,recap_table[,'better_model'] )
# Render table
kable(post_numb_block_across_years_table, format = 'latex', digits = 3)


# ------------------------------------------------------------------------
# Jittered scatterplot: Probability of top-block membership by season
# ------------------------------------------------------------------------
prob_assignment_across_years <- prob_assignment_across_years %>%
  left_join(num_blocks_season, by = 'season')

# Convert 'season' to a factor for nicer plotting
prob_assignment_across_years$season <- factor(
  prob_assignment_across_years$season,
  levels = unique(prob_assignment_across_years$season)
)


p_top_across_time = ggplot(prob_assignment_across_years, aes(x = season, y = Cluster_1)) +
  geom_jitter(aes(color = factor(num_blocks)), width = 0.2, alpha = 0.6, size = 3) +
  labs(
    x     = "Season",
    y     = "P(Top Block)",
    color = "Nº of blocks"
  ) +
  scale_color_manual(values = c(
    "3" = "#2E8B57",  # verde Wimbledon
    "4" = "#D2691E",  # terra rossa
    "5" = "#1E90FF"   # blu cemento
  )) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.minor = element_blank()
)

ggsave(plot = p_top_across_time,filename = "./images/Ptop_across_time1.png",
       width = 9, height=5)

cluster_stats <- prob_assignment_across_years %>%
  group_by(season) %>%
  mutate(Cluster_1 = Cluster_1/sum(Cluster_1))%>%
  summarize(
    shannon_cluster1 = shannon_entropy(Cluster_1),
    .groups = "drop"
  )



entropy_plot<- entropy_per_season %>%
  mutate(season_start = as.numeric(sub("/.*", "", season))) %>%  # estrai 2000 da "2000/2001"
  ggplot(aes(x = season_start, y = mean_entropy)) +
  geom_ribbon(aes(ymin = quantile005, ymax = quantile095), 
              fill = "#2E8B57", alpha = 0.5) +     # verde Wimbledon
  geom_line(color = "#FFD700", size = 2) +         # giallo pallina
  geom_point(color = "#FFD700", size = 4) +
  scale_x_continuous(breaks = unique(as.numeric(sub("/.*", "", entropy_per_season$season)))) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Season",
    y = "Shannon Entropy"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("./images/entropy_plot1.png",plot = entropy_plot,
       width = 13, height=5)


# ------------------------------------------------------------------------
# Bar plot: Estimated number of players in the top block by season
# ------------------------------------------------------------------------
num_block_plot = prob_assignment_across_years %>%
  pivot_longer(
    cols       = -c(season, pl_name, num_blocks,eos_ranking),
    names_to   = "cluster",
    values_to  = "prob"
  ) %>%
  group_by(season, pl_name) %>%
  # For each player, choose the cluster with highest prob
  reframe(ass_cluster = cluster[ which.max(prob) ]) %>%
  ungroup() %>%
  group_by(season) %>%
  count(ass_cluster) %>%
  # We only want the count of top-block membership
  filter(ass_cluster == "Cluster_1") %>%
  left_join(num_blocks_season, by = 'season') %>%
  ggplot(aes(x = season, y = n, fill = num_blocks)) +
  geom_col() +
  scale_fill_manual(values = c(
    "3" = "#2E8B57",  # verde Wimbledon
    "4" = "#D2691E",  # terra rossa
    "5" = "#1E90FF"   # blu cemento
  ))+
  theme_minimal() +
  labs(x     = "Season",
    y     = "Nº of Players in Top Block",
    fill  = "Nº blocks"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )


ggsave(file = "./images/num_block_plot1.png",num_block_plot,
       width=9,height = 5)

# ------------------------------------------------------------------------
# Player trajectory across the seasons
# Compare P(top block) vs (100 - end-of-season ranking), for example
# ------------------------------------------------------------------------
pl_selected <- c("Rafael Nadal", "Roger Federer", "Novak Djokovic", "Andy Murray")

prob_assignment_across_years %>%
  filter(pl_name %in% pl_selected) %>%
  ggplot(aes(x = season, y = Cluster_1 * 100, group = pl_name, color = pl_name)) +
  geom_line() +
  theme_minimal() +
  labs(
    y     = "P(Top Block) (%)",
    x     = "Season",
    color = "Player"
  ) +
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~pl_name)




# --------------------------------------------
# Focused Analysis on a Single Season (BT-SBM)
# --------------------------------------------

# Setup
first_year <- 1999
yr <- 18  # season index (e.g. 2000/01 + 18 = 2018/2019)
model <- 'GN'
season_label <- paste0(first_year + yr, "/", first_year + yr + 1)
print(season_label)

# Load posterior draws
res_i <- readRDS(paste0("results/augmented_multiple_seasons", model, ".rds"))[[yr]]
x_samples <- res_i$x_samples
lambda_samples <- res_i$lambda_samples

# External data for this season
tennis_years = readRDS("./data/2000_2022_data.rds")
Y_ij <- tennis_years[[yr]]$Y_ij
N_ij <- tennis_years[[yr]]$N_ij
pl_df <- tennis_years[[yr]]$players_df

# ------------------------------
# Helper: relabel draws by λ rank
# ------------------------------
relabel_by_lambda <- function(x_draws, lambda_draws) {
  K <- ncol(lambda_draws)
  out_x <- x_draws
  out_lambda <- lambda_draws
  
  for (i in seq_len(nrow(lambda_draws))) {
    ord <- order(lambda_draws[i, ], decreasing = TRUE)
    out_lambda[i, ] <- lambda_draws[i, ord]
    renamer <- match(seq_len(K), ord)
    out_x[i, ] <- renamer[x_draws[i, ]]
  }
  list(x = out_x, lambda = out_lambda)
}

#define custom palette
custom_palette <- c(
  "0" = "white",       # missing data
  "1" = "#CDEB8B",     
  "2" = "#78AB46",     
  "3" = "#FFD700",    
  "4" = "#FF8C00",    
  "5" = "#00441B"      
)

# ------------------------------
# Subset draws with target K=5
# ------------------------------
k_target <- 5
keep <- apply(x_samples, 1, \(z) length(unique(z)) == k_target)
x_k5 <- x_samples[keep, , drop = FALSE]
lambda_k5 <- lambda_samples[keep, , drop = FALSE]

relab <- relabel_by_lambda(x_k5, lambda_k5)
x_k5 <- relab$x
lambda_k5 <- relab$lambda

# ------------------------------
# Build λ matrix by player ID
# ------------------------------
idx_mat <- cbind(rep(seq_len(nrow(x_k5)), each = ncol(x_k5)),
                 as.vector(x_k5))
lambda_player <- matrix(lambda_k5[idx_mat],
                        nrow = nrow(x_k5),
                        ncol = ncol(x_k5),
                        byrow = FALSE)
colnames(lambda_player) <- colnames(x_k5)


# ------------------------------
# λ uncertainty plot
# ------------------------------
log_lp <- log10(lambda_player)
hpd90 <- t(apply(log_lp, 2, \(v) HPDinterval(as.mcmc(v), prob = 0.90)))

player_summ <- data.frame(
  player  = colnames(Y_ij),
  mean    = 10^(colMeans(log_lp)),
  low90   = 10^(hpd90[, 1]),
  up90    = 10^(hpd90[, 2]),
  cluster = inference_helper(x_samples, lambda_samples)$partition_expected
)

plot_lambda <- player_summ |>
  arrange(cluster, desc(mean)) |>
  mutate(
    player = sapply(str_to_title(gsub("_", " ", player)), clean_players_names),
  )|>
  mutate(
    player = factor(player, levels = rev(player)),
    surname = gsub(".*[ _]", "", player)
  ) |>
  ggplot(aes(x = mean, y = player, colour = factor(cluster))) +
  geom_pointrange(aes(xmin = low90, xmax = up90), size = 0.4, fatten = 0.6) +
  scale_x_log10() +
  scale_colour_manual(values = custom_palette) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(x = expression(lambda~"(posterior mean, log"[10]*" scale)"),
       y = NULL,
       color = "Cluster") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# ggsave("images/lambda_uncertainty.png", plot_lambda, height = 5, width = 5)

# ------------------------------
# Heatmap of Block-Ordered Adjacency
# ------------------------------

inf_i <- inference_helper(res_i$x_samples, lambda_samples = res_i$lambda_samples)
x_relabel <- inf_i$x_samples_relabel
lambdas_reordered <- inf_i$lambda_samples_relabel


# Add cluster info to players
partition_minVI <- inference_helper(x_relabel, lambdas_reordered)$partition_expected
df_cl <- data.frame(
  players = rownames(Y_ij),
  cl = partition_minVI,
  marginal_win = rowSums(Y_ij)
)

# Summarize block compositions
block_players <- aggregate(players ~ cl, data = df_cl, FUN = \(x) paste(x, collapse = ", "))
block_players$Block_Size <- sapply(strsplit(block_players$players, ", "), length)
kable(df_cl, format = 'latex', row.names = FALSE)
kable(block_players[, c("cl", "Block_Size", "players")], format = 'latex')

# Prepare long format for adjacency matrix
Y_long <- melt(Y_ij)
colnames(Y_long) <- c("Winner", "Loser", "Win_Count")
Y_long$Matches_Count <- melt(N_ij)$value

K <- 5

Y_long_plot <- Y_long %>%
  mutate(perc_success = Win_Count / Matches_Count) %>%
  left_join(df_cl, by = c("Loser" = "players")) %>%
  rename(row_cl = cl, marginal_win_row = marginal_win) %>%
  left_join(df_cl, by = c("Winner" = "players")) %>%
  rename(col_cl = cl, marginal_win_col = marginal_win) %>%
  mutate(
    Winner = sapply(str_to_title(gsub("_", " ", Winner)), clean_players_names),
    Loser  = sapply(str_to_title(gsub("_", " ", Loser)),  clean_players_names)
  )

# Reorder factors by cluster and marginal wins
Y_long_plot <- Y_long_plot %>%
  mutate(
    Winner = factor(Winner, levels = unique(Winner[order(col_cl, -marginal_win_col,decreasing = T)])),
    Loser  = factor(Loser,  levels = unique(Loser[order(row_cl, -marginal_win_row)])),
    col_cl = factor(col_cl, ordered = TRUE)
  )
# 1. Get block boundaries for the x-axis (i.e. losers) using row_cl:
v_lines_list <- Y_long_plot %>%
  group_by(row_cl) %>%
  summarize(x_break = max(as.numeric(Loser)), .groups = "drop") %>%
  pull(x_break)

# Remove the last boundary so we don't draw a line at the extreme edge:
v_lines_list <- v_lines_list[-length(v_lines_list)]

# 2. Get block boundaries for the y-axis (i.e. winners) using col_cl:
h_lines_list <- Y_long_plot %>%
  group_by(col_cl) %>%
  summarize(y_break = min(as.numeric(Winner)), .groups = "drop") %>%
  pull(y_break)

# Remove the last boundary similarly:
h_lines_list <- h_lines_list[-length(h_lines_list)]

geom_adjacency_fixed <- ggplot(Y_long_plot, aes(x = Loser, y = Winner)) +
  geom_tile(aes(fill = perc_success), color = 'grey40') +
  scale_fill_gradient(low = "#FFFFCC", high = "#006400", na.value = "#009680") +
  geom_ysidecol(aes(color = factor(col_cl))) +
  scale_color_manual(values = custom_palette) +
  geom_vline(xintercept = unlist(v_lines_list) + 0.5, color = 'black', linewidth = 0.3) +
  geom_hline(yintercept = unlist(h_lines_list) - 0.5, color = 'black', linewidth = 0.3) +
  labs(
    x = "Players (ordered by block)",
    y = "Players (ordered by block)",
    fill = "% victories",
    color = "Block"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "left"
  ) +
  theme_ggside_void() +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
  coord_fixed(ratio = 1)

geom_adjacency_fixed

# ------------------------------
# Heatmap: Assignment Probabilities
# ------------------------------

count_cl = function(x){
  length(unique(x))
}

x_samples       <- res_i$x_samples
lambda_samples  <- res_i$lambda_samples

unique_count = apply(x_samples,1,count_cl)
unique_count <- apply(x_samples, 1, count_cl)
x_samples_sub <- x_samples[which(unique_count == 5), ]
lambda_samples_sub <- lambda_samples[which(unique_count == 5), ]
inf_i <- inference_helper(x_samples_sub, lambda_samples_sub)

block_prob <- inf_i$player_block_assignment_probs[, 1:5]
block_prob <- cbind(block_prob, pl_name = rownames(Y_ij))

assignment_probs_long <- block_prob %>%
  pivot_longer(cols = 1:5, names_to = "Cluster", values_to = "prob") %>%
  mutate(
    Cluster = gsub(x = Cluster, replacement = " ", pattern = "_"),
    pl_name = gsub(x = pl_name, replacement = " ", pattern = "_")
  )

# Cluster with highest assignment probability
max_prob_clusters <- assignment_probs_long %>%
  group_by(pl_name) %>%
  summarize(Cl_ass = Cluster[which.max(prob)], .groups = "drop")

# Marginal win percentages
marg_pro_win <- data.frame(
  pl_name = rownames(Y_ij),
  marg_pro_win  = rowSums(Y_ij),
  marg_pro_loss = colSums(Y_ij)
) %>%
  mutate(
    pct_win = marg_pro_win / (marg_pro_loss + marg_pro_win),
    pl_name = gsub(x = pl_name, replacement = " ", pattern = "_")
  )

assignment_probs_long_plot <- assignment_probs_long %>%
  left_join(max_prob_clusters, by = "pl_name") %>%
  left_join(marg_pro_win, by = "pl_name") %>%
  ungroup()|>
  mutate(
    pl_name = sapply(str_to_title(gsub("_", " ", pl_name)), clean_players_names),
  )|>
  mutate(pl_name = factor(pl_name,
                          levels = unique(pl_name[order(Cl_ass, -marg_pro_win,decreasing = T)]),
                          ordered = TRUE))

ass_prob_plot <- ggplot(assignment_probs_long_plot) +
  geom_tile(aes(x = Cluster, y = pl_name, fill = prob)) +
  scale_fill_gradient(low = "#FFFFCC", high = "#006400", na.value = "#009680") +
  labs(x = "", y = "", fill = "Assign. Prob.") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2))

 ass_prob_plot
 
 ggsave(filename = "./images/ass_prob_plot.png",ass_prob_plot, height = 8, width = 10)
 ggsave(filename = "./images/geom_adjacency_fixed.png",geom_adjacency_fixed, height = 8, width = 10)
 ggsave(filename = "./images/lambda_uncertainty.png",plot_lambda, height = 6, width = 7)
 
 
