# -------------------------------
# Single seasons Analysis
# -------------------------------

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

ggsave(filename = "./images/ass_prob_plot1.png",ass_prob_plot, height = 8, width = 10)
ggsave(filename = "./images/geom_adjacency_fixed1.png",geom_adjacency_fixed, height = 8, width = 10)
ggsave(filename = "./images/lambda_uncertainty1.png",plot_lambda, height = 6, width = 7)


