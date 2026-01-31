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
library(BTSBM)

# -------------------------------
# Single seasons Analysis
# -------------------------------
#change this
current_wd = "/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model/"
setwd(current_wd)

# ensure images dir exists
if (!dir.exists("./images")) dir.create("./images", recursive = TRUE)


# Setup
first_year <- 1999
yr <- 18  # season index (e.g. 2000/01 + 18 = 2018/2019)
model <- 'GN'
season_label <- paste0(first_year + yr, "/", first_year + yr + 1)
print(season_label)

# Load the output of gibbs_bt_sbm()
fit <- readRDS(paste0("./raw_output_ext/MCMC_raw_output_ext.rds"))[[yr]]
x_samples <- fit$x_samples
lambda_samples <- fit$lambda_samples

# External data for this season
tennis_years <- readRDS("./data/kallog_2000_2026_top105_data.rds")

w_ij <- tennis_years[[yr]]$Y_ij #pairwise success matrix
pl_df <- tennis_years[[yr]]$players_df #info about players

#Relabeling posterior draws
inf_i <- relabel_by_lambda(x_samples, lambda_samples = lambda_samples)

BTSBM::pretty_table_K_distribution(inf_i)

count_K = function(x) length(unique(x))
coda::HPDinterval(mcmc(apply(x_samples,1,count_K)))

#Reordered heatmap
reordered_heatmap <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$minVI_partition,
                                          clean_fun = function(x) as.character(x))  # remove outer padding)
reordered_heatmap

#Uncertainty over the assignment
ass_prob_plot <- plot_assignment_probabilities(inf_i,w_ij = w_ij,max_n_clust = 4,
                                               clean_fun = function(x) as.character(x))
ass_prob_plot

#Posterior Lambdas, and their uncertainty
plot_lambda <- plot_lambda_uncertainty(inf_i,w_ij = w_ij,max_n_clust = 3,
                                      clean_fun = function(x) as.character(x))
plot_lambda

ggsave(filename = "./images/reordered_heatmap.png",reordered_heatmap, height = 8, width = 10)
ggsave(filename = "./images/ass_prob_plot.png",ass_prob_plot, height = 8, width = 10)
ggsave(filename = "./images/lambda_uncertainty2.png",plot_lambda, height = 6, width = 7)


#-------------------------------
# Saving also the credible balls

#vertical upper bound (coarsest partition)
reordered_heatmap_vert_ub <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$credible_ball_upper_partition)+theme( plot.margin = unit(c(0, 0, 0, 0), "pt"))  # remove outer padding
#vertical lower bound (finest partition)
reordered_heatmap_vert_lb <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$credible_ball_lower_partition)+theme( plot.margin = unit(c(0, 0, 0, 0), "pt"))
#Horizontal bound (complexity unaware)
reordered_heatmap_horiz <- plot_block_adjacency(fit = inf_i,w_ij = w_ij,x_hat = inf_i$credible_ball_horiz_partition)+theme( plot.margin = unit(c(0, 0, 0, 0), "pt"))


ggsave(filename = "./images/reordered_heatmap_v_ub.png",reordered_heatmap_vert_ub, height = 8, width = 10)
ggsave(filename = "./images/reordered_heatmap_v_lb.png",reordered_heatmap_vert_lb, height = 8, width = 10)
ggsave(filename = "./images/reordered_heatmap_horiz.png",reordered_heatmap_horiz, height = 8, width = 10)






library(cowplot)
library(grid)  # unit()

# Remove *everything* axis/legend-ish
strip_xy_legend <- theme(
  legend.position = "none",
  axis.title.y   = element_blank(),
  axis.title.x   = element_blank(),
  axis.text.y    = element_blank(),
  axis.text.x    = element_blank(),
  axis.ticks.y   = element_blank(),
  axis.ticks.x   = element_blank(),
  strip.text.y   = element_blank(),
  strip.text.x   = element_blank(),
  strip.background = element_blank()
)

pad_lr1 <- theme(plot.margin = unit(c(0, 1, 0, 1), "pt"))

# ---- Plots ----
# Point estimate (rename to whatever your PE partition is)
reordered_heatmap_pe <- plot_block_adjacency(
  fit = inf_i, w_ij = w_ij, x_hat = inf_i$point_estimate_partition  # <<< adjust if different
) + strip_xy_legend + pad_lr1 + guides(fill = "none", colour = "none")

reordered_heatmap_vert_ub <- plot_block_adjacency(
  fit = inf_i, w_ij = w_ij, x_hat = inf_i$credible_ball_upper_partition
) + strip_xy_legend + pad_lr1 + guides(fill = "none", colour = "none")

reordered_heatmap_vert_lb <- plot_block_adjacency(
  fit = inf_i, w_ij = w_ij, x_hat = inf_i$credible_ball_lower_partition
) + strip_xy_legend + pad_lr1 + guides(fill = "none", colour = "none")

reordered_heatmap_horiz <- plot_block_adjacency(
  fit = inf_i, w_ij = w_ij, x_hat = inf_i$credible_ball_horiz_partition
) + strip_xy_legend + pad_lr1 + guides(fill = "none", colour = "none")

# ---- Helper: bottom caption ----
add_bottom_title <- function(p, label, fontsize = 11, pad_bottom = 0.08) {
  ggdraw() +
    draw_plot(p, x = 0, y = pad_bottom, width = 1, height = 1 - pad_bottom) +
    draw_label(label, x = 0.5, y = pad_bottom * 0.5, vjust = 1, hjust = 0.5,
               fontface = "bold", size = fontsize)
}

p_pe  <- add_bottom_title(reordered_heatmap_pe,   "Point estimate")
p_ub  <- add_bottom_title(reordered_heatmap_vert_ub,  "Vertical upper bound")
p_lb  <- add_bottom_title(reordered_heatmap_vert_lb,  "Vertical lower bound")
p_hor <- add_bottom_title(reordered_heatmap_horiz,    "Horizontal bound")

# ---- Spacing controls ----
# Horizontal gap (fraction of row width) and vertical gap (fraction of column height)
hgap <- 0.04
vgap <- 0.06

pad_lr <- theme(plot.margin = unit(c(0, 4, 0, 4), "pt"))  # left/right = 4pt
p_ub2  <- p_ub  + pad_lr
p_lb2  <- p_lb  + pad_lr
p_hor2 <- p_hor + pad_lr

combined_heatmaps <- plot_grid(p_ub2, p_lb2, p_hor2, nrow = 1, align = "hv", axis = "tb")
combined_heatmaps
combined_heatmaps

combined_heatmaps
ggsave('./images/combined_1.png', combined_heatmaps)




# ============================================================
# Pairwise posterior:  Pr( item i ranked LOWER than item j )
# i lower than j  <=>  i is WEAKER than j  <=>  x_i > x_j
# (assuming relabel_by_lambda() makes label 1 = strongest tier)
# ============================================================

# --- Grab relabelled MCMC labels ---
# (Adjust the field name if your relabel_by_lambda() output differs)
x_post <- inf_i$x_samples_relabel
n <- ncol(x_post)
nrow(x_post)
# --- Player names + ordering by end-of-season ranking (rank 1 on top) ---
player_names <- pl_df$player_slug
if (is.null(player_names) || length(player_names) != n) {
  player_names <- rownames(Y_ij)
}
if (is.null(player_names)) player_names <- rownames(w_ij)

ord <- order(pl_df$last_rank, na.last = TRUE)   # smaller rank = better
player_ord <- player_names[ord]

# --- Compute P[i,j] = Pr( x_i < x_j | data ) over the full chain (no filtering on K) ---
P <- matrix(NA_real_, n, n)
for (i in seq_len(n)) {
  # vectorised: compares x_post[,i] to every column of x_post
  P[i, ] <- colMeans(x_post[, i] < x_post)
}
diag(P) <- NA_real_

# reorder rows/cols according to ranking order
P <- P[ord, ord]

# --- Long format for ggplot ---
P_long <- reshape2::melt(P, varnames = c("i", "j"), value.name = "prob") |>
  dplyr::mutate(
    i_name = factor(player_ord[i], levels = rev(player_ord)),  # rank 1 at top
    j_name = factor(player_ord[j], levels = player_ord)
  )

# --- Heatmap ---
prob_heatmap <- ggplot(P_long, aes(x = j_name, y = i_name, fill = prob)) +
  geom_tile(colour = "grey80", linewidth = 0.1) +
  coord_fixed() +
  scale_fill_gradient(
    limits = c(0, 1),
    low  = "#F7FCF5",   # very light green
    high = "#006D2C",   # dark green
    na.value = "grey85",
    name = "Posterior Pr(i above j)"  # <- tweak wording here
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_colourbar(title.position = "top", direction = "horizontal"))

prob_heatmap

ggsave("./images/prob_i_below_j.png", prob_heatmap, width = 10, height = 8)


# ============================================================
# Quantitative pairwise heatmap:  E[ log(lambda_i) - log(lambda_j) | data ]
# where lambda_i := lambda_{x_i} at each iteration (player inherits block strength)
# ============================================================

# relabelled draws
x_post <- inf_i$x_samples_relabel
lam_post <- inf_i$lambda_samples_relabel  # matrix: iterations x K_max (or list) — see note below


T <- nrow(x_post)
n <- ncol(x_post)

# player ordering (rank 1 on top)
player_names <- pl_df$player_slug
ord <- order(pl_df$last_rank, na.last = TRUE)
player_ord <- player_names[ord]

# ---- build player-level log strength per iteration: log_lambda_player[t,i] ----
# assumes lam_post[t, k] exists for all k that appear in x_post[t,]
log_lambda_player <- matrix(NA_real_, nrow = T, ncol = n)

for (t in seq_len(T)) {
  log_lambda_player[t, ] <- log(lam_post[t, x_post[t, ]])
}

# ---- pairwise posterior mean of differences ----
# D[i,j] = E[ logλ_i - logλ_j ]
D <- matrix(0, n, n)
for (i in seq_len(n)) {
  D[i, ] <- colMeans(log_lambda_player[, i] - log_lambda_player)
}
diag(D) <- NA_real_

# reorder by ranking
D <- D[ord, ord]

# long for ggplot
D_long <- reshape2::melt(D, varnames = c("i", "j"), value.name = "dloglambda") |>
  dplyr::mutate(
    i_name = factor(player_ord[i], levels = rev(player_ord)),
    j_name = factor(player_ord[j], levels = player_ord)
  )

# heatmap (diverging palette)
dlog_heatmap <- ggplot(D_long, aes(x = j_name, y = i_name, fill = dloglambda)) +
  geom_tile(colour = "grey80", linewidth = 0.1) +
  coord_fixed() +
  scale_fill_gradient2(
    low = "red",
    high = "#006D2C",
    na.value = "grey85",
    name = expression(E[log~tilde(lambda)[i] - log~tilde(lambda)[j]])
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_colourbar(title.position = "top", direction = "horizontal"))

dlog_heatmap
ggsave("./images/heatmap_Eloglambda_diff.png", dlog_heatmap, width = 10, height = 8)




library(dplyr)
library(tidyr)
library(knitr)
library(tibble)
library(kableExtra)

# ----------------------------
# Inputs
# ----------------------------
x_hat <- inf_i$minVI_partition
Khat  <- length(unique(x_hat))

lambda_draws <- inf_i$lambda_samples_relabel
# If lambda_draws is a list of draws, bind to a matrix
if (is.list(lambda_draws)) lambda_draws <- do.call(rbind, lambda_draws)

stopifnot(ncol(lambda_draws) == length(x_hat))

# If your stored parameter is actually log-skill (often called alpha),
# uncomment the next line:
# lambda_draws <- exp(lambda_draws)

# ----------------------------
# Block lambdas per MCMC draw
# ----------------------------
block_ids <- sort(unique(x_hat))
players_in_block <- lapply(block_ids, function(k) which(x_hat == k))

# lambda_block: (n_draws x Khat), each column is a block-level lambda per draw
lambda_block <- sapply(players_in_block, function(idx) rowMeans(lambda_draws[, idx, drop = FALSE]))
colnames(lambda_block) <- paste0("B", block_ids)

# ----------------------------
# Posterior summaries for p_rs = lambda_r / (lambda_r + lambda_s)
# ----------------------------
pair_grid <- expand.grid(r = block_ids, s = block_ids) |>
  dplyr::filter(r != s) |>
  dplyr::arrange(r, s)


library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(tibble)

make_block_matrix <- function(long_tbl, block_ids, value_col, caption) {
  block_ids <- sort(block_ids)
  
  wide <- long_tbl %>%
    mutate(
      r = factor(r, levels = block_ids),
      s = factor(s, levels = block_ids)
    ) %>%
    select(r, s, {{ value_col }}) %>%
    tidyr::pivot_wider(names_from = s, values_from = {{ value_col }}) %>%
    arrange(r)
  
  df <- as.data.frame(wide)
  r_lab <- as.character(df$r)
  df$r <- NULL
  
  # enforce column order = row order (critical bit)
  df <- df[, as.character(block_ids), drop = FALSE]
  
  rownames(df) <- paste0("B", r_lab)
  colnames(df) <- paste0("B", colnames(df))
  
  diag(df) <- NA_real_
  
  out <- tibble::rownames_to_column(df, var = "Group r")
  
  knitr::kable(
    out,
    format   = "latex",
    booktabs = TRUE,
    digits   = 3,
    na       = "",          # makes diagonal blank instead of "NA"
    caption  = caption
  ) %>%
    kableExtra::kable_styling(latex_options = c("hold_position"))
}


p_model_summ <- pair_grid |>
  rowwise() |>
  mutate(
    p_draw = list({
      lr <- lambda_block[, match(r, block_ids)]
      ls <- lambda_block[, match(s, block_ids)]
      lr / (lr + ls)
    }),
    p_mean = mean(p_draw[[1]]),
    p_sd   = sd(p_draw[[1]]),
    p_q025 = quantile(p_draw[[1]], 0.025),
    p_q975 = quantile(p_draw[[1]], 0.975)
  ) |>
  ungroup() |>
  select(r, s, p_mean, p_sd, p_q025, p_q975)


library(dplyr)
library(tidyr)
library(knitr)
library(tibble)
library(kableExtra)

# w_ij <- tennis_years[[yr]]$Y_ij
# If you already have N_ij, use it; otherwise:
n_ij <- w_ij + t(w_ij)

block_ids <- sort(unique(x_hat))

emp_summ <- expand.grid(r = block_ids, s = block_ids) |>
  dplyr::filter(r != s) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    W_rs = sum(w_ij[x_hat == r, x_hat == s, drop = FALSE]),
    N_rs = sum(n_ij[x_hat == r, x_hat == s, drop = FALSE]),
    p_emp = ifelse(N_rs > 0, W_rs / N_rs, NA_real_)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(r, s)


block_ids <- sort(unique(x_hat))

make_block_matrix(
  long_tbl  = p_model_summ,
  block_ids = block_ids,
  value_col = p_mean,
  caption   = "Model-based win probabilities between inferred blocks: E[\\lambda_r/(\\lambda_r+\\lambda_s)\\mid data]"
)

make_block_matrix(
  long_tbl  = emp_summ,
  block_ids = block_ids,
  value_col = p_emp,
  caption   = "Empirical win rates between inferred blocks: W_{r\\to s}/N_{rs}"
)
