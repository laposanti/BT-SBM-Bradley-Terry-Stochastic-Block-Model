# -------------------------------
# Multiple seasons Analysis
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

# Set your project root 
#setwd('/.../current folder')

tennis_years <- readRDS("./data/ATP_2000_2025_SN_extended.rds")
# Load full MCMC results across seasons
res_list <- readRDS("raw_output_ext/MCMC_raw_output_ext.rds")

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

  inf_i <- relabel_by_lambda(x_samples, lambda_samples)
  x_relabeled <- inf_i$x_samples_relabel
  lambdas_reordered <- inf_i$lambda_samples_relabel
  T_iter <- nrow(lambdas_reordered)
  n_players <- ncol(lambdas_reordered)

  # Player-level skill computation and entropy per iteration
  pl_lambda <- matrix(0, nrow = T_iter, ncol = n_players)
  entropy_container <- numeric(T_iter)

  for (i in seq_len(T_iter)) {
    x_cur <- x_relabeled[i, ]
    p1 <- mean(x_cur == 1)
    p <- c(p1, 1 - p1)
    p <- p[p > 0]                    # ignore empty categories

    if (length(p) <= 1) {
      entropy_container[i] <- 0      # degenerate: all mass in one bin
    } else {
      H <- -sum(p * log(p))
      entropy_container[i] <- H / log(length(p))  # here length(p) == 2
    }


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
  
  w_ij <- tennis_years[[yr]]$Y_ij #pairwise success matrix
  pl_df <- tennis_years[[yr]]$players_df #info about players
  
  block_assignment_i <- inf_i$item_cluster_assignment_probs
  block_assignment_i$season <- season_label
  block_assignment_i$pl_name  <- pl_df$player_label
  block_assignment_i$eos_ranking <- pl_df$last_rank

  # Clean underscores
  block_assignment_i$pl_name <- pl_df$player_label

  prob_assignment_across_years <- rbind(prob_assignment_across_years, block_assignment_i)

  post_numb_blocks_i <- inf_i$block_count_distribution
  post_numb_blocks_i$season <- season_label
  post_numb_block_across_years <- rbind(post_numb_block_across_years, post_numb_blocks_i)

  message("Processed season: ", season_label)
}


# Generate plots and tables

# ------------------------------------------------------------------------
# Probability on the number of clusters for each season
# Table in kable that has rows = season, columns = #blocks,
# and cells = probability of that #blocks
# ------------------------------------------------------------------------

post_numb_block_across_years_wide <- post_numb_block_across_years %>%
  filter(num_blocks < 9) %>%
  dplyr::select(-count)%>%
  pivot_wider(
    names_from  = num_blocks,
    values_from = prob
  )


# Render table
latex_table <- kable(post_numb_block_across_years_wide, format = "latex", digits = 3, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "hold_position"))



# ------------------------------------------------------------------------
# Jittered scatterplot: Probability of top-block membership by season
# ------------------------------------------------------------------------
num_blocks_season <- post_numb_block_across_years %>%
  dplyr::group_by(season) %>%
  dplyr::summarise(
    num_blocks = num_blocks[which.max(prob)],
    .groups = "drop"
  )

prob_assignment_across_years <- prob_assignment_across_years %>%
  left_join(num_blocks_season, by = 'season') #num_blocks_season does not exist


# Convert 'season' to a factor for nicer plotting
prob_assignment_across_years$season <- factor(
  prob_assignment_across_years$season,
  levels = unique(prob_assignment_across_years$season)
)

p_top_across_time = ggplot(prob_assignment_across_years, aes(x = season, y = Cluster_1)) +
  geom_jitter(aes(color = factor(num_blocks.y)), width = 0.2, alpha = 0.6, size = 3) +
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

# ------------------------------------------------------------------------
# Entropy plot: indicator of the dispersion of the top-block membership
# ------------------------------------------------------------------------

#functions that computes the Shannon entropy
shannon_entropy <- function(x, base = exp(1), na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  # Use only occupied categories to match H / log(K_occ)
  x <- x[x > 0]
  s <- sum(x)
  if (s <= 0) return(NA_real_)   # undefined if no mass
  p <- x / s
  H <- -sum(p * log(p, base = base))
  K <- length(p)
  if (K <= 1) return(0)          # degenerate: one occupied block → zero entropy
  H / log(K, base = base)
}

cluster_stats <- prob_assignment_across_years %>%
  group_by(season) %>%
  mutate(Cluster_1 = Cluster_1/sum(Cluster_1))%>%
  summarize(
    shannon_cluster1 = shannon_entropy(Cluster_1),
    .groups = "drop"
  )

#Plot
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


# ------------------------------------------------------------------------
# Bar plot: Estimated number of players in the top block by season
# ------------------------------------------------------------------------
num_block_plot = num_blocks_season%>%
  mutate(num_blocks = as.factor(num_blocks)) %>%
  ggplot(aes(x = season, y = num_blocks, fill = num_blocks)) +
  geom_col() +
  scale_fill_manual(values = c(
    "2" = "#FFD700",  # giallo pallina
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



# ------------------------------------------------------------------------
# Extra: Player trajectory across the seasons
# Compare P(top block) vs (100 - end-of-season ranking), for example
# ------------------------------------------------------------------------
pl_selected <- c("Nadal R.", "Federer R.", "Djokovic N.", "Murray A.", "Alcaraz C.", "Sinner J.")

prob_assignment_across_years %>%
  filter(pl_label %in% pl_selected) %>%
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

#-----
# Saving the outputs
#-----


# Save the .tex file in the tables folder 
writeLines(latex_table, "./tables/post_numb_block_across_years_table1.tex")

#Save the plots in the images folder
ggsave(filename = "./images/entropy_plot.pdf",
       plot = entropy_plot,
       width = 13, height=5)

ggsave(filename = "./images/Ptop_across_time.pdf",
       plot = p_top_across_time,
       width = 9, height=5)


ggsave(filename = "./images/num_block_plot.pdf",
       plot = num_block_plot,
       width=9,height = 5)


