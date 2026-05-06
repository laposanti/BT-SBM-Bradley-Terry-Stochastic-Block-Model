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
library(BTSBM)
# Set your project root 
#setwd('/.../current folder')
theme_btsbm <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size * 0.9),
      legend.text  = ggplot2::element_text(size = base_size * 0.85),
      plot.title   = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6))
    )
}

# Extendable palette for discrete #blocks (no hardcoding)
btsbm_block_palette <- function(K, palette = NULL) {
  K <- as.integer(K)
  if (!is.finite(K) || K < 1) stop("K must be a positive integer")

  if (!is.null(palette) && length(palette) >= K) return(palette[seq_len(K)])
  if (!is.null(palette) && length(palette) >= 2) return(grDevices::colorRampPalette(palette)(K))

  # Base R qualitative palette, works for many categories
  grDevices::hcl.colors(K, palette = "Dark 3")
}

theme_season_x <- function() {
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
}

season_start_year <- function(season_label) {
  as.integer(sub("/.*", "", season_label))
}

tennis_surface_palette <- function(levels_vec) {
  lv <- sort(unique(as.character(levels_vec)))
  lv <- lv[is.finite(as.numeric(lv))]
  if (length(lv) == 0L) return(setNames(character(0), character(0)))

  # Default palette for any K
  pal <- setNames(btsbm_block_palette(length(lv)), lv)

  # Prefer classic tennis-surface colors where applicable
  overrides <- c(
    "3" = "#2E8B57",  # green (grass)
    "4" = "#D2691E",  # clay
    "5" = "#1E90FF"   # hard court
  )
  for (nm in names(overrides)) {
    if (nm %in% names(pal)) pal[[nm]] <- overrides[[nm]]
  }
  pal
}

if (!exists("relabel_by_lambda", mode = "function")) {
  relabel_by_lambda <- function(x_samples, lambda_samples) {
    stopifnot(is.matrix(x_samples))
    S <- nrow(x_samples)
    N <- ncol(x_samples)
    is_list_format <- is.list(lambda_samples)
    get_lambda_vec <- function(iter) {
      if (is_list_format) {
        v <- lambda_samples[[iter]]
        if (!is.numeric(v)) {
          stop("lambda_samples[[iter]] must be numeric.")
        }
        v
      } else {
        lambda_samples[iter, ]
      }
    }
    x_relabeled <- matrix(NA_integer_, S, N)
    lambda_per_item <- matrix(NA_real_, S, N)
    cluster_lambda_ordered <- vector("list", S)
    n_clusters_each_iter <- integer(S)
    top_block_count_per_iter <- integer(S)
    for (iter in seq_len(S)) {
      xi <- as.integer(x_samples[iter, ])
      occ_raw <- sort(unique(xi))
      K <- length(occ_raw)
      lam_vec_full <- get_lambda_vec(iter)
      lam_occ <- rep(NA_real_, length(occ_raw))
      ok_idx <- occ_raw <= length(lam_vec_full)
      lam_occ[ok_idx] <- lam_vec_full[occ_raw[ok_idx]]
      ord <- order(lam_occ, decreasing = TRUE, na.last = TRUE)
      occ_ord <- occ_raw[ord]
      lam_ord <- lam_occ[ord]
      if (anyNA(lam_ord)) {
        lam_ord[is.na(lam_ord)] <- .Machine$double.xmin
      }
      # Fast relabel: map each raw label in xi to its rank in occ_ord
      xi_new <- match(xi, occ_ord)
      x_relabeled[iter, ] <- xi_new
      lambda_per_item[iter, ] <- lam_ord[xi_new]
      cluster_lambda_ordered[[iter]] <- lam_ord
      n_clusters_each_iter[iter] <- K
      top_block_count_per_iter[iter] <- sum(xi_new == 1L)
    }

    psm <- mcclust::comp.psm(x_samples)
    partition_binder <- mcclust.ext::minbinder.ext(psm, cls.draw = x_samples,
                                                   method = "all")$cl[1, ]
    partition_minVI <- mcclust.ext::minVI(psm, cls.draw = x_samples,
                                          method = "all")$cl[1, ]
    x_ball <- mcclust.ext::credibleball(c.star = partition_minVI,
                                        cls.draw = x_samples, c.dist = "VI")

    relabel_partition_by_item_mean_lambda <- function(z, lambda_item_mean) {
      stopifnot(length(z) == length(lambda_item_mean))
      z <- as.integer(z)
      labs <- sort(unique(z))
      cl_means <- vapply(labs, function(k) mean(lambda_item_mean[z == k], na.rm = TRUE), numeric(1))
      ord <- order(cl_means, decreasing = TRUE)
      new_ids <- seq_along(labs)
      names(new_ids) <- labs[ord]
      z_new <- new_ids[as.character(z)]
      as.integer(z_new)
    }
    lambda_item_mean <- colMeans(lambda_per_item, na.rm = TRUE)
    partition_minVI <- relabel_partition_by_item_mean_lambda(partition_minVI, lambda_item_mean)
    partition_binder <- relabel_partition_by_item_mean_lambda(partition_binder, lambda_item_mean)

    get_part <- function(obj, name1, name2) {
      if (!is.null(obj[[name1]])) obj[[name1]] else obj[[name2]]
    }
    c_lower_raw <- get_part(x_ball, "c.lower", "c.lowervert")
    c_upper_raw <- get_part(x_ball, "c.upper", "c.uppervert")
    c_horiz_raw <- x_ball$c.horiz

    pick_row <- function(obj, centre) {
      if (is.vector(obj) && length(obj) == N) return(as.integer(obj))
      if (is.matrix(obj) && ncol(obj) == N) {
        d <- apply(obj, 1, function(z) mcclust::vi.dist(as.integer(z), as.integer(centre)))
        return(as.integer(obj[which.min(d), ]))
      }
      stop("Unexpected credibleball partition format.")
    }

    c_lower_vec <- pick_row(c_lower_raw, partition_minVI)
    c_upper_vec <- pick_row(c_upper_raw, partition_minVI)
    c_horiz_vec <- pick_row(c_horiz_raw, partition_minVI)

    c_lower_rl <- relabel_partition_by_item_mean_lambda(c_lower_vec, lambda_item_mean)
    c_upper_rl <- relabel_partition_by_item_mean_lambda(c_upper_vec, lambda_item_mean)
    c_horiz_rl <- relabel_partition_by_item_mean_lambda(c_horiz_vec, lambda_item_mean)

    K_VI_upper <- length(unique(c_upper_rl))
    K_VI_lower <- length(unique(c_lower_rl))
    K_VI_horiz <- length(unique(c_lower_rl))
    # Item-cluster assignment probabilities
    # Previous code did colMeans(x_relabeled == k) for k=1..N => O(S*N^2).
    # This version is O(S*N) using per-item tabulation.
    Kmax <- N
    assignment_probs <- matrix(0, nrow = N, ncol = Kmax)
    for (j in seq_len(N)) {
      v <- x_relabeled[, j]
      v <- v[!is.na(v)]
      if (length(v) > 0L) {
        assignment_probs[j, ] <- tabulate(v, nbins = Kmax) / length(v)
      }
    }
    colnames(assignment_probs) <- paste0("Cluster_", seq_len(Kmax))
    rownames(assignment_probs) <- paste0("Item_", seq_len(N))
    assignment_probs_df <- as.data.frame(assignment_probs)
    bc_tab <- table(n_clusters_each_iter)
    block_count_df <- data.frame(
      num_blocks = as.integer(names(bc_tab)),
      count = as.vector(bc_tab),
      prob = as.vector(bc_tab) / sum(bc_tab)
    )

    tbc_tab <- table(top_block_count_per_iter)
    top_block_size_df <- data.frame(
      top_block_size = as.integer(names(tbc_tab)),
      count = as.vector(tbc_tab),
      prob = as.vector(tbc_tab) / sum(tbc_tab)
    )

    list(
      x_samples_relabel = x_relabeled,
      lambda_samples_relabel = lambda_per_item,
      cluster_lambda_ordered = cluster_lambda_ordered,
      co_clustering = psm,
      minVI_partition = partition_minVI,
      partition_binder = partition_binder,
      n_clusters_each_iter = n_clusters_each_iter,
      block_count_distribution = block_count_df,
      item_cluster_assignment_probs = assignment_probs_df,
      avg_top_block_count = mean(top_block_count_per_iter),
      top_block_count_per_iter = top_block_count_per_iter,
      top_block_size_distribution = top_block_size_df,
      credible_ball_lower_partition = c_lower_rl,
      credible_ball_upper_partition = c_upper_rl,
      credible_ball_horiz_partition = c_horiz_rl,
      K_VI_lower = K_VI_lower,
      K_VI_upper = K_VI_upper,
      K_VI_horiz = K_VI_horiz
    )
  }
}
setwd('/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM-Bradley-Terry-Stochastic-Block-Model/')
tennis_years <- readRDS("./data/ATP_2000_2025_SN_extended.rds")
# Load full MCMC results across seasons
res_list <- readRDS("raw_output_ext/MCMC_raw_output_ext.rds")

first_year <- 1999

# Players to track over time (must match players_df naming)
pl_selected <- c("Nadal R.", "Federer R.", "Djokovic N.", "Murray A.", "Alcaraz C.", "Sinner J.")

# Initialize containers (use lists + bind_rows at end for speed)
top_block_counts_list <- vector("list", length(res_list))
prob_assignment_list  <- vector("list", length(res_list))
post_numb_block_list  <- vector("list", length(res_list))
avg_strength_list     <- vector("list", length(res_list))
entropy_list          <- vector("list", length(res_list))

# Point-estimate TOP-block size per season (for `num_block_plot` bar chart)
top_block_size_point_list <- vector("list", length(res_list))

# For LaTeX table: store top-block player names for 2012–2020
top_block_players_2012_2020 <- list()

# For player trajectories: store posterior mean + 95% CI of P(top block)
ptop_selected_across_years <- data.frame(
  season = character(),
  player = character(),
  p_mean = numeric(),
  p_low = numeric(),
  p_high = numeric(),
  stringsAsFactors = FALSE
)

# -------------------------------
# Main Loop Over Seasons
# -------------------------------
for (yr in seq_along(res_list)) {
  season_label <- paste0(first_year + yr, "/", first_year + yr + 1)
  season_start <- first_year + yr
  res_i <- res_list[[yr]]
  x_samples <- res_i$x_samples
  lambda_samples <- res_i$lambda_samples

  inf_i <- relabel_by_lambda(x_samples, lambda_samples)
  x_relabeled <- inf_i$x_samples_relabel
  lambdas_reordered <- inf_i$lambda_samples_relabel
  T_iter <- nrow(lambdas_reordered)
  n_players <- ncol(lambdas_reordered)

  # Entropy per iteration (vectorized): p1 = fraction in top block
  p1_vec <- rowMeans(x_relabeled == 1L, na.rm = TRUE)
  entropy_container <- rep(NA_real_, T_iter)
  ok <- is.finite(p1_vec)
  deg <- ok & (p1_vec <= 0 | p1_vec >= 1)
  mid <- ok & !deg
  entropy_container[deg] <- 0
  entropy_container[mid] <- -(
    p1_vec[mid] * log(p1_vec[mid]) + (1 - p1_vec[mid]) * log(1 - p1_vec[mid])
  ) / log(2)

  # Player-level skill computation (same output, simpler mapping)
  pl_lambda <- matrix(0, nrow = T_iter, ncol = n_players)
  for (i in seq_len(T_iter)) {
    lambda_cur <- lambdas_reordered[i, ]
    u <- unique(lambda_cur)
    denom <- exp(mean(log(u)))
    u_norm <- u / denom
    pl_lambda[i, ] <- u_norm[match(lambda_cur, u)]
  }

  HPD_entropy <- HPDinterval(as.mcmc(entropy_container))
  entropy_list[[yr]] <- data.frame(
    season = season_label,
    ci_low = HPD_entropy[1],
    ci_high = HPD_entropy[2],
    mean_entropy = mean(entropy_container, na.rm = TRUE)
  )

  avg_strength_list[[yr]] <- data.frame(
    season = rep(season_label, n_players),
    mean_str = apply(pl_lambda, 2, median, na.rm = TRUE),
    lower_quantile = apply(pl_lambda, 2, quantile, probs = 0.025, na.rm = TRUE),
    upper_quantile = apply(pl_lambda, 2, quantile, probs = 0.975, na.rm = TRUE)
  )

  top_block_counts_list[[yr]] <- data.frame(
    season = season_label,
    avg_top_block_cnt = inf_i$avg_top_block_count
  )

  w_ij <- tennis_years[[yr]]$Y_ij #pairwise success matrix
  pl_df <- tennis_years[[yr]]$players_df #info about players

  # Store P(top block) trajectory for selected players (with 95% Beta credible interval)
  nm_track <- pl_df$player_label
  if (is.null(nm_track) || length(nm_track) != ncol(x_relabeled)) nm_track <- pl_df$player_slug
  if (is.null(nm_track) || length(nm_track) != ncol(x_relabeled)) nm_track <- rownames(w_ij)
  if (!is.null(nm_track) && length(nm_track) == ncol(x_relabeled)) {
    nm_track_clean <- gsub("_", " ", as.character(nm_track))
    pl_selected_clean <- gsub("_", " ", as.character(pl_selected))
    idx_sel <- which(tolower(nm_track_clean) %in% tolower(pl_selected_clean))
    if (length(idx_sel) > 0) {
      for (j in idx_sel) {
        s <- sum(x_relabeled[, j] == 1L, na.rm = TRUE)
        Tj <- sum(!is.na(x_relabeled[, j]))
        if (Tj > 0) {
          a <- 1 + s
          b <- 1 + (Tj - s)
          p_mean <- s / Tj
          p_low <- stats::qbeta(0.025, a, b)
          p_high <- stats::qbeta(0.975, a, b)
          ptop_selected_across_years <- rbind(
            ptop_selected_across_years,
            data.frame(
              season = season_label,
              player = as.character(nm_track_clean[j]),
              p_mean = p_mean,
              p_low = p_low,
              p_high = p_high,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  # ------------------------------------------------------------
  # Print players in the strongest block (label 1) for 2012–2020
  # ------------------------------------------------------------
  # Point estimate partition (used for barplot + optional 2012–2020 listing)
  x_hat <- NULL
  if (!is.null(inf_i$minVI_partition)) {
    x_hat <- inf_i$minVI_partition
  } else if (!is.null(inf_i$partition_binder)) {
    x_hat <- inf_i$partition_binder
  } else if (!is.null(inf_i$x_samples_relabel)) {
    x_hat <- apply(inf_i$x_samples_relabel, 2, function(v) {
      tab <- table(v)
      as.integer(names(tab)[which.max(tab)])
    })
  }

  top_block_size_hat <- NA_real_
  if (!is.null(x_hat)) {
    top_block_size_hat <- sum(as.integer(x_hat) == 1L, na.rm = TRUE)
  }

  if (season_start >= 2012 && season_start <= 2020) {
    if (!is.null(x_hat)) {
      top_idx <- which(as.integer(x_hat) == 1L)
      if (length(top_idx) == 0L) {
        message("Top block is empty for season ", season_label)
      } else {
        nm <- pl_df$player_label
        if (is.null(nm) || length(nm) != ncol(x_samples)) nm <- pl_df$player_slug
        if (is.null(nm) || length(nm) != ncol(x_samples)) nm <- rownames(w_ij)
        if (is.null(nm) || length(nm) != ncol(x_samples)) nm <- paste0("Player_", seq_len(ncol(x_samples)))

        ord_top <- order(pl_df$last_rank[top_idx], na.last = TRUE)
        top_names <- nm[top_idx][ord_top]

        # store for LaTeX table
        top_block_players_2012_2020[[season_label]] <- as.character(top_names)

        cat("\n", season_label, " — Top block (n=", length(top_idx), "):\n", sep = "")
        cat(paste0("  ", top_names, collapse = "\n"), "\n", sep = "")
      }
    } else {
      message("Cannot determine x_hat partition for season ", season_label)
    }
  }
  
  block_assignment_i <- inf_i$item_cluster_assignment_probs
  block_assignment_i$season <- season_label
  block_assignment_i$pl_name  <- pl_df$player_label
  block_assignment_i$eos_ranking <- pl_df$last_rank

  # Add a season-level summary for colouring: modal number of blocks (posterior mode)
  modal_K <- NA_integer_
  if (!is.null(inf_i$block_count_distribution) && nrow(inf_i$block_count_distribution) > 0) {
    modal_K <- inf_i$block_count_distribution$num_blocks[which.max(inf_i$block_count_distribution$prob)]
  } else if (!is.null(inf_i$n_clusters_each_iter)) {
    modal_K <- as.integer(names(which.max(table(inf_i$n_clusters_each_iter))))
  }
  block_assignment_i$num_blocks <- modal_K

  prob_assignment_list[[yr]] <- block_assignment_i

  post_numb_blocks_i <- inf_i$block_count_distribution
  post_numb_blocks_i$season <- season_label
  post_numb_block_list[[yr]] <- post_numb_blocks_i

  top_block_size_point_list[[yr]] <- data.frame(
    season = season_label,
    top_block_size_hat = top_block_size_hat,
    num_blocks = modal_K,
    stringsAsFactors = FALSE
  )

  message("Processed season: ", season_label)
}

# Bind outputs (fast)
top_block_counts_across_years <- dplyr::bind_rows(top_block_counts_list)
prob_assignment_across_years  <- dplyr::bind_rows(prob_assignment_list)
post_numb_block_across_years  <- dplyr::bind_rows(post_numb_block_list)
avg_strength_each_player      <- dplyr::bind_rows(avg_strength_list)
entropy_per_season            <- dplyr::bind_rows(entropy_list)
top_block_size_point_across_years <- dplyr::bind_rows(top_block_size_point_list)

# ------------------------------------------------------------------------
# LaTeX table: Top-block players by season (2012–2020)
# Columns = seasons; rows = player names (padded with blanks)
# ------------------------------------------------------------------------
season_labels_target <- paste0(2012:2020, "/", 2013:2021)
season_labels_present <- season_labels_target[season_labels_target %in% names(top_block_players_2012_2020)]

if (length(season_labels_present) > 0) {
  max_len <- max(vapply(season_labels_present, function(s) length(top_block_players_2012_2020[[s]]), integer(1)))
  top_block_table <- as.data.frame(lapply(season_labels_present, function(s) {
    v <- top_block_players_2012_2020[[s]]
    c(v, rep("", max_len - length(v)))
  }), stringsAsFactors = FALSE)
  names(top_block_table) <- season_labels_present

  if (!dir.exists("./tables")) dir.create("./tables", recursive = TRUE)

  top_block_tex <- knitr::kable(
    top_block_table,
    format = "latex",
    booktabs = TRUE,
    escape = TRUE,
    caption = "Players in inferred top block (cluster 1) by season (minVI partition)"
  ) %>%
    kableExtra::kable_styling(latex_options = c("hold_position", "striped"), full_width = FALSE)

  writeLines(top_block_tex, "./tables/top_block_players_2012_2020.tex")
} else {
  message("No top-block players stored for 2012–2020; table not written.")
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
# (palette now extends to any number of blocks)
# ------------------------------------------------------------------------
prob_assignment_across_years$season_start <- season_start_year(prob_assignment_across_years$season)
num_blocks_vals <- sort(unique(stats::na.omit(prob_assignment_across_years$num_blocks)))
col_map <- tennis_surface_palette(num_blocks_vals)

p_top_across_time = ggplot(prob_assignment_across_years, aes(x = factor(season_start), y = Cluster_1)) +
  geom_jitter(aes(color = factor(num_blocks)), width = 0.2, alpha = 0.6, size = 3) +
  labs(
    x     = "Season",
    y     = "P(Top Block)",
    color = "Nº of blocks"
  ) +
  scale_color_manual(values = col_map) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------
# Entropy plot (needed for saving)
# ------------------------------------------------------------------------
entropy_per_season$season_start <- season_start_year(entropy_per_season$season)
entropy_plot <- ggplot(entropy_per_season, aes(x = factor(season_start))) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, group = 1),
              fill = "grey60", alpha = 0.35, inherit.aes = TRUE) +
  geom_line(aes(y = mean_entropy, group = 1), color = "grey10", linewidth = 0.7) +
  geom_point(aes(y = mean_entropy), color = "grey10", size = 1.6) +
  labs(x = "Season", y = "Entropy") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )
# ------------------------------------------------------------------------
# Bar plot: Estimated number of players in the top block by season
# (Barplot: point estimate partition size per season)
# ------------------------------------------------------------------------
top_block_size_point_across_years$season <- factor(
  top_block_size_point_across_years$season,
  levels = unique(prob_assignment_across_years$season)  # keep same ordering as other plots
)
top_block_size_point_across_years$season_start <- season_start_year(top_block_size_point_across_years$season)

fill_map <- tennis_surface_palette(stats::na.omit(top_block_size_point_across_years$num_blocks))

num_block_plot <- ggplot(
  top_block_size_point_across_years,
  aes(x = factor(season_start), y = top_block_size_hat, fill = factor(num_blocks))
) +
  geom_col(alpha = 0.9, color = "grey10", linewidth = 0.35, width = 0.85) +
  scale_fill_manual(values = fill_map, na.value = "grey75", name = "Nº of blocks") +
  labs(
    x = "Season",
    y = "Nº of Players in Top Block"
  ) +
  coord_cartesian(ylim = c(0, 25)) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------
# Extra: Player trajectory across the seasons
# Compare P(top block) vs (100 - end-of-season ranking), for example
# ------------------------------------------------------------------------
ptop_selected_plot <- ptop_selected_across_years %>%
  {
    if (nrow(.) == 0) {
      ggplot() +
        theme_void() +
        labs(title = "No selected players found in this dataset")
    } else {
      . %>%
        mutate(season_start = season_start_year(season)) %>%
        mutate(season_start = factor(season_start, levels = sort(unique(season_start)))) %>%
        mutate(
          p_mean = 100 * p_mean,
          p_low = 100 * p_low,
          p_high = 100 * p_high
        ) %>%
        ggplot(aes(x = season_start, group = player)) +
        geom_ribbon(aes(ymin = p_low, ymax = p_high, fill = player), alpha = 0.18, colour = NA) +
        geom_line(aes(y = p_mean, colour = player), linewidth = 1) +
        geom_point(aes(y = p_mean, colour = player), size = 1.7) +
        theme_minimal() +
        labs(
          y = "P(Top Block) (%)",
          x = "Season",
          colour = "Player",
          fill = "Player"
        ) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.minor = element_blank()
        ) +
        facet_wrap(~player)
    }
  }

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

ggsave(filename = "./images/Ptop_selected_players.pdf",
  plot = ptop_selected_plot,
  width = 11, height = 6)


