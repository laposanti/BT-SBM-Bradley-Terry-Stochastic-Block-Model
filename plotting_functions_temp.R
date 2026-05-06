block_winprob_table_from_lambdas = function(
    lambda_item,
    z,
    lambda_stat = c("mean", "median"),
    N_ij = NULL,
    include_diag = FALSE,
    diag_value = NA_real_
) {
  if (is.null(z)) stop("`z` must be provided (block memberships).")
  z <- as.integer(z)
  
  if (is.list(lambda_item)) lambda_item <- do.call(rbind, lambda_item)
  if (is.matrix(lambda_item)) {
    lambda_stat <- match.arg(lambda_stat)
    if (lambda_stat == "mean") {
      lambda_point <- colMeans(lambda_item, na.rm = TRUE)
    } else {
      lambda_point <- apply(lambda_item, 2, stats::median, na.rm = TRUE)
    }
  } else {
    lambda_point <- as.numeric(lambda_item)
  }
  
  n <- length(lambda_point)
  if (length(z) != n) stop("`z` must have the same length as `lambda_item` (players).")
  if (any(!is.finite(lambda_point))) stop("`lambda_item` contains non-finite values.")
  if (any(lambda_point < 0)) stop("`lambda_item` must be non-negative.")
  
  if (!is.null(N_ij)) {
    if (!is.matrix(N_ij) || nrow(N_ij) != n || ncol(N_ij) != n) {
      stop("`N_ij` must be an (n_players x n_players) matrix aligned with `lambda_item`.")
    }
  }
  
  block_labels <- sort(unique(z))
  K <- length(block_labels)
  zK <- match(z, block_labels)
  block_sizes <- tabulate(zK, nbins = K)
  
  P_block <- matrix(diag_value, nrow = K, ncol = K)
  N_pairs <- matrix(0, nrow = K, ncol = K)
  dimnames(P_block) <- list(paste0("Block_", block_labels), paste0("Block_", block_labels))
  dimnames(N_pairs) <- dimnames(P_block)
  
  for (i in seq_len(K)) {
    idx_i <- which(zK == i)
    li <- lambda_point[idx_i]
    
    for (j in seq_len(K)) {
      if (i == j && !isTRUE(include_diag)) next
      
      idx_j <- which(zK == j)
      lj <- lambda_point[idx_j]
      
      # Probabilities for all ordered pairs (a in i, b in j)
      denom <- outer(li, lj, `+`)
      num <- matrix(li, nrow = length(li), ncol = length(lj))
      Pij <- num / denom
      
      if (i == j) {
        diag(Pij) <- NA_real_  # exclude self-pairs
      }
      
      if (is.null(N_ij)) {
        P_block[i, j] <- mean(Pij, na.rm = TRUE)
        N_pairs[i, j] <- sum(is.finite(Pij))
      } else {
        W <- N_ij[idx_i, idx_j, drop = FALSE]
        if (i == j) diag(W) <- 0
        wsum <- sum(W, na.rm = TRUE)
        if (wsum <= 0) {
          P_block[i, j] <- NA_real_
          N_pairs[i, j] <- 0
        } else {
          P_block[i, j] <- sum(Pij * W, na.rm = TRUE) / wsum
          N_pairs[i, j] <- wsum
        }
      }
    }
  }
  
  list(
    P_block = P_block,
    N_pairs = N_pairs,
    lambda_point = lambda_point,
    block_sizes = stats::setNames(block_sizes, paste0("Block_", block_labels)),
    block_labels = block_labels
  )
}
#' Adjacency heatmap (wins per matches)
#'
#' Exploratory adjacency heatmap showing wins for each ordered pair of players.
#' Optionally restrict to a subset of players and/or order rows/cols by an
#' external metric such as year-end ranking.
#'
#' @param Y_ij Integer matrix of wins. Entry (i, j) is wins of i over j.
#' @param N_ij Integer matrix of matches. If `NULL`, computed as `Y_ij + t(Y_ij)`.
#' @param names Optional character vector of player names to include. Can be
#'   either raw matrix dimnames (slugs) or display names (after cleaning).
#'   If `NULL`, plots all players.
#' @param last_rank Optional numeric vector of year-end ranks. Can be named
#'   (names are player identifiers) or an unnamed vector aligned with the players
#'   being plotted.
#' @param order_by Ordering criterion. `"last_rank"` orders by ascending rank
#'   (rank 1 at the top). `"none"` keeps matrix order. `"auto"` uses
#'   `"last_rank"` if ranking info is provided, otherwise `"none"`.
#' @param players_df Optional data.frame with columns `player_slug` and `last_rank`
#'   (e.g., season-level metadata). Used when `order_by = "last_rank"`.
#' @param clean_fun Name-cleaning function applied to display labels.
#' @param max_wins Maximum win count shown explicitly in the legend; larger
#'   counts are bucketed into `paste0(max_wins, "+")`.
#' @param show_margins If `TRUE` and package `ggside` is available, draws a
#'   side bar with total matches per player.
#' @return A `ggplot` object.
#'
# Figure 1 plotting function ---------
exploratory_adjacency <- function(
    Y_ij,
    N_ij = NULL,
    names = NULL,
    last_rank = NULL,
    order_by = c("auto", "none", "last_rank"),
    players_df = NULL,
    players_id_col = NULL,
    players_rank_col = NULL,
    labels = NULL,
    labels_key_col = NULL,
    labels_value_col = NULL,
    clean_display_labels = NULL,
    clean_fun = clean_players_names,
    max_wins = 4L,
    show_margins = TRUE,
    legend_title = "Number of wins",
    no_match_label = "No match played",
    no_match_color = "white",         # <- important for B/W
    mark_unplayed = TRUE,             # <- redundant cue
    unplayed_mark = "\u00D7",         # "×"
    unplayed_mark_size = 2.2,
    bw_preview = TRUE                 # <- prints an extra B/W check plot
){
  order_by <- match.arg(order_by)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  if (!is.matrix(Y_ij)) stop("`Y_ij` must be a matrix.")
  if (is.null(N_ij)) N_ij <- Y_ij + t(Y_ij)
  if (!is.matrix(N_ij)) stop("`N_ij` must be a matrix (or NULL).")
  if (nrow(Y_ij) != nrow(N_ij) || ncol(Y_ij) != ncol(N_ij)) {
    stop("`Y_ij` and `N_ij` must have the same dimensions.")
  }
  
  # Ensure dimnames exist / are aligned
  if (is.null(rownames(Y_ij)) && !is.null(colnames(Y_ij))) rownames(Y_ij) <- colnames(Y_ij)
  if (is.null(colnames(Y_ij)) && !is.null(rownames(Y_ij))) colnames(Y_ij) <- rownames(Y_ij)
  if (is.null(rownames(N_ij)) && !is.null(colnames(N_ij))) rownames(N_ij) <- colnames(N_ij)
  if (is.null(colnames(N_ij)) && !is.null(rownames(N_ij))) colnames(N_ij) <- rownames(N_ij)
  
  if (is.null(rownames(Y_ij)) || is.null(colnames(Y_ij))) {
    ids <- paste0("Item_", seq_len(nrow(Y_ij)))
    rownames(Y_ij) <- ids; colnames(Y_ij) <- ids
    rownames(N_ij) <- ids; colnames(N_ij) <- ids
  }
  
  to_title_safe <- function(x) {
    x <- gsub("[-_]", " ", x)
    if (requireNamespace("stringr", quietly = TRUE)) return(stringr::str_to_title(x))
    tools::toTitleCase(tolower(x))
  }
  clean_name_vec <- function(x) {
    x2 <- to_title_safe(x)
    if (!is.null(clean_fun)) x2 <- vapply(x2, clean_fun, character(1))
    x2
  }
  
  ids <- rownames(Y_ij)
  id_to_disp <- stats::setNames(clean_name_vec(ids), ids)
  
  # Optional: override axis labels via a user mapping or via players_df
  infer_labels_from_players_df <- function(df, ids, key_col = NULL, value_col = NULL) {
    if (is.null(df) || !is.data.frame(df) || length(ids) < 1) return(NULL)
    
    if (is.null(key_col)) {
      candidates <- c("player_slug", "player_id", "player", "player_label")
      candidates <- candidates[candidates %in% colnames(df)]
      if (length(candidates) == 0) return(NULL)
      overlaps <- vapply(candidates, function(cc) sum(as.character(df[[cc]]) %in% ids, na.rm = TRUE), numeric(1))
      key_col <- candidates[which.max(overlaps)]
      if (is.na(key_col) || overlaps[which.max(overlaps)] == 0) return(NULL)
    }
    
    if (is.null(value_col)) {
      candidates <- c("player_label", "label", "player", key_col)
      value_col <- candidates[candidates %in% colnames(df)][1]
      if (is.na(value_col) || is.null(value_col)) return(NULL)
    }
    
    key <- as.character(df[[key_col]])
    val <- as.character(df[[value_col]])
    ok <- !is.na(key) & !is.na(val) & nzchar(key)
    if (!any(ok)) return(NULL)
    
    map <- tapply(val[ok], key[ok], function(v) v[which(!is.na(v) & nzchar(v))[1]])
    map <- map[!is.na(map) & nzchar(map)]
    if (length(map) == 0) return(NULL)
    
    list(map = map, key_col = key_col, value_col = value_col)
  }
  
  labels_map <- NULL
  labels_from <- NULL
  
  if (!is.null(labels)) {
    if (is.character(labels) && !is.null(names(labels))) {
      labels_map <- labels; labels_from <- "labels"
    } else if (is.character(labels) && length(labels) == length(ids)) {
      labels_map <- stats::setNames(labels, ids); labels_from <- "labels"
    } else {
      stop("`labels` must be a named character vector (names are matrix ids), or a character vector aligned with the plotted ids.")
    }
  } else if (!is.null(players_df)) {
    inferred <- infer_labels_from_players_df(players_df, ids, key_col = labels_key_col, value_col = labels_value_col)
    if (!is.null(inferred)) {
      labels_map <- inferred$map
      labels_from <- "players_df"
      labels_key_col <- inferred$key_col
      labels_value_col <- inferred$value_col
    }
  }
  
  if (!is.null(labels_map)) {
    if (is.null(clean_display_labels)) {
      clean_display_labels <- !(isTRUE(labels_from == "players_df") &&
                                  !is.null(labels_value_col) &&
                                  grepl("label", labels_value_col, ignore.case = TRUE))
    }
    labels_use <- labels_map
    if (isTRUE(clean_display_labels)) {
      labels_use <- clean_name_vec(labels_use)
      names(labels_use) <- names(labels_map)
    }
    overlap <- intersect(ids, names(labels_use))
    id_to_disp[overlap] <- unname(labels_use[overlap])
  }
  
  # Subset by `names` if requested
  if (!is.null(names)) {
    if (!is.character(names) || length(names) < 1) stop("`names` must be a non-empty character vector (or NULL).")
    names_clean <- clean_name_vec(names)
    keep <- (id_to_disp %in% names_clean) | (ids %in% names)
    if (!any(keep)) stop("None of the requested `names` match the matrix dimnames (after cleaning).")
    Y_ij <- Y_ij[keep, keep, drop = FALSE]
    N_ij <- N_ij[keep, keep, drop = FALSE]
    ids <- rownames(Y_ij)
    id_to_disp <- id_to_disp[ids]
  }
  
  # Resolve ranking information if requested
  ranks_by_disp <- NULL
  ranks_by_id <- NULL
  
  if (!is.null(players_df)) {
    if (!is.data.frame(players_df)) stop("`players_df` must be a data.frame (or NULL).")
    
    if (is.null(players_id_col)) {
      candidates <- c("player_slug", "player_id", "player", "player_label")
      players_id_col <- candidates[candidates %in% colnames(players_df)][1]
    }
    if (is.null(players_rank_col)) {
      candidates <- c("last_rank", "rank", "year_end_rank", "yearend_rank")
      players_rank_col <- candidates[candidates %in% colnames(players_df)][1]
    }
    if (is.na(players_id_col) || is.null(players_id_col)) stop("Could not infer `players_id_col` from `players_df`.")
    if (is.na(players_rank_col) || is.null(players_rank_col)) stop("Could not infer `players_rank_col` from `players_df`.")
    
    pl_key_raw <- as.character(players_df[[players_id_col]])
    pl_disp <- clean_name_vec(pl_key_raw)
    pl_rank <- as.numeric(players_df[[players_rank_col]])
    ok <- is.finite(pl_rank)
    
    if (any(ok)) {
      ranks_by_disp <- tapply(pl_rank[ok], pl_disp[ok], function(v) min(v, na.rm = TRUE))
      ranks_by_disp <- as.numeric(ranks_by_disp)
      names(ranks_by_disp) <- names(tapply(pl_rank[ok], pl_disp[ok], function(v) min(v, na.rm = TRUE)))
      
      if (any(pl_key_raw %in% ids, na.rm = TRUE)) {
        ranks_by_id <- tapply(pl_rank[ok], pl_key_raw[ok], function(v) min(v, na.rm = TRUE))
        ranks_by_id <- as.numeric(ranks_by_id)
        names(ranks_by_id) <- names(tapply(pl_rank[ok], pl_key_raw[ok], function(v) min(v, na.rm = TRUE)))
      }
    }
  }
  
  if (!is.null(last_rank)) {
    if (!is.numeric(last_rank)) stop("`last_rank` must be numeric (or NULL).")
    if (!is.null(names(last_rank))) {
      rk_names_disp <- clean_name_vec(names(last_rank))
      ranks_by_disp <- tapply(as.numeric(last_rank), rk_names_disp, function(v) min(v, na.rm = TRUE))
      if (any(names(last_rank) %in% ids, na.rm = TRUE)) {
        ranks_by_id <- tapply(as.numeric(last_rank), names(last_rank), function(v) min(v, na.rm = TRUE))
      }
    } else if (length(last_rank) == length(ids)) {
      ranks_by_id <- stats::setNames(as.numeric(last_rank), ids)
    } else {
      stop("`last_rank` must be either a named vector, or have length equal to the number of players plotted.")
    }
  }
  
  if (order_by == "auto") {
    order_by <- if (!is.null(ranks_by_id) || !is.null(ranks_by_disp)) "last_rank" else "none"
  }
  if (order_by == "last_rank" && is.null(ranks_by_id) && is.null(ranks_by_disp)) {
    stop("`order_by = 'last_rank'` requires `last_rank` or `players_df`.")
  }
  
  ids_order <- ids
  if (order_by == "last_rank") {
    rk <- NULL
    if (!is.null(ranks_by_id)) rk <- ranks_by_id[ids]
    if (is.null(rk) || all(is.na(rk))) {
      disp_now <- unname(id_to_disp[ids])
      rk <- ranks_by_disp[disp_now]
    }
    ranked <- !is.na(rk)
    disp_now <- unname(id_to_disp[ids])
    ord <- order(rk[ranked], disp_now[ranked], ids[ranked])
    ids_order <- c(ids[ranked][ord], ids[!ranked])
  }
  
  # Long format
  Y_long <- as.data.frame(as.table(Y_ij), stringsAsFactors = FALSE)
  colnames(Y_long) <- c("Winner_id", "Loser_id", "Count")
  N_long <- as.data.frame(as.table(N_ij), stringsAsFactors = FALSE)
  colnames(N_long) <- c("Winner_id", "Loser_id", "N_matches")
  
  Y_long <- Y_long |>
    dplyr::left_join(N_long, by = c("Winner_id", "Loser_id")) |>
    dplyr::mutate(
      Winner = unname(id_to_disp[Winner_id]),
      Loser  = unname(id_to_disp[Loser_id]),
      Played = !is.na(N_matches) & N_matches > 0
    )
  
  max_wins <- as.integer(max_wins)
  if (is.na(max_wins) || max_wins < 0) stop("`max_wins` must be a non-negative integer.")
  win_levels <- as.character(0:max_wins)
  
  Y_long <- Y_long |>
    dplyr::mutate(
      Count_capped = as.character(pmin(Count, max_wins)),
      Count_disp   = ifelse(!Played, no_match_label, Count_capped)
    ) |>
    dplyr::mutate(
      Winner_id = factor(Winner_id, levels = ids_order),
      Loser_id  = factor(Loser_id,  levels = ids_order),
      Count_disp = factor(Count_disp, levels = c(no_match_label, win_levels), ordered = TRUE)
    )
  
  marginal_df <- Y_long |>
    dplyr::group_by(Winner_id) |>
    dplyr::summarise(Total_matches = sum(N_matches, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(Winner_id = factor(Winner_id, levels = levels(Y_long$Winner_id)))
  
  # Colour-blind friendly WITH GREEN (viridis "D")
  # Key B/W fix: do NOT allow the lightest played colour to be too close to white.
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    ramp <- viridisLite::viridis(length(win_levels), option = "D", begin = 0.25, end = 0.95)
  } else {
    ramp <- grDevices::colorRampPalette(c("#2c7fb8", "#41ab5d", "#fddc6c"))(length(win_levels))
  }
  names(ramp) <- win_levels
  
  fill_values <- c(setNames(no_match_color, no_match_label), ramp)
  
  x_lab <- function(x) unname(id_to_disp[x])
  
  base_plot <- ggplot2::ggplot(Y_long, ggplot2::aes(x = Loser_id, y = Winner_id)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = Count_disp),
      colour = "grey70", linewidth = 0.25
    ) +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_y_discrete(
      limits = rev(levels(Y_long$Winner_id)),
      labels = x_lab,
      guide  = ggplot2::guide_axis(n.dodge = 2)
    ) +
    ggplot2::scale_x_discrete(
      labels = x_lab,
      guide  = ggplot2::guide_axis(n.dodge = 2)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, fill = legend_title) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key.width  = grid::unit(1.2, "cm"),
      legend.key.height = grid::unit(0.45, "cm")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        title.position = "top",
        direction = "horizontal",
        nrow = 1,
        byrow = TRUE
      )
    )
  
  # Redundant cue for unplayed: add a small "×" overlay (works in B/W and photocopies)
  if (isTRUE(mark_unplayed)) {
    unplayed_df <- Y_long |> dplyr::filter(!Played)
    base_plot <- base_plot +
      ggplot2::geom_text(
        data = unplayed_df,
        ggplot2::aes(label = unplayed_mark),
        colour = "grey55",
        size = unplayed_mark_size / ggplot2::.pt,
        show.legend = FALSE
      )
  }
  
  if (isTRUE(show_margins)) {
    if (requireNamespace("ggside", quietly = TRUE)) {
      base_plot <- base_plot +
        ggside::geom_ysidecol(
          data = marginal_df,
          ggplot2::aes(y = Winner_id, x = Total_matches),
          fill = "grey40",
          width = 0.8,
          show.legend = FALSE
        ) +
        ggside::scale_ycolor_continuous(name = "Total matches")
    } else {
      warning("Package 'ggside' not installed: drawing heatmap without side marginal totals.", call. = FALSE)
    }
  }
  
  print(base_plot)
  
  # ---- B/W preview: convert your chosen colours to greys by luminance ----
  # This is an explicit "will this survive printing?" check.
  if (isTRUE(bw_preview)) {
    hex_to_grey <- function(hex) {
      rgb <- grDevices::col2rgb(hex) / 255
      # sRGB -> linear
      lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
      # relative luminance (WCAG)
      L <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
      grDevices::rgb(L, L, L)
    }
    fill_bw <- vapply(fill_values, hex_to_grey, character(1))
    
    bw_plot <- base_plot +
      ggplot2::scale_fill_manual(values = fill_bw, drop = FALSE) +
      ggplot2::labs(fill = paste0(legend_title, " (B/W preview)"))
    
    print(bw_plot)
    
    return(list(color_plot = base_plot, bw_plot = bw_plot))
  }
  
  return(base_plot)
}








#' Figure 3 plotting function
#' Block-ordered adjacency heatmap (BT-SBM)
#'
#' Plot the adjacency (wins / matches) heatmap with rows/cols ordered by
#' the hard partition (\code{fit$estimates$x_hat}) and, within blocks,
#' by marginal wins. Draws block boundary lines and (optionally) a side color strip
#' for block IDs if \pkg{ggside} is installed.
#' Optionally relabels blocks so that Block 1 has the largest mean player
#' strength \eqn{\lambda} (mean computed across players in the block).
#'
#' @param fit Output list from \code{gibbs_BT_SBM()}.
#' @param w_ij Integer matrix of wins (same players & order used in \code{fit}).
#' @param x_hat partition point estimate. One n-length integer vector.
#' @param clean_fun Optional function to prettify player names. Default: identity.
#' @param palette Named colors for blocks (as character vector). Defaults to a
#'   Wimbledon-ish palette.
#' @param fill_low,fill_high Colors for the heatmap pixels gradient low/high.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' p <- plot_block_adjacency(fit, w_ij)
#' }
#' @export
#'
plot_block_adjacency <- function(
    fit,
    w_ij,
    x_hat = NULL,
  relabel_blocks = c("none", "avg_lambda"),
  lambda_hat = NULL,
  relabel_decreasing = TRUE,
    clean_fun = clean_players_names,
    players_df = NULL,
    players_id_col = NULL,
    labels = NULL,
    labels_key_col = NULL,
    labels_value_col = NULL,
    clean_display_labels = NULL,
    order_ids = NULL,
    palette = NULL,
    fill_low = "#FFFFCC",
    fill_high = "#006400",
    no_match_label = "No match played",
    no_match_color = "white",      # <- crucial for B/W printing
    mark_unplayed = TRUE,          # <- redundant cue, no colour needed
    unplayed_mark = "\u00D7",      # "×"
    unplayed_mark_size = 2.2,
    bw_preview = TRUE              # <- print a B/W check plot too
){
  stopifnot(is.matrix(w_ij))
  relabel_blocks <- match.arg(relabel_blocks)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.")
  if (!requireNamespace("ggside", quietly = TRUE)) stop("Package 'ggside' is required.")
  if (!requireNamespace("reshape2", quietly = TRUE)) stop("Package 'reshape2' is required.")
  
  N_ij <- w_ij + t(w_ij)
  if (is.null(x_hat)) x_hat <- fit$minVI_partition
  stopifnot(length(x_hat) == nrow(w_ij))
  
  # Input names
  pl_ids <- rownames(w_ij)
  if (is.null(pl_ids)) pl_ids <- paste0("Item_", seq_len(nrow(w_ij)))
  colnames(w_ij) <- pl_ids
  rownames(w_ij) <- pl_ids
  colnames(N_ij) <- pl_ids
  rownames(N_ij) <- pl_ids

  # ---- Optional: relabel blocks by mean player lambda ----
  # This does NOT change within-block ordering (still based on marginal wins);
  # it only renumbers blocks so that 1 is the strongest block on average.
  if (identical(relabel_blocks, "avg_lambda")) {
    lambda_item_mean <- NULL

    # Prefer explicit lambda_hat if provided.
    if (!is.null(lambda_hat)) {
      if (is.list(lambda_hat)) lambda_hat <- do.call(cbind, lambda_hat)
      if (is.matrix(lambda_hat) || is.data.frame(lambda_hat)) {
        lambda_item_mean <- colMeans(as.matrix(lambda_hat), na.rm = TRUE)
        if (!is.null(colnames(lambda_hat))) names(lambda_item_mean) <- colnames(lambda_hat)
      } else {
        lambda_item_mean <- as.numeric(lambda_hat)
        if (!is.null(names(lambda_hat))) names(lambda_item_mean) <- names(lambda_hat)
      }
    } else if (!is.null(fit$lambda_samples_relabel)) {
      lam <- fit$lambda_samples_relabel
      if (is.list(lam)) lam <- do.call(rbind, lam)
      lam <- as.matrix(lam)
      lambda_item_mean <- colMeans(lam, na.rm = TRUE)
      if (!is.null(colnames(lam))) names(lambda_item_mean) <- colnames(lam)
    } else if (!is.null(fit$lambda_hat)) {
      lambda_item_mean <- as.numeric(fit$lambda_hat)
      if (!is.null(names(fit$lambda_hat))) names(lambda_item_mean) <- names(fit$lambda_hat)
    } else if (!is.null(fit$lambda_samples)) {
      lam <- fit$lambda_samples
      if (is.list(lam)) lam <- do.call(rbind, lam)
      lam <- as.matrix(lam)
      lambda_item_mean <- colMeans(lam, na.rm = TRUE)
      if (!is.null(colnames(lam))) names(lambda_item_mean) <- colnames(lam)
    }

    if (is.null(lambda_item_mean)) {
      warning("relabel_blocks='avg_lambda' requested but no player lambda found in `lambda_hat`, fit$lambda_samples_relabel, fit$lambda_hat, or fit$lambda_samples. Using original x_hat labels.", call. = FALSE)
    } else {
      # Align lambda vector to the plotted ids.
      if (!is.null(names(lambda_item_mean))) {
        lambda_item_mean <- lambda_item_mean[pl_ids]
      }
      if (length(lambda_item_mean) != length(pl_ids)) {
        stop("`lambda_hat` (or lambda in `fit`) must be length n_items (or an iters x n_items matrix).")
      }

      labs_chr <- unique(as.character(x_hat))
      labs_chr <- labs_chr[!is.na(labs_chr)]
      cl_means <- vapply(labs_chr, function(k) {
        mean(lambda_item_mean[as.character(x_hat) == k], na.rm = TRUE)
      }, numeric(1))

      ord <- order(cl_means, decreasing = isTRUE(relabel_decreasing), labs_chr)
      map <- setNames(seq_along(ord), labs_chr[ord])
      x_hat <- as.integer(map[as.character(x_hat)])
    }
  }
  
  make_palette <- function(n, palette_name = "Dark 3") {
    if (is.function(grDevices::hcl.colors)) return(grDevices::hcl.colors(n, palette = palette_name))
    grDevices::rainbow(n)
  }
  
  K <- length(unique(x_hat))
  if (is.null(palette)) palette <- make_palette(K)
  stopifnot(length(palette) == K)
  levs <- levels(factor(x_hat))
  if (!is.null(names(palette)) && all(levs %in% names(palette))) {
    palette <- palette[levs]
  }
  if (is.null(names(palette)) || !all(levs %in% names(palette))) names(palette) <- levs
  
  # --- Display labels (axis text) ---
  to_title_safe <- function(x) {
    x <- gsub("[-_]", " ", x)
    if (requireNamespace("stringr", quietly = TRUE)) return(stringr::str_to_title(x))
    tools::toTitleCase(tolower(x))
  }
  
  clean_name_vec <- function(x) {
    x2 <- to_title_safe(x)
    if (!is.null(clean_fun)) x2 <- vapply(x2, clean_fun, character(1))
    x2
  }
  
  infer_labels_from_players_df <- function(df, ids, key_col = NULL, value_col = NULL) {
    if (is.null(df) || !is.data.frame(df) || length(ids) < 1) return(NULL)
    
    if (is.null(key_col)) {
      candidates <- c("player_slug", "player_id", "player", "player_label")
      candidates <- candidates[candidates %in% colnames(df)]
      if (length(candidates) == 0) return(NULL)
      overlaps <- vapply(candidates, function(cc) sum(as.character(df[[cc]]) %in% ids, na.rm = TRUE), numeric(1))
      key_col <- candidates[which.max(overlaps)]
      if (is.na(key_col) || overlaps[which.max(overlaps)] == 0) return(NULL)
    }
    
    if (is.null(value_col)) {
      candidates <- c("player_label", "label", "player", key_col)
      value_col <- candidates[candidates %in% colnames(df)][1]
      if (is.na(value_col) || is.null(value_col)) return(NULL)
    }
    
    key <- as.character(df[[key_col]])
    val <- as.character(df[[value_col]])
    ok <- !is.na(key) & !is.na(val) & nzchar(key)
    if (!any(ok)) return(NULL)
    
    map <- tapply(val[ok], key[ok], function(v) v[which(!is.na(v) & nzchar(v))[1]])
    map <- map[!is.na(map) & nzchar(map)]
    if (length(map) == 0) return(NULL)
    list(map = map, key_col = key_col, value_col = value_col)
  }
  
  id_to_disp <- stats::setNames(clean_name_vec(pl_ids), pl_ids)
  labels_map <- NULL
  labels_from <- NULL
  
  if (!is.null(labels)) {
    if (is.character(labels) && !is.null(names(labels))) {
      labels_map <- labels; labels_from <- "labels"
    } else if (is.character(labels) && length(labels) == length(pl_ids)) {
      labels_map <- stats::setNames(labels, pl_ids); labels_from <- "labels"
    } else {
      stop("`labels` must be a named character vector (names are matrix ids), or a character vector aligned with the plotted ids.")
    }
  } else {
    inferred <- infer_labels_from_players_df(players_df, pl_ids, key_col = labels_key_col, value_col = labels_value_col)
    if (!is.null(inferred)) {
      labels_map <- inferred$map
      labels_from <- "players_df"
      labels_key_col <- inferred$key_col
      labels_value_col <- inferred$value_col
    }
  }
  
  if (!is.null(labels_map)) {
    if (is.null(clean_display_labels)) {
      clean_display_labels <- !(isTRUE(labels_from == "players_df") &&
                                  !is.null(labels_value_col) &&
                                  grepl("label", labels_value_col, ignore.case = TRUE))
    }
    labels_use <- labels_map
    if (isTRUE(clean_display_labels)) {
      labels_use <- clean_name_vec(labels_use)
      names(labels_use) <- names(labels_map)
    }
    overlap <- intersect(pl_ids, names(labels_use))
    id_to_disp[overlap] <- unname(labels_use[overlap])
  }
  
  # Add cluster info to players
  df_cl <- data.frame(
    players = pl_ids,
    cl = x_hat,
    marginal_win = rowSums(w_ij),
    marginal_matches = rowSums(N_ij)
  )
  
  df_cl <- df_cl |>
    dplyr::mutate(
      marginal_win_prob = dplyr::if_else(
        is.finite(marginal_matches) & marginal_matches > 0,
        marginal_win / marginal_matches,
        NA_real_
      )
    )
  
  # --- Common ordering across plots ---
  if (is.null(order_ids)) {
    df_ord <- df_cl |>
      dplyr::arrange(cl, dplyr::desc(marginal_win_prob), dplyr::desc(marginal_win), players)
    order_ids <- df_ord$players
  }
  
  order_ids <- as.character(order_ids)
  order_ids <- order_ids[order_ids %in% pl_ids]
  order_ids <- unique(order_ids)
  if (length(order_ids) != length(pl_ids)) {
    order_ids <- c(order_ids, setdiff(pl_ids, order_ids))
  }
  
  # Prepare long format for adjacency matrix
  Y_long <- reshape2::melt(w_ij)
  colnames(Y_long) <- c("Winner_id", "Loser_id", "Win_Count")
  Y_long$Matches_Count <- reshape2::melt(N_ij)$value
  
  Y_long_plot <- Y_long |>
    dplyr::mutate(
      # IMPORTANT: define unplayed cleanly
      Played = is.finite(Matches_Count) & Matches_Count > 0,
      perc_success = dplyr::if_else(Played, Win_Count / Matches_Count, NA_real_)
    ) |>
    dplyr::left_join(df_cl, by = c("Loser_id" = "players")) |>
    dplyr::rename(row_cl = cl, marginal_win_row = marginal_win) |>
    dplyr::left_join(df_cl, by = c("Winner_id" = "players")) |>
    dplyr::rename(col_cl = cl, marginal_win_col = marginal_win) |>
    dplyr::mutate(
      Winner = unname(id_to_disp[Winner_id]),
      Loser  = unname(id_to_disp[Loser_id])
    ) |>
    dplyr::mutate(
      Winner_id = factor(Winner_id, levels = order_ids),
      Loser_id  = factor(Loser_id,  levels = order_ids),
      col_cl = factor(col_cl, ordered = TRUE),
      row_cl = factor(row_cl, ordered = TRUE)
    )
  
  # Block boundaries aligned to order_ids
  ord_cl <- x_hat[match(order_ids, pl_ids)]
  ord_cl_chr <- as.character(ord_cl)
  ord_cl_chr[is.na(ord_cl_chr)] <- "__NA__"
  run_ends <- cumsum(rle(ord_cl_chr)$lengths)
  breaks <- run_ends[-length(run_ends)]
  v_lines_list <- breaks
  h_lines_list <- length(order_ids) - breaks
  
  x_lab <- function(x) unname(id_to_disp[x])
  
  # Heatmap fill: colour-blind friendly WITH GREEN (viridis "D")
  # Key B/W fix: NA is white and the low end is kept away from white via begin=0.25.
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    fill_scale <- ggplot2::scale_fill_gradientn(
      colours = viridisLite::viridis(256, option = "D", begin = 0.25, end = 0.95),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = no_match_color
    )
  } else {
    # fallback: keep some green but ensure low isn't too close to white
    fill_scale <- ggplot2::scale_fill_gradientn(
      colours = c("#2b8cbe", "#41ab5d", "#fddc6c"),
      values = scales::rescale(c(0, 0.5, 1)),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = no_match_color
    )
  }
  
  p <- ggplot2::ggplot(Y_long_plot, ggplot2::aes(x = Loser_id, y = Winner_id)) +
    ggplot2::geom_tile(ggplot2::aes(fill = perc_success), colour = "grey70", linewidth = 0.25) +
    fill_scale +
    ggside::geom_ysidecol(ggplot2::aes(color = factor(col_cl))) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::geom_vline(xintercept = unlist(v_lines_list) + 0.5, color = "black", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = unlist(h_lines_list) + 0.5, color = "black", linewidth = 0.3) +
    ggplot2::labs(
      x = "Items (ordered by block)",
      y = "Items (ordered by block)",
      fill = "% success",
      color = "Block"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "left"
    ) +
    ggside::theme_ggside_void() +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          fill = unname(palette),
          colour = unname(palette),
          linewidth = 0.3
        )
      )
    ) +
    ggplot2::scale_y_discrete(limits = rev(order_ids), labels = x_lab, guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::scale_x_discrete(labels = x_lab, guide = ggplot2::guide_axis(n.dodge = 2)) +
    ggplot2::coord_fixed(ratio = 1)
  
  # Redundant cue: mark unplayed tiles with × so B/W printing is unambiguous
  if (isTRUE(mark_unplayed)) {
    unplayed_df <- Y_long_plot |>
      dplyr::filter(!Played)
    p <- p +
      ggplot2::geom_text(
        data = unplayed_df,
        ggplot2::aes(label = unplayed_mark),
        colour = "grey55",
        size = unplayed_mark_size / ggplot2::.pt,
        show.legend = FALSE
      )
  }
  
  print(p)
  
  # ---- B/W preview: convert both heatmap scale and block palette to greys ----
  if (isTRUE(bw_preview)) {
    hex_to_grey <- function(hex) {
      rgb <- grDevices::col2rgb(hex) / 255
      lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
      L <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
      grDevices::rgb(L, L, L)
    }
    
    pal_bw <- vapply(unname(palette), hex_to_grey, character(1))
    names(pal_bw) <- names(palette)
    
    # For the fill: explicit grey gradient, NA stays white
    p_bw <- p +
      ggplot2::scale_fill_gradient(
        low = "grey85",
        high = "grey10",
        limits = c(0, 1),
        oob = scales::squish,
        na.value = no_match_color
      ) +
      ggplot2::scale_color_manual(values = pal_bw) +
      ggplot2::labs(fill = "% success (B/W preview)", color = "Block (B/W preview)")
    
    print(p_bw)
    
  }
  
  return(p)
}


#' Figure 4 plotting function
#' Assignment probability heatmap (BT-SBM)
#'
#' Heatmap of item-wise posterior assignment probabilities for clusters
#' (relabeled so that Cluster 1 is the top block by decreasing \eqn{\lambda}).
#' Items are ordered by their most probable cluster and marginal wins.
#'
#' @param fit Output list from \code{gibbs_BT_SBM()} (must include \code{relabeled$assign_prob}).
#' @param w_ij Optional wins matrix to compute marginal wins for ordering and annotation.
#'   If \code{NULL}, items are ordered by most-probable cluster only.
#' @param players_df Optional data.frame for label lookup (e.g. season metadata).
#' @param players_id_col Optional column name in \code{players_df} that matches item IDs.
#' @param labels Optional character vector of display labels. Either named by item ID, or aligned with the items.
#' @param labels_key_col Optional column name in \code{players_df} to match item IDs.
#' @param labels_value_col Optional column name in \code{players_df} containing display labels.
#' @param clean_display_labels If \code{TRUE}, applies \code{clean_fun} to display labels.
#' @param clean_fun Optional function to prettify names. Default: identity.
#' @param k_show Optional integer number of clusters to show (defaults to all columns in \code{assign_prob}).
#' @param max_n_clust Where to filter the mcmc x_t. If not specified we use the modal K
#' @param fill_low,fill_high Colors for the heatmap gradient low/high.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' p <- plot_assignment_probabilities(fit, w_ij)
#' }
#' @export
#'
plot_assignment_probabilities <- function(
    fit,
    w_ij = NULL,
    max_n_clust = NULL,
    players_df = NULL,
    players_id_col = NULL,
    labels = NULL,
    labels_key_col = NULL,
    labels_value_col = NULL,
    clean_display_labels = NULL,
    clean_fun = clean_players_names,
    x_hat = NULL,
    order_ids = NULL,
    k_show = NULL,                # currently unused, kept for API compat
    fill_low = "#FFFFCC",
    fill_high = "#006400"
){
  # --- Item ids (stable, used internally) ---
  if (!is.null(w_ij)) {
    item_ids <- rownames(w_ij)
    if (is.null(item_ids)) item_ids <- paste0("Item_", seq_len(nrow(w_ij)))
  } else if (!is.null(rownames(fit$item_cluster_assignment_probs))) {
    item_ids <- rownames(fit$item_cluster_assignment_probs)
  } else {
    item_ids <- paste0("Item_", seq_len(nrow(fit$item_cluster_assignment_probs)))
  }
  
  # --- Display labels mapping (id -> label), consistent with other plot helpers ---
  to_title_safe <- function(x) {
    x <- gsub("[-_]", " ", x)
    if (requireNamespace("stringr", quietly = TRUE)) {
      return(stringr::str_to_title(x))
    }
    tools::toTitleCase(tolower(x))
  }
  clean_name_vec <- function(x) {
    x2 <- to_title_safe(x)
    if (!is.null(clean_fun)) x2 <- sapply(x2, clean_fun)
    x2
  }
  
  infer_labels_from_players_df <- function(df, ids, key_col = NULL, value_col = NULL) {
    if (is.null(df) || !is.data.frame(df) || length(ids) < 1) return(NULL)
    
    if (is.null(key_col)) {
      candidates <- c("player_slug", "player_id", "player", "player_label")
      candidates <- candidates[candidates %in% colnames(df)]
      if (length(candidates) == 0) return(NULL)
      overlaps <- vapply(candidates, function(cc) {
        sum(as.character(df[[cc]]) %in% ids, na.rm = TRUE)
      }, numeric(1))
      key_col <- candidates[which.max(overlaps)]
      if (is.na(key_col) || overlaps[which.max(overlaps)] == 0) return(NULL)
    }
    
    if (is.null(value_col)) {
      candidates <- c("player_label", "label", "player", key_col)
      value_col <- candidates[candidates %in% colnames(df)][1]
      if (is.na(value_col) || is.null(value_col)) return(NULL)
    }
    
    key <- as.character(df[[key_col]])
    val <- as.character(df[[value_col]])
    ok <- !is.na(key) & !is.na(val) & nzchar(key)
    if (!any(ok)) return(NULL)
    
    map <- tapply(val[ok], key[ok], function(v) v[which(!is.na(v) & nzchar(v))[1]])
    map <- map[!is.na(map) & nzchar(map)]
    if (length(map) == 0) return(NULL)
    list(map = map, key_col = key_col, value_col = value_col)
  }
  
  id_to_disp <- stats::setNames(clean_name_vec(item_ids), item_ids)
  labels_map <- NULL
  labels_from <- NULL
  
  if (!is.null(labels)) {
    if (is.character(labels) && !is.null(names(labels))) {
      labels_map <- labels
      labels_from <- "labels"
    } else if (is.character(labels) && length(labels) == length(item_ids)) {
      labels_map <- stats::setNames(labels, item_ids)
      labels_from <- "labels"
    } else {
      stop("`labels` must be either a named character vector (names are item ids), or a character vector aligned with the items.")
    }
  } else if (!is.null(players_df)) {
    if (is.null(labels_key_col) && !is.null(players_id_col)) labels_key_col <- players_id_col
    inferred <- infer_labels_from_players_df(players_df, item_ids, key_col = labels_key_col, value_col = labels_value_col)
    if (!is.null(inferred)) {
      labels_map <- inferred$map
      labels_from <- "players_df"
      labels_key_col <- inferred$key_col
      labels_value_col <- inferred$value_col
    }
  }
  
  if (!is.null(labels_map)) {
    if (is.null(clean_display_labels)) {
      clean_display_labels <- !(isTRUE(labels_from == "players_df") && !is.null(labels_value_col) && grepl("label", labels_value_col, ignore.case = TRUE))
    }
    labels_use <- labels_map
    if (isTRUE(clean_display_labels)) {
      labels_use <- clean_name_vec(labels_use)
      names(labels_use) <- names(labels_map)
    }
    overlap <- intersect(item_ids, names(labels_use))
    id_to_disp[overlap] <- unname(labels_use[overlap])
  }
  
  # --- Optional: common ordering based on a point-estimate partition ---
  if (is.null(x_hat)) {
    if (!is.null(fit$minVI_partition)) {
      x_hat <- fit$minVI_partition
    } else if (!is.null(fit$partition_binder)) {
      x_hat <- fit$partition_binder
    }
  }
  if (!is.null(x_hat) && length(x_hat) == length(item_ids)) {
    if (is.null(order_ids)) {
      if (!is.null(w_ij)) {
        N_ij <- w_ij + t(w_ij)
        marg_win <- rowSums(w_ij[item_ids, item_ids, drop = FALSE], na.rm = TRUE)
        marg_matches <- rowSums(N_ij[item_ids, item_ids, drop = FALSE], na.rm = TRUE)
        marg_win_prob <- ifelse(is.finite(marg_matches) & marg_matches > 0, marg_win / marg_matches, NA_real_)
      } else {
        marg_win <- rep(NA_real_, length(item_ids))
        marg_win_prob <- rep(NA_real_, length(item_ids))
      }
      ord_df <- data.frame(
        players = item_ids,
        cl = as.integer(x_hat),
        marginal_win_prob = marg_win_prob,
        marginal_win = marg_win,
        stringsAsFactors = FALSE
      )
      ord_df <- ord_df %>%
        dplyr::arrange(cl, dplyr::desc(marginal_win_prob), dplyr::desc(marginal_win), players)
      order_ids <- ord_df$players
    }
  }
  
  # --- Determine modal K if not provided ---
  x_samples <- fit$x_samples
  lambda_samples <- fit$lambda_samples
  
  unique_count <- apply(x_samples, 1, function(z) length(unique(z)))
  modal_K <- as.numeric(names(which.max(table(unique_count))))
  if (is.null(max_n_clust)) max_n_clust <- modal_K
  
  # Optionally subset samples (kept for future use / consistency)
  x_samples_sub <- x_samples[which(unique_count == max_n_clust), , drop = FALSE]
  lambda_samples_sub <- lambda_samples[which(unique_count == max_n_clust)]
  
  # --- Take first max_n_clust columns and renormalise per row ---
  if (is.null(fit$item_cluster_assignment_probs)) {
    stop("fit$item_cluster_assignment_probs is NULL.")
  }
  if (ncol(fit$item_cluster_assignment_probs) < max_n_clust) {
    stop("Requested max_n_clust exceeds available columns in item_cluster_assignment_probs.")
  }
  
  block_prob_mat <- as.matrix(fit$item_cluster_assignment_probs[, seq_len(max_n_clust), drop = FALSE])
  rownames(block_prob_mat) <- item_ids
  
  row_mass <- rowSums(block_prob_mat, na.rm = TRUE)
  zero_mass_idx <- which(!is.finite(row_mass) | row_mass <= 0)
  
  if (length(zero_mass_idx) > 0) {
    msg <- sprintf(
      "Dropping %d item(s) with zero probability mass within the shown %d clusters: %s",
      length(zero_mass_idx), max_n_clust, paste(rownames(block_prob_mat)[zero_mass_idx], collapse = ", ")
    )
    message(msg)
    block_prob_mat <- block_prob_mat[-zero_mass_idx, , drop = FALSE]
    item_ids <- item_ids[-zero_mass_idx]
    row_mass <- row_mass[-zero_mass_idx]
    if (nrow(block_prob_mat) == 0) stop("All items had zero mass in the filtered clusters.")
  }
  
  # Renormalise to sum to 1 across the shown clusters
  norm_block_prob_mat <- sweep(block_prob_mat, 1, row_mass, FUN = "/")
  # Numerical cleanup
  norm_block_prob_mat[!is.finite(norm_block_prob_mat)] <- 0
  
  # Optional: sanity check (can be commented out)
  # max(abs(rowSums(norm_block_prob_mat) - 1))
  
  # --- Long format for plotting ---
  block_prob <- as.data.frame(norm_block_prob_mat, check.names = FALSE)
  block_prob <- cbind(block_prob, pl_id = rownames(norm_block_prob_mat))
  
  cluster_cols <- names(block_prob)[seq_len(max_n_clust)]
  
  assignment_probs_long <- block_prob %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(cluster_cols),
      names_to = "Cluster",
      values_to = "prob"
    ) %>%
    dplyr::mutate(
      Cluster = gsub(x = Cluster, replacement = " ", pattern = "_")
    )
  
  # Cluster with highest assignment probability (post-renormalisation)
  max_prob_clusters <- assignment_probs_long %>%
    dplyr::group_by(pl_id) %>%
    dplyr::summarize(Cl_ass = Cluster[which.max(prob)], .groups = "drop")
  
  # --- Marginal win percentages (if w_ij is provided) ---
  if (!is.null(w_ij)) {
    if (is.null(rownames(w_ij))) stop("w_ij must have rownames to match item ids.")
    if (is.null(colnames(w_ij))) colnames(w_ij) <- rownames(w_ij)
    missing_ids <- setdiff(rownames(norm_block_prob_mat), rownames(w_ij))
    if (length(missing_ids) > 0L) {
      stop("w_ij is missing some item ids used in assignment probs: ", paste(missing_ids, collapse = ", "))
    }
    marg_pro_win <- data.frame(
      pl_id = rownames(norm_block_prob_mat),
      marg_pro_win  = rowSums(w_ij[rownames(norm_block_prob_mat), , drop = FALSE], na.rm = TRUE),
      marg_pro_loss = colSums(w_ij[, rownames(norm_block_prob_mat), drop = FALSE], na.rm = TRUE)
    ) %>%
      dplyr::mutate(
        pct_win = marg_pro_win / (marg_pro_loss + marg_pro_win)
      )
  } else {
    marg_pro_win <- data.frame(pl_id = assignment_probs_long$pl_id |> unique(),
                               marg_pro_win = NA_real_,
                               marg_pro_loss = NA_real_,
                               pct_win = NA_real_)
  }
  
  assignment_probs_long_plot <- assignment_probs_long %>%
    dplyr::left_join(max_prob_clusters, by = "pl_id") %>%
    dplyr::left_join(marg_pro_win, by = "pl_id") %>%
    dplyr::ungroup() |>
    dplyr::mutate(
      pl_name = unname(id_to_disp[pl_id])
    ) |>
    dplyr::mutate(
      pl_id = {
        if (!is.null(order_ids)) {
          factor(pl_id, levels = unique(order_ids), ordered = TRUE)
        } else {
          factor(pl_id, levels = unique(pl_id[order(Cl_ass, -marg_pro_win, decreasing = TRUE)]), ordered = TRUE)
        }
      }
    )
  
  # --- Plot ---
  ggplot2::ggplot(assignment_probs_long_plot) +
    ggplot2::geom_tile(ggplot2::aes(x = Cluster, y = pl_id, fill = prob)) +
    ggplot2::scale_fill_gradient(low = fill_low, high = fill_high, na.value = "#009680") +
    ggplot2::labs(x = "", y = "", fill = "Assign. Prob.") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_discrete(
      limits = if (!is.null(order_ids)) rev(unique(order_ids)) else waiver(),
      labels = function(x) unname(id_to_disp[x]),
      guide = ggplot2::guide_axis(n.dodge = 2)
    )
}

#' Figure 5 plotting function
#' Lambda uncertainty plot (per player)
#'
#' Forest plot of eqn{\lambda} with uncertainty intervals,
#' using relabeled draws. Points colored by the hard partition (\code{fit$estimates$x_hat}).
#'
#' @param fit Output from \code{gibbs_BT_SBM()} with \code{opt_lambda$lambda_item} computed
#'   (set \code{keep_lambda=TRUE} when sampling).
#' @param w_ij Optional wins matrix to compute marginal wins for ordering and annotation.
#' @param log_base Base for the x-axis logarithm (10 or e). Defaults to 10.
#' @param prob Interval probability for HPD (e.g., 0.90).
#' @param max_n_clust Where to filter the mcmc x_t. If not specified we use the modal K
#' @param palette Named colors for clusters.
#' @param clean_fun Optional player-name cleaner. Default: identity.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' # fit with keep_lambda=TRUE
#' p <- plot_lambda_uncertainty(fit, prob = 0.90)
#' }
#' @export
#'
plot_lambda_uncertainty <- function(fit,
                                    w_ij = NULL,
                                    players_df = NULL,
                                    players_id_col = NULL,
                                    labels = NULL,
                                    labels_key_col = NULL,
                                    labels_value_col = NULL,
                                    clean_display_labels = NULL,
                                    order_ids = NULL,
                                    log_base = exp(1),
                                    max_n_clust = NULL,
                                    prob = 0.9,
                                    palette = NULL,
                                    clean_fun = function(x) x,
                                    x_hat = NULL,
                                    filter_lambdas = TRUE,
                                    conditional = TRUE,
                                    ...) {
  if (is.null(fit$lambda_samples_relabel)) {
    stop("Expected fit$lambda_samples_relabel (player-level λ draws). Did you run relabel_by_lambda()?")
  }
  if (is.null(fit$x_samples_relabel)) {
    stop("Expected fit$x_samples_relabel (relabelled partition draws). Did you run relabel_by_lambda()?")
  }
  
  lambda_item <- fit$lambda_samples_relabel
  x_samples <- fit$x_samples_relabel
  
  if (is.list(lambda_item)) {
    lambda_item <- do.call(rbind, lambda_item)
  }
  lambda_item <- as.matrix(lambda_item)
  x_samples <- as.matrix(x_samples)
  
  S <- nrow(lambda_item)
  n <- ncol(lambda_item)
  if (nrow(x_samples) != S || ncol(x_samples) != n) {
    stop("x_samples_relabel and lambda_samples_relabel must have the same (iters x n_items) shape.")
  }
  if (S < 2) stop("Need at least 2 posterior draws for uncertainty plot")
  
  # ---- Optional: filter draws to a target K ----
  keep_idx <- seq_len(S)
  if (!is.null(max_n_clust)) {
    K_target <- as.integer(max_n_clust)
    if (!is.finite(K_target) || K_target < 1) stop("max_n_clust must be a positive integer (target K).")
    
    K_each <- apply(x_samples, 1, function(z) length(unique(z)))
    keep_idx <- which(K_each == K_target)
    if (length(keep_idx) == 0L) {
      warning(paste0("No iterations with exactly K=", K_target, ". Using all iterations."))
      keep_idx <- seq_len(S)
    }
  }
  
  x_samples_use <- x_samples[keep_idx, , drop = FALSE]
  
  # Filter lambdas ONLY when conditional=TRUE (because we index lambdas using x iterations).
  if (isTRUE(conditional)) {
    if (!isTRUE(filter_lambdas)) {
      warning("conditional=TRUE requires matching x/lambda iterations; overriding filter_lambdas=TRUE.")
    }
    lambda_use <- lambda_item[keep_idx, , drop = FALSE]
  } else {
    # Unconditional summaries: keep full lambda chain (no filtering),
    # even if we filtered x_samples to build x_hat / ordering.
    lambda_use <- lambda_item
  }
  
  # ---- Choose/compute the plotting partition x_hat ----
  if (is.null(x_hat)) {
    if (!is.null(max_n_clust) && length(keep_idx) < S) {
      psm_filt <- mcclust::comp.psm(x_samples_use)
      cl_obj <- mcclust.ext::minbinder.ext(psm_filt, cls.draw = x_samples_use, method = "draws")$cl
      if (is.matrix(cl_obj)) {
        x_hat <- as.integer(cl_obj[1, ])
      } else {
        x_hat <- as.integer(cl_obj)
      }
    } else if (!is.null(fit$minVI_partition)) {
      x_hat <- as.integer(fit$minVI_partition)
    } else if (!is.null(fit$partition_binder)) {
      x_hat <- as.integer(fit$partition_binder)
    } else {
      x_hat <- apply(x_samples_use, 2, function(v) {
        tab <- table(v)
        as.integer(names(tab)[which.max(tab)])
      })
    }
  }
  x_hat <- as.integer(x_hat)
  if (length(x_hat) != n) stop("x_hat must have length n_items")
  
  # ---- Player ids + display labels ----
  item_ids <- NULL
  if (!is.null(w_ij) && length(dim(w_ij)) == 2 && !is.null(rownames(w_ij))) {
    item_ids <- rownames(w_ij)
  }
  if (is.null(item_ids)) item_ids <- colnames(lambda_item)
  if (is.null(item_ids)) item_ids <- paste0("Item_", seq_len(n))
  if (length(item_ids) != n) item_ids <- paste0("Item_", seq_len(n))
  
  to_title_safe <- function(x) {
    x <- gsub("[-_]", " ", x)
    if (requireNamespace("stringr", quietly = TRUE)) {
      return(stringr::str_to_title(x))
    }
    tools::toTitleCase(tolower(x))
  }
  clean_name_vec <- function(x) {
    x2 <- to_title_safe(x)
    if (!is.null(clean_fun)) x2 <- sapply(x2, clean_fun)
    x2
  }
  infer_labels_from_players_df <- function(df, ids, key_col = NULL, value_col = NULL) {
    if (is.null(df) || !is.data.frame(df) || length(ids) < 1) return(NULL)
    
    if (is.null(key_col)) {
      candidates <- c("player_slug", "player_id", "player", "player_label")
      candidates <- candidates[candidates %in% colnames(df)]
      if (length(candidates) == 0) return(NULL)
      overlaps <- vapply(candidates, function(cc) {
        sum(as.character(df[[cc]]) %in% ids, na.rm = TRUE)
      }, numeric(1))
      key_col <- candidates[which.max(overlaps)]
      if (is.na(key_col) || overlaps[which.max(overlaps)] == 0) return(NULL)
    }
    
    if (is.null(value_col)) {
      candidates <- c("player_label", "label", "player", key_col)
      value_col <- candidates[candidates %in% colnames(df)][1]
      if (is.na(value_col) || is.null(value_col)) return(NULL)
    }
    
    key <- as.character(df[[key_col]])
    val <- as.character(df[[value_col]])
    ok <- !is.na(key) & !is.na(val) & nzchar(key)
    if (!any(ok)) return(NULL)
    
    map <- tapply(val[ok], key[ok], function(v) v[which(!is.na(v) & nzchar(v))[1]])
    map <- map[!is.na(map) & nzchar(map)]
    if (length(map) == 0) return(NULL)
    list(map = map, key_col = key_col, value_col = value_col)
  }
  
  id_to_disp <- stats::setNames(clean_name_vec(item_ids), item_ids)
  labels_map <- NULL
  labels_from <- NULL
  
  if (!is.null(labels)) {
    if (is.character(labels) && !is.null(names(labels))) {
      labels_map <- labels
      labels_from <- "labels"
    } else if (is.character(labels) && length(labels) == length(item_ids)) {
      labels_map <- stats::setNames(labels, item_ids)
      labels_from <- "labels"
    } else {
      stop("`labels` must be either a named character vector (names are item ids), or a character vector aligned with the items.")
    }
  } else if (!is.null(players_df)) {
    if (is.null(labels_key_col) && !is.null(players_id_col)) labels_key_col <- players_id_col
    inferred <- infer_labels_from_players_df(players_df, item_ids, key_col = labels_key_col, value_col = labels_value_col)
    if (!is.null(inferred)) {
      labels_map <- inferred$map
      labels_from <- "players_df"
      labels_key_col <- inferred$key_col
      labels_value_col <- inferred$value_col
    }
  }
  
  if (!is.null(labels_map)) {
    if (is.null(clean_display_labels)) {
      clean_display_labels <- !(isTRUE(labels_from == "players_df") && !is.null(labels_value_col) && grepl("label", labels_value_col, ignore.case = TRUE))
    }
    labels_use <- labels_map
    if (isTRUE(clean_display_labels)) {
      labels_use <- clean_name_vec(labels_use)
      names(labels_use) <- names(labels_map)
    }
    overlap <- intersect(item_ids, names(labels_use))
    id_to_disp[overlap] <- unname(labels_use[overlap])
  }
  
  # Common ordering across plots: x_hat cluster then within-cluster win rate (if w_ij provided)
  if (is.null(order_ids)) {
    if (!is.null(w_ij) && !is.null(rownames(w_ij))) {
      N_ij <- w_ij + t(w_ij)
      marg_win <- rowSums(w_ij[item_ids, item_ids, drop = FALSE], na.rm = TRUE)
      marg_matches <- rowSums(N_ij[item_ids, item_ids, drop = FALSE], na.rm = TRUE)
      marg_win_prob <- ifelse(is.finite(marg_matches) & marg_matches > 0, marg_win / marg_matches, NA_real_)
    } else {
      marg_win <- rep(NA_real_, n)
      marg_win_prob <- rep(NA_real_, n)
    }
    ord_df <- data.frame(
      players = item_ids,
      cl = as.integer(x_hat),
      marginal_win_prob = marg_win_prob,
      marginal_win = marg_win,
      stringsAsFactors = FALSE
    )
    ord_df <- ord_df %>%
      dplyr::arrange(cl, dplyr::desc(marginal_win_prob), dplyr::desc(marginal_win), players)
    order_ids <- ord_df$players
  }
  order_ids <- as.character(order_ids)
  order_ids <- order_ids[order_ids %in% item_ids]
  order_ids <- unique(order_ids)
  if (length(order_ids) != length(item_ids)) {
    order_ids <- c(order_ids, setdiff(item_ids, order_ids))
  }
  
  # ---- Summaries on log-scale (HPD), back-transform for plotting ----
  eps <- .Machine$double.eps
  logb <- function(x) log(pmax(x, eps), base = log_base)
  powb <- function(x) log_base^x
  
  # Conditional player strength: summarise λ_i only on iterations where the
  # player is assigned to their plotted tier x_hat[i].
  mean_hat <- rep(NA_real_, n)
  low_hat <- rep(NA_real_, n)
  up_hat <- rep(NA_real_, n)
  
  if (isTRUE(conditional)) {
    for (i in seq_len(n)) {
      idx_i <- which(x_samples_use[, i] == x_hat[i])
      if (length(idx_i) < 2L) {
        # Fallback to marginal for this player only.
        v <- lambda_use[, i]
      } else {
        v <- lambda_use[idx_i, i]
      }
      vlog <- logb(v)
      mean_hat[i] <- powb(mean(vlog, na.rm = TRUE))
      h <- coda::HPDinterval(coda::as.mcmc(vlog), prob = prob)[1, ]
      low_hat[i] <- powb(h[1])
      up_hat[i] <- powb(h[2])
    }
  } else {
    log_lp <- apply(lambda_use, 2, logb)
    if (is.vector(log_lp)) log_lp <- matrix(log_lp, ncol = 1)
    hpd_mat <- t(apply(log_lp, 2, function(v) {
      coda::HPDinterval(coda::as.mcmc(v), prob = prob)[1, ]
    }))
    mean_hat <- powb(colMeans(log_lp))
    low_hat <- powb(hpd_mat[, 1])
    up_hat <- powb(hpd_mat[, 2])
  }
  
  # ---- Relabel clusters so 1 = strongest ----
  lambda_mean_item <- mean_hat
  labs <- sort(unique(x_hat))
  cl_means <- vapply(labs, function(k) mean(lambda_mean_item[x_hat == k], na.rm = TRUE), numeric(1))
  ord <- order(cl_means, decreasing = TRUE)
  map <- setNames(seq_along(ord), as.character(labs[ord]))
  cluster_ord <- as.integer(map[as.character(x_hat)])
  K_plot <- length(unique(cluster_ord))
  if (!is.null(max_n_clust) && K_plot != as.integer(max_n_clust)) {
    warning(paste0("Requested K=", as.integer(max_n_clust), " but x_hat has K=", K_plot, "."))
  }
  
  player_summ <- data.frame(
    player_id = item_ids,
    player = unname(id_to_disp[item_ids]),
    mean = mean_hat,
    low = low_hat,
    up = up_hat,
    cluster = cluster_ord,
    stringsAsFactors = FALSE
  )
  
  # Use common ordering for y-axis; show strongest at top via reversed limits.
  player_summ$player_id <- factor(player_summ$player_id, levels = order_ids, ordered = TRUE)
  player_summ$cluster <- factor(player_summ$cluster, levels = seq_len(K_plot), ordered = TRUE)
  
  # ---- palette handling ----
  # Use a high-contrast qualitative palette by default; fall back gracefully
  # if btsbm_block_palette is not available in the current session.
  if (exists("btsbm_block_palette", mode = "function")) {
    pal <- btsbm_block_palette(K_plot, palette)
  } else {
    if (!is.null(palette) && length(palette) >= K_plot) {
      pal <- palette[seq_len(K_plot)]
    } else if (!is.null(palette) && length(palette) >= 2) {
      pal <- grDevices::colorRampPalette(palette)(K_plot)
    } else {
      pal <- grDevices::hcl.colors(K_plot, palette = "Dark 3")
    }
  }
  names(pal) <- as.character(seq_len(K_plot))
  
  if (isTRUE(conditional)) {
    p <- ggplot2::ggplot(player_summ, ggplot2::aes(x = mean, y = player_id, colour = cluster)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = low, xmax = up),
        width = 0, linewidth = 0.4, alpha = 0.7, orientation = "y"
      ) +
      ggplot2::geom_point(size = 1.7, alpha = 0.9) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, alpha = 0.7) +
      ggplot2::scale_x_continuous(
        trans = scales::log_trans(base = log_base),
        labels = scales::label_number(accuracy = 0.01),
        breaks = scales::breaks_log(base = log_base)
      ) +
      ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = "Cluster\n(1=strongest)") +
      ggplot2::scale_y_discrete(
        limits = rev(order_ids),
        labels = function(x) unname(id_to_disp[x]),
        guide = ggplot2::guide_axis(n.dodge = 2)
      ) +
      ggplot2::labs(x = bquote(lambda ~ " (posterior mean; log scale)"), y = NULL) +
      theme_btsbm() +
      ggplot2::theme(legend.position = "right")
  } else {
    p <- ggplot2::ggplot(player_summ, ggplot2::aes(x = mean, y = player_id)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = low, xmax = up),
        width = 0, linewidth = 0.4, alpha = 0.75, orientation = "y",
        colour = "grey35"
      ) +
      ggplot2::geom_point(size = 1.7, alpha = 0.95, colour = "black") +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, alpha = 0.7) +
      ggplot2::scale_x_continuous(
        trans = scales::log_trans(base = log_base),
        labels = scales::label_number(accuracy = 0.01),
        breaks = scales::breaks_log(base = log_base)
      ) +
      ggplot2::scale_y_discrete(
        limits = rev(order_ids),
        labels = function(x) unname(id_to_disp[x]),
        guide = ggplot2::guide_axis(n.dodge = 2)
      ) +
      ggplot2::labs(x = bquote(lambda ~ " (posterior mean; log scale)"), y = NULL) +
      theme_btsbm() +
      ggplot2::theme(legend.position = "none")
  }
  
  attr(p, "player_summary") <- dplyr::mutate(
    dplyr::arrange(as.data.frame(player_summ), cluster, dplyr::desc(mean)),
    cluster = as.integer(as.character(cluster))
  )
  p
}







# Relabel clusters by decreasing lambda each draw (helps with label switching).
# Also computes representative partitions (minVI, Binder) and credible balls.
relabel_by_lambda <- function(x_samples, lambda_samples) {
  stopifnot(is.matrix(x_samples))
  S <- nrow(x_samples)
  N <- ncol(x_samples)

  if (!requireNamespace("mcclust", quietly = TRUE)) {
    stop("Package 'mcclust' is required for relabel_by_lambda().")
  }
  if (!requireNamespace("mcclust.ext", quietly = TRUE)) {
    stop("Package 'mcclust.ext' is required for relabel_by_lambda().")
  }

  is_list_format <- is.list(lambda_samples)
  get_lambda_vec <- function(iter) {
    if (is_list_format) {
      v <- lambda_samples[[iter]]
      if (!is.numeric(v)) stop("lambda_samples[[iter]] must be numeric.")
      return(v)
    }
    lambda_samples[iter, ]
  }

  x_relabeled <- matrix(NA_integer_, S, N)
  lambda_per_item <- matrix(NA_real_, S, N)
  cluster_lambda_ordered <- vector("list", S)
  n_clusters_each_iter <- integer(S)
  top_block_count_per_iter <- integer(S)

  for (iter in seq_len(S)) {
    xi <- as.integer(x_samples[iter, ])
    occ_raw <- sort(unique(xi))
    xi_seq <- match(xi, occ_raw)
    K <- max(xi_seq)

    lam_vec_full <- get_lambda_vec(iter)
    lam_occ <- rep(NA_real_, length(occ_raw))
    ok_idx <- occ_raw <= length(lam_vec_full)
    lam_occ[ok_idx] <- lam_vec_full[occ_raw[ok_idx]]

    ord <- order(lam_occ, decreasing = TRUE, na.last = TRUE)
    occ_ord <- occ_raw[ord]
    lam_ord <- lam_occ[ord]
    if (anyNA(lam_ord)) lam_ord[is.na(lam_ord)] <- .Machine$double.xmin

    raw_to_ord_id <- integer(max(occ_ord))
    raw_to_ord_id[occ_ord] <- seq_len(K)
    xi_new <- raw_to_ord_id[occ_raw[xi_seq]]

    x_relabeled[iter, ] <- xi_new
    lambda_per_item[iter, ] <- lam_ord[xi_new]
    cluster_lambda_ordered[[iter]] <- lam_ord
    n_clusters_each_iter[iter] <- K
    top_block_count_per_iter[iter] <- sum(xi_new == 1L)
  }

  psm <- mcclust::comp.psm(x_samples)

  partition_binder <- mcclust.ext::minbinder.ext(
    psm, cls.draw = x_samples, method = "all"
  )$cl[1, ]
  partition_minVI <- mcclust.ext::minVI(
    psm, cls.draw = x_samples, method = "all"
  )$cl[1, ]

  x_ball <- mcclust.ext::credibleball(
    c.star = partition_minVI, cls.draw = x_samples, c.dist = "VI"
  )

  relabel_partition_by_item_mean_lambda <- function(z, lambda_item_mean) {
    stopifnot(length(z) == length(lambda_item_mean))
    z <- as.integer(z)
    labs <- sort(unique(z))
    cl_means <- vapply(
      labs,
      function(k) mean(lambda_item_mean[z == k], na.rm = TRUE),
      numeric(1)
    )
    ord <- order(cl_means, decreasing = TRUE)
    new_ids <- seq_along(labs)
    names(new_ids) <- labs[ord]
    as.integer(new_ids[as.character(z)])
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
      d <- apply(obj, 1, function(z) mcclust::vi.dist(cl1 = as.integer(z), cl2 = as.integer(centre)))
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

  Kmax <- N
  assignment_probs <- matrix(0, nrow = N, ncol = Kmax)
  for (k in seq_len(Kmax)) {
    assignment_probs[, k] <- colMeans(x_relabeled == k, na.rm = TRUE)
  }
  colnames(assignment_probs) <- paste0("Cluster_", seq_len(Kmax))
  rownames(assignment_probs) <- paste0("Item_", seq_len(N))

  list(
    x_samples_relabel = x_relabeled,
    lambda_samples_relabel = lambda_per_item,
    cluster_lambda_ordered = cluster_lambda_ordered,
    co_clustering = psm,
    minVI_partition = partition_minVI,
    partition_binder = partition_binder,
    n_clusters_each_iter = n_clusters_each_iter,
    block_count_distribution = block_count_df,
    item_cluster_assignment_probs = as.data.frame(assignment_probs),
    avg_top_block_count = mean(top_block_count_per_iter),
    top_block_count_per_iter = top_block_count_per_iter,
    top_block_size_distribution = top_block_size_df,
    credible_ball_lower_partition = c_lower_rl,
    credible_ball_upper_partition = c_upper_rl,
    credible_ball_horiz_partition = c_horiz_rl
  )
}




# ---- Rank by posterior expected wins, summarised via midranks + position probabilities ----
# ---- Modified: posterior rank distribution via midranks of lambda draws ----
compute_expected_wins_rank_posterior <- function(lambda_item,
                                                 w_ij = NULL,          # kept for compatibility, ignored
                                                 player_names = NULL,
                                                 burn_in = 0,
                                                 thin = 1,
                                                 ci = 0.95,
                                                 ties.method = "average",
                                                 rank_grid_by = 0.5) {
  
  if (!is.null(w_ij)) {
    warning("w_ij is ignored: ranking is computed by sorting lambda draws only (no schedule weighting).")
  }
  
  if (is.list(lambda_item)) lambda_item <- do.call(rbind, lambda_item)
  lambda_item <- as.matrix(lambda_item)
  
  S <- nrow(lambda_item)
  n <- ncol(lambda_item)
  
  # burn-in + thinning
  B <- burn_in
  if (B > 0 && B < 1) B <- floor(S * B)
  B <- as.integer(B)
  thin <- as.integer(thin)
  
  keep <- seq.int(B + 1L, S, by = thin)
  lam_use <- lambda_item[keep, , drop = FALSE]
  Tkeep <- nrow(lam_use)
  
  # rank samples: rank 1 = best (largest lambda)
  rank_samples <- matrix(NA_real_, nrow = Tkeep, ncol = n)
  for (t in seq_len(Tkeep)) {
    lam <- lam_use[t, ]
    rank_samples[t, ] <- rank(-lam, ties.method = ties.method)
  }
  
  # player names
  if (is.null(player_names) || length(player_names) != n) player_names <- colnames(lambda_item)
  if (is.null(player_names) || length(player_names) != n) player_names <- paste0("Player_", seq_len(n))
  colnames(rank_samples) <- as.character(player_names)
  
  # posterior summaries of (mid)rank
  alpha <- (1 - ci) / 2
  rank_med  <- apply(rank_samples, 2, stats::median, na.rm = TRUE)
  rank_low  <- apply(rank_samples, 2, stats::quantile, probs = alpha, na.rm = TRUE)
  rank_high <- apply(rank_samples, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE)
  rank_mean <- colMeans(rank_samples, na.rm = TRUE)
  
  summary <- data.frame(
    player = as.character(player_names),
    rank_mean   = as.numeric(rank_mean),
    rank_median = as.numeric(rank_med),
    rank_low    = as.numeric(rank_low),
    rank_high   = as.numeric(rank_high),
    stringsAsFactors = FALSE
  )
  
  # posterior distribution over midrank grid (needed because midranks can be 2.5 etc.)
  grid <- seq(1, n, by = rank_grid_by)
  
  rank_prob <- t(vapply(seq_len(n), function(j) {
    idx <- match(rank_samples[, j], grid)
    tabulate(idx, nbins = length(grid)) / Tkeep
  }, numeric(length(grid))))
  
  colnames(rank_prob) <- paste0("R", grid)
  rownames(rank_prob) <- as.character(player_names)
  
  list(
    rank_samples = rank_samples,   # Tkeep x n, numeric midranks
    rank_prob    = rank_prob,      # n x length(grid), pmf over midranks
    summary      = summary,
    keep         = keep,
    grid         = grid,
    ties_method  = ties.method
  )
}



plot_rank_intervals <- function(x,
                                max_players = NULL,
                                w_ij = NULL,
                                player_names = NULL,
                                burn_in = 0,
                                thin = 1,
                                ci = 0.95,
                                ties.method = "average",
                                rank_grid_by = 0.5) {
  
  # x can be one of:
  #  - inf_i output from relabel_by_lambda() (has $lambda_samples_relabel)
  #  - lambda draws (matrix/list) -> compute_expected_wins_rank_posterior() is called
  #  - a result list returned by compute_expected_wins_rank_posterior() (has $summary)
  #  - a precomputed summary data.frame (must contain player + rank_* or midrank_* cols)
  
  if (is.list(x) && !is.null(x$summary) && is.data.frame(x$summary)) {
    rank_summary <- x$summary
  } else if (is.data.frame(x)) {
    rank_summary <- x
  } else if (is.list(x) && !is.null(x$lambda_samples_relabel)) {
    # Accept relabel_by_lambda() output directly.
    rank_res <- compute_expected_wins_rank_posterior(
      lambda_item   = x$lambda_samples_relabel,
      w_ij          = w_ij,
      player_names  = player_names,
      burn_in       = burn_in,
      thin          = thin,
      ci            = ci,
      ties.method   = ties.method,
      rank_grid_by  = rank_grid_by
    )
    rank_summary <- rank_res$summary
  } else {
    rank_res <- compute_expected_wins_rank_posterior(
      lambda_item   = x,
      w_ij          = w_ij,
      player_names  = player_names,
      burn_in       = burn_in,
      thin          = thin,
      ci            = ci,
      ties.method   = ties.method,
      rank_grid_by  = rank_grid_by
    )
    rank_summary <- rank_res$summary
  }
  
  df <- rank_summary
  
  # Support both historical naming conventions.
  if (all(c("rank_mean", "rank_low", "rank_high") %in% names(df))) {
    mean_col <- "rank_mean"; low_col <- "rank_low"; high_col <- "rank_high"
  } else if (all(c("midrank_mean", "midrank_low", "midrank_high") %in% names(df))) {
    mean_col <- "midrank_mean"; low_col <- "midrank_low"; high_col <- "midrank_high"
  } else if (all(c("midrank_median", "midrank_low", "midrank_high") %in% names(df))) {
    mean_col <- "midrank_median"; low_col <- "midrank_low"; high_col <- "midrank_high"
  } else if (all(c("rank_median", "rank_low", "rank_high") %in% names(df))) {
    mean_col <- "rank_median"; low_col <- "rank_low"; high_col <- "rank_high"
  } else {
    stop(
      "rank_summary must contain either rank_mean/rank_low/rank_high (preferred) or ",
      "midrank_* columns. Found: ",
      paste(names(df), collapse = ", ")
    )
  }
  
  if (!"player" %in% names(df)) {
    stop("rank_summary must contain a 'player' column.")
  }
  
  df <- df[order(df[[mean_col]], df[[high_col]], df[[low_col]]), , drop = FALSE]
  if (!is.null(max_players)) df <- utils::head(df, max_players)
  df$player <- factor(df$player, levels = rev(df$player))
  
  ggplot(df, aes(x = !!rlang::sym(mean_col), y = player)) +
    geom_errorbar(aes(xmin = !!rlang::sym(low_col), xmax = !!rlang::sym(high_col)),
                  width = 0, linewidth = 0.4, alpha = 0.8) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      x = "Posterior rank (1 = best)",
      y = NULL,
      caption = "Points: posterior mean (or median) rank. Bars: credible interval."
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
}



