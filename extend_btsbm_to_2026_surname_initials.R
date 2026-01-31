#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(readr)
  library(purrr)
  library(BTSBM)
})

# =========================================================
# User inputs
# =========================================================
matches_path <- "/Users/lapo_santi/Downloads/atp_tennis.csv"   # <-- change me
top_n <- 105                                                  # players per season
out_rds <- "./data/ATP_2000_2026_SN_extended.rds"             # <-- change me

# =========================================================
# Helpers
# =========================================================

# Make labels unique in a readable way: "Name", "Name (2)", "Name (3)", ...
make_unique_labels <- function(x) {
  idx <- ave(seq_along(x), x, FUN = seq_along)
  ifelse(idx == 1, x, paste0(x, " (", idx, ")"))
}

# Convert BTSBM "Given_Surname" / "Given Given_Surname" -> "Surname I." / "Surname I.I."
btsbm_to_surname_initials <- function(x) {
  vapply(x, function(name) {
    parts <- str_split(name, "_", simplify = TRUE)
    parts <- parts[parts != ""]
    if (length(parts) < 2) return(NA_character_)

    surname <- parts[length(parts)]
    given_chunks <- parts[-length(parts)]

    given_tokens <- unlist(str_split(given_chunks, "[[:space:]-]+"))
    given_tokens <- given_tokens[given_tokens != ""]

    initials <- paste0(substr(given_tokens, 1, 1), collapse = ".")
    paste0(surname, " ", initials, ".")
  }, character(1))
}

# Build a BTSBM-like season object from a "yearly_data[[yr]]" object
season_from_yearly_data <- function(season_obj) {
  stopifnot(!is.null(season_obj$Y_ij), !is.null(season_obj$N_ij), !is.null(season_obj$players_df))

  Y <- season_obj$Y_ij
  N <- season_obj$N_ij

  labs <- rownames(Y)
  stopifnot(identical(labs, colnames(Y)))
  stopifnot(identical(labs, rownames(N)))
  stopifnot(identical(labs, colnames(N)))

  df <- season_obj$players_df

  # Ensure player_label exists and matches matrix order
  if (!("player_label" %in% names(df))) {
    if ("player" %in% names(df)) {
      df <- df %>% mutate(player_label = player)
    } else {
      stop("players_df missing player_label and player columns.")
    }
  }

  df <- df %>%
    mutate(.ord = match(player_label, labs)) %>%
    arrange(.ord)

  if (any(is.na(df$.ord))) {
    bad <- df %>% filter(is.na(.ord)) %>% dplyr::select(any_of(c("player", "player_label"))) %>% head(10)
    stop("Some players_df labels do not match matrix rownames. Example rows:\n",
         paste(capture.output(print(bad)), collapse = "\n"))
  }

  # BTSBM-like players_df (age/height unknown for added seasons)
  players_df_out <- df %>%
    transmute(
      player       = seq_along(labs),
      worst_rank   = as.integer(worst_rank),
      median_rank  = as.numeric(median_rank),
      last_rank    = as.integer(last_rank),
      age_year     = NA_real_,
      ht_year      = NA_real_,
      player_slug  = player_label,
      player_label = player_label
    )

  list(Y_ij = Y, N_ij = N, players_df = players_df_out)
}

# Per-season duplicate check on matrix labels
dup_report_by_season <- function(ATP_obj) {
  imap_dfr(ATP_obj, function(season, yr) {
    labs <- rownames(season$Y_ij)
    dups <- unique(labs[duplicated(labs)])
    if (length(dups) == 0) return(NULL)
    tibble(year = yr, duplicate_label = sort(dups))
  })
}

# =========================================================
# 1) Read + clean Kallog match DB
# =========================================================
raw <- read_csv(matches_path, show_col_types = FALSE)

needed <- c("Tournament","Date","Player_1","Player_2","Winner","Rank_1","Rank_2","Pts_1","Pts_2")
stopifnot(all(needed %in% names(raw)))

matches_df <- raw %>%
  mutate(
    Date = as.Date(Date),
    match_year = year(Date),
    Player_1 = str_squish(Player_1),
    Player_2 = str_squish(Player_2),
    Winner   = str_squish(Winner),
    Rank_1   = na_if(as.integer(Rank_1), -1),
    Rank_2   = na_if(as.integer(Rank_2), -1),
    Pts_1    = na_if(as.integer(Pts_1), -1),
    Pts_2    = na_if(as.integer(Pts_2), -1)
  ) %>%
  mutate(
    winner_name = Winner,
    loser_name  = case_when(
      Winner == Player_1 ~ Player_2,
      Winner == Player_2 ~ Player_1,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Date), !is.na(winner_name), !is.na(loser_name))

# =========================================================
# 2) Extract weekly rankings from matches
# =========================================================
rankings <- bind_rows(
  matches_df %>% transmute(ranking_date = Date, player = Player_1, rank = Rank_1, points = Pts_1),
  matches_df %>% transmute(ranking_date = Date, player = Player_2, rank = Rank_2, points = Pts_2)
) %>%
  filter(!is.na(rank), rank > 0, !is.na(player), player != "") %>%
  mutate(
    ranking_year = year(ranking_date),
    ranking_week = isoweek(ranking_date),
    week_start   = floor_date(ranking_date, unit = "week", week_start = 1)
  ) %>%
  group_by(player, ranking_year, ranking_week, week_start) %>%
  summarise(
    rank   = median(rank, na.rm = TRUE),
    points = median(points, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(ranking_date = week_start)

# =========================================================
# 3) Build Kallog player key using "Surname I." format
#     IMPORTANT: do NOT slugify as label; keep human-readable style.
# =========================================================
player_key <- tibble(player = sort(unique(rankings$player))) %>%
  mutate(
    player_label = str_squish(player)
  ) %>%
  mutate(player_label = make_unique_labels(player_label)) %>%
  dplyr::select(player, player_label)

# =========================================================
# 4) Build yearly Y_ij / N_ij from Kallog (calendar-year seasons)
# =========================================================
start_year <- min(matches_df$match_year, na.rm = TRUE)
end_year   <- max(matches_df$match_year, na.rm = TRUE)
all_years  <- start_year:end_year

yearly_data <- list()

for (yr in all_years) {

  filtered_matches <- matches_df %>% filter(match_year == yr)
  if (nrow(filtered_matches) == 0) next

  rankings_this_year <- rankings %>% filter(ranking_year == yr)
  if (nrow(rankings_this_year) == 0) next

  top_players <- rankings_this_year %>%
    arrange(player, ranking_date) %>%
    group_by(player) %>%
    summarise(
      worst_rank  = max(rank, na.rm = TRUE),
      median_rank = median(rank, na.rm = TRUE),
      last_rank   = rank[which.max(ranking_date)],
      last_date   = max(ranking_date),
      .groups = "drop"
    ) %>%
    filter(!is.na(last_rank)) %>%
    arrange(last_rank) %>%
    slice_head(n = top_n)

  top_players_info <- top_players %>%
    left_join(player_key, by = "player") %>%
    filter(!is.na(player_label))

  if (nrow(top_players_info) == 0) next

  top_names <- top_players_info$player

  filtered_matches_top <- filtered_matches %>%
    dplyr::select(Date, winner_name, loser_name) %>%
    filter(winner_name %in% top_names, loser_name %in% top_names) %>%
    left_join(top_players_info %>% dplyr::select(player, player_label), by = c("winner_name" = "player")) %>%
    rename(winner_label = player_label) %>%
    left_join(top_players_info %>% dplyr::select(player, player_label), by = c("loser_name" = "player")) %>%
    rename(loser_label = player_label) %>%
    filter(!is.na(winner_label), !is.na(loser_label))

  if (nrow(filtered_matches_top) == 0) next

  # Ensure labels are unique (within-year) and stable
  labels_raw <- top_players_info$player_label
  labels <- make_unique_labels(labels_raw)
  top_players_info <- top_players_info %>% mutate(player_label = labels)

  n <- length(labels)
  Y <- matrix(0L, nrow = n, ncol = n, dimnames = list(labels, labels))

  win_counts <- filtered_matches_top %>% dplyr::count(winner_label, loser_label, name = "wins")
  for (k in seq_len(nrow(win_counts))) {
    Y[win_counts$winner_label[k], win_counts$loser_label[k]] <- win_counts$wins[k]
  }

  N <- Y + t(Y)
  diag(N) <- 0L

  yearly_data[[as.character(yr)]] <- list(
    Y_ij       = Y,
    N_ij       = N,
    players_df = top_players_info
  )
}

# =========================================================
# 5) Normalise BTSBM 2000–2022 to "Surname I." and make labels unique per season
# =========================================================
ATP_2000_2022_SN <- BTSBM::ATP_2000_2022

for (yr in names(ATP_2000_2022_SN)) {
  obj <- ATP_2000_2022_SN[[yr]]

  old_labels <- rownames(obj$Y_ij)
  new_labels <- btsbm_to_surname_initials(old_labels)
  new_labels <- make_unique_labels(new_labels)

  dimnames(obj$Y_ij) <- list(new_labels, new_labels)
  dimnames(obj$N_ij) <- list(new_labels, new_labels)

  obj$players_df <- obj$players_df %>%
    mutate(player_label = btsbm_to_surname_initials(player_slug)) %>%
    mutate(player_label = make_unique_labels(player_label))

  ATP_2000_2022_SN[[yr]] <- obj
}

# =========================================================
# 6) Extend with new seasons beyond BTSBM (e.g., 2023+)
# =========================================================
ATP_extended_SN <- ATP_2000_2022_SN

new_years <- setdiff(names(yearly_data), names(ATP_2000_2022_SN))
new_years <- new_years[1:3]
new_years <- new_years[order(as.integer(new_years))]

for (yr in new_years) {
  ATP_extended_SN[[yr]] <- season_from_yearly_data(yearly_data[[yr]])
}

# =========================================================
# 7) Checks + save
# =========================================================
dup_matrix <- dup_report_by_season(ATP_extended_SN)
if (nrow(dup_matrix) > 0) {
  message("WARNING: duplicate labels detected in some seasons (after normalisation):")
  print(dup_matrix)
}

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(ATP_extended_SN, out_rds)

message("Saved extended dataset to: ", out_rds)
message("Seasons available: ", paste(tail(names(ATP_extended_SN), 8), collapse = ", "))

# quick sanity: show format examples
if ("2021" %in% names(ATP_extended_SN)) {
  message("Example 2021 labels: ", paste(head(ATP_extended_SN$`2021`$players_df$player_label, 5), collapse = " | "))
}
last_year <- as.character(max(as.integer(names(ATP_extended_SN))))
message("Example ", last_year, " labels: ", paste(head(ATP_extended_SN[[last_year]]$players_df$player_label, 5), collapse = " | "))
