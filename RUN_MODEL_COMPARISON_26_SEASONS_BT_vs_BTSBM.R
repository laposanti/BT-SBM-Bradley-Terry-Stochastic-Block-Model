
# If run via `Rscript`, set working directory to this script's folder
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- sub("^--file=", "", file_arg)
  # Rscript can encode spaces as ~+~ in command args
  script_path <- gsub("~\\+~", " ", script_path)
  script_dir <- dirname(normalizePath(script_path))
  setwd(script_dir)
}

if (!dir.exists("./results")) dir.create("./results", recursive = TRUE)
if (!dir.exists("./images")) dir.create("./images", recursive = TRUE)
if (!dir.exists("./tables")) dir.create("./tables", recursive = TRUE)

results_csv <- "./results/model_comparison_26_seasons_bt_vs_btsbm.csv"
out_plot    <- "./images/DELPD_plot1.png"
out_plot_pdf <- "./images/DELPD_plot1.pdf"
out_table_by_season <- "./tables/model_comparison_by_season.tex"
out_table_summary   <- "./tables/model_comparison_summary.tex"

# Run mode:
# - "all": re-run all available seasons and overwrite results_csv
# - "append_missing_last_n": compute missing seasons among last N and append
run_mode <- Sys.getenv("RUN_MODE", unset = "all")
only_last_n <- as.integer(Sys.getenv("ONLY_LAST_N", unset = "3"))

# Also allow `Rscript RUN_MODEL_COMPARISON_26_SEASONS_BT_vs_BTSBM.R append` to switch modes
args_trailing <- commandArgs(trailingOnly = TRUE)
if (length(args_trailing) >= 1L && nzchar(args_trailing[1])) {
  if (tolower(args_trailing[1]) %in% c("append", "append_missing_last_n")) {
    run_mode <- "append_missing_last_n"
  } else if (tolower(args_trailing[1]) %in% c("all", "rerun", "full")) {
    run_mode <- "all"
  }
}
######################################################################
## 0.  Packages  –––––––––––––––––––––––––––––––––––––––––––––––––– ##
######################################################################
pkgs <- c("loo", "purrr", "dplyr", "tibble", "ggplot2", "tidyr","BTSBM")

suppressPackageStartupMessages(lapply(pkgs, require, character.only = TRUE))

set.seed(2025)       # reproducibility

######################################################################
## 1.  Load the data  ––––––––––––––––––––––––––––––––––––––––––––– ##
######################################################################
tennis_years <- readRDS("./data/ATP_2000_2025_SN_extended.rds")



######################################################################
## 3.  Main loop over seasons  –––––––––––––––––––––––––––––––––––– ##
######################################################################
available_years <- sort(as.integer(names(tennis_years)))

existing_results <- NULL
if (file.exists(results_csv)) {
  existing_results <- tryCatch(utils::read.csv(results_csv), error = function(e) NULL)
}
existing_years <- integer(0)
if (!is.null(existing_results) && ("season" %in% names(existing_results))) {
  existing_years <- sort(unique(as.integer(existing_results$season)))
}

years_to_run <- integer(0)
if (identical(run_mode, "append_missing_last_n")) {
  if (!is.finite(only_last_n) || only_last_n < 1L) only_last_n <- 3L
  target_years <- tail(available_years, only_last_n)
  years_to_run <- sort(setdiff(target_years, existing_years))
  if (length(years_to_run) == 0L) {
    message("No new seasons to compute (append mode); using existing results in ", results_csv)
    results_new <- tibble::tibble()
  } else {
    message("Append mode: computing missing seasons: ", paste(years_to_run, collapse = ", "))
    results_new <- purrr::map_dfr(years_to_run, function(yr) {
      
      dat <- tennis_years[[as.character(yr)]]
      Y   <- dat$Y_ij       # wins
      N   <- dat$N_ij       # matches played
      
      ############################################################
      ## 3a.  Fit both models
      ############################################################
      
      fit_simple  <- gibbs_bt_simple(w_ij = Y,
                                     T_iter = 10000, T_burn = 1500,
                                     verbose = T)
      
      fit_cluster <- gibbs_bt_sbm(w_ij = Y, T_iter = 10000, T_burn = 1500,
                                 prior = "GN", a = 1.5,
                                 gamma_GN = 0.8,
                                 verbose = T)
      
      ############################################################
      ## 3b.  Build point‑wise log‑likelihood matrices
      ############################################################
      
      lambda_simple <- fit_simple$lambda_samples
      if (is.list(lambda_simple)) lambda_simple <- do.call(rbind, lambda_simple)
      ll_simple <- make_bt_simple_loo(w_ij = Y, lambda_samples = lambda_simple)
      
      lambda_cluster <- fit_cluster$lambda_samples
      if (is.list(lambda_cluster)) lambda_cluster <- do.call(rbind, lambda_cluster)
      ll_cluster <- make_bt_cluster_loo(
        w_ij = Y,
        lambda_samples = lambda_cluster,
        x_samples = fit_cluster$x_samples
      )
      
      ############################################################
      ## 3c.  LOO, ELPD differences, stacking weights
      ############################################################
      loo_s <- loo::loo(ll_simple$ll)
      loo_c <- loo::loo(ll_cluster$ll)
      
      # elpd difference (cluster ‑ simple) and its SE
      elpd_diff <- loo_c$estimates["elpd_loo","Estimate"] -
        loo_s$estimates["elpd_loo","Estimate"]
      
      se_diff   <- sqrt(loo_c$estimates["elpd_loo","SE"]^2 +
                          loo_s$estimates["elpd_loo","SE"]^2)/2
      
      ## --- WAIC --------------------------------------------------------------
      waic_s <- loo::waic(ll_simple$ll)
      waic_c <- loo::waic(ll_cluster$ll)
      
      delta_waic   <- waic_s$estimate["waic", "Estimate"] -
        waic_c$estimate["waic", "Estimate"]      # >0 ⇒ cluster better
      se_delta_w   <- sqrt(sum(waic_s$estimate["waic", "SE"]^2,
                               waic_c$estimate["waic", "SE"]^2))/2
      
      ############################################################
      ## 3e.  Collect and return one row
      ############################################################
      to_be_returned = tibble::tibble(
        season              = as.integer(yr),
        elpd_loo_simple     = loo_s$estimates["elpd_loo","Estimate"],
        elpd_loo_cluster    = loo_c$estimates["elpd_loo","Estimate"],
        elpd_diff           = elpd_diff,
        se_diff             = se_diff,
        waic_s              = waic_s$estimate["waic", "Estimate"],
        waic_c              = waic_c$estimate["waic", "Estimate"],
        waic_diff           = delta_waic,
        waic_se_diff        = se_delta_w
      )
      print(to_be_returned[1,'elpd_diff'])
      to_be_returned
    })
  }
} else {
  run_mode <- "all"
  years_to_run <- available_years
  message("Full re-run mode: computing all seasons: ", min(years_to_run), "-", max(years_to_run),
          " (n=", length(years_to_run), ")")
  results_new <- purrr::map_dfr(years_to_run, function(yr) {
  
  dat <- tennis_years[[as.character(yr)]]
  Y   <- dat$Y_ij       # wins
  N   <- dat$N_ij       # matches played
  
  ############################################################
  ## 3a.  Fit both models
  ############################################################
  
  fit_simple  <- gibbs_bt_simple(w_ij = Y,
                                 T_iter = 10000, T_burn = 1500,
                                 verbose = T)
  
  fit_cluster <- gibbs_bt_sbm(w_ij = Y, T_iter = 10000, T_burn = 1500,
                              prior = "GN", a = 1.5,
                              gamma_GN = 0.8,
                              verbose = T)
  
  ############################################################
  ## 3b.  Build point‑wise log‑likelihood matrices
  ############################################################
  
  
  lambda_simple <- fit_simple$lambda_samples
  if (is.list(lambda_simple)) lambda_simple <- do.call(rbind, lambda_simple)
  ll_simple <- make_bt_simple_loo(w_ij = Y, lambda_samples = lambda_simple)

  lambda_cluster <- fit_cluster$lambda_samples
  if (is.list(lambda_cluster)) lambda_cluster <- do.call(rbind, lambda_cluster)
  ll_cluster <- make_bt_cluster_loo(
    w_ij = Y,
    lambda_samples = lambda_cluster,
    x_samples = fit_cluster$x_samples
  )
  
  ############################################################
  ## 3c.  LOO, ELPD differences, stacking weights
  ############################################################
  loo_s <- loo::loo(ll_simple$ll)
  loo_c <- loo::loo(ll_cluster$ll)
  
  # elpd difference (cluster ‑ simple) and its SE
  elpd_diff <- loo_c$estimates["elpd_loo","Estimate"] -
    loo_s$estimates["elpd_loo","Estimate"]
  
  se_diff   <- sqrt(loo_c$estimates["elpd_loo","SE"]^2 +
                      loo_s$estimates["elpd_loo","SE"]^2)/2
  
  
  ## --- WAIC --------------------------------------------------------------
  waic_s <- loo::waic(ll_simple$ll)
  waic_c <- loo::waic(ll_cluster$ll)
  
  delta_waic   <- waic_s$estimate["waic", "Estimate"] -
    waic_c$estimate["waic", "Estimate"]      # >0 ⇒ cluster better
  se_delta_w   <- sqrt(sum(waic_s$estimate["waic", "SE"]^2,
                           waic_c$estimate["waic", "SE"]^2))/2
  
  ############################################################
  ## 3e.  Collect and return one row
  ############################################################
  to_be_returned = tibble::tibble(
    season              = as.integer(yr),
    elpd_loo_simple     = loo_s$estimates["elpd_loo","Estimate"],
    elpd_loo_cluster    = loo_c$estimates["elpd_loo","Estimate"],
    elpd_diff           = elpd_diff,
    se_diff             = se_diff,
    waic_s              = waic_s$estimate["waic", "Estimate"],
    waic_c              = waic_c$estimate["waic", "Estimate"],
    waic_diff           = delta_waic,
    waic_se_diff        = se_delta_w
  )
  print(to_be_returned[1,'elpd_diff'])
  to_be_returned
})
}

results_all <- NULL
if (identical(run_mode, "append_missing_last_n")) {
  results_all <- dplyr::bind_rows(
    if (!is.null(existing_results)) existing_results else NULL,
    results_new
  )
} else {
  results_all <- results_new
}
results_all <- results_all %>%
  dplyr::mutate(season = as.integer(.data$season)) %>%
  dplyr::arrange(.data$season) %>%
  dplyr::distinct(.data$season, .keep_all = TRUE)

######################################################################
## 4.  Quick plots  ––––––––––––––––––––––––––––––––––––– ##
######################################################################

utils::write.csv(x = results_all, results_csv, row.names = FALSE)

results_all <- utils::read.csv(results_csv)
## 4a.  ΔELPD (cluster – simple)
ggplot(results_all, aes(season, elpd_diff)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  labs(title = "LOO comparison by season",
       x     = "Season (year)",
       y     = expression(italic(ELPD)[cluster] - italic(ELPD)[simple])) +
  theme_minimal()

## 4b.  BIC difference  (simple – cluster)
# ggplot(results, aes(season, bic_diff)) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_line() +
#   geom_point() +
#   labs(title = "BIC comparison by season",
#        x     = "Season (year)",
#        y     = expression(italic(BIC)[simple] - italic(BIC)[cluster])) +
#   theme_minimal()
results_plot_df <- results_all %>%
  dplyr::mutate(
    season_f = factor(.data$season, levels = sort(unique(.data$season))),
    ci_low = .data$elpd_diff - .data$se_diff,
    ci_high = .data$elpd_diff +  .data$se_diff
  )

DELPD_plot <- ggplot(results_plot_df, aes(x = season_f)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey35") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, group = 1),
              fill = "grey60", alpha = 0.35, inherit.aes = TRUE) +
  geom_line(aes(y = elpd_diff, group = 1), color = "grey10", linewidth = 0.7) +
  geom_point(aes(y = elpd_diff), color = "grey10", size = 1.6) +
  labs(x = "Season", y = "Delta ELPD (cluster - simple)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ggplot(results, aes(season, bic_diff)) + 
#   geom_col() + geom_hline(yintercept=0, lty=2) + 
#   labs(y="ΔBIC (simple − cluster)")


ggsave(filename = out_plot, plot = DELPD_plot, width = 13, height = 5)
ggsave(filename = out_plot_pdf, plot = DELPD_plot, width = 13, height = 5)


by_season_tab <- results_all %>%
  dplyr::mutate(
    ci_low = .data$elpd_diff - 1.96 * .data$se_diff,
    ci_high = .data$elpd_diff + 1.96 * .data$se_diff
  ) %>%
  dplyr::select(season, elpd_diff, se_diff, ci_low, ci_high, waic_diff, waic_se_diff)

tab1 <- by_season_tab %>%
  kableExtra::kable(format = "latex", digits = 2, booktabs = TRUE,
                    caption = "Model comparison by season: Delta ELPD (cluster - simple) with +/- 1.96*SE interval") %>%
  kableExtra::kable_styling(latex_options = c("striped", "hold_position"), full_width = FALSE)
writeLines(tab1, out_table_by_season)

tab2 <- results_all %>% summarise(across(c(elpd_diff), 
                             list(min=min, median=median, mean=mean, max=max)), 
                      prop_ELPD_gt_SE = mean(elpd_diff>se_diff))%>%
  kableExtra::kable(format = 'latex',digits = 2, booktabs = TRUE) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"), full_width = FALSE)
writeLines(tab2, out_table_summary)


sum(results_all$elpd_diff > results_all$se_diff/2) / nrow(results_all)
median(results_all$se_diff/2)
