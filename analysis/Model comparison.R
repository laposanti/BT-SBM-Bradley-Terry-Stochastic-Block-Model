
setwd("/Users/lapo_santi/Desktop/Nial/Bterry/BT-SBM/")
######################################################################
## 0.  Packages  –––––––––––––––––––––––––––––––––––––––––––––––––– ##
######################################################################
pkgs <- c("loo", "purrr", "dplyr", "tibble", "ggplot2", "tidyr")

lapply(pkgs, require, character.only = TRUE)

set.seed(2025)       # reproducibility

######################################################################
## 1.  Load the data  ––––––––––––––––––––––––––––––––––––––––––––– ##
######################################################################
tennis_years <- readRDS("./data/2000_2022_data.rds")

source("./main/mcmc_bt_sbm.R")
source("main/utils.R")
######################################################################
## 2.  Helper: a thin BIC wrapper  –––––––––––––––––––––––––––––––– ##
######################################################################
bic_from_ll <- function(loglik, k, n) {
  # BIC = k*log(n) - 2*logLik
  k*log(n) - 2*loglik
}


######################################################################
## 3.  Main loop over seasons  –––––––––––––––––––––––––––––––––––– ##
######################################################################
results <- purrr::map_dfr(names(tennis_years), function(yr) {
  
  dat <- tennis_years[[yr]]
  Y   <- dat$Y_ij       # wins
  N   <- dat$N_ij       # matches played
  
  ############################################################
  ## 3a.  Fit both models
  ############################################################
  fit_simple  <- gibbs_bt_simple(Y, N,
                                 n_iter = 4000, burnin = 1000,
                                 verbose = T)
  
  fit_cluster <- gibbs_bt_sbm(Y, N,
                                  n_iter = 4000, burnin = 1000,
                                  prior = "GN", a = 0.5,b = 1,
                                  gamma_GN = 0.8,
                                  verbose = T)
  
  ############################################################
  ## 3b.  Build point‑wise log‑likelihood matrices
  ############################################################
  ll_simple  <- make_bt_simple_loo(Y, N,  fit_simple$lambda_samples)
  ll_cluster <- make_bt_cluster_loo(Y, N,
                                    fit_cluster$lambda_samples,
                                    fit_cluster$x_samples)
  
  ############################################################
  ## 3c.  LOO, ELPD differences, stacking weights
  ############################################################
  loo_s <- loo::loo(ll_simple$ll)
  loo_c <- loo::loo(ll_cluster$ll)
  
  # elpd difference (cluster ‑ simple) and its SE
  elpd_diff <- loo_c$estimates["elpd_loo","Estimate"] -
    loo_s$estimates["elpd_loo","Estimate"]
  
  se_diff   <- sqrt(loo_c$estimates["elpd_loo","SE"]^2 +
                      loo_s$estimates["elpd_loo","SE"]^2)
  
  
  ## --- WAIC --------------------------------------------------------------
  waic_s <- loo::waic(ll_simple$ll)
  waic_c <- loo::waic(ll_cluster$ll)
  
  delta_waic   <- waic_s$estimate["waic", "Estimate"] -
    waic_c$estimate["waic", "Estimate"]      # >0 ⇒ cluster better
  se_delta_w   <- sqrt(sum(waic_s$estimate["waic", "SE"]^2,
                           waic_c$estimate["waic", "SE"]^2))
  
  ############################################################
  ## 3d.  Pragmatic BIC
  ############################################################
  # k   – effective #parameters (crude: #λ’s)
  # n   – #observed comparison pairs
  k_simple   <- ncol(Y)                              # one λ per player
  k_cluster  <- max(fit_cluster$x_samples, na.rm = TRUE) # one λ per *cluster*
  n_pairs    <- sum(upper.tri(N) & N > 0)
  
  # Monte‑Carlo estimate of log‑likelihood at posterior mean
  llk_simple   <- sum(colMeans(ll_simple$ll))
  llk_cluster  <- sum(colMeans(ll_cluster$ll))
  
  bic_simple  <- bic_from_ll(llk_simple,  k_simple,  n_pairs)
  bic_cluster <- bic_from_ll(llk_cluster, k_cluster, n_pairs)
  
  ############################################################
  ## 3e.  Collect and return one row
  ############################################################
  to_be_returned = tibble::tibble(
    season              = as.integer(yr),
    elpd_loo_simple     = loo_s$estimates["elpd_loo","Estimate"],
    elpd_loo_cluster    = loo_c$estimates["elpd_loo","Estimate"],
    elpd_diff           = elpd_diff,
    se_diff             = se_diff,
    bic_simple          = bic_simple,
    bic_cluster         = bic_cluster,
    bic_diff            = bic_simple - bic_cluster,
    waic_s              = waic_s$estimate["waic", "Estimate"],
    waic_c              = waic_c$estimate["waic", "Estimate"],
    waic_diff           = delta_waic,
    waic_se_diff        = se_delta_w
  )
  print(to_be_returned[1,'elpd_diff'])
  to_be_returned
})

######################################################################
## 4.  Quick‑and‑dirty plots  ––––––––––––––––––––––––––––––––––––– ##
######################################################################

write.csv("/Users/lapo_santi/Desktop/results.csv")

## 4a.  ΔELPD (cluster – simple)
ggplot(results, aes(season, elpd_diff)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  labs(title = "LOO comparison by season",
       x     = "Season (year)",
       y     = expression(italic(ELPD)[cluster] - italic(ELPD)[simple])) +
  theme_minimal()

## 4b.  BIC difference  (simple – cluster)
ggplot(results, aes(season, bic_diff)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  labs(title = "BIC comparison by season",
       x     = "Season (year)",
       y     = expression(italic(BIC)[simple] - italic(BIC)[cluster])) +
  theme_minimal()


DELPD_plot = ggplot(results, aes(season, elpd_diff)) + 
  geom_hline(yintercept=0, lty=2) + 
  geom_ribbon(aes(ymin=elpd_diff-se_diff/2, ymax=elpd_diff+se_diff/2), fill="#2E8B57",alpha=.5) + 
  geom_line(color = "#FFD700", linewidth = 2) +         # tennis ball yellow
  geom_point(color = "#FFD700", size = 4) + 
  labs(y="ΔELPD (cluster − simple)")+
  theme_minimal()

ggplot(results, aes(season, bic_diff)) + 
  geom_col() + geom_hline(yintercept=0, lty=2) + 
  labs(y="ΔBIC (simple − cluster)")



ggsave(filename = "./images/DELPD_plot.png", DELPD_plot, width = 13, height=5)


results %>% summarise(across(c(elpd_diff), 
                             list(min=min, median=median, mean=mean, max=max)), 
                      prop_ELPD_gt_SE = mean(elpd_diff>se_diff))%>%
  kableExtra::kable(format = 'latex',digits = 2)


sum(results$elpd_diff>results$se_diff/2)/nrow(results)
median(results$se_diff/2)
