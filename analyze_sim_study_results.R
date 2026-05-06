# Minimal post-processing for simulation study:
# 1) BT-SBM recovery table for K*=3,5,7
# 2) One-run baseline table: BT vs RCBTL vs BT-SBM
# 3) ARI plot
#
# This script does not run simulation. It only reads existing CSV outputs.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) == 1L) {
  script_path <- sub("^--file=", "", file_arg)
  script_path <- gsub("~\\+~", " ", script_path)
  script_dir <- dirname(normalizePath(script_path))
  setwd(script_dir)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(knitr)
  library(kableExtra)
  library(ggplot2)
})

res_dir <- "results"
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)
if (!dir.exists("images")) dir.create("images", recursive = TRUE)

sim_csv <- file.path(res_dir, "simulation_comparison_all_models.csv")
rcbtl_csv <- file.path(res_dir, "rcbtl_baseline_design1_K3_seed123", "comparison_pearce_ereshova_design1_K3_iters10000_runs1.csv")

missing_inputs <- c(sim_csv, rcbtl_csv)[!file.exists(c(sim_csv, rcbtl_csv))]
if (length(missing_inputs) > 0) {
  stop("Missing required input files: ", paste(missing_inputs, collapse = ", "))
}

sim <- readr::read_csv(sim_csv, show_col_types = FALSE)
base <- readr::read_csv(rcbtl_csv, show_col_types = FALSE)

if (!"mae_rel" %in% names(sim)) {
  sim$mae_rel <- NA_real_
}

# -----------------------------
# Table 1: BT-SBM recovery by K*=3,5,7
# -----------------------------
recovery <- sim %>%
  filter(model == "BT-SBM", K_true %in% c(3, 5, 7)) %>%
  group_by(K_true) %>%
  summarise(
    MAE_rel = mean(mae_rel, na.rm = TRUE),
    Spearman = mean(spearman, na.rm = TRUE),
    ARI = mean(ari_minVI, na.rm = TRUE),
    VI = mean(vi_minVI, na.rm = TRUE),
    K_hat = round(mean(K_median, na.rm = TRUE), 0),
    K_ci_low = round(mean(HPD_low, na.rm = TRUE), 0),
    K_ci_high = round(mean(HPD_high, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  mutate(
    `K*` = K_true,
    `K_hat (95% CI)` = paste0(K_hat, " [", K_ci_low, ",", K_ci_high, "]")
  ) %>%
  select(`K*`, MAE_rel, Spearman, ARI, VI, `K_hat (95% CI)`) %>%
  arrange(`K*`) %>%
  mutate(across(c(MAE_rel, Spearman, ARI, VI), ~ round(.x, 3)))

recovery_out <- recovery %>%
  mutate(MAE_rel = ifelse(is.na(MAE_rel), "--", sprintf("%.3f", MAE_rel)))

write_csv(recovery_out, file.path(res_dir, "table_btsbm_recovery_k357.csv"))

tab_recovery <- knitr::kable(
  recovery_out,
  format = "latex",
  booktabs = TRUE,
  caption = "BT--SBM recovery by target number of clusters K* (averaged across designs).",
  align = "cccccc"
) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"))

cat(tab_recovery, file = file.path(res_dir, "table_btsbm_recovery_k357.tex"))

# -----------------------------
# Table 2: One-run baseline comparison (BT / RCBTL / BT-SBM)
# -----------------------------
model_order <- c("BT", "BT-SBM", "RCBTL")

baseline <- base %>%
  filter(model %in% model_order) %>%
  mutate(
    model = factor(model, levels = model_order),
    K_hat = ifelse(is.na(K_mean), "--", as.character(round(K_mean, 0))),
    K_CI = ifelse(is.na(K_ci_low), "--", paste0("[", round(K_ci_low, 0), ",", round(K_ci_high, 0), "]")),
    `K_hat (95% CI)` = ifelse(K_hat == "--", "--", paste(K_hat, K_CI)),
    Spearman = round(spearman, 3),
    ARI = ifelse(is.na(ari), "--", sprintf("%.3f", ari)),
    VI = ifelse(is.na(vi), "--", sprintf("%.3f", vi)),
    `Time (s)` = round(elapsed, 1)
  ) %>%
  transmute(
    Model = as.character(model),
    `K_hat (95% CI)` = `K_hat (95% CI)`,
    Spearman = Spearman,
    ARI = ARI,
    VI = VI,
    `Time (s)` = `Time (s)`
  ) %>%
  arrange(factor(Model, levels = model_order))

write_csv(baseline, file.path(res_dir, "table_rcbtl_baseline_one_run.csv"))

tab_baseline <- knitr::kable(
  baseline,
  format = "latex",
  booktabs = TRUE,
  caption = "Baseline comparison on one synthetic dataset (design1, K*=3, seed=123).",
  align = "lccccc"
) %>%
  kableExtra::kable_styling(latex_options = c("hold_position"))

cat(tab_baseline, file = file.path(res_dir, "table_rcbtl_baseline_one_run.tex"))

# -----------------------------
# ARI plot (BT-SBM recovery only)
# -----------------------------
ari_plot <- ggplot(recovery, aes(x = factor(`K*`), y = ARI, group = 1)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  theme_bw(base_size = 11) +
  labs(
    title = "BT-SBM ARI recovery by K*",
    x = "K*",
    y = "ARI"
  )

ggsave(
  filename = file.path("images", "ARI_plot.png"),
  plot = ari_plot,
  width = 6,
  height = 4,
  dpi = 300
)

message("Wrote:")
message("- ", file.path(res_dir, "table_btsbm_recovery_k357.csv"))
message("- ", file.path(res_dir, "table_btsbm_recovery_k357.tex"))
message("- ", file.path(res_dir, "table_rcbtl_baseline_one_run.csv"))
message("- ", file.path(res_dir, "table_rcbtl_baseline_one_run.tex"))
message("- images/ARI_plot.png")
