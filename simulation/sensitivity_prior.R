# ---------------------------
# BT–SBM simulation + P(new)
# ---------------------------
set.seed(123)

# Logistic win-prob from lambdas
lambda_to_theta <- function(lam) {
  outer(lam, lam, function(a, b) a / (a + b))
}

simulate_bt_sbm <- function(n = 60, K_true = 6, ratio = 2.5,
                            mean_matches = 20, p_edge = 0.3) {
  # Balanced labels (approx)
  z <- rep(seq_len(K_true), length.out = n)
  z <- sample(z)
  # Geometric block strengths, GM = 1
  lam_blocks <- ratio^(0:(K_true - 1))
  lam_blocks <- lam_blocks / exp(mean(log(lam_blocks)))
  lam_item <- lam_blocks[z]
  # Match topology
  N <- matrix(0L, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    if (runif(1) < p_edge) {
      nij <- rpois(1, mean_matches)
      N[i, j] <- N[j, i] <- nij
    }
  }
  # Outcomes
  W <- matrix(0L, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    nij <- N[i, j]
    if (nij > 0L) {
      p <- lam_item[i] / (lam_item[i] + lam_item[j])
      wij <- rbinom(1, nij, p)
      W[i, j] <- wij
      W[j, i] <- nij - wij
    }
  }
  # Gamma augmentation draws Z; only row sums matter
  Z <- matrix(0, n, n)
  Z_row <- numeric(n); w_row <- rowSums(W)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    nij <- N[i, j]
    if (nij > 0L) {
      rate <- lam_item[i] + lam_item[j]
      zdraw <- rgamma(1, shape = nij, rate = rate) # scale=1/rate
      Z[i, j] <- Z[j, i] <- zdraw
      Z_row[i] <- Z_row[i] + zdraw
      Z_row[j] <- Z_row[j] + zdraw
    }
  }
  list(W = W, N = N, Z = Z, Z_row = Z_row, w_row = w_row,
       z = z, lam_item = lam_item, lam_blocks = lam_blocks)
}

# Anchored collapsed new score (increment rel. to empty cluster)
log_new_anchored <- function(wi, Zi, a, b) {
  lgamma(a + wi) - (a + wi) * log(b + Zi)
}

# One-step P(new) vs best existing competitor (hybrid rule)
p_new_for_item <- function(wi, Zi, lam_exist, a, optionA = F) {
  if (optionA) {
    b <- exp(digamma(a))  # scale alignment
    # center occupied lambdas to GM = 1
    lam <- lam_exist / exp(mean(log(lam_exist)))
  } else {
    b <- 1; lam <- lam_exist
  }
  # plug-in scores (existing); take best competitor
  logp_exist <- wi * log(lam) - lam * Zi
  best_exist <- max(logp_exist)
  # anchored collapsed new score
  logp_new <- log_new_anchored(wi, Zi, a, b)
  # 2-class softmax prob
  m <- max(best_exist, logp_new)
  p_new <- exp(logp_new - m) / (exp(best_exist - m) + exp(logp_new - m))
  c(p_new = p_new, b_eff = b)
}

# ----- Run a sensitivity sweep over 'a' -----
library(dplyr)

sim <- simulate_bt_sbm(n = 60, K_true = 6, ratio = 2.5,
                       mean_matches = 20, p_edge = 0.3)

lam_exist <- sim$lam_blocks  # existing blocks' lambdas (GM=1 by construction)

a_grid <- c(0.5, 1, 2, 3, 4, 6, 8, 12)

res <- lapply(a_grid, function(a) {
  out <- sapply(seq_along(sim$w_row), function(i) {
    p_new_for_item(sim$w_row[i], sim$Z_row[i], lam_exist, a, optionA = F)
  })
  tibble(a = a,
         p_new_mean   = mean(out["p_new", ]),
         p_new_median = median(out["p_new", ]),
         p_new_q90    = quantile(out["p_new", ], 0.9),
         b_eff        = as.numeric(out["b_eff", 1]))
}) %>% bind_rows()

print(res)

# Optional: plot
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  res_long <- tidyr::pivot_longer(res, cols = starts_with("p_new_"),
                                  names_to = "stat", values_to = "value")
  plot1 = ggplot(res_long, aes(a, value, shape = stat,color = stat)) +
    geom_point() + geom_line() +
    labs(title = "Effect of 'a' on P(new) with b = 1",
         x = "a", y = "P(new)") +
    theme_minimal()
  plot1
  ggsave(filename = "./images/sensitivity_b1.png",plot1)
}




library(dplyr)
library(tidyr)
library(ggplot2)

# Build an item-level frame with observed wi and Zi
items <- tibble(
  i   = seq_along(sim$w_row),
  wi  = as.integer(sim$w_row),
  Zi  = sim$Z_row
) %>%
  mutate(
    Z_bin_id = ntile(Zi, 3),
    Z_bin    = factor(Z_bin_id, labels = c("low z_i", "mid z_i", "high z_i"))
  )

# Compute P(new) for each item across the a-grid
res_item <- lapply(a_grid, function(a) {
  items %>%
    rowwise() %>%
    mutate(p_new = p_new_for_item(wi, Zi, lam_exist, a, optionA = T)[["p_new"]]) %>%
    ungroup() %>%
    mutate(a = a)
}) %>% bind_rows()

# OPTIONAL: pick representative wi values within each Z panel to avoid spaghetti
# Here we take 6 quantiles of wi per Z_bin and keep unique integer values
rep_w_by_z <- res_item %>%
  distinct(Z_bin, wi) %>%
  group_by(Z_bin) %>%
  summarize(
    wi_rep = unique(as.integer(quantile(wi, probs = c(.05, .20, .4, .6, .8, .95), type = 1))),
    .groups = "drop"
  ) %>%
  unnest_longer(wi_rep) %>%
  distinct(Z_bin, wi = wi_rep)

# Aggregate across items that share the same (Z_bin, wi) at each a
plot_df <- res_item %>%
  inner_join(rep_w_by_z, by = c("Z_bin", "wi")) %>%
  group_by(Z_bin, wi, a) %>%
  summarize(p_new = mean(p_new), .groups = "drop") %>%
  mutate(wi_f = factor(wi))  # for a discrete legend

# Plot: lines by wi, facets by Z_bin
plot2 <- ggplot(plot_df, aes(x = a, y = p_new, color = wi_f, group = wi_f)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  facet_wrap(~ Z_bin, ncol = 3) +
  labs(
    x = "a",
    y = "P(new)",
    color = expression(w[i])
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(plot2)
# ggsave("./images/pnew_lines_by_wi_faceted_by_zibpsia.png", plot2, width = 9, height = 5.5, dpi = 300)



# Load required package
library(Deriv)

# Define the function to invert the trigamma numerically
inv_trigamma <- function(tau2) {
  f <- function(a) trigamma(a) - tau2
  uniroot(f, interval = c(1e-8, 100))$root
}

# Values of tau
tau <- c(0.25, 0.33, 0.50, 0.70, 1.00)

# Compute corresponding a values by solving ψ₁(a) = τ²
a <- sapply(tau, function(t) inv_trigamma(t^2))

# Create the table
tab <- data.frame(
  a = seq(6,1)
)%>%
  mutate(tau = round(sqrt(trigamma(a)),2))%>%
  mutate(sd_lambda = round(sqrt(a)/exp(digamma(a)),2))

kable(tab, format = "latex", digits = 3, booktabs = TRUE) %>%
  kable_styling(latex_options = c("striped", "hold_position"))






