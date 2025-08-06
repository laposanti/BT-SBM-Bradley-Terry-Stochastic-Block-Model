library(ggplot2)
library(dplyr)

# Parameters
  # fixed latent variable sum
b <- 1          # fixed Gamma rate

# Function to compute log marginal contribution for new cluster
new_cluster_integral_log <- function(a_current, w_i_val, Z_i_val, b) {
  log_top <- (a_current * log(b)) + lgamma(a_current + w_i_val)
  log_bot <- lgamma(a_current) + (a_current + w_i_val) * log(b + Z_i_val)
  log_top - log_bot
}

# Grid of a values and w_i values
a_grid <- seq(0.1, 5, length.out = 300)
w_i_vals <- c(0)
Z_i_vals <- c(0)


# Create dataframe
df_plot <- expand.grid(a = a_grid, w_i = w_i_vals, Z = Z_i_vals) %>%
  rowwise() %>%
  mutate(log_prob = new_cluster_integral_log(a, w_i, Z, b)) %>%
  mutate(Z_label = paste0("Z[i] == ", Z)) %>%  # <- This will be parsed
  ungroup()

# Plot
ggplot(df_plot, aes(x = a, y = log_prob, color = factor(w_i))) +
  geom_line(size = 1.1) +
  scale_color_brewer(palette = "Set1", name = expression(w[i])) +
  labs(
    title = "Sensitivity of New-Cluster Probability to Prior Shape a",
    subtitle = paste0("Rate parameter b = ", b),
    x = expression(a),
    y = expression(log ~ p(W ~ "|" ~ x[i] == "new"))
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom") +
  facet_wrap(~Z_label, labeller = label_parsed)
