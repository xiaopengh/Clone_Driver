intersected_dt <- as.data.table(intersected_gr)

intersected_dt[, coverage := (t_alt_count + t_ref_count)]
intersected_dt <- intersected_dt[complete.cases(intersected_dt)]
intersected_dt[, binom_pval := mapply(function(x, n, p_hyp) {
  binom.test(x, n, p = p_hyp, alternative = "less")$p.value
}, x = t_alt_count, n = coverage, p_hyp = purity * 0.5)]

# add is_sysnonymous, is_subclonal column based on Variant_Classification
intersected_dt[, is_subclonal := (binom_pval < 0.05)]
intersected_dt[, name := NULL] # Remove 'name' column if it exists

setnames(intersected_dt, old = c("name.1"), new = c("name"))

clonal_summary_binom <- intersected_dt[, .(
  clonal_count = sum(!is_subclonal),
  subclonal_count = sum(is_subclonal)
), by = c("feature", "name")]


# Fit models
x <- clonal_summary_binom$clonal_count
y <- clonal_summary_binom$subclonal_count

# Negative binomial regression
negbin_model <- glm.nb(y ~ x, control = glm.control(maxit = 200))

# Calculate response metrics
clonal_summary_binom[, predicted := predict(negbin_model, type = "response")]

clonal_summary_binom[, `:=`(
  left_prob = pnbinom(y, size = negbin_model$theta, mu = predicted, lower.tail = TRUE),
  right_prob = pnbinom(y, size = negbin_model$theta, mu = predicted, lower.tail = FALSE)
)]

clonal_summary_binom[, p_val := 2 * pmin(left_prob, right_prob)]

# Set a normal alpha 
alpha <- 0.05

{
  theta <- negbin_model$theta
  mu_vec <- clonal_summary_binom$predicted
  
  upper_q <- qnbinom(alpha / 2, size = theta, mu = mu_vec, lower.tail = FALSE)
  lower_q <- qnbinom(alpha / 2, size = theta, mu = mu_vec, lower.tail = TRUE)
  
  clonal_summary_binom[, `:=`(
    upper_quantile = upper_q,
    lower_quantile = lower_q
  )]
  
  rm(theta, mu_vec, upper_q, lower_q) ; gc()
}

# ==============================================================================


g <- ggplot(data = clonal_summary_binom)

# negbin regression with clonal to subclonal plot
ord <- order(x)
g_reg <- g + stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count), bins = 200, geom = "tile", alpha = 0.8) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  # geom_line(mapping = aes(x = x[ord], y = Predicted_Poisson[ord], color = "Poisson"),
  #           linewidth = 1) +
  geom_line(mapping = aes(x = x[ord], y = predicted[ord], color = "Negative Binomial"),
            linewidth = 1) +
  xlim(0, 1000) + ylim(0, 350) +
  scale_color_manual(name = "Model", values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

g_reg

# negbin regression with clonal to subclonal plot and quantiles of upper and lower levels  
ord <- order(x)
g_classification <- g + stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count), 
                                                     bins = 200, geom = "tile", alpha = 0.8) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis",  # or "plasma", "inferno", "magma"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  geom_line(mapping = aes(x = x[ord], y = lower_quantile[ord], color = "Quantiles"),
            linewidth = 0.8) +
  geom_line(mapping = aes(x = x[ord], y = upper_quantile[ord], color = "Quantiles"),
            linewidth = 0.8) +
  geom_line(mapping = aes(x = x[ord], y = predicted[ord], color = "Negative Binomial"),
            linewidth = 1) +
  # geom_point(data = clonal_summary_binom[p_val <= alpha],
  #            mapping = aes(x = clonal_count, y = subclonal_count, shape = "Outliers by quantiles"), 
  #            size = 3, alpha = 0.3) +
  # geom_ribbon(aes(x = x[ord], ymin = lower_quantile[ord], ymax = upper_quantile[ord]), fill = "lightgreen", alpha = 0.2) +
  xlim(0, 1000) + ylim(0, 350) +
  scale_shape_manual(values = c("Outliers by quantiles" = 1)) +
  scale_color_manual(values = c("Quantiles" = "lightgreen", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

g_classification
