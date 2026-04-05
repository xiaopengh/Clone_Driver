# Produce clonal summary for poisson test
intersected_dt <- as.data.table(intersected_gr)
intersected_dt[, coverage := (t_alt_count + t_ref_count)]
# intersected_dt <- intersected_dt[complete.cases(intersected_dt)]

intersected_dt[, binom_pval := mapply(function(x, n, p_hyp) {
  binom.test(x, n, p = p_hyp, alternative = "less")$p.value
}, x = t_alt_count, n = coverage, p_hyp = purity * 0.5)]

# add is_subclonal column based on Variant_Classification
intersected_dt[, is_subclonal := (binom_pval < 0.05)]
intersected_dt[, name := NULL] # Remove 'name' column if it exists

setnames(intersected_dt, old = c("name.1"), new = c("name"))

clonal_summary_binom <- intersected_dt[, .(
  clonal_count = sum(!is_subclonal),
  subclonal_count = sum(is_subclonal)
), by = c("feature", "name")]

# Estimate a ratio of overall clonal to subclonal counts for poisson test
ratio <- sum(clonal_summary_binom$subclonal_count) / sum(clonal_summary_binom$clonal_count)

# Calculate p-values for clonal and subclonal using exact conditional binomial test
clonal_summary_binom[, pval_poissonTest := mapply(function(x, y) {
  poisson.test(c(x, y), r = ratio)$p.value
}, subclonal_count, clonal_count)]

# Multiple testing correction
clonal_summary_binom[, pval_poissonTest_adj := p.adjust(pval_poissonTest, method = "fdr")]

# Get outliers
outliers <- clonal_summary_binom[pval_poissonTest< 0.05]
outliers_adj <- clonal_summary_binom[pval_poissonTest_adj < 0.05]

# ==============================================================================

ggplot(data = clonal_summary_binom) + 
  stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count), bins = 200) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  xlim(0, 600) + ylim(0, 200) +
  # scale_color_manual(name = "Model", values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))
  
ggplot(data = outliers) + 
  stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count), bins = 200) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  xlim(0, 600) + ylim(0, 200) +
  # scale_color_manual(name = "Model", values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

ggplot(data = outliers_adj) + 
  stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count), bins = 200) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  xlim(0, 600) + ylim(0, 200) +
  # scale_color_manual(name = "Model", values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

# ==============================================================================

ggplot(data = clonal_summary_binom) + 
  stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count),
              alpha = 0.2, bins = 200) +
  stat_bin_2d(data = outliers, 
              mapping = aes(x = clonal_count, y = subclonal_count), 
              alpha = 0.6, bins = 200) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis", direction = -1, name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  xlim(0, 600) + ylim(0, 200) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))


ggplot(data = clonal_summary_binom) + 
  stat_bin_2d(mapping = aes(x = clonal_count, y = subclonal_count),
              alpha = 1, bins = 200) +
  stat_bin_2d(data = outliers, 
              mapping = aes(x = clonal_count, y = subclonal_count), 
              alpha = 0.6, bins = 200) +
  geom_point(data = outliers_adj, mapping = aes(x = clonal_count, y = subclonal_count),
             shape = 3, color = "blue") +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis", direction = -1, name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  xlim(0, 600) + ylim(0, 200) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))






