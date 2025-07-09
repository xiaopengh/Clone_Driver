library(ggplot2)
library(dplyr) # For pipline operators 
library(pROC)  # For ROC curves
library(data.table) 
library(patchwork)

reverse_offset <- function(response, offset_term) {
  exp(log(response) - offset_term)
  }

plots <- list()

plots$g <- ggplot(data = clonal_summary_binom)

# dN/dS clonal to subclonal plot
ord.dnds <- order(v)
plots$g_dnds <- plots$g + stat_bin_2d(mapping = aes(x = dN_to_dS_clonal, y = dN_to_dS_subclonal), bins = 200) +
  scale_fill_binned(type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  geom_point(x = 1, y = 1, color = "red", size = 3, shape = 17) + # Corrected line
  geom_line(mapping = aes(x = v[ord.dnds], color = "Poisson",
                          y = reverse_offset(Predicted_Poisson_dNdS, offset_term)[ord.dnds]),
            linewidth = 1) +
  geom_line(mapping = aes(x = v[ord.dnds], color = "Negative Binomial",
                          y = reverse_offset(Predicted_NegBin_dNdS, offset_term)[ord.dnds]),
            linewidth = 1) +
  xlim(0, 10) + ylim(0, 10) +
  scale_color_manual(name = "Model", values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

plots$g_dnds

# negbin regression with clonal to subclonal plot
ord <- order(x)
plots$g_reg <- plots$g + stat_bin_2d(mapping = aes(x = Clonal_Count, y = Subclonal_Count), bins = 200, geom = "tile", alpha = 0.8) +
  scale_fill_binned(breaks = c(0, 5, 10, 20, 50, 100),
                    type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  geom_line(mapping = aes(x = x[ord], y = Predicted_Poisson[ord], color = "Poisson"),
            linewidth = 1) +
  geom_line(mapping = aes(x = x[ord], y = Predicted_NegBin[ord], color = "Negative Binomial"),
            linewidth = 1) +
  xlim(0, 600) + ylim(0, 600) +
  scale_color_manual(name = "Model", values = c("Poisson" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

plots$g_reg

plots$g_dnds | plots$g_reg

# Residuals histograms
plots$h1 <- plots$g + geom_histogram(aes(x = Pearson_Residual, y = after_stat(density)), 
                                     bins = 200, fill = "lightblue", color = "white", lwd = 0.2) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(title = "Distribution of Pearson Residuals", x = "Residual Value", y = "Frequency") +
  theme_minimal()

plots$h2 <- plots$g + geom_histogram(aes(x = Deviance_Residual, y = after_stat(density)), 
                                     bins = 200, fill = "lightgreen", color = "white", lwd = 0.2) +
  coord_cartesian(xlim = c(-5, 5)) +
  labs(title = "Distribution of Deviance Residuals", x = "Residual Value", y = "Frequency") + 
  theme_minimal()

plots$h1 | plots$h2

# Plot the ROC curves
plot(roc_pearson,
     main = "ROC Curve for Cancer Gene Prediction",
     col = "blue",
     lwd = 2,
     print.auc = TRUE,
     legacy.axes = TRUE)

plot(roc_deviance,
     add = TRUE,         
     col = "red",        
     lwd = 2,
     print.auc = TRUE,   
     print.auc.y = 0.4)  # Adjust the vertical position of the second AUC value

plot(roc_p_val,
     add = TRUE,         
     col = "green",      
     lwd = 2,
     print.auc = TRUE,   
     print.auc.y = 0.3)  # Adjust the vertical position of the third AUC value

# plots$p_roc_pearson <- ggroc(roc_pearson, color = "blue", lwd = 1.2, legacy.axes = TRUE) +
#   geom_segment(mapping = aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", lwd = 0.8) +
#   labs(title = "Pearson", x = "specificity", subtitle = paste0("AUC: ", round(auc(roc_pearson), 3)))
# 
# plots$p_roc_deviance <- ggroc(roc_deviance, color = "blue", lwd = 1.2, legacy.axes = TRUE) +
#   geom_segment(mapping = aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", lwd = 0.8) +
#   labs(title = "Deviance", x = "specificity", subtitle = paste0("AUC: ", round(auc(roc_deviance), 3)))

plots$p_roc <- ggroc(roc_list, aes = "color", lwd = 0.5, legacy.axes = TRUE) +
  geom_segment(mapping = aes(x = 0, y = 0, xend = 1, yend = 1), color = "lightgray", lwd = 0.8) + 
  labs(title = "ROC Curve Comparison")

plots$F1_thres <- ggplot(data = all_coords) + 
  geom_line(mapping = aes(x = threshold, y = F1_score), lwd = 0.8, color = "red") + 
  coord_cartesian(xlim = c(0, 1))

plots$precision_recall <- ggplot(data = all_coords) + 
  geom_line(mapping = aes(x = recall, y = precision), lwd = 0.8, color = "blue") + 
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision")

plots$precision <- ggplot(data = all_coords) + 
  geom_line(mapping = aes(x = threshold, y = precision), lwd = 0.8, color = "blue") + 
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Precision vs Threshold", x = "Threshold", y = "Precision")

# negbin regression with clonal to subclonal plot and quantiles of upper and lower levels  
ord <- order(x)
plots$negbin_classification <- plots$g + stat_bin_2d(mapping = aes(x = Clonal_Count, y = Subclonal_Count), 
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
  geom_line(mapping = aes(x = x[ord], y = Predicted_NegBin[ord], color = "Negative Binomial"),
            linewidth = 1) +
  geom_point(data = outliers, mapping = aes(x = Clonal_Count, y = Subclonal_Count), shape = 1, size = 3, alpha = 0.3) +
  # geom_ribbon(aes(x = x[ord], ymin = lower_quantile[ord], ymax = upper_quantile[ord]), fill = "lightgreen", alpha = 0.2) +
  xlim(0, 600) + ylim(0, 600) +
  scale_color_manual(values = c("Quantiles" = "lightgreen", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

plots$negbin_classification | plots$g_reg | plots$g_dnds

# negbin regression with dN/dS clonal to subclonal plot and quantiles of upper and lower levels
ord.dnds <- order(v)
plots$negbin_classification_dnds <- plots$g + stat_bin_2d(mapping = aes(x = dN_to_dS_clonal, y = dN_to_dS_subclonal), bins = 200) +
  scale_fill_binned(type = "viridis",  # "gradient", "viridis"
                    direction = -1,
                    name = "Count",
                    guide = guide_coloursteps(show.limits = TRUE)) +
  geom_point(x = 1, y = 1, color = "cyan", size = 3, shape = 17) + # Corrected line
  geom_line(mapping = aes(x = v[ord.dnds], y = reverse_offset(lower_quantile_dNdS, offset_term)[ord.dnds], color = "Quantiles"),
            linewidth = 0.8) +
  geom_line(mapping = aes(x = v[ord.dnds], y = reverse_offset(upper_quantile_dNdS, offset_term)[ord.dnds], color = "Quantiles"),
            linewidth = 0.8) +
  geom_line(mapping = aes(x = v[ord.dnds], color = "Negative Binomial",
                          y = reverse_offset(Predicted_NegBin_dNdS, offset_term)[ord.dnds]),
            linewidth = 1) +
  geom_point(data = outliers_dNdS, mapping = aes(x = dN_to_dS_clonal, y = dN_to_dS_subclonal), shape = 1, size = 3, alpha = 0.3) +
  xlim(0, 10) + ylim(0, 10) +
  scale_color_manual(name = "Model", values = c("Quantiles" = "blue", "Negative Binomial" = "red")) +
  theme_light() + guides(fill = guide_coloursteps(order = 1), color = guide_legend(order = 2))

plots$negbin_classification_dnds | plots$negbin_classification

plots$test <- plots$g + 
  # geom_line(mapping = aes(x = v[ord.dnds], y = lower_quantile_dNdS[ord.dnds], color = "Quantiles"),
  #           linewidth = 0.8) +
  # geom_line(mapping = aes(x = v[ord.dnds], y = upper_quantile_dNdS[ord.dnds], color = "Quantiles"),
  #           linewidth = 0.8) +
  geom_line(mapping = aes(x = v[ord.dnds], color = "Predicted Raw",
                          y = Predicted_NegBin_dNdS[ord.dnds] / 10), linewidth = 0.05) +
  geom_line(mapping = aes(x = v[ord.dnds], color = "Predicted Trans",
                          y = reverse_offset(Predicted_NegBin_dNdS, offset_term)[ord.dnds]),
            linewidth = 1) +
  scale_color_manual(name = "Model", values = c("Quantiles" = "blue", "Predicted Raw" = "#ffb732", "Predicted Trans" = "#ff3153")) + 
  ylim(0, 50)

plots$test

plots$negbin_classification_true_val <- plots$negbin_classification +
  geom_point(data = clonal_summary_binom[Cancer_Gene == TRUE], 
             mapping = aes(x = Clonal_Count, y = Subclonal_Count), 
             color = "red", shape = 3, size = 3) +
  annotate("text", 
           x = 200, y = 400,  # choose appropriate coordinates
           label = test_res, 
           hjust = 0, size = 5, color = "black", fontface = "italic")

plots$negbin_classification_true_val

plots$negbin_classification_dnds_true_val <- plots$negbin_classification_dnds +
  geom_point(data = clonal_summary_binom[Cancer_Gene == TRUE], 
             mapping = aes(x = dN_to_dS_clonal, y = dN_to_dS_subclonal), 
             color = "red", shape = 3, size = 3) +
  annotate("text", 
           x = 3.75, y = 10,  # choose appropriate coordinates
           label = test_res_dnds, 
           hjust = 0, size = 5, color = "black", fontface = "italic")

plots$negbin_classification_dnds_true_val

# =======================================================================
# Save plots to files
ggsave("plots/g_reg.svg", plot = plots$g_reg, width = 12, height = 9)
ggsave("plots/g_dnds.svg", plot = plots$g_dnds, width = 12, height = 9)
ggsave("plots/negbin_classification.svg", plot = plots$negbin_classification, width = 12, height = 9)
ggsave("plots/negbin_classification_dnds.svg", plot = plots$negbin_classification_dnds, width = 12, height = 9)
ggsave("plots/p_roc.svg", plot = plots$p_roc, width = 12, height = 9)
ggsave("plots/negbin_classification_true_val.svg", plot = plots$negbin_classification_true_val, width = 12, height = 9)
ggsave("plots/negbin_classification_dnds_true_val.svg", plot = plots$negbin_classification_dnds_true_val, width = 12, height = 9)
# =======================================================================