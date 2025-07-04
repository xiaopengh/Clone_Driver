library(data.table)
library(MASS)  # For negative binomial
library(pROC)  # For ROC curves
library(dplyr) # For pipe operators

# Fit model with an offset 
v <- clonal_summary_binom$dN_to_dS_clonal
w <- clonal_summary_binom$obs_dN_subclonal
offset_term <- log(clonal_summary_binom$non_syn_count * (clonal_summary_binom$obs_dS_subclonal / clonal_summary_binom$syn_count) )

# dN_dS models 
dNdS_models <- list()

# Poisson regression
dNdS_models$poisson <- glm(w ~ v + offset(offset_term), family = poisson(link = "log"))

# Negative binomial regression
dNdS_models$negbin <- glm.nb(w ~ v + offset(offset_term), control = glm.control(maxit = 200))

# Calculate response metrics for dN_dS models 
clonal_summary_binom[, `:=`(
  Predicted_Poisson_dNdS = predict(dNdS_models$poisson, type = "response"),
  Predicted_NegBin_dNdS = predict(dNdS_models$negbin, type = "response")
)]

clonal_summary_binom[, `:=`(
  Left_Prob_dNdS = pnbinom(y, size = dNdS_models$negbin$theta, mu = Predicted_NegBin_dNdS, lower.tail = TRUE),
  Right_Prob_dNdS = pnbinom(y, size = dNdS_models$negbin$theta, mu = Predicted_NegBin_dNdS, lower.tail = FALSE)
)]

clonal_summary_binom[, `:=`(
  p_val_dNdS = 2 * pmin(Left_Prob_dNdS, Right_Prob_dNdS),
  Pearson_Residual_dNdS = residuals(dNdS_models$negbin, type = "pearson"),
  Deviance_Residual_dNdS = residuals(dNdS_models$negbin, type = "deviance")
)]

# ROC analysis for cancer gene validation using dN dS method
roc_list_dNdS <- list(
  `Pearson Residuals (AUC: 0.XX)` = roc(clonal_summary_binom$Cancer_Gene, abs(clonal_summary_binom$Pearson_Residual_dNdS)),
  `Deviance Residuals (AUC: 0.YY)` = roc(clonal_summary_binom$Cancer_Gene, abs(clonal_summary_binom$Deviance_Residual_dNdS)),
  `p-value (AUC: 0.ZZ)` = roc(clonal_summary_binom$Cancer_Gene, abs(clonal_summary_binom$p_val_dNdS))
)

# Populate the actual AUC values into the list names
names(roc_list_dNdS)[1] <- paste0("Pearson Residuals (AUC: ", round(auc(roc_list_dNdS$`Pearson Residuals (AUC: 0.XX)`), 3), ")")
names(roc_list_dNdS)[2] <- paste0("Deviance Residuals (AUC: ", round(auc(roc_list_dNdS$`Deviance Residuals (AUC: 0.YY)`), 3), ")")
names(roc_list_dNdS)[3] <- paste0("p-value (AUC: ", round(auc(roc_list_dNdS$`p-value (AUC: 0.ZZ)`), 3), ")")

# Set an empirical alpha
alpha_dNdS <- 0.05

# Calculate upper and lower quantiles based on selected alpha
clonal_summary_binom[, `:=`(
  upper_quantile_dNdS = qnbinom(alpha_dNdS / 2, size = dNdS_models$negbin$theta, mu = Predicted_NegBin_dNdS, lower.tail = FALSE),
  lower_quantile_dNdS = qnbinom(alpha_dNdS / 2, size = dNdS_models$negbin$theta, mu = Predicted_NegBin_dNdS, lower.tail = TRUE)
)]

# Get outliers and their intersections
cgc_summary_binom <- clonal_summary_binom[Hugo_Symbol %in% cgc_data[["Gene Symbol"]]]
outliers <- clonal_summary_binom[p_val < 0.05]
outliers_dNdS <- clonal_summary_binom[p_val_dNdS < 0.05]
# outliers_intersect <- outliers_dNdS[Transcript_ID %in% outliers$Transcript_ID]
outliers_intersect <- cgc_summary_binom[Transcript_ID %in% outliers$Transcript_ID]
outliers_intersect_dNdS <- cgc_summary_binom[Transcript_ID %in% outliers_dNdS$Transcript_ID]

# Left and right tail p value for hypergeometric test, using the number of outliers in dN/dS method as K
# phyper(
#   nrow(outliers_intersect), 
#   nrow(outliers_dNdS),
#   nrow(clonal_summary_binom) - nrow(outliers_dNdS),
#   nrow(outliers),
#   lower.tail = TRUE # ALTER: Gene find in our method too little 
# )
# phyper(
#   nrow(outliers_intersect), # k_0 
#   nrow(outliers_dNdS), # K
#   nrow(clonal_summary_binom) - nrow(outliers_dNdS), # N - K
#   nrow(outliers), # n 
#   lower.tail = FALSE # ALTER: Gene find in our method is too much
# )


# Right tail p values for hypergeometric test, using the number ture cancer genes as K, testing for our method and dN/dS method
phyper(
  nrow(outliers_intersect), 
  nrow(cgc_summary_binom),
  nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
  nrow(outliers),
  lower.tail = FALSE # ALTER: Gene find in our method is much more than finding cancer genes by random chance
)

phyper(
  nrow(outliers_intersect_dNdS),
  nrow(cgc_summary_binom),
  nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
  nrow(outliers_dNdS),
  lower.tail = FALSE
) 
