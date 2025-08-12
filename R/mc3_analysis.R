library(data.table)
library(MASS)  # For negative binomial
library(pROC)  # For ROC curves
library(dplyr) # For pipe operators

# Create summary_dt
summary_dt <- rbindlist(summaries, use.names = TRUE, fill = TRUE, idcol = "Chromosome")

# Read mc3 dataset which contains the observed mutations.
mc3_raw <- fread("../data/mc3.maf", sep = "\t")
mc3 <- mc3_raw[, .(Hugo_Symbol, t_alt_count, t_ref_count, Segment_Mean, Purity, Variant_Classification, Transcript_ID)]
rm(mc3_raw)

# Perform stat tests on mc3 dataset and build the summary table by txid
mc3_summary_dt <- build_mc3_summary_bytx(mc3)

clonal_summary_binom <- merge.data.table(mc3_summary_dt[, c("Transcript_ID", "Hugo_Symbol")], 
                                         summary_dt, by.x = "Transcript_ID", by.y = "maf_txid") %>% na.omit()

clonal_summary_binom <- clonal_summary_binom[, c("Clonal_Count", "Subclonal_Count") :=
                                               .(obs_dS_clonal + obs_dN_clonal, 
                                                 obs_dS_subclonal + obs_dN_subclonal)]

# Fit models
x <- clonal_summary_binom$Clonal_Count
y <- clonal_summary_binom$Subclonal_Count

# Poisson regression
poisson_model <- glm(y ~ x, family = poisson(link = "log"))

# Negative binomial regression
negbin_model <- glm.nb(y ~ x, control = glm.control(maxit = 200))

# Load cancer gene data
cgc_data <- fread("../data/cgc.tsv", sep = "\t")
tsv_hugo <- cgc_data$`Gene Symbol`

# Add cancer gene indicator
clonal_summary_binom[, Cancer_Gene := Hugo_Symbol %in% tsv_hugo]

# Model diagnostics and overdispersion testing
cat("=== Model Diagnostics ===\n")
cat("Poisson Model AIC:", AIC(poisson_model), "\n")
cat("Negative Binomial Model AIC:", AIC(negbin_model), "\n")
cat("Theta (dispersion parameter):", negbin_model$theta, "\n")

# Calculate response metrics
clonal_summary_binom[, `:=`(
  Predicted_Poisson = predict(poisson_model, type = "response"),
  Predicted_NegBin = predict(negbin_model, type = "response")
)]

clonal_summary_binom[, `:=`(
  Left_Prob = pnbinom(y, size = negbin_model$theta, mu = Predicted_NegBin, lower.tail = TRUE),
  Right_Prob = pnbinom(y, size = negbin_model$theta, mu = Predicted_NegBin, lower.tail = FALSE)
)]

clonal_summary_binom[, `:=`(
  p_val = 2 * pmin(Left_Prob, Right_Prob),
  Pearson_Residual = residuals(negbin_model, type = "pearson"),
  Deviance_Residual = residuals(negbin_model, type = "deviance")
)]

# ROC analysis for cancer gene prediction
roc_pearson <- roc(clonal_summary_binom$Cancer_Gene, abs(clonal_summary_binom$Pearson_Residual))
roc_deviance <- roc(clonal_summary_binom$Cancer_Gene, abs(clonal_summary_binom$Deviance_Residual))
roc_p_val <- roc(clonal_summary_binom$Cancer_Gene, abs(clonal_summary_binom$p_val))

roc_list <- list(
  `Pearson Residuals (AUC: 0.XX)` = roc_pearson,
  `Deviance Residuals (AUC: 0.YY)` = roc_deviance,
  `p-value (AUC: 0.ZZ)` = roc_p_val
)

# Populate the actual AUC values into the list names
names(roc_list)[1] <- paste0("Pearson Residuals (AUC: ", round(auc(roc_pearson), 3), ")")
names(roc_list)[2] <- paste0("Deviance Residuals (AUC: ", round(auc(roc_deviance), 3), ")")
names(roc_list)[3] <- paste0("p-value (AUC: ", round(auc(roc_p_val), 3), ")")

# Choosing a threshold p_val
all_coords <- coords(roc_p_val, seq(min(roc_pearson$predictor), max(roc_pearson$predictor), 0.01),
                     ret = c("threshold", "sensitivity", "specificity", "precision", "recall"))

all_coords$F1_score <- 2 * (all_coords$precision * all_coords$recall) / (all_coords$precision + all_coords$recall)

# Set an optimal alpha base on maximum precision
alpha <- all_coords[which.max(all_coords$precision), ]$threshold
# Set a normal alpha 
alpha <- 0.05

# Calculate upper and lower quantiles based on selected alpha
clonal_summary_binom[, `:=`(
  upper_quantile = qnbinom(alpha / 2, size = negbin_model$theta, mu = Predicted_NegBin, lower.tail = FALSE),
  lower_quantile = qnbinom(alpha / 2, size = negbin_model$theta, mu = Predicted_NegBin, lower.tail = TRUE)
)]