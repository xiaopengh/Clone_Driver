# Get outliers and their intersections
cgc_summary_binom <- clonal_summary_binom[Cancer_Gene == TRUE]
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
test_res <- paste0(
  "Hypergeometric p = ", 
  phyper(
    nrow(outliers_intersect), 
    nrow(cgc_summary_binom),
    nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
    nrow(outliers),
    lower.tail = FALSE # ALTER: Gene find in our method is much more than finding cancer genes by random chance
  ) %>% signif(3)
)

test_res_dnds <- paste0(
  "Hypergeometric p = ",
  phyper(
    nrow(outliers_intersect_dNdS),
    nrow(cgc_summary_binom),
    nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
    nrow(outliers_dNdS),
    lower.tail = FALSE
  ) %>% signif(3)
)

# Calculate p-values for clonal and subclonal using exact conditional binomial test
clonal_summary_binom[, `:=`(
  p_val_poissonTest_clonal = poisson.test(c(obs_dN_clonal, obs_dS_clonal),
                                          c(non_syn_count, syn_count))$p.value,
  p_val_poissonTest_subclonal = poisson.test(c(obs_dN_subclonal, obs_dS_subclonal),
                                             c(non_syn_count, syn_count))$p.value
), by = 1:nrow(clonal_summary_binom)]

# Label significance based on p-values from previous tests
clonal_summary_binom[, significance := fifelse(p_val_poissonTest_clonal < 0.05 & p_val_poissonTest_subclonal < 0.05, "Both",
                                       fifelse(p_val_poissonTest_clonal < 0.05, "Clonal",
                                       fifelse(p_val_poissonTest_subclonal < 0.05, "Subclonal", NA_character_)))]

# Adjust p-values for multiple testing for p-values from negative binomial regression using Benjamini-Hochberg method
clonal_summary_binom[, `:=`(
  p_val_adjBH_negbintest = p.adjust(p_val, method = "BH")
)]

# Adjust p-values for multiple testing with Holm method
clonal_summary_binom[, `:=`(
  p_val_adjHolm_poissonTest_clonal = p.adjust(p_val_poissonTest_clonal, method = "holm"),
  p_val_adjHolm_poissonTest_subclonal = p.adjust(p_val_poissonTest_subclonal, method = "holm")
)]

# Combine adjusted p-values for clonal and subclonal tests by choosing the maximum (Conservative approach)
clonal_summary_binom[, p_val_dNdS_adjHolm := 
  max(p_val_adjHolm_poissonTest_clonal, p_val_adjHolm_poissonTest_subclonal), 
  by = 1:nrow(clonal_summary_binom)
]

# Combine adjusted p-values for clonal and subclonal tests by Fisher's method (Less conservative)
clonal_summary_binom[, p_val_dNdS_adjHolm_Fisher := 
  pchisq(
    -2 * (log(p_val_poissonTest_clonal) + log(p_val_poissonTest_subclonal)),
    df = 4,
    lower.tail = FALSE
  ), 
  by = 1:nrow(clonal_summary_binom)
]

outliers_adj <- clonal_summary_binom[p_val_adjBH_negbintest < 0.05] # Our method outliers
outliers_dNdS_adj <- clonal_summary_binom[p_val_dNdS_adjHolm < 0.05] # dN/dS method max fusion outliers
outliers_dNdS_adj_Fisher <- clonal_summary_binom[p_val_dNdS_adjHolm_Fisher < 0.05] # dN/dS method Fisher fusion outliers

# outliers_intersect <- outliers_dNdS[Transcript_ID %in% outliers$Transcript_ID]
outliers_intersect_adj <- cgc_summary_binom[Transcript_ID %in% outliers_adj$Transcript_ID]
outliers_intersect_dNdS_adj <- cgc_summary_binom[Transcript_ID %in% outliers_dNdS_adj$Transcript_ID]
outliers_intersect_dNdS_adj_Fisher <- cgc_summary_binom[Transcript_ID %in% outliers_dNdS_adj_Fisher$Transcript_ID]
outliers_intersect_our_dnds <- outliers_adj[Transcript_ID %in% outliers_dNdS_adj$Transcript_ID]
outliers_intersect_our_dnds_f <- outliers_adj[Transcript_ID %in% outliers_dNdS_adj_Fisher$Transcript_ID] 

test_res_adj <- paste0(
  "Hypergeometric p (adjusted for multiple testing) = ", 
  phyper(
    nrow(outliers_intersect_adj), 
    nrow(cgc_summary_binom),
    nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
    nrow(outliers_adj),
    lower.tail = FALSE
  ) %>% signif(3)
)

test_res_dnds_adj <- paste0(
  "Hypergeometric p (adjusted for multiple testing) = ",
  phyper(
    nrow(outliers_intersect_dNdS_adj),
    nrow(cgc_summary_binom),
    nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
    nrow(outliers_dNdS_adj),
    lower.tail = FALSE
  ) %>% signif(3)
)

test_res_dnds_adj_f <- paste0(
  "Hypergeometric p (adjusted for multiple testing Fisher) = ",
  phyper(
    nrow(outliers_intersect_dNdS_adj_Fisher),
    nrow(cgc_summary_binom),
    nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
    nrow(outliers_dNdS_adj_Fisher),
    lower.tail = FALSE
  ) %>% signif(3)
)

test_res_dnds_ours <- paste0(
  "Hypergeometric p (our method vs dnds) = ",
  phyper(
    nrow(outliers_intersect_our_dnds), # k_0
    nrow(outliers_intersect_dNdS_adj), # K
    nrow(clonal_summary_binom) - nrow(outliers_intersect_dNdS_adj), # N - K
    nrow(outliers_adj), # n
    lower.tail = FALSE
  ) %>% signif(3)
)

test_res_dnds_ours_f <- paste0(
  "Hypergeometric p (our method vs dnds Fisher) = ",
  phyper(
    nrow(outliers_intersect_our_dnds_f), # k_0
    nrow(outliers_dNdS_adj_Fisher), # K
    nrow(clonal_summary_binom) - nrow(outliers_dNdS_adj_Fisher), # N - K
    nrow(outliers_adj), # n
    lower.tail = FALSE
  ) %>% signif(3)
)
# test_res_dnds_adj <- paste0(
#   "Hypergeometric p (adjusted) = ",
#   phyper(
#     nrow(outliers_intersect_dNdS),
#     nrow(cgc_summary_binom),
#     nrow(clonal_summary_binom) - nrow(cgc_summary_binom),
#     nrow(outliers_dNdS),
#     lower.tail = FALSE
#   ) %>% signif(3)
# )