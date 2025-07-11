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

clonal_summary_binom[, `:=`(
  p_val_poissonTest_clonal = poisson.test(c(obs_dN_clonal, obs_dS_clonal),
                                          c(non_syn_count, syn_count))$p.value,
  p_val_poissonTest_subclonal = poisson.test(c(obs_dN_subclonal, obs_dS_subclonal),
                                             c(non_syn_count, syn_count))$p.value
), by = 1:nrow(clonal_summary_binom)]