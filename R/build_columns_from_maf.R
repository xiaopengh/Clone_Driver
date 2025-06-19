# dataset import
mc3_raw <- fread("../data/mc3.maf", sep = "\t")
mc3 <- mc3_raw[, .(Hugo_Symbol, t_alt_count, t_ref_count, Segment_Mean, Purity, Variant_Classification, Transcript_ID)]
rm(mc3_raw)
gc()

# add coverage, pval and restraint data to CN == 2
mc3 <- mc3[(Segment_Mean >= -0.3) & (Segment_Mean <= 0.3)]
mc3[, coverage := (t_alt_count + t_ref_count)]
mc3 <- mc3[complete.cases(mc3)]
mc3[, binom_pval := mapply(function(x, n, p_hyp) {
  binom.test(x, n, p = p_hyp, alternative = "less")$p.value
  }, x = t_alt_count, n = coverage, p_hyp = Purity * 0.5)]

# add is_sysnonymous, is_subclonal column based on Variant_Classification
mc3[, is_subclonal := (binom_pval < 0.05)] 
mc3[, is_synonymous := (Variant_Classification == "Silent")]

# split transcript_id in summary obtained from pmutdt from its raw id and version numbers
split_matrix <- stringr::str_split_fixed(summary_txid$transcript_id, pattern = "\\.", n = 2)
summary_txid[, maf_txid := split_matrix[, 1]]
summary_txid[, txid_version := split_matrix[, 2]]
rm(split_matrix)

# aggregate by txid and add columns to summary_txid
mc3_summary_bytx <- mc3[, .(
  obs_dN_clonal = sum(!is_subclonal & !is_synonymous),
  obs_dS_clonal = sum(!is_subclonal & is_synonymous),
  obs_dN_subclonal = sum(is_subclonal & !is_synonymous),
  obs_dS_subclonal = sum(is_subclonal & is_synonymous)
), by = Transcript_ID]

# summary_txid has 'maf_txid' but mc3_summary_bytx has 'Transcript_ID'
summary_txid[mc3_summary_bytx, 
             c("obs_dN_clonal", "obs_dS_clonal", "obs_dN_subclonal", "obs_dS_subclonal") := 
               .(i.obs_dN_clonal, i.obs_dS_clonal, i.obs_dN_subclonal, i.obs_dS_subclonal),
             on = c("maf_txid" = "Transcript_ID")]

# set NAs to 0 because it means 0 obs 
summary_txid[is.na(summary_txid)] <- 0

save.image()