length(transcript_seqs)
transcript_seqs[[1]][1:3]
transcript_seqs[[1]][1:100]
for (i in 1:20) {
  transcript_seqs[[i]][1:100]
}

rep(
  xscat(
    siteRanges$codingRef[siteRanges$frame == 1],
    siteRanges$codingRef[siteRanges$frame == 2],
    siteRanges$codingRef[siteRanges$frame == 3]
  ),
  each = 3
)

for (id_val in siteRanges$transcript_id %>% unique()) {
  
  subset_granges <- siteRanges[siteRanges$transcript_id == id_val]
  l <- length(subset_granges)
  
  # cat("ID:", id_val, " Length:", l, "Multiple of 3?", l%%3 == 0, "\n")
  if (l%%3 != 0) cat("ID:", id_val, " Length:", l, "Modulus 3", l%%3, "\n")
}

seq3 <- siteRanges$codingRef[siteRanges$frame == 3][1:12]
seq2 <- siteRanges$codingRef[siteRanges$frame == 2][1:10]
seq1 <- siteRanges$codingRef[siteRanges$frame == 1][1:10]
xscat( seq1, seq2, seq3)
?xscat

rm(seq1, seq2, seq3)

sr <- siteRanges[1:100]

rm(sr)

multiple3 <- sapply(siteRanges, function(gr) length(gr) %% 3) == 0

(sapply(siteRanges[multiple3], function(gr) length(gr) %% 3) == 0) %>% sum()

sum(multiple3)
sum(!multiple3)

sr <- siteRanges[1:100]

pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"), A = c("C", "G", "T"), G = c("A", "C", "T"))[as.character(sr$codingRef)]

sr$codingRef

str(pmuts)

pmutRanges$mut
pmutRanges

?Biostrings::subseq()

rm(gr, sr)

pmutRanges

library(ggplot2)

df <- data.frame(ratio = ratio_bytxid)
sd_factor <- 0.94
ggplot(df, aes(x = ratio)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "blue", alpha = 0.3, col = "gray31") + 
  stat_function(fun = dnorm, 
                args = list(mean = mean(df$ratio, na.rm = TRUE), 
                            sd = sd_factor * sd(df$ratio, na.rm = TRUE)),
                color = "navy",
                linewidth = 1
                )

rm(df, sd_factor)


gtf_dt$type %>% unique()

ls()

?fread


g <- ggplot(data = mc3_raw) 



# ==============================================================================

mc3_raw <- fread("../data/mc3.maf", sep = "\t")

# mc3_raw[, .(count_of_rows = .N), by = Transcript_ID]
# 
# library(ggplot2)
# ggplot(data = mc3_raw[, .(count_of_rows = .N), by = Transcript_ID]) +
#   geom_histogram(mapping = aes(x = count_of_rows), bins = 50)

# add is_sysnonymous column based on Variant_Classification 
mc3_raw[, is_synonymous := (Variant_Classification == "Silent")]

# Split transcript_id in summary obtained from pmutdt from its raw id and version numbers
split_matrix <- stringr::str_split_fixed(summary_txid$transcript_id, pattern = "\\.", n = 2)
summary_txid[, maf_txid := split_matrix[, 1]]
summary_txid[, txid_version := split_matrix[, 2]]
rm(split_matrix)








