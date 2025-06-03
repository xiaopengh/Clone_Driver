library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(txdbmaker)

# ---- Load GTF and extract CDS regions ----
txdb <- makeTxDbFromGFF("../data/gencode.v48.annotation.gtf", format = "gtf")

cds_by_tx <- cdsBy(txdb, by = "tx", use.names = TRUE)

# ---- Load genome FASTA ----
fasta_file <- "gencode.vXX.genome.fa"
genome <- readDNAStringSet(fasta_file)
names(genome) <- sub(" .*", "", names(genome))

# ---- STEP 3: Extract and translate coding sequences ----
get_cds_seq <- function(cds_gr, genome) {
  seq <- getSeq(genome, cds_gr)
  if (as.character(strand(cds_gr)[1]) == "-") seq <- reverseComplement(seq)
  return(unlist(seq))
}

# Function to count synonymous vs nonsynonymous variants per codon
count_syn_nonsyn <- function(dna_seq) {
  if (length(dna_seq) < 3 || (length(dna_seq) %% 3) != 0) return(NULL)
  codons <- as.character(subseq(dna_seq, start = seq(1, length(dna_seq) - 2, by = 3), width = 3))
  codons <- codons[codons %in% names(GENETIC_CODE)]  # Only valid codons
  
  aa_seq <- GENETIC_CODE[codons]
  syn_count <- 0
  nonsyn_count <- 0
  
  for (i in seq_along(codons)) {
    original_codon <- codons[i]
    original_aa <- aa_seq[i]
    
    for (pos in 1:3) {
      for (base in c("A", "C", "G", "T")) {
        mutated_codon <- original_codon
        substr(mutated_codon, pos, pos) <- base
        if (mutated_codon == original_codon || !(mutated_codon %in% names(GENETIC_CODE)))
          next
        new_aa <- GENETIC_CODE[mutated_codon]
        if (new_aa == original_aa) {
          syn_count <- syn_count + 1
        } else {
          nonsyn_count <- nonsyn_count + 1
        }
      }
    }
  }
  
  return(c(synonymous = syn_count, nonsynonymous = nonsyn_count))
}

# ---- STEP 4: Loop over transcripts and summarize ----
results <- lapply(names(cds_by_tx), function(tx_id) {
  cds_gr <- cds_by_tx[[tx_id]]
  if (length(cds_gr) == 0) return(NULL)
  
  seq <- get_cds_seq(cds_gr, genome)
  if (length(seq) == 0) return(NULL)
  
  count_syn_nonsyn(seq)
})

# ---- STEP 5: Aggregate and report ----
results <- do.call(rbind, results)
total_synonymous <- sum(results[, "synonymous"], na.rm = TRUE)
total_nonsynonymous <- sum(results[, "nonsynonymous"], na.rm = TRUE)

cat("Total synonymous mutations:", total_synonymous, "\n")
cat("Total non-synonymous mutations:", total_nonsynonymous, "\n")
cat("Ratio (nonsyn/syn):", total_nonsynonymous / total_synonymous, "\n")
