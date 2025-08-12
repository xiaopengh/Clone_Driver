library(data.table)
library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)
library(progress)
source("helpers.R")

# ==============================================================================

# Data cleaning and preparation for calculating synonymous and non-synonymous mutations 

# ==============================================================================

# Define the file paths
gtf_file_path <- "../data/gencode.v19.annotation.gtf"
fasta_file_path <- "../data/GRCh37.p13.genome.fa"

# Read GTF file
message("Reading GTF file...")
gtf_data <- tryCatch({
  rtracklayer::import(gtf_file_path)
}, error = function(e) {
  stop("Error reading GTF file: ", e$message)
})
message("GTF file loaded.")

# transform gtf_data to data.table format
gtf_dt <- as.data.table(as.data.frame(gtf_data))
rm(gtf_data)

# Read transcript FASTA file
message("Reading FASTA file...")
chromosomes <- tryCatch({
  Biostrings::readDNAStringSet(fasta_file_path)
}, error = function(e) {
  stop("Error reading FASTA file: ", e$message)
})
message("FASTA file loaded.")

# Clean transcript IDs FASTA headers - GTF transcript_id format
# GENCODE FASTA headers are often like ">ENST0000012345.6|ENSG00000XYZ.1|..."
original_fasta_names <- names(chromosomes)
cleaned_fasta_names <- sapply(strsplit(original_fasta_names, "\\|"), `[`, 1)
cleaned_fasta_names <- sub(" .*", "", cleaned_fasta_names) # Remove anything after a space if ID is first
names(chromosomes) <- cleaned_fasta_names
message("FASTA sequence names cleaned.")

rm(original_fasta_names, cleaned_fasta_names)

# Read mc3 dataset which contains the observed mutations.
mc3_raw <- fread("../data/mc3.maf", sep = "\t")
mc3 <- mc3_raw[, .(Hugo_Symbol, t_alt_count, t_ref_count, Segment_Mean, Purity, Variant_Classification, Transcript_ID)]
rm(mc3_raw)

# Perform stat tests on mc3 dataset
mc3_summary_bytx <- build_mc3_summary_bytx(mc3)

# Free up memory
gc()

# ==============================================================================

# Main loop over all 25 chromosomes (chr1 to chr22 + X Y M)

# ==============================================================================

all_names <- names(chromosomes)[1:25]
total_chromosomes <- length(all_names)
summaries <- vector("list", total_chromosomes)
names(summaries) <- all_names

for (i in 1:total_chromosomes){
  
  chr_name <- all_names[i]
  
  cat(sprintf("Processing %d/%d, chromosome name : %s\n", i, total_chromosomes, chr_name))
  
  chromosome <- chromosomes[i]
  
  # For each coding site, determine all possible single-nucleotide mutations.
  siteRanges <- extract_coding_sites(gtf_dt)
    
  # Calculate the number of synonymous and non-synonymous mutation opportunities for each transcript.
  pmutdt <- determine_SNP_effects(siteRanges)
    
  # Build summary data.table of theoretical mutation types opportunities and observed mutation frequency by txid
  summary_txid <- compute_dN_dS_metrics(pmutdt)
  
  # Write summary_txid into summaries list
  summaries[[chr_name]] <- summary_txid
  
}

rm(chromosome, chromosomes, gtf_dt, mc3, mc3_summary_bytx, pmutdt, siteRanges, summary_txid,
   all_names, chr_name, i, total_chromosomes)