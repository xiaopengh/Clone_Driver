library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(txdbmaker)

# Define the file paths
gtf_file_path <- "../data/gencode.v48.annotation.gtf"
fasra_file_path <- "../data/gencode.v48.transcripts.fa"

# Read GTF file
message("Reading GTF file...")
gtf_data <- tryCatch({
  rtracklayer::import(gtf_file_path)
}, error = function(e) {
  stop("Error reading GTF file: ", e$message)
})
message("GTF file loaded.")


































































































































