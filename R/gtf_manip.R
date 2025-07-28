library(GenomicRanges)
library(BiocParallel) # for parallel processing
param <- MulticoreParam(workers = 16)  # Adjust the number of cores

# Define gtf file path
gtf_file_path <- "../data/gencode.v48.chr_patch_hapl_scaff.annotation.gtf"

# Read GTF file
gtf_data <- tryCatch({
  rtracklayer::import(gtf_file_path)
}, error = function(e) {
  stop("Error reading GTF file: ", e$message)
})

exons <- gtf_data[gtf_data$type == "exon"]
utrs <- gtf_data[gtf_data$type == "UTR"]

exon_list <- split(exons, exons$transcript_id)
utr_list  <- split(utrs, utrs$transcript_id)

nonutr_exon_list <- bplapply(names(exon_list), function(tx_id) {
  ex <- exon_list[[tx_id]]
  ut <- utr_list[[tx_id]]
  if (is.null(ut)) return(ex)
  setdiff(ex, ut)
}, BPPARAM = param)

nonutr_exons <- do.call(c, nonutr_exon_list)
nonutr_exons$type <- "nonUTR_exon"


