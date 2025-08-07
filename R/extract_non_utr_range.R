library(data.table)
library(GenomicRanges)
library(BiocParallel) # for parallel processing
source("helpers_utr.R")

# # Detect number of logical cores
# n_cores <- parallel::detectCores()
# bp_param <- SnowParam(workers = n_cores, type = "SOCK", progressbar = TRUE) # Adjust the number of cores
# register(bp_param)

# Define gtf file path
gtf_file_path <- "../data/gencode.v48.chr_patch_hapl_scaff.annotation.gtf"

# Read GTF file
gtf_data <- tryCatch({
  rtracklayer::import(gtf_file_path)
}, error = function(e) {
  stop("Error reading GTF file: ", e$message)
})

{ # Annotate 5'UTR and 3'UTR (using a function that only accepts DT)
  gtf_dt <- as.data.table(gtf_data)
  mrnaGtf <- gtf_dt[(gtf_dt$type %in% c("CDS", "UTR"))]
  utr5Gtf <- get5utr(mrnaGtf)
  utr3Gtf <- get3utr(mrnaGtf)
  utr5Gtf$type = "5'UTR"
  utr3Gtf$type = "3'UTR"
  utr_dt <- rbindlist(list(utr3Gtf, utr5Gtf))
  rm(utr5Gtf, utr3Gtf) ; gc()
  utr_gr <- GRanges(
    seqnames = utr_dt$seqnames,
    ranges = IRanges(start = utr_dt$start, end = utr_dt$end),
    strand = utr_dt$strand
  )
  mcols(utr_gr) <- utr_dt[, 6:27] # 6:27 is the meta data column
}

utrs <- utr_gr[utr_gr$gene_type == "protein_coding"]

genes <- gtf_data[gtf_data$type == "gene" & gtf_data$gene_type == "protein_coding"]

exons <- gtf_data[gtf_data$type == "exon" & gtf_data$gene_type == "protein_coding"]

substract_regions <- function(gr1, gr2){ 
  # Produce a final data table 
  ov <- findOverlaps(gr1, gr2, ignore.strand = FALSE)
  igr <- pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])
  gr3 <- setdiff(gr1, igr)
  # mov <- findOverlaps(non_utr, gr1, ignore.strand = TRUE) # to get metadata
  # mov <- mov[!duplicated(queryHits(mov))]
  # mcols(non_utr) <- mcols(gr1)[subjectHits(mov), ]
  return(gr3)
}

non_utr <- substract_regions(exons, utrs)
introns <- substract_regions(genes, exons)

output_dt <- rbind(
  data.table(
    seqname = as.character(seqnames(non_utr)),
    start = start(non_utr),
    end = end(non_utr),
    feature = "cds"
  ),
  data.table(
    seqname = as.character(seqnames(introns)),
    start = start(introns),
    end = end(introns),
    feature = "intron"
  )
)

fwrite(output_dt, "../data/output_dt.tsv")
