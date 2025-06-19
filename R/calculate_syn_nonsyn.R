library(data.table)
library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)


# ==============================================================================

# Use gencode release 19

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

gtf_dt <- as.data.table(as.data.frame(gtf_data))

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
# We only want the "ENST0000012345.6" part.
original_fasta_names <- names(chromosomes)
cleaned_fasta_names <- sapply(strsplit(original_fasta_names, "\\|"), `[`, 1)
cleaned_fasta_names <- sub(" .*", "", cleaned_fasta_names) # Remove anything after a space if ID is first
names(chromosomes) <- cleaned_fasta_names
message("FASTA sequence names cleaned.")

rm(original_fasta_names, cleaned_fasta_names)

chromosome <- chromosomes[1]

# ==============================================================================

cdsgtf <- gtf_dt[gtf_dt$type %in% c("CDS", "stop_codon")] # verify

cdsgtf <- cdsgtf[cdsgtf$seqnames == "chr1"]

cdsgtf <- cdsgtf[order(start)]

gtfRanges <- GenomicRanges::GRanges(cdsgtf$seqnames, IRanges::IRanges(cdsgtf$start, cdsgtf$end))

siteRanges <- unlist(GenomicRanges::slidingWindows(gtfRanges, 1))

siteRanges$transcript_id <- rep(cdsgtf$transcript_id, GenomicRanges::width(gtfRanges))

siteRanges$codingStrand <- rep(cdsgtf$strand, GenomicRanges::width(gtfRanges))

rm(gtfRanges, cdsgtf, gtf_data)

# split sites by transcript, if they are in negative strand reverse order so that everything is in transcription order
siteRanges <- split(siteRanges, siteRanges$transcript_id)

# get only the transcripts whose length is a multiple of 3 
multiple3 <- sapply(siteRanges, function(gr) length(gr) %% 3) == 0
siteRanges <- siteRanges[multiple3]

# if they are in negative strand reverse order so that everything is in transcription order
reverse <- sapply(siteRanges, function(gr) gr$codingStrand[1]) == "-"
siteRanges[reverse] <- S4Vectors::endoapply(siteRanges[reverse], rev)

# get relative position of each site within its transcript CDS
relpos <- lapply(IRanges::elementNROWS(siteRanges), function(l) 1:l)
siteRanges <- unlist(siteRanges)
siteRanges$relpos <- unlist(relpos)
rm(relpos)

# get nucleotides on + strand
siteRanges$ref <- chromosome[siteRanges]

# get nucleotides on - strand
siteRanges$codingRef <- siteRanges$ref
reverse <- siteRanges$codingStrand == "-"
siteRanges$codingRef[reverse] <- Biostrings::reverseComplement(siteRanges$codingRef[reverse])

# ?determine pyrimidine strand and get pyrimidine
# siteRanges$pyrimidineStrand <- ifelse(as.character(siteRanges$ref) %in% c("C", "T"), "+", "-")
# siteRanges$pyrimidine <- siteRanges$ref
# reverse <- siteRanges$pyrimidineStrand == "-"
# siteRanges$pyrimidine[reverse] <- Biostrings::reverseComplement(siteRanges$pyrimidine[reverse])

# get the codon and aminoacid that each site influences
siteRanges$frame <- ((siteRanges$relpos + 2L) %% 3L) + 1L # gives 1 if in frame, 2 for one base out of frame and 3 for two bases out of frame
siteRanges$codon <- rep(
  xscat(
    siteRanges$codingRef[siteRanges$frame == 1],
    siteRanges$codingRef[siteRanges$frame == 2],
    siteRanges$codingRef[siteRanges$frame == 3]
  ),
  each = 3
)
siteRanges$wtAA <- Biostrings::translate(siteRanges$codon, no.init.codon = TRUE)

# get all possible mutations in the examined sites
pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"), A = c("C", "G", "T"), G = c("A", "C", "T"))[as.character(siteRanges$codingRef)]
pmutRanges <- siteRanges[rep(1:length(siteRanges), each = 3)]
# pmutRanges$pyrimidineMut <- Biostrings::DNAStringSet(unlist(pmuts))
# reverse <- rep(reverse, each = 3)
# pmutRanges$mut <- pmutRanges$pyrimidineMut
# pmutRanges$mut[reverse] <- Biostrings::reverseComplement(pmutRanges$mut[reverse])
pmutRanges$mut <- Biostrings::DNAStringSet(unlist(pmuts))
rm(pmuts)

# # get the mutation in the coding strand according to pyrimidine strand
# reverse <- pmutRanges$pyrimidineStrand != pmutRanges$codingStrand
# pmutRanges$codingMut <- pmutRanges$pyrimidineMut
# pmutRanges$codingMut[reverse] <- Biostrings::reverseComplement(pmutRanges$codingMut[reverse])
# 
# # get the mutated codon and mutated aminoacid according to mutation in coding strand
# pmutRanges$mutCodon <- pmutRanges$codon
# Biostrings::subseq(pmutRanges$mutCodon, pmutRanges$frame, pmutRanges$frame) <- pmutRanges$codingMut
# pmutRanges$mutAA <- Biostrings::translate(pmutRanges$mutCodon, no.init.codon = TRUE)

# get the mutated codon and mutated aminoacid according to mutation in coding strand
pmutRanges$mutCodon <- pmutRanges$codon
Biostrings::subseq(pmutRanges$mutCodon, pmutRanges$frame, pmutRanges$frame) <- pmutRanges$mut
pmutRanges$mutAA <- Biostrings::translate(pmutRanges$mutCodon, no.init.codon = TRUE)

# determine the type of mutation
pmutRanges$type <- rep("missense", length(siteRanges))
pmutRanges$type[pmutRanges$wtAA == pmutRanges$mutAA] <- "syn"
pmutRanges$type[pmutRanges$wtAA != "*" & pmutRanges$mutAA == "*"] <- "nonsense"
pmutRanges$type[pmutRanges$wtAA == "*" & pmutRanges$mutAA != "*"] <- "nonstop"
rm(siteRanges)

# convert to data.table
names(pmutRanges) <- 1:length(pmutRanges)
pmutdt <- data.table::data.table(as.data.frame(pmutRanges))
pmutdt[, c("width", "strand", "end") := NULL]
names(pmutdt)[2] <- "position"
rm(pmutRanges)

# regroup data table by transcript_id then calculate dn/ds by transcript_id
# pmutdt_bytxid <- split(pmutdt, by = "transcript_id")
# dnds_bytxid <- split(pmutdt, by = "transcript_id") %>%
#   sapply(
#     function(dt) {
#       # Number of non-synonymous mutations (or not "syn")
#       non_syn_count <- sum(dt$type != "syn")
#       
#       # Number of synonymous mutations ("syn")
#       syn_count <- sum(dt$type == "syn")
#       
#       return()
#             
#       # Handle cases where syn_count might be zero to avoid division by zero
#       if (syn_count > 0) {
#         return(non_syn_count / syn_count)
#       } else {
#         # Decide how to handle this: NA, 0, or Inf depending on your interpretation
#         return(NA) # Or Inf if you consider it an infinite ratio
#       }
#     }
#   )

summary_txid <- pmutdt[, .(
  non_syn_count = sum(type != "syn"),
  syn_count = sum(type == "syn")
), by = transcript_id]

# Remove temp data for one chromosome loop
rm(chromosome, gtf_dt, multiple3, reverse)
