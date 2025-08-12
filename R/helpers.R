# ==============================================================================

# Helper functions to be executed in main.R

# ==============================================================================

extract_coding_sites <- function(gtf_dt, show_progress = TRUE) {
  
  if (show_progress) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Package 'progress' is required for progress bar. Install with: install.packages('progress')")
    }
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent :current/:total (:elapsed) :step",
      total = 10,
      clear = FALSE,
      width = 80
    )
  }
  
  message("Extracting coding sites from GTF data...")
  
  # Step 1: Filter for CDS and stop_codon
  if (show_progress) pb$tick(tokens = list(step = "Filtering CDS entries"))
  cdsgtf <- gtf_dt[gtf_dt$type %in% c("CDS", "stop_codon")]
  
  # Step 2: Filter by chromosome
  if (show_progress) pb$tick(tokens = list(step = "Filtering by chromosome"))
  if (exists("chr_name")) {
    cdsgtf <- cdsgtf[cdsgtf$seqnames == chr_name]
  } else {
    stop("chr_name doesn't exist in scope")
  }
  
  # Step 3: Order by start position
  if (show_progress) pb$tick(tokens = list(step = "Ordering by position"))
  cdsgtf <- cdsgtf[order(start)]
  
  # Step 4: Create genomic ranges
  if (show_progress) pb$tick(tokens = list(step = "Creating genomic ranges"))
  gtfRanges <- GenomicRanges::GRanges(cdsgtf$seqnames, IRanges::IRanges(cdsgtf$start, cdsgtf$end))
  siteRanges <- unlist(GenomicRanges::slidingWindows(gtfRanges, 1))
  
  # Step 5: Add transcript and strand information
  if (show_progress) pb$tick(tokens = list(step = "Adding transcript info"))
  siteRanges$transcript_id <- rep(cdsgtf$transcript_id, GenomicRanges::width(gtfRanges))
  siteRanges$codingStrand <- rep(cdsgtf$strand, GenomicRanges::width(gtfRanges))
  
  # Step 6: Split by transcript and filter for multiples of 3
  if (show_progress) pb$tick(tokens = list(step = "Splitting by transcript"))
  siteRanges <- split(siteRanges, siteRanges$transcript_id)
  multiple3 <- sapply(siteRanges, function(gr) length(gr) %% 3) == 0
  siteRanges <- siteRanges[multiple3]
  
  # Step 7: Reverse negative strand transcripts
  if (show_progress) pb$tick(tokens = list(step = "Processing strand orientation"))
  reverse <- sapply(siteRanges, function(gr) gr$codingStrand[1]) == "-"
  siteRanges[reverse] <- S4Vectors::endoapply(siteRanges[reverse], rev)
  
  # Step 8: Calculate relative positions
  if (show_progress) pb$tick(tokens = list(step = "Calculating positions"))
  relpos <- lapply(IRanges::elementNROWS(siteRanges), function(l) 1:l)
  siteRanges <- unlist(siteRanges)
  siteRanges$relpos <- unlist(relpos)
  rm(relpos)
  
  # Step 9: Extract reference nucleotides
  if (show_progress) pb$tick(tokens = list(step = "Extracting nucleotides"))
  if (exists("chromosome")) {
    siteRanges$ref <- chromosome[siteRanges]
  } else {
    stop("chromosome doesn't exist in scope")
  }
  
  # Step 10: Get coding strand nucleotides
  if (show_progress) pb$tick(tokens = list(step = "Processing coding strand"))
  siteRanges$codingRef <- siteRanges$ref
  reverse <- siteRanges$codingStrand == "-"
  siteRanges$codingRef[reverse] <- Biostrings::reverseComplement(siteRanges$codingRef[reverse])
  
  message("Extracted coding sites from GTF data.")
  
  return(siteRanges)
}

# extract_coding_sites <- function(gtf_dt){
#   
#   message("Extracting coding sites from GTF data...")
# 
#   cdsgtf <- gtf_dt[gtf_dt$type %in% c("CDS", "stop_codon")]
#   
#   if (exists("chr_name")) {
#     cdsgtf <- cdsgtf[cdsgtf$seqnames == chr_name]
#   } else {
#       stop("chr_name doesn't exist in scope")
#     }
#   
#   cdsgtf <- cdsgtf[order(start)]
#   
#   gtfRanges <- GenomicRanges::GRanges(cdsgtf$seqnames, IRanges::IRanges(cdsgtf$start, cdsgtf$end))
#   
#   siteRanges <- unlist(GenomicRanges::slidingWindows(gtfRanges, 1))
#   
#   siteRanges$transcript_id <- rep(cdsgtf$transcript_id, GenomicRanges::width(gtfRanges))
#   
#   siteRanges$codingStrand <- rep(cdsgtf$strand, GenomicRanges::width(gtfRanges))
#   
#   # split sites by transcript, if they are in negative strand reverse order so that everything is in transcription order
#   siteRanges <- split(siteRanges, siteRanges$transcript_id)
#   
#   # get only the transcripts whose length is a multiple of 3 
#   multiple3 <- sapply(siteRanges, function(gr) length(gr) %% 3) == 0
#   siteRanges <- siteRanges[multiple3]
#   
#   # if they are in negative strand reverse order so that everything is in transcription order
#   reverse <- sapply(siteRanges, function(gr) gr$codingStrand[1]) == "-"
#   siteRanges[reverse] <- S4Vectors::endoapply(siteRanges[reverse], rev)
#   
#   # get relative position of each site within its transcript CDS
#   relpos <- lapply(IRanges::elementNROWS(siteRanges), function(l) 1:l)
#   siteRanges <- unlist(siteRanges)
#   siteRanges$relpos <- unlist(relpos)
#   rm(relpos)
#   
#   # get nucleotides on + strand
#   if (exists("chromosome")) {
#     siteRanges$ref <- chromosome[siteRanges]
#   } else {
#       stop("chromosome doesn't exist in scope")
#     }
#   
#   # get nucleotides on - strand
#   siteRanges$codingRef <- siteRanges$ref
#   reverse <- siteRanges$codingStrand == "-"
#   siteRanges$codingRef[reverse] <- Biostrings::reverseComplement(siteRanges$codingRef[reverse])
#   
#   message("Extracted coding sites from GTF data.")
#   
#   return(siteRanges)
# }

# ==============================================================================

determine_SNP_effects <- function(siteRanges, show_progress = TRUE) {
  
  if (show_progress) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Package 'progress' is required for progress bar. Install with: install.packages('progress')")
    }
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent :current/:total (:elapsed) :step",
      total = 7,
      clear = FALSE,
      width = 80
    )
  }
  
  message("Determining SNP effects on coding sites...")
  
  # Step 1: Calculate frame
  if (show_progress) pb$tick(tokens = list(step = "Calculating reading frames"))
  siteRanges$frame <- ((siteRanges$relpos + 2L) %% 3L) + 1L
  
  # Step 2: Determine codons
  if (show_progress) pb$tick(tokens = list(step = "Determining codons"))
  siteRanges$codon <- rep(
    xscat(
      siteRanges$codingRef[siteRanges$frame == 1],
      siteRanges$codingRef[siteRanges$frame == 2],
      siteRanges$codingRef[siteRanges$frame == 3]
    ),
    each = 3
  )
  
  # Step 3: Translate to amino acids
  if (show_progress) pb$tick(tokens = list(step = "Translating to amino acids"))
  siteRanges$wtAA <- Biostrings::translate(siteRanges$codon, no.init.codon = TRUE)
  
  # Step 4: Generate possible mutations
  if (show_progress) pb$tick(tokens = list(step = "Generating possible mutations"))
  pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"), A = c("C", "G", "T"), G = c("A", "C", "T"))[as.character(siteRanges$codingRef)]
  pmutRanges <- siteRanges[rep(1:length(siteRanges), each = 3)]
  pmutRanges$mut <- Biostrings::DNAStringSet(unlist(pmuts))
  
  # Step 5: Calculate mutated codons and amino acids
  if (show_progress) pb$tick(tokens = list(step = "Calculating mutated sequences"))
  pmutRanges$mutCodon <- pmutRanges$codon
  Biostrings::subseq(pmutRanges$mutCodon, pmutRanges$frame, pmutRanges$frame) <- pmutRanges$mut
  pmutRanges$mutAA <- Biostrings::translate(pmutRanges$mutCodon, no.init.codon = TRUE)
  
  # Step 6: Classify mutation types
  if (show_progress) pb$tick(tokens = list(step = "Classifying mutation types"))
  pmutRanges$type <- rep("missense", length(pmutRanges))
  pmutRanges$type[pmutRanges$wtAA == pmutRanges$mutAA] <- "syn"
  pmutRanges$type[pmutRanges$wtAA != "*" & pmutRanges$mutAA == "*"] <- "nonsense"
  pmutRanges$type[pmutRanges$wtAA == "*" & pmutRanges$mutAA != "*"] <- "nonstop"
  
  # Step 7: Convert to data.table
  if (show_progress) pb$tick(tokens = list(step = "Converting to data.table"))
  names(pmutRanges) <- 1:length(pmutRanges)
  pmutdt <- data.table::data.table(as.data.frame(pmutRanges))
  pmutdt[, c("width", "strand", "end") := NULL]
  names(pmutdt)[2] <- "position"
  
  message("Determined SNP effects on coding sites.")
  
  return(pmutdt)
}

# determine_SNP_effects <- function(siteRanges){
#   
#   message("Determining SNP effects on coding sites...")
#   
#   # get the codon and aminoacid that each site influences
#   siteRanges$frame <- ((siteRanges$relpos + 2L) %% 3L) + 1L # gives 1 if in frame, 2 for one base out of frame etc
#   siteRanges$codon <- rep(
#     xscat(
#       siteRanges$codingRef[siteRanges$frame == 1],
#       siteRanges$codingRef[siteRanges$frame == 2],
#       siteRanges$codingRef[siteRanges$frame == 3]
#     ),
#     each = 3
#   )
#   siteRanges$wtAA <- Biostrings::translate(siteRanges$codon, no.init.codon = TRUE)
#   
#   # get all possible mutations in the examined sites
#   pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"), A = c("C", "G", "T"), G = c("A", "C", "T"))[as.character(siteRanges$codingRef)]
#   pmutRanges <- siteRanges[rep(1:length(siteRanges), each = 3)]
#   pmutRanges$mut <- Biostrings::DNAStringSet(unlist(pmuts))
#   
#   # get the mutated codon and mutated aminoacid according to mutation in coding strand
#   pmutRanges$mutCodon <- pmutRanges$codon
#   Biostrings::subseq(pmutRanges$mutCodon, pmutRanges$frame, pmutRanges$frame) <- pmutRanges$mut
#   pmutRanges$mutAA <- Biostrings::translate(pmutRanges$mutCodon, no.init.codon = TRUE)
#   
#   # determine the type of mutation
#   pmutRanges$type <- rep("missense", length(siteRanges))
#   pmutRanges$type[pmutRanges$wtAA == pmutRanges$mutAA] <- "syn"
#   pmutRanges$type[pmutRanges$wtAA != "*" & pmutRanges$mutAA == "*"] <- "nonsense"
#   pmutRanges$type[pmutRanges$wtAA == "*" & pmutRanges$mutAA != "*"] <- "nonstop"
#   
#   # convert to data.table
#   names(pmutRanges) <- 1:length(pmutRanges)
#   pmutdt <- data.table::data.table(as.data.frame(pmutRanges))
#   pmutdt[, c("width", "strand", "end") := NULL]
#   names(pmutdt)[2] <- "position"
#   
#   message("Determined SNP effects on coding sites.")
#   
#   return(pmutdt)
# 
# }

# ==============================================================================

build_mc3_summary_bytx <- function(mc3) {
  
  message("Building mc3 summary by transcript...")
  
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
  
  # aggregate mc3 by txid
  mc3_summary_bytx <- mc3[, .(
    obs_dN_clonal = sum(!is_subclonal & !is_synonymous),
    obs_dS_clonal = sum(!is_subclonal & is_synonymous),
    obs_dN_subclonal = sum(is_subclonal & !is_synonymous),
    obs_dS_subclonal = sum(is_subclonal & is_synonymous)
  ), by = .(Transcript_ID, Hugo_Symbol)]
  
  message("Built mc3 summary by transcript.")
  
  return(mc3_summary_bytx)
}

# ==============================================================================

compute_dN_dS_metrics <- function(pmutdt) {
  
  message("Computing dN/dS metrics from pmutdt...")
  
  # Init of summary txid
  summary_txid <- pmutdt[, .(
    non_syn_count = sum(type != "syn"),
    syn_count = sum(type == "syn")
  ), by = transcript_id]
  
  # split transcript_id in summary obtained from pmutdt from its raw id and version numbers
  split_matrix <- stringr::str_split_fixed(summary_txid$transcript_id, pattern = "\\.", n = 2)
  summary_txid[, maf_txid := split_matrix[, 1]]
  summary_txid[, txid_version := split_matrix[, 2]]
  
  
  if (exists("mc3_summary_bytx")) {
    # summary_txid has 'maf_txid' but mc3_summary_bytx has 'Transcript_ID'
    summary_txid[mc3_summary_bytx, 
                 c("obs_dN_clonal", "obs_dS_clonal", "obs_dN_subclonal", "obs_dS_subclonal") := 
                   .(i.obs_dN_clonal, i.obs_dS_clonal, i.obs_dN_subclonal, i.obs_dS_subclonal),
                 on = c("maf_txid" = "Transcript_ID")]
  } else {
    stop("mc3_summary_bytx does not exist in scope. Please run build_mc3_summary_bytx() first.")
  }
  
  
  # set NAs to 0 because it means 0 obs 
  summary_txid[is.na(summary_txid)] <- 0
  
  # calculate dN/dS ratios for clonal and subclonal mutations
  summary_txid[, dN_to_dS_clonal := ifelse( 
    (non_syn_count == 0 | syn_count == 0 | obs_dS_clonal == 0),
    NA_real_, (obs_dN_clonal / non_syn_count) / (obs_dS_clonal / syn_count) )]
  
  summary_txid[, dN_to_dS_subclonal := ifelse(
    (non_syn_count == 0 | syn_count == 0 | obs_dS_subclonal == 0),
    NA_real_, (obs_dN_subclonal / non_syn_count) / (obs_dS_subclonal / syn_count) )]
  
  message("Computed dN/dS metrics from pmutdt.")
  
  return(summary_txid)
}











