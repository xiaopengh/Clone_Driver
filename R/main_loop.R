library(rtracklayer)
library(Biostrings)
library(dplyr)
library(stringr)
library(GenomicFeatures) # For TxDb and related functions

# --- 1. Define File Paths and Parameters ---
gtf_file_path <- "../data/gencode.v48.annotation.gtf" # User provided
fasta_file_path <- "../data/gencode.v48.transcripts.fa" # User provided

# Standard DNA bases
DNA_BASES <- c("A", "C", "G", "T")

# --- 2. Load and Preprocess Data ---

# Read GTF file (raw import for chromosome list and basic CDS feature check)
message("Reading GTF file (for initial processing)...")
gtf_data_raw <- tryCatch({
  rtracklayer::import(gtf_file_path)
}, error = function(e) {
  stop("Error reading GTF file: ", e$message)
})
message("GTF file (raw) loaded.")

# Create TxDb from GTF (this is the primary source for CDS/UTR info)
message("Creating Transcript Database (TxDb) from GTF... This may take a few minutes.")
txdb <- tryCatch({
  txdbmaker::makeTxDbFromGFF(gtf_file_path, format="gtf")
}, error = function(e) {
  stop("Error creating TxDb object: ", e$message)
})
message("TxDb created.")

# Get 5' UTRs and CDS regions by transcript from TxDb
message("Extracting 5'UTR and CDS information from TxDb...")
five_utr_by_tx <- tryCatch({
  fiveUTRsByTranscript(txdb, use.names=TRUE) # GRangesList
}, error = function(e) {
  message("Warning: Could not extract 5'UTRs by transcript. Assuming no 5'UTRs. Error: ", e$message)
  GRangesList() # Return empty GRangesList to allow script to proceed, but CDS extraction might be affected
})

cds_grlist_by_tx <- tryCatch({
  cdsBy(txdb, by="tx", use.names=TRUE) # GRangesList of CDS genomic coordinates
}, error = function(e) {
  stop("Error extracting CDS by transcript from TxDb: ", e$message)
})
message("5'UTR and CDS information extracted.")

# Filter for features of type CDS from the raw GTF import
# This is mainly to get a list of chromosomes and transcripts that *have* CDS annotations
# for the main looping logic.
cds_features_gtf <- gtf_data_raw[gtf_data_raw$type == "CDS"]
if (!"transcript_id" %in% names(mcols(cds_features_gtf))) {
  # Attempt to find transcript_id in other potential attribute columns if GTF is non-standard
  if ("transcript" %in% names(mcols(cds_features_gtf))) {
    cds_features_gtf$transcript_id <- cds_features_gtf$transcript
  } else {
    stop("GTF metadata (from direct import) does not contain 'transcript_id' or 'transcript' column.")
  }
}
cds_features_gtf <- cds_features_gtf[!is.na(cds_features_gtf$transcript_id)]
message(paste("Found", length(cds_features_gtf), "CDS features in raw GTF import (used for chromosome/transcript iteration list)."))


# Read transcript FASTA file
message("Reading FASTA file (full transcripts expected)...")
transcript_seqs_full <- tryCatch({
  Biostrings::readDNAStringSet(fasta_file_path)
}, error = function(e) {
  stop("Error reading FASTA file: ", e$message)
})
message("FASTA file loaded.")

# Clean transcript IDs from FASTA headers
original_fasta_names <- names(transcript_seqs_full)
# Handles cases like ">ENST...|ENSG...|..." or ">ENST... gene=..."
cleaned_fasta_names <- sapply(strsplit(original_fasta_names, "\\|"), `[`, 1) # Get part before first pipe
cleaned_fasta_names <- sub(" .*", "", cleaned_fasta_names) # Remove anything after a space
# Remove potential ">" if it wasn't handled by strsplit/sub alone (e.g. ">ENST...")
cleaned_fasta_names <- sub(">", "", cleaned_fasta_names)
names(transcript_seqs_full) <- cleaned_fasta_names
message("FASTA sequence names cleaned.")

# --- 3. Get Chromosome List ---
# Using chromosomes present in the CDS features from the raw GTF import
chromosomes <- unique(as.character(seqnames(cds_features_gtf)))
if (length(chromosomes) == 0) {
  stop("No chromosomes found with CDS features in the GTF file. Cannot proceed.")
}
message(paste("Chromosomes to process (based on GTF CDS features):", paste(chromosomes, collapse=", ")))

# ==============================================================================

# --- 1. Get Chromosome List ---
# Using chromosomes present in the CDS features from the raw GTF import
chromosomes <- unique(as.character(seqnames(cds_features_gtf)))
if (length(chromosomes) == 0) {
  stop("No chromosomes found with CDS features in the GTF file. Cannot proceed.")
}
message(paste("Chromosomes to process (based on GTF CDS features):", paste(chromosomes, collapse=", ")))

# --- 2. Main Processing Loop ---
all_mutation_results <- list()

for (chr_name in chromosomes) {
  message(paste("\nProcessing chromosome:", chr_name))
  
  # Filter CDS features from raw GTF for the current chromosome to get transcript list
  cds_on_chr_gtf <- cds_features_gtf[seqnames(cds_features_gtf) == chr_name]
  if (length(cds_on_chr_gtf) == 0) {
    message(paste("No CDS features found in raw GTF on chromosome", chr_name, ". Skipping."))
    next
  }
  
  # Get unique transcript IDs on this chromosome that have CDS features in the GTF
  transcripts_on_chr <- unique(cds_on_chr_gtf$transcript_id)
  if (length(transcripts_on_chr) == 0) {
    message(paste("No transcripts with CDS features found on chromosome", chr_name, "after filtering raw GTF. Skipping."))
    next
  }
  
  message(paste("Found", length(transcripts_on_chr), "transcripts with CDS annotations on", chr_name, "to process."))
  
  
  chromosome_results <- list()
  
  current_idx <- 1
  
  for (current_transcript_id in transcripts_on_chr) {
    
    message(paste("Processing transcript:", current_transcript_id, current_idx, "/", length(transcripts_on_chr)))
    
    # Check if the full transcript sequence exists in the FASTA data
    if (!current_transcript_id %in% names(transcript_seqs_full)) {
      # warning(paste("Transcript", current_transcript_id, "not found in FASTA file. Skipping."))
      next
    }
    full_transcript_dna <- transcript_seqs_full[[current_transcript_id]]
    if (length(full_transcript_dna) == 0) {
      # warning(paste("Transcript", current_transcript_id, "is empty in FASTA. Skipping."))
      next
    }
    
    # Get CDS genomic ranges for this transcript from TxDb
    cds_genomic_ranges_for_tx <- cds_grlist_by_tx[[current_transcript_id]]
    if (is.null(cds_genomic_ranges_for_tx) || length(cds_genomic_ranges_for_tx) == 0) {
      # warning(paste("No CDS genomic ranges found in TxDb for transcript", current_transcript_id, ". Skipping."))
      next
    }
    len_cds_from_txdb <- sum(width(cds_genomic_ranges_for_tx))
    if (len_cds_from_txdb == 0 || len_cds_from_txdb %% 3 != 0) {
      # warning(paste("Total CDS length from TxDb is 0 or not a multiple of 3 for transcript", current_transcript_id, "(Length:", len_cds_from_txdb,"). Skipping."))
      next
    }
    
    # Get 5' UTR ranges for this transcript from TxDb
    utr5_ranges_for_tx <- five_utr_by_tx[[current_transcript_id]] # This can be NULL if no 5'UTR
    len_5utr_from_txdb <- if (!is.null(utr5_ranges_for_tx)) sum(width(utr5_ranges_for_tx)) else 0
    
    # Define coding_dna_sequence based on 5'UTR length and CDS length from TxDb
    cds_start_in_full_tx_seq <- len_5utr_from_txdb + 1
    cds_end_in_full_tx_seq <- len_5utr_from_txdb + len_cds_from_txdb
    
    if (cds_end_in_full_tx_seq > length(full_transcript_dna)) {
      #   warning(paste("Calculated CDS end (", cds_end_in_full_tx_seq,
      #                 ") exceeds full transcript length (", length(full_transcript_dna),
      #                 ") for transcript", current_transcript_id, ". Check UTR/CDS annotations. Skipping."))
      next
    }
    if (cds_start_in_full_tx_seq > cds_end_in_full_tx_seq) { # Should not happen if len_cds_from_txdb > 0
      #    warning(paste("Calculated CDS start (", cds_start_in_full_tx_seq,
      #                 ") is after CDS end (", cds_end_in_full_tx_seq,
      #                 ") for transcript", current_transcript_id, ". Skipping."))
      next
    }
    
    
    coding_dna_sequence <- subseq(full_transcript_dna, start = cds_start_in_full_tx_seq, end = cds_end_in_full_tx_seq)
    
    # Final check on the extracted CDS sequence
    if (length(coding_dna_sequence) == 0 || length(coding_dna_sequence) %% 3 != 0) {
      #   warning(paste("Extracted CDS for transcript", current_transcript_id,
      #                 "has length not divisible by 3 or is empty after subseq. Length:",
      #                 length(coding_dna_sequence), ". Original CDS length from TxDb:", len_cds_from_txdb, ". Skipping."))
      next
    }
    
    transcript_mutation_data <- list()
    num_codons <- length(coding_dna_sequence) / 3
    
    for (codon_idx in 1:num_codons) {
      codon_start_pos_in_cds <- (codon_idx - 1) * 3 + 1
      original_codon_dna <- subseq(coding_dna_sequence, start = codon_start_pos_in_cds, end = codon_start_pos_in_cds + 2)
      
      # Skip codon if it contains 'N'
      if (grepl("N", as.character(original_codon_dna), ignore.case = TRUE)) {
        # warning(paste("Original codon contains 'N':", as.character(original_codon_dna), "for T:", current_transcript_id, "Codon:", codon_idx, ". Skipping codon."))
        next
      }
      
      original_aa <- tryCatch({
        # Use if.fuzzy.codon="error" as Ns should be caught above
        as.character(Biostrings::translate(original_codon_dna, if.fuzzy.codon = "error"))
      }, error = function(e) {
        # warning(paste("Could not translate original codon:", as.character(original_codon_dna), "for T:", current_transcript_id, ". Error:", e$message, ". Skipping codon."))
        NA_character_
      })
      
      if (is.na(original_aa)){
        next
      }
      
      for (pos_in_codon in 1:3) {
        original_base <- as.character(subseq(original_codon_dna, pos_in_codon, pos_in_codon))
        # Ensure original base is not 'N' (already handled by grepl above for the whole codon)
        # If somehow an 'N' got through or if we want to be extra careful at base level:
        if (original_base == "N") next 
        
        alternative_bases <- DNA_BASES[DNA_BASES != original_base]
        
        for (alt_base in alternative_bases) {
          # Create mutated codon
          mutated_codon_dna_mut_char <- strsplit(as.character(original_codon_dna), "")[[1]]
          mutated_codon_dna_mut_char[pos_in_codon] <- alt_base
          mutated_codon_dna_mut <- DNAString(paste0(mutated_codon_dna_mut_char, collapse=""))
          
          mutated_aa <- tryCatch({
            as.character(Biostrings::translate(mutated_codon_dna_mut, if.fuzzy.codon = "error"))
          }, error = function(e) {
            # warning(paste("Could not translate mutated codon:", as.character(mutated_codon_dna_mut), "for T:", current_transcript_id, ". Error:", e$message, ". Skipping mutation."))
            NA_character_
          })
          
          if (is.na(mutated_aa)){
            next
          }
          
          mutation_type <- ifelse(original_aa == mutated_aa, "synonymous", "non-synonymous")
          
          transcript_mutation_data[[length(transcript_mutation_data) + 1]] <- data.frame(
            chromosome = chr_name,
            transcript_id = current_transcript_id,
            codon_index_in_cds = codon_idx,
            original_codon = as.character(original_codon_dna),
            mutated_codon = as.character(mutated_codon_dna_mut),
            mutation_nucleotide_change = paste0(original_base, pos_in_codon, "->", alt_base),
            original_aa = original_aa,
            mutated_aa = mutated_aa,
            mutation_type = mutation_type,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(transcript_mutation_data) > 0) {
      chromosome_results[[current_transcript_id]] <- dplyr::bind_rows(transcript_mutation_data)
    }
    current_idx <- current_idx + 1
    } # End transcript loop
  
  if (length(chromosome_results) > 0) {
    all_mutation_results[[chr_name]] <- dplyr::bind_rows(chromosome_results)
  }
  message(paste("\nFinished processing chromosome:", chr_name))
  
} # End chromosome loop

# --- 3. Combine All Results ---
message("\nCombining all results...")
final_results_df <- dplyr::bind_rows(all_mutation_results)

if (nrow(final_results_df) > 0) {
  message(paste("Analysis complete. Total potential mutations analyzed:", nrow(final_results_df)))
  # Example: View or save results
  # print(head(final_results_df))
  # write.csv(final_results_df, "codon_mutation_analysis_results_corrected.csv", row.names = FALSE)
} else {
  message("Analysis complete. No results generated (this might indicate issues with input files, filtering, or annotations).")
}