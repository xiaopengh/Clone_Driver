library(Biostrings)
library(rtracklayer)
library(dplyr)
library(stringr)
library(progress)



# Define the file paths
gtf_file_path <- "../data/gencode.v48.annotation.gtf"
fasta_file_path <- "../data/gencode.v48.transcripts.fa"

# Standard DNA bases
DNA_BASES <- c("A", "C", "G", "T")

# Read GTF file
message("Reading GTF file...")
gtf_data <- tryCatch({
  rtracklayer::import(gtf_file_path)
}, error = function(e) {
  stop("Error reading GTF file: ", e$message)
})
message("GTF file loaded.")

# Read transcript FASTA file
message("Reading FASTA file...")
transcript_seqs <- tryCatch({
  Biostrings::readDNAStringSet(fasta_file_path)
}, error = function(e) {
  stop("Error reading FASTA file: ", e$message)
})
message("FASTA file loaded.")

# Filter for CDS features and ensure essential columns are present
cds_features <- gtf_data[gtf_data$type == "CDS"]
# Remove NA transcript_ids if any
cds_features <- cds_features[!is.na(cds_features$transcript_id)]
message(paste("Found", length(cds_features), "CDS features."))

# Clean transcript IDs FASTA headers - GTF transcript_id format
# GENCODE FASTA headers are often like ">ENST0000012345.6|ENSG00000XYZ.1|..."
# We want the "ENST0000012345.6" part.
original_fasta_names <- names(transcript_seqs)
cleaned_fasta_names <- sapply(strsplit(original_fasta_names, "\\|"), `[`, 1)
cleaned_fasta_names <- sub(" .*", "", cleaned_fasta_names) # Remove anything after a space if ID is first
names(transcript_seqs) <- cleaned_fasta_names
message("FASTA sequence names cleaned.")

# ==============================================================================

transcript_seqs <- transcript_seqs[width(transcript_seqs) %% 3 == 0] # Ensure all sequences are divisible by 3

# Using chromosomes present in the CDS features
chromosomes <- unique(as.character(seqnames(cds_features)))
message(paste("Chromosomes to process:", paste(chromosomes, collapse=", ")))

# --- Main Processing Loop ---
all_mutation_results <- list()

for (chr_name in chromosomes) {
  message(paste("\nProcessing chromosome:", chr_name))

  # Filter CDS features for the current chromosome
  cds_on_chr <- cds_features[seqnames(cds_features) == chr_name]
  if (length(cds_on_chr) == 0) {
    message(paste("No CDS features found on chromosome", chr_name, ". Skipping."))
    next
  }

  # Get unique transcript IDs on this chromosome that have CDS features
  transcripts_on_chr <- unique(cds_on_chr$transcript_id)
  if (length(transcripts_on_chr) == 0) {
    message(paste("No transcripts with CDS features found on chromosome", chr_name, "after filtering. Skipping."))
    next
  }

  message(paste("Found", length(transcripts_on_chr), "transcripts with CDS on", chr_name))

  # Initialize progress bar for transcripts on this chromosome
  pb <- progress_bar$new(
    format = paste("  Chromosome", chr_name, "[:bar] :percent ETA: :eta (:current/:total transcripts)"),
    total = length(transcripts_on_chr),
    width = 80
  )

  chromosome_results <- list()

  for (current_transcript_id in transcripts_on_chr) {
    pb$tick() # Update progress bar

    # Check if the transcript sequence exists in the FASTA data
    if (!current_transcript_id %in% names(transcript_seqs)) {
      warning(paste("Transcript", current_transcript_id, "not found in FASTA file. Skipping."))
      next
    }

    coding_dna_sequence <- transcript_seqs[[current_transcript_id]]

    if (length(coding_dna_sequence) == 0 || length(coding_dna_sequence) %% 3 != 0) {
      warning(paste("Transcript", current_transcript_id, "has CDS length not divisible by 3 or is empty. Length:", length(coding_dna_sequence), ". Skipping."))
      next
    }

    transcript_mutation_data <- list()
    num_codons <- length(coding_dna_sequence) / 3

    for (codon_idx in 1:num_codons) {
      codon_start_pos_in_cds <- (codon_idx - 1) * 3 + 1
      original_codon_dna <- subseq(coding_dna_sequence, start = codon_start_pos_in_cds, end = codon_start_pos_in_cds + 2)

      # Translate original codon, handle potential errors (e.g., Ns in codon)
      original_aa <- tryCatch({
        as.character(Biostrings::translate(original_codon_dna, if.fuzzy.codon = "solve"))
      }, error = function(e) { NA_character_ })

      if (is.na(original_aa)) {
        # warning(paste("Could not translate original codon:", as.character(original_codon_dna), "for transcript", current_transcript_id, ". Skipping codon."))
        next
      }

      # Generate 9 point mutations for this codon
      for (pos_in_codon in 1:3) {
        original_base <- as.character(subseq(original_codon_dna, pos_in_codon, pos_in_codon))
        alternative_bases <- DNA_BASES[DNA_BASES != original_base]

        for (alt_base in alternative_bases) {
          mutated_codon_dna <- original_codon_dna
          mutated_codon_char_vec <- strsplit(as.character(mutated_codon_dna), "")[[1]]
          mutated_codon_char_vec[pos_in_codon] <- alt_base
          mutated_codon_dna_mut <- DNAString(paste0(mutated_codon_char_vec, collapse = ""))

          mutated_aa <- tryCatch({
            as.character(Biostrings::translate(mutated_codon_dna_mut, if.fuzzy.codon = "solve"))
          }, error = function(e) { NA_character_ })

          if (is.na(mutated_aa)) {
            # warning(paste("Could not translate mutated codon:", as.character(mutated_codon_dna_mut), "for transcript", current_transcript_id, ". Skipping mutation."))
            next
          }

          mutation_type <- ifelse(original_aa == mutated_aa, "synonymous", "non-synonymous")

          # Store result
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
  } # End transcript loop

  if (length(chromosome_results) > 0) {
    all_mutation_results[[chr_name]] <- dplyr::bind_rows(chromosome_results)
  }
  message(paste("\nFinished processing chromosome:", chr_name))

} # End chromosome loop