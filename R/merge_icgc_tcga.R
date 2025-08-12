# ==============================================================================
# This script recquire
# ==============================================================================
# This suppose to do the data cleaning that extends the current analysis
wg_data <- fread("../data/ProducedData/icgcPCAWG_clonality.maf.gz")
chunk <- fread("../data/ProducedData/tcgaPCAWG_clonality.maf.gz")

wg_data <- wg_data[, .(Hugo_Symbol, Chromosome, Start_position, End_position, Strand, 
                       Tumor_Sample_Barcode, t_alt_count, t_ref_count, Purity, Variant_Classification, CNA, CCF)]
chunk <- chunk[, .(Hugo_Symbol, Chromosome, Start_position, End_position, Strand, 
                   Tumor_Sample_Barcode, t_alt_count, t_ref_count, Purity, Variant_Classification, CNA, CCF)]

wg_data <- rbindlist(list(wg_data, chunk), fill = TRUE)

rm(chunk) ; gc()

purity_data <- fread("../data/consensus.20170218.purity.ploidy.txt")

# Ensure the matching column names
setnames(purity_data, "samplename", "Tumor_Sample_Barcode")

# Merge purity data with wgdata
merged_dt <- merge.data.table(wg_data, purity_data, by = "Tumor_Sample_Barcode", all.x = TRUE)

rm(wg_data, purity_data) ; gc()

# Omit rows with whole-genome duplication and wgd_uncertainty
merged_dt <- merged_dt[wgd_status == "no_wgd" & wgd_uncertain == FALSE]

# Filter out rows based on purity and CCF
merged_dt <- merged_dt[
  (CNA == 2) & 
  (Purity >= 0.5) &
  (CCF >= 0.2)
  ]

merged_dt[purity != Purity]
merged_dt[, Purity := NULL] ; gc()

# Basically it changes 1 to chr1
old <- merged_dt$Chromosome
new <- sub("^MT$", "chrM", sub("^([0-9]+|X|Y)$", "chr\\1", old))
# Apply the renaming
merged_dt$Chromosome <- new # mapping from old to new
rm(old, new) ; gc()

setnames(merged_dt, 
         old = c("Chromosome", "Hugo_Symbol", "Start_position", "End_position", "Strand"), 
         new = c("seqname", "name", "start", "end", "strand"))

