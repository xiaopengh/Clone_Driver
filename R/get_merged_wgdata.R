# This suppose to do the data cleaning that extends the current analysis
load("wgdata.RData")
purity_data <- fread("../data/consensus.20170218.purity.ploidy.txt")

# Ensure the matching column names
setnames(purity_data, "samplename", "Tumor_Sample_Barcode")

# Merge purity data with wgdata
merged_dt <- merge.data.table(wgdata, purity_data, by = "Tumor_Sample_Barcode", all.x = TRUE)

# Omit rows with whole-genome duplication and wgd_uncertainty
merged_dt <- merged_dt[wgd_status == "no_wgd" & wgd_uncertain == FALSE]

# Filter out rows based on purity and CCF
merged_dt <- merged_dt[
  (CNA == 2) & 
  (Purity >= 0.5) &
  (CCF >= 0.2)
  ]

rm(wgdata, purity_data)
gc()
