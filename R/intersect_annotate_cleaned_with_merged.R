library(data.table)
library(GenomicRanges)

cleaned_dt <- fread("../data/ProducedData/cleanAnnotations_hg19lift.tsv.gz")

cleaned_gr <- GRanges(
  seqnames = cleaned_dt$seqname,
  ranges   = IRanges(start = cleaned_dt$start, end = cleaned_dt$end),
  strand   = "*",
  # metadata columns go here
  cleaned_dt[, !c("seqname", "start", "end")]
)

merged_gr <- GRanges(
  seqnames = merged_dt$seqname,
  ranges   = IRanges(start = merged_dt$start, end = merged_dt$end),
  strand   = "*",
  # metadata columns go here
  merged_dt[, !c("seqname", "start", "end", "strand")]
)

hits <- findOverlaps(merged_gr, cleaned_gr, ignore.strand = TRUE)
hits <- hits[!duplicated(queryHits(hits))]  # Remove duplicates

intersected_gr <- GRanges(
  seqnames = seqnames(merged_gr)[queryHits(hits)],
  ranges   = ranges(merged_gr)[queryHits(hits)],
  strand   = strand(merged_gr)[queryHits(hits)],
  # metadata columns from merged_gr
  mcols(merged_gr)[queryHits(hits), ]
)

# add the mcols from cleaned_gr
mcols(intersected_gr) <- cbind(
  mcols(intersected_gr),
  mcols(cleaned_gr)[subjectHits(hits), ]
)

rm(hits) ; gc()
