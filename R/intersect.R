if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("rtracklayer") # For reading GTF files

library(GenomicRanges)
library(rtracklayer)

# 1. Read the BED file of sequenced regions
exome_regions <- import("../data/gencode.v19.basic.exome.bed")

# 2. Read the GTF file of gene annotations
gene_annotations <- import("../data/gencode.v19.annotation.gtf")

# Filter for only gene features
genes <- gene_annotations[gene_annotations$type == "gene"]

summary(genes)
length(genes)
names(mcols(genes)) # Use mcols() to access metadata columns
seqinfo(genes)
str(genes)

# 3. Perform the intersection
overlapping_genes <- intersect(genes, exome_regions)

# 'overlapping_genes' is now a GRanges object containing the regions of genes
# that overlap with the exome regions.

# 4. Extract the gene names (assuming the GTF has a 'gene_name' attribute)
overlapping_gene_names <- unique(overlapping_genes$gene_name)

# Print the number of overlapping genes
cat("Number of genes overlapping with exome regions:", length(overlapping_gene_names), "\n")

# Get the genes that fall INSIDE the exome regions (i.e., completely contained),
overlaps_within <- findOverlaps(genes, exome_regions, type = "within")
genes_within_exome <- genes[queryHits(overlaps_within)]
unique_genes_within <- unique(genes_within_exome$gene_name)
cat("Number of genes completely within exome regions:", length(unique_genes_within), "\n")

# Get the genes that have ANY overlap (as intersect):
overlaps_any <- findOverlaps(genes, exome_regions)
overlapping_genes_any <- genes[queryHits(overlaps_any)]
unique_overlapping_genes_any_names <- unique(overlapping_genes_any$gene_name)
cat("Number of genes with any overlap with exome regions (using findOverlaps):", length(unique_overlapping_genes_any_names), "\n")
