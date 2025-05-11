library(GenomicRanges)
library(rtracklayer)
library(maftools)

# 1. Read the BED file of sequenced regions
exome_regions <- import("../data/gencode.v19.basic.exome.bed")

# 2. Read the GTF file of gene annotations
gene_annotations <- import("../data/gencode.v19.annotation.gtf")

# 3. Read the MAF file for mutated genes 
maf <- read.maf(maf = "../data/mc3.maf")

# Filter for only gene features
genes <- gene_annotations[gene_annotations$type == "gene"]

# Get mutated genes for maf file
mutated_gene_names <- getGeneSummary(maf)$Hugo_Symbol %>% unique()
mutated_gene_names

length(exome_regions)
length(gene_annotations)
length(genes)
length(gene_annotations) - length(genes)

# Define new UCSC-style sequence levels
old <- seqlevels(exome_regions)
new <- sub("^MT$", "chrM", sub("^([0-9]+|X|Y)$", "chr\\1", old))

# Apply the renaming
names(new) <- old  # mapping from old to new
exome_regions <- renameSeqlevels(exome_regions, new)

# Sequence levels
exome_seqlevels <- seqlevels(exome_regions)
gene_seqlevels <- seqlevels(gene_annotations)

# Find sequence levels in exome_regions but not in gene_annotations
setdiff(exome_seqlevels, gene_seqlevels)

# Find sequence levels in gene_annotations but not in exome_regions
setdiff(gene_seqlevels, exome_seqlevels)

# Find the common sequence levels
intersect(exome_seqlevels, gene_seqlevels)

# Find overlaps keyword any
hits <- findOverlaps(gene_annotations, exome_regions)
hits

# Do we need pintersect the hits? 
sequenced_genes <- gene_annotations[queryHits(hits)] %>% unique() # Keeping only one individual
gene_names <- mcols(sequenced_genes)$gene_name  # or "gene_id" depending on the GTF

# Genes that were sequenced but not mutated
zero_mutation_genes <- sequenced_genes[!(gene_names %in% mutated_gene_names)]
zero_mutation_gene_names <- mcols(zero_mutation_genes)$gene_name %>% unique()
zero_mutation_gene_names

# Save a csv file for this
# Convert to data frame with one column named "gene_name"
df <- data.frame(hugo_symbol = zero_mutation_gene_names)

# Write to CSV (no row names)
write.csv(df, file = "../data/zero_mutation_genes.csv", row.names = FALSE)

