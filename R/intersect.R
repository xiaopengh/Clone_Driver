library(GenomicRanges)
library(rtracklayer)
library(maftools)
library(data.table)
library(dplyr)

# ==============================================================================
# Prepare working data  
# ==============================================================================

# 1. Read the BED file of sequenced regions
exome_regions <- import("../data/gencode.v19.basic.exome.bed")


# 2. Read the GTF file of gene annotations
gene_annotations <- import("../data/gencode.v19.annotation.gtf")

# Count table
tab <- table(gene_annotations$type)
# Basic barplot
barplot(tab,
        las = 2,        # make axis labels horizontal
        col = "steelblue",
        border = "white")
tab <- table(gene_annotations$gene_name)
barplot(tab)

# Extract only individuals with type "gene"
genes <- gene_annotations[gene_annotations$type == "gene"]

# Barplot of 'gene_type' categories
tab <- table(genes$gene_type)
barplot(tab, col="lightgreen", las=2)

# Extract only individuals with type "gene" 
# also filter by gene_type == "protein_coding" & level %in% c('1', '2')
protein_genes <- genes[genes$gene_type == "protein_coding" & genes$level %in% c('1', '2')]

# 3. Read the MAF file for mutated genes
maf_data <- fread(
  input = "../data/mc3.maf",
  sep = "\t",                # MAF is tab-delimited
  header = TRUE,             # First non-comment line is the header
  fill = TRUE,               # Fill missing columns if any
  quote = "",                # Disable quoting to avoid parsing errors
  na.strings = c("NA", ".", "")  # Common NA representations
)
# str(maf_data)
# maf_obj <- read.maf(maf = maf_data)
# plotmafSummary(maf)
# mafSummary(maf)

# Get mutated genes for maf file
# mutated_gene_names <- getGeneSummary(maf)$Hugo_Symbol
mutated_gene_names <- unique(maf_data$Hugo_Symbol) 
mutated_gene_names

# Define new UCSC-style sequence levels
# ==============================================================================
# Basically it changes 1 to chr1
old <- seqlevels(exome_regions)
new <- sub("^MT$", "chrM", sub("^([0-9]+|X|Y)$", "chr\\1", old))
# Apply the renaming
names(new) <- old  # mapping from old to new
exome_regions <- renameSeqlevels(exome_regions, new)
rm(old, new)

# Sequence levels verification
exome_seqlevels <- seqlevels(exome_regions)
gene_seqlevels <- seqlevels(gene_annotations)
setdiff(exome_seqlevels, gene_seqlevels) # expected 0
setdiff(gene_seqlevels, exome_seqlevels) # expected 0
intersect(exome_seqlevels, gene_seqlevels) # expected chr1 ... chrY
rm(exome_seqlevels, gene_seqlevels)

# ==============================================================================
# This extract hugo symbols for all feature types(gene, exon, etc) in gtf file 
# ==============================================================================
# Find overlaps keyword any
hits <- findOverlaps(gene_annotations, exome_regions)
length(hits)

# Get sequenced genes during the experiment
sequenced_genes <- gene_annotations[queryHits(hits)]
exome_hit <- exome_regions[subjectHits(hits)]
gnew <- pintersect(sequenced_genes, exome_hit)


gene_names <- mcols(gnew)$gene_name %>% unique()
length(gene_names)

# Get the gene_names that doesn't appear in maf file 
zero_mutation_gene_names <- gene_names[!(gene_names %in% mutated_gene_names)]
length(zero_mutation_gene_names)

# Save a csv file for this
# Convert to data frame with one column named "hugo_symbol"
df <- data.frame(hugo_symbol = zero_mutation_gene_names)

# Write to CSV (no row names)
write.csv(df, file = "../data/zero_mutation_genes.csv", row.names = FALSE)

# ==============================================================================
# This only extract hugo symbol of individuals with type 'gene' in gtf
# ==============================================================================
# Find overlaps keyword any
hits <- findOverlaps(genes, exome_regions)
length(hits)

# Get sequenced genes during the experiment
sequenced_genes <- genes[queryHits(hits)]
exome_hit <- exome_regions[subjectHits(hits)]
gnew <- pintersect(sequenced_genes, exome_hit)

gene_names <- mcols(gnew)$gene_name %>% unique()
length(gene_names)

# Get the gene_names that doesn't appear in maf file 
zero_mutation_gene_names <- gene_names[!(gene_names %in% mutated_gene_names)]
length(zero_mutation_gene_names)

# Save a csv file for this
# Convert to data frame with one column named "hugo_symbol"
df <- data.frame(hugo_symbol = zero_mutation_gene_names)

# Write to CSV (no row names)
write.csv(df, file = "../data/zero_mutation_genes.csv", row.names = FALSE)

# ==============================================================================
# extract hugo symbol of individuals with type 'gene' in gtf
# + filtered by gene_type == "protein_coding" & level %in% c('1', '2')
# ==============================================================================
# Find overlaps keyword any
hits <- findOverlaps(protein_genes, exome_regions)
length(hits)

# Get sequenced genes during the experiment
sequenced_genes <- protein_genes[queryHits(hits)]
exome_hit <- exome_regions[subjectHits(hits)]
gnew <- pintersect(sequenced_genes, exome_hit)

gene_names <- mcols(gnew)$gene_name %>% unique()
length(gene_names)

# Get the gene_names that doesn't appear in maf file 
zero_mutation_gene_names <- gene_names[!(gene_names %in% mutated_gene_names)]
length(zero_mutation_gene_names)

# Save a csv file for this
# Convert to data frame with one column named "hugo_symbol"
df <- data.frame(hugo_symbol = zero_mutation_gene_names)

# Write to CSV (no row names)
write.csv(df, file = "../data/zero_mutation_genes.csv", row.names = FALSE)
