library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(circlize)
library(ggbio)
library(GenomeInfoDb)
library(RColorBrewer)
library(viridis)

getwd()

# 1. Read the BED file of sequenced regions
exome_regions <- import("../data/gencode.v19.basic.exome.bed")

# 2. Read the GTF file of gene annotations
gene_annotations <- import("../data/gencode.v19.annotation.gtf")

# Filter for only gene features
genes <- gene_annotations[gene_annotations$type == "gene"]

# Create a function to save plots
save_plot <- function(plot, filename, width = 10, height = 6) {
  ggsave(filename, plot, width = width, height = height, dpi = 300)
}

# 1. Distribution of gene types
gene_type_dist <- ggplot(as.data.frame(genes), aes(x = gene_type)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Gene Types",
       x = "Gene Type",
       y = "Count")
gene_type_dist

# 2. Gene length distribution
gene_lengths <- width(genes)
gene_length_dist <- ggplot(data.frame(length = gene_lengths), aes(x = length)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "Distribution of Gene Lengths",
       x = "Gene Length (log scale)",
       y = "Count")
gene_length_dist

# 3. Chromosome coverage
chr_coverage <- ggplot(as.data.frame(genes), aes(x = seqnames)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Coverage by Chromosome",
       x = "Chromosome",
       y = "Number of Genes")
chr_coverage

# 4. Strand distribution
strand_dist <- ggplot(as.data.frame(genes), aes(x = strand)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Distribution of Genes by Strand",
       x = "Strand",
       y = "Count")
strand_dist

# 5. Gene status distribution
status_dist <- ggplot(as.data.frame(genes), aes(x = gene_status)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Gene Status",
       x = "Gene Status",
       y = "Count")
status_dist

# 6. Create a summary statistics table
summary_stats <- data.frame(
  Total_Genes = length(genes),
  Unique_Gene_Types = length(unique(genes$gene_type)),
  Unique_Chromosomes = length(unique(seqnames(genes))),
  Mean_Gene_Length = mean(width(genes)),
  Median_Gene_Length = median(width(genes))
)
summary_stats

# 7. Create a genomic coverage plot for a subset of chromosomes
# Select first 5 chromosomes for visualization
top_chromosomes <- names(sort(table(seqnames(genes)), decreasing = TRUE))
genes_subset <- genes[seqnames(genes) %in% top_chromosomes]

coverage_plot <- ggplot(as.data.frame(genes_subset),
                        aes(x = start, y = seqnames, color = strand)) +
  geom_segment(aes(xend = end, yend = seqnames), linewidth = 5) +
  theme_minimal() +
  labs(title = "Gene Coverage on Human Chromosomes",
       x = "Genomic Position",
       y = "Chromosome")
coverage_plot

# 8. Create a multi-panel scientific visualization

# Prepare the genomic data
genes_gr <- genes
seqlevelsStyle(genes_gr) <- "UCSC"  # Convert to UCSC style (e.g., "chr1")

# Calculate gene density per chromosome
chr_density <- data.frame(
  chromosome = seqnames(genes_gr),
  position = start(genes_gr),
  gene_type = genes_gr$gene_type,
  strand = as.character(strand(genes_gr))
) %>%
  group_by(chromosome) %>%
  mutate(
    position_scaled = scale(position),
    bin = cut(position, breaks = 100)
  )

# Main visualization with multiple panels
scientific_plot <- ggplot(chr_density) +
  # Panel A: Chromosome-wise gene density
  geom_violin(aes(x = as.factor(chromosome), y = position),
              fill = "lightblue", alpha = 0.6) +
  facet_grid(. ~ "A: Gene Density Distribution") +

  # Panel B: Gene type distribution
  geom_bar(aes(x = as.factor(chromosome), fill = gene_type),
           position = "fill") +
  facet_grid(. ~ "B: Gene Type Composition") +

  # Panel C: Strand bias
  geom_bar(aes(x = as.factor(chromosome), fill = strand),
           position = "fill") +
  facet_grid(. ~ "C: Strand Distribution") +

  # Common theme elements
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  scale_fill_viridis_d() +
  labs(
    title = "Multi-dimensional Analysis of Gene Distribution",
    x = "Chromosome",
    y = "Position/Proportion",
    fill = "Category"
  )
scientific_plot

# Create statistical summary visualization
gene_stats <- chr_density %>%
  group_by(chromosome) %>%
  summarise(
    gene_count = n(),
    type_diversity = n_distinct(gene_type),
    strand_ratio = sum(strand == "+") / n()
  )

stats_plot <- ggplot(gene_stats) +
  geom_point(aes(x = gene_count, y = type_diversity,
                 size = strand_ratio, color = chromosome)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal() +
  labs(
    title = "Chromosome Characteristics Overview",
    x = "Number of Genes",
    y = "Gene Type Diversity",
    size = "Strand Ratio (+/-)",
    color = "Chromosome"
  )
stats_plot

# Generate numerical summaries
summary_stats <- chr_density %>%
  group_by(chromosome) %>%
  summarise(
    n_genes = n(),
    n_types = n_distinct(gene_type),
    density = n() / (max(position) - min(position)),
    strand_balance = sum(strand == "+") / n()
  ) %>%
  arrange(desc(n_genes))

# Print summary statistics
print(summary_stats)

# Print summary information
cat("\nGene Annotation Summary:\n")
cat("------------------------\n")
cat("Total number of genes:", length(genes), "\n")
cat("Number of unique gene types:", length(unique(genes$gene_type)), "\n")
cat("Number of chromosomes:", length(unique(seqnames(genes))), "\n")
cat("Mean gene length:", mean(width(genes)), "bp\n")
cat("Median gene length:", median(width(genes)), "bp\n")