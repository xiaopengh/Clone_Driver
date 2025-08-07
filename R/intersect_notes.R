# install.packages(c("ggplot2", "dplyr", "tidyr", "ggpubr", "circlize", "RColorBrewer", "viridis"))
# install.packages("BiocManager")
# BiocManager::install(c("GenomicRanges", "rtracklayer", "ggbio", "GenomeInfoDb"))


library(GenomicRanges)
gr1 <- GRanges(c("chr1 1", "chr1 1", "chr1 1"), IRanges(c(20001, 2, 20), c(20060, 10, 30)))
gr2 <- GRanges(c("chr1", "chr1", "chr1"), IRanges(c(5, 7, 25), c(6, 8, 26)))
gr1
gr2
seqlevels(gr1)
ov <- findOverlaps(gr1, gr2)
ov
queryHits(ov)
subjectHits(ov)
gnew <- pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])
gnew
gdiff <- setdiff(gr1, gnew)
gr2 <- GRanges(c("chr1", "chr1"), IRanges(c(5, 25), c(12, 50)))
gr1
gr2
ov <- findOverlaps(gr1, gr2)
ov
gnew <- pintersect(gr1[queryHits(ov)], gr2[subjectHits(ov)])
gnew
gr2 <- GRanges(c("chr1", "chr1"), IRanges(c(5, 25), c(12, 50)), c("*", "*"), gene_symbol = c("TP53", "SOX1"))
gr2
gnew
gnew$gene_symbol <- gr2$gene_symbol[subjectHits(ov)]
gnew

library(Gviz)
library(GenomicRanges)

# Extract a small representative subset (first 5 for more illustrative example)
gr <- gtf_data[1]

# Define your sequence as a DNAStringSet
seq <- DNAStringSet("ATGGTGACTG")

# Create a SequenceTrack correctly
seqTrack <- SequenceTrack(sequence = seq, chromosome = "chr1")

# Plot the track clearly marking positions
plotTracks(seqTrack, from = 1, to = width(seq))


# Create tracks
atrack <- AnnotationTrack(
  gr,
  name = "gtf"
  # fill = "steelblue",
  # col.line = "darkblue",
  # rotation.title = 0,
  # background.title = "gray90"
)

gtrack <- GenomeAxisTrack()

plotTracks(list(gtrack, atrack))



itrack <- IdeogramTrack(
  genome = "hg19",
  chromosome = as.character(seqnames(gr))[1]
)


# Codon mod by site
codon <- "ACT"

mut_map <- list(
  C = c("A", "G", "T"),
  T = c("A", "C", "G"),
  A = c("C", "G", "T"),
  G = c("A", "C", "T")
)

codon_split <- strsplit(codon, "")[[1]]

# Generate flat named vector
mutants <- unlist(lapply(1:3, function(pos) {
  ref <- codon_split[pos]
  alts <- mut_map[[ref]]
  
  codons <- sapply(alts, function(alt) {
    new_codon <- codon_split
    new_codon[pos] <- alt
    paste0(new_codon, collapse = "")
  })
  
  # Name all elements with the ref nt only
  names(codons) <- rep(ref, length(codons))
  codons
}))

# Output
mutants

library(ggplot2)
library(tidyr)
library(dplyr)

# Convert to data frame
df_mutants <- data.frame(
  codon = mutants,
  ref = names(mutants),
  stringsAsFactors = FALSE
)

# Split each codon into columns
df_split <- df_mutants %>%
  mutate(row = row_number()) %>%
  mutate(nt_list = strsplit(codon, "")) %>%
  unnest_wider(nt_list, names_sep = "") %>%
  pivot_longer(cols = starts_with("nt_list"), names_to = "pos", values_to = "nt") %>%
  mutate(pos = as.integer(gsub("nt_list", "", pos)))



library(ggplot2)
library(patchwork)

seq <- strsplit("ACT", "")[[1]]
df <- data.frame(
  pos = 1:length(seq),
  nt = seq
)

# Soft color palette mimicking your image
nt_colors <- c(
  A = "#FFCCCC",  # light red
  T = "#FFE5B4",  # light orange
  G = "#CCFFFF",  # light cyan
  C = "#CCFFCC"   # light green
)

ggplot(df_split, aes(x = pos, y = -row, label = nt, fill = nt)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 1) +
  geom_text(size = 5, fontface = "bold") +
  scale_fill_manual(values = nt_colors) +
  scale_x_continuous(breaks = 1:3, limits = c(0.5, 10.5), labels = c("", "", "")) +
  # scale_y_continuous(breaks = -1:-9, labels = names(mutants)) +
  coord_fixed(ratio = 0.9) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, hjust = 1),
    legend.position = "none"
  )

# --- add this to create background tiles for y-axis labels
label_bg_df <- df_split %>%
  group_by(row, ref) %>%
  summarise(y = -unique(row), .groups = "drop") %>%
  mutate(x = 0)  # position just before x-axis for visual

# Plot
ggplot() +
  # Background rectangles behind y-axis labels
  geom_tile(data = label_bg_df, aes(x = x, y = y, fill = ref), width = 0.8, height = 1) +
  
  # Codon tiles
  geom_tile(data = df_split, aes(x = pos, y = -row, fill = nt), width = 1, height = 1, color = "white", linewidth = 1) +
  geom_text(data = df_split, aes(x = pos, y = -row, label = nt), size = 5, fontface = "bold") +
  
  # Y-axis labels (over tile backgrounds)
  geom_text(data = label_bg_df, aes(x = x, y = y, label = ref), hjust = 1.2, size = 5, fontface = "bold", color = "black") +
  
  # Scales
  scale_fill_manual(values = nt_colors) +
  scale_x_continuous(breaks = 1:3, limits = c(-1.5, 10.5), labels = c("", "", "")) +
  scale_y_continuous(breaks = -1:-9, labels = NULL) +
  
  # Layout
  coord_fixed(ratio = 0.9) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )


ggplot(df, aes(x = pos, y = 0, label = nt, fill = nt)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 1) +
  geom_text(size = 5, fontface = "bold") +
  scale_fill_manual(values = nt_colors) +
  scale_x_continuous(
    breaks = 1 : length(seq),
    limits = c(0.5, 10 + 0.5),  # ensure space for edge tiles
    labels = c("", "", "")
  ) +
  coord_fixed(ratio = 0.9) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )


ggplot(df, aes(x = pos, y = 0, label = nt, fill = nt)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 1) +  # No border, shorter tiles
  geom_text(size = 5, fontface = "bold") +
  scale_fill_manual(values = nt_colors) +
  scale_x_continuous(breaks = 1 : length(seq), labels = c("", "", "")) +
  coord_fixed(ratio = 0.90) +  # makes it horizontal and flat
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

seq <- strsplit("ATGATGATGGTGGTGGTGACTACTACTGCG", "")[[1]]
triplet_id <- rep(1:ceiling(length(seq)/3), each = 3)[1:length(seq)]

# Custom x position: 1 unit per base + extra space every 3rd base
gap_width <- 0.3  # adjust this to control the inter-triplet gap

# Within each codon: offset 0, 1, 2; between codons: add gap
offset_in_triplet <- rep(0:2, length.out = length(seq))
extra_spacing <- (triplet_id - 1) * gap_width
xpos <- (triplet_id - 1) * 3 + offset_in_triplet + extra_spacing + 1

df <- data.frame(
  pos = xpos,
  nt = seq
)

# Your colors
nt_colors <- c(A = "#FFCCCC", T = "#FFE5B4", G = "#CCFFFF", C = "#CCFFCC")

ggplot(df, aes(x = pos, y = 0, label = nt, fill = nt)) +
  geom_tile(width = 1, height = 1, color = "white", linewidth = 1) +
  geom_text(size = 5, fontface = "bold") +
  scale_fill_manual(values = nt_colors) +
  scale_x_continuous(breaks = df$pos, labels = seq_along(seq)) +
  coord_fixed(ratio = 0.90) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )
  



