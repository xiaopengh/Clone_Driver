# install.packages(c("ggplot2", "dplyr", "tidyr", "ggpubr", "circlize", "RColorBrewer", "viridis"))
# install.packages("BiocManager")
# BiocManager::install(c("GenomicRanges", "rtracklayer", "ggbio", "GenomeInfoDb"))


library(GenomicRanges)
gr1 <- GRanges(c("chr1", "chr1", "chr1"), IRanges(c(1, 2, 20), c(10, 10, 30)))
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