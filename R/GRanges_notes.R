library(GenomicRanges)

gtf_dt <- as.data.table(gtf_data)
str(gtf_dt)

# Access functions for GRange object

# Basic information about the GRange List data
length(gtf_data)
head(gtf_data)
seqnames(gtf_data) # Different sequence of chromosome level
ranges(gtf_data) # Range of the data without metadata columns

# More details about the GRange List data
start(gtf_data) # Start position of the range
end(gtf_data) # End position of the range
strand(gtf_data) # Strand information
names(gtf_data) # Names of the ranges typically null 

# metadata columns (the columns that are not part of the range(columns 1 : 5))
mcols(gtf_data) # Metadata columns typeof(...) = DataFrame


gtf_data[mcols(gtf_data)$type == "gene"]
gtf_data[mcols(gtf_data)$type != "gene"]


mcols(gtf_data)$type %>% table()


