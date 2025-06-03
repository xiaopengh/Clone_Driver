# First, install Bioconductor (if not already installed)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
    "GenomicRanges",
    "rtracklayer",
    "ggbio",
    "GenomeInfoDb"
))

# Install CRAN packages
install.packages(c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "ggpubr",
    "circlize",
    "RColorBrewer",
    "viridis"
))

# Verify installations by loading all packages
lapply(c(
    "GenomicRanges",
    "rtracklayer",
    "ggplot2",
    "dplyr",
    "tidyr",
    "ggpubr",
    "circlize",
    "ggbio",
    "GenomeInfoDb",
    "RColorBrewer",
    "viridis"
), function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        stop("Package ", pkg, " not installed!")
    } else {
        message(pkg, " loaded successfully")
    }
})
