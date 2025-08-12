# Mutation Rate-Free Natural Selection Inference in Cancer

A computational framework for detecting cancer driver genes using clonal expansion patterns, without requiring complex mutational rate models.

## Project Overview

This project implements a novel R-based genomic analysis approach for cancer driver gene detection, building upon initial explorations in Python. The main composition of this project is the R implementation, which provides a comprehensive pipeline for mutation rate-free natural selection inference in cancer.

### The Problem

Traditional cancer driver detection methods rely on complex mutation rate models that account for:
- Trinucleotide contexts and mutational signatures
- Regional mutation rate covariates
- Local hypermutation patterns

These approaches are prone to false positives, especially in non-coding regions where modeling mutational covariates is particularly challenging.

### Our Solution

We propose a **mutation rate-free** approach that leverages the temporal dynamics of cancer evolution:

1. **Clonal vs Subclonal Classification**: Using binomial tests on variant allele frequencies (VAF) adjusted for tumor purity
2. **Cross-tumor Analysis**: Aggregating clonal/subclonal mutation counts across multiple tumors for each genomic element
3. **Statistical Modeling**: Fitting negative binomial models to detect deviations from neutral evolution
4. **Driver Detection**: Identifying genomic elements with significantly skewed clonal:subclonal ratios

### Key Innovation

Driver mutations tend to occur early in tumor development (clonally) because tumors need these alterations to become cancerous. This creates a detectable bias toward clonal mutations in driver genes when analyzed across multiple independent tumors.

## Project Structure

```
Clone_Driver/
├── R/                                              # Main R analysis pipeline
│   ├── hyp_test.R                                  # ⭐ Key statistical analysis (Poisson & quantile tests)
│   ├── apply_our_method_hyp_test.R                # ⭐ Core outlier detection method
│   ├── helpers.R                                   # Essential helper functions
│   ├── intersect_BED_GTF.R                        # Maps mutations to genes (isolated workflow)
│   ├── produce_mutation_theo_chances_per_txid.R   # Generates theoretical mutation chances
│   ├── mc3_analysis.R                             # Processes MC3 mutation data
│   ├── syn_nonsyn_analysis.R                      # Synonymous vs non-synonymous analysis
│   ├── plotting.R                                 # Result visualization (coupled with hyp_test.R)
│   ├── helpers_utr.R                              # UTR extraction functions
│   ├── extract_non_utr_range.R                    # Produces output_dt dataset
│   ├── merge_icgc_tcga.R                          # Merges ICGC and TCGA data
│   ├── intersect_annotate_cleaned_with_merged.R   # Intersects & annotates datasets
│   ├── extract_merged.R                           # Extracts merged data
│   └── data/                                       # Processed R data files
├── src/                                            # Python trial implementation
│   ├── models.py                                   # Trial statistical models
│   ├── data_handlers.py                           # Data processing utilities
│   └── visualizations.py                          # Plotting functions
├── data/                                           # Input datasets
│   ├── mc3.maf                                     # MC3 mutation data
│   ├── gencode.v19.annotation.gtf                 # Gene annotations
│   ├── cgc.tsv                                     # Cancer Gene Census
│   └── consensus.*.txt                             # Tumor purity data
├── notebook_v1.ipynb                              # Python exploration notebook
└── dashboard.py                                    # Interactive Streamlit dashboard
```

## Methodology

### 1. Data Processing Pipeline

**Input Data:**
- TCGA MC3 mutation calls (MAF format)
- ICGC PCAWG mutation data
- GENCODE v19 gene annotations
- Tumor purity estimates
- Cancer Gene Census (CGC) for validation

**Processing Steps:**
1. **Clonality Classification** (`R/helpers.R:build_mc3_summary_bytx`)
   - Filter mutations with copy number ≈ 2 (diploid regions)
   - Calculate VAF adjusted by tumor purity
   - Binomial test: H₀: VAF ~ Binom(coverage, 0.5 × purity)
   - Mutations with p < 0.05 classified as subclonal

2. **Genomic Element Aggregation**
   - Sum clonal/subclonal counts across tumors for each gene
   - Create 2D histogram of (clonal_count, subclonal_count)

### 2. Statistical Modeling

**Negative Binomial Regression** (`R/apply_new_method.R`)
```r
# Model: subclonal_count ~ clonal_count
negbin_model <- glm.nb(subclonal_count ~ clonal_count)
```

**Why Negative Binomial?**
- Accounts for overdispersion in count data
- Better fit than Poisson model (confirmed by likelihood ratio tests)
- Provides robust confidence intervals for outlier detection

**Alternative Models Tested:**
- Canonical log link: `log(μ) = β₀ + β₁x`
- Identity link: `μ = β₀ + β₁x` (Python implementation)

### 3. Statistical Methods for Outlier Detection

Our method identifies significant genes or regions by applying two complementary statistical tests:

#### I. Poisson Test

**Goal:** Detect whether the observed ratio of mutation types (e.g., clonal vs subclonal counts, or non-synonymous vs synonymous mutations) deviates significantly from the expected ratio under a Poisson process model.

**In hyp_test.R:**
1. Uses exact conditional Poisson tests for both clonal and subclonal mutation categories, comparing observed counts against expectations derived from genome-wide mutation patterns
2. Applies multiple testing corrections (Benjamini–Hochberg and Holm) and combines results using both conservative (max) and less conservative (Fisher's method) approaches

**In apply_our_method_hyp_test.R:**
1. Extends the Poisson test beyond MC3 coding data to broader genomic regions (including non-coding sites)
2. Calculates a global clonal-to-subclonal mutation ratio across the dataset
3. Tests each region's clonal/subclonal counts against this global ratio using the exact Poisson test
4. Applies FDR correction and identifies significant regions as outliers

#### II. Quantile-Based Test

**Goal:** Detect genes or regions with mutation burdens that fall in the extreme tails of the empirical distribution, without assuming a specific parametric form.

**In hyp_test.R:**
- Computes the expected distribution of mutation counts under neutral evolution assumptions
- Defines significance thresholds based on selected quantiles (e.g., upper tail 5%)
- Classifies genes or regions exceeding these thresholds as outliers

**Key Difference:**
- The **Poisson test** is a parametric test that models mutation events as independent and occurring at a constant rate, allowing explicit p-value computation and adjustment for multiple testing
- The **quantile-based test** is a non-parametric approach that flags statistical extremes purely based on the observed distribution, making fewer assumptions but potentially less powerful for small-effect signals

### 4. Validation Methods

- Hypergeometric test for enrichment in Cancer Gene Census
- ROC curve analysis using known cancer genes
- Comparison with dN/dS methods

## Key Results

### Method Performance
- **Hypergeometric p-value**: 1.08e-05 (significant enrichment in CGC genes)
- **ROC AUC**: ~0.52-0.53 for cancer gene classification
- **Superior to dN/dS**: Better performance than traditional selection ratio methods

### Biological Insights
- Early driver events show strong clonal bias
- Late driver events may appear subclonal
- Method works for both coding and non-coding regions
- Reduced false positives compared to mutation rate-dependent approaches

## Usage

### R Pipeline (Main Implementation)

The R implementation provides the complete analysis pipeline. Execute scripts in the following sequence:

#### Step 1: Data Preparation
```r
# Intersect BED mutation data with GTF annotations (isolated workflow)
source("R/intersect_BED_GTF.R")

# Merge ICGC and TCGA datasets
source("R/merge_icgc_tcga.R")

# Extract and annotate non-UTR regions
source("R/extract_non_utr_range.R")  # Produces output_dt

# Intersect, annotate, and clean datasets
source("R/intersect_annotate_cleaned_with_merged.R")  # Uses cleaned_dt
```

#### Step 2: Mutation Analysis
```r
# Load helper functions (contains extract_coding_sites, determine_SNP_effects, compute_dNdS_metrics)
source("R/helpers.R")

# Generate theoretical mutation chances per transcript
source("R/produce_mutation_theo_chances_per_txid.R")

# Process MC3 mutation data
source("R/mc3_analysis.R")

# Perform synonymous vs non-synonymous analysis
source("R/syn_nonsyn_analysis.R")
```

#### Step 3: Statistical Testing (Core Analysis)
```r
# ⭐ MOST IMPORTANT: Statistical hypothesis testing
source("R/hyp_test.R")                     # Poisson & quantile-based tests
source("R/apply_our_method_hyp_test.R")    # Core outlier detection method

# Generate visualization plots
source("R/plotting.R")                     # Coupled with hyp_test.R
```

#### Dependencies & Notes
- **External datasets required**: MC3 mutation data, GENCODE GTF annotations, tumor purity data, Cancer Gene Census
- **Key helper scripts**: Functions in `helpers.R` are essential for multiple steps
- **Data flow**: `extract_non_utr_range.R` produces `output_dt` → annotated by Mike with features (cds, enhancer_K562, lncrna_exon, mrna_intron, pri_mirna, promoter, utr3_exon, utr5_exon) → results in `cleaned_dt` used by `intersect_annotate_cleaned_with_merged.R`
- **Most critical analysis**: Steps in `hyp_test.R` and `apply_our_method_hyp_test.R` represent the core methodological contribution

### Python Implementation

The Python code provides a trial implementation of the new method, allowing exploration of identity link in negative binomial regression and includes an interactive Streamlit dashboard:

```bash
# Install dependencies
pip install -r requirements.txt

# Run Jupyter notebook for method exploration
jupyter notebook notebook_v1.ipynb

# Launch interactive dashboard
streamlit run dashboard.py
```

### Key Functions

**R Functions (helpers.R):**
- `extract_coding_sites()`: Extract CDS regions from GTF annotations  
- `determine_SNP_effects()`: Calculate synonymous/non-synonymous mutation potential
- `compute_dNdS_metrics()`: Calculate dN/dS selection ratios
- Used in `produce_mutation_theo_chances_per_txid.R`

**Other Key Functions:**
- **hyp_test.R**: Poisson and quantile-based statistical tests
- **apply_our_method_hyp_test.R**: Core outlier detection for broader genomic regions
- **plotting.R**: Visualization functions (coupled with hyp_test.R)

**Python Functions:**
- `fit_negbin_model()`: Trial implementation of negative binomial regression with identity link exploration
- `calculate_summary_stats()`: Compute descriptive statistics for dashboard
- `negbin_quantiles()`: Calculate confidence intervals for outlier detection

## Implementation Notes

### R vs Python

- **R pipeline**: Complete implementation of the published method, handles genomic data processing, performs chromosome-by-chromosome analysis
- **Python pipeline**: Trial implementation for exploring the new method, specifically testing identity link in negative binomial regression, includes interactive dashboard for data exploration

The R folder contains the continuation and extension of the initial Python trial, implementing the full genomic pipeline for mutation rate-free cancer driver detection.

## Dependencies

**R Packages:**
```r
library(data.table)
library(rtracklayer)
library(Biostrings)
library(MASS)           # For glm.nb()
library(GenomicRanges)
library(dplyr)
library(stringr)
library(progress)
```

**Python Packages:**
```python
import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats
import matplotlib.pyplot as plt
import streamlit as st
import pickle
```

## Data Requirements

1. **MC3 MAF file**: TCGA mutation calls with VAF and coverage information
2. **GENCODE GTF**: Gene annotations (v19) for genomic element definitions  
3. **Tumor purity data**: Consensus purity estimates for VAF adjustment
4. **Reference genome**: GRCh37/GRCh38 for sequence context analysis
5. **Cancer Gene Census**: For method validation and benchmarking

## Workflow Overview

1. **Data Preprocessing** (`main.R`)
   - Load GTF annotations and reference genome
   - Process MC3 mutation data
   - Extract coding sites for each chromosome

2. **Mutation Analysis** (`helpers.R`)
   - Classify mutations as synonymous/non-synonymous
   - Determine clonal vs subclonal status using binomial tests
   - Aggregate counts across tumors by genomic element

3. **Statistical Modeling** (`apply_new_method.R`)
   - Fit negative binomial regression models
   - Calculate confidence intervals and identify outliers
   - Perform multiple testing correction

4. **Validation** (`hyp_test.R`, `syn_nonsyn_analysis.R`)
   - Compare with Cancer Gene Census
   - ROC curve analysis
   - Benchmark against dN/dS methods

## Citation

If you use this method, please cite:

```
Zhang, X. (2025). Mutation rate-free natural selection inference in cancer. 
Centre for Genomic Regulation - Weghorn Lab.
```

## Contact

**Author**: Xiaopeng Zhang  
**Supervisor**: Miguel A. Cortés G.  
**Institution**: Centre for Genomic Regulation - Weghorn Lab

For questions about implementation or methodology, please refer to the research presentation and code documentation.

## License

This project is available for academic and research use. Please contact the authors for commercial applications.