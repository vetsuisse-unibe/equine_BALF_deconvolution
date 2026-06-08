# BALF Deconvolution Analysis

Computational deconvolution of equine bronchoalveolar lavage fluid (BALF) bulk RNA-seq using BayesPrism, with validation against matched single-cell RNA sequencing (scRNA-seq).

## Overview

This repository contains the analysis scripts for the manuscript:

> *[Manuscript title]*

The pipeline:
1. Deconvolves bulk RNA-seq into cell-type proportions using BayesPrism with a matched scRNA-seq reference
2. Validates estimated proportions against experimental scRNA-seq (cell counts, mRNA-corrected, and mRNA-weighted)
3. Compares cell-type-specific differential expression (SEA vs CTL) between deconvoluted and experimental scRNA-seq data

## Repository Structure

```
├── scripts/
│   ├── 01_BayesPrism_deconvolution.R          # Main deconvolution and validation
│   ├── 02_DEG_comparison.R                     # DEG comparison: scRNA-seq vs deconvolution
│   ├── 03_Figure4_with_cytology.R              # Figure 4: composition with cytology overlay
│   ├── 04_generate_DEG_tables.R                # Generate filtered DEG summary tables
│   ├── 05_extract_per_cell_umi.R               # Per-cell UMI extraction (input to 06/07)
│   ├── 06_mRNA_condition_specific_scaling.R    # mRNA correction: CTL vs SEA sensitivity (Suppl. Table S6)
│   ├── 07_mRNA_LOSO.R                          # mRNA correction: leave-one-sample-out sensitivity (Suppl. Table S6)
│   └── 08_cibersortx_comparison.R              # Second-method comparison: CIBERSORTx vs BayesPrism
├── data/
│   ├── sample_id_mapping.csv                   # Sequencing ID to sample ID mapping
│   ├── per_cell_umi.csv                        # Per-cell UMI + cell type + sample + condition (output of script 05)
│   ├── validation_real_bulk_matched.csv        # Raw BayesPrism theta + scRNA-seq cell-count truth
│   ├── validation_real_bulk_all_methods.csv    # BayesPrism (raw/corrected) + cell-count and mRNA-weighted truth
│   └── CIBERSORTx_Job3_Results.csv             # CIBERSORTx relative-mode fractions (15 bulk samples)
├── results/                                    # Created by scripts 06-08
└── README.md
```


### Sample ID Mapping

Bulk RNA-seq samples were sequenced with internal IDs (e.g., S_3113) and renamed to anonymized IDs (e.g., A30) for analysis. The `A` prefix is stripped in scripts to match scRNA-seq sample numbering. See `data/sample_id_mapping.csv` for the full mapping.

Of 16 bulk RNA-seq samples, 11 had matched scRNA-seq data (6 SEA, 5 CTL). One control sample (A61) was excluded due to insufficient RNA yield.

## Running the Analysis

### Prerequisites

R version 4.4.1 with the following packages:

```r
# Bioconductor
BiocManager::install(c("DESeq2", "BayesPrism"))

# CRAN
install.packages(c("Seurat", "tidyverse", "ggplot2", "cowplot",
                    "gridExtra", "ggrepel", "pheatmap", "readxl",
                    "VennDiagram", "ggbreak"))
```

### Execution Order

```bash
# Step 1: Run BayesPrism deconvolution and validation (generates Figures 1-4, Tables 1-6)
# Note: This step takes ~1 hour due to Gibbs sampling
Rscript scripts/01_BayesPrism_deconvolution.R

# Step 2: DEG comparison between scRNA-seq and deconvolution (generates Figure 5, Tables 7-8)
Rscript scripts/02_DEG_comparison.R

# Step 3: Figure 4 with cytology overlay (requires DCC_cyto.xlsx)
Rscript scripts/03_Figure4_with_cytology.R

# Step 4: Generate filtered DEG summary tables (requires output from Step 2)
Rscript scripts/04_generate_DEG_tables.R

# Step 5: (optional) Regenerate data/per_cell_umi.csv from the scRNA-seq object.
# Requires sc22_all_seed.rds (~5 GB; see data accession in the manuscript).
# A pre-computed data/per_cell_umi.csv is already provided, so Steps 6-7 can be run directly.
Rscript scripts/05_extract_per_cell_umi.R

# Step 6: mRNA-content correction sensitivity - condition-specific (CTL vs SEA) scaling (Suppl. Table S6)
Rscript scripts/06_mRNA_condition_specific_scaling.R

# Step 7: mRNA-content correction sensitivity - leave-one-sample-out cross-validation (Suppl. Table S6)
Rscript scripts/07_mRNA_LOSO.R

# Step 8: Second-method comparison - CIBERSORTx vs BayesPrism against both references
Rscript scripts/08_cibersortx_comparison.R
```

### Configuration

All scripts use a `data_dir` variable at the top (defaults to `"data"`). Update this path if your input files are located elsewhere.

## Session Info

```
R version 4.4.1 (2024-06-14)
Platform: x86_64-apple-darwin20
Running under: macOS 15.6

Key packages:
- BayesPrism 2.2.2
- Seurat 5.3.0
- DESeq2 1.44.1
- tidyverse 2.0.0
- ggplot2 4.0.0
```

See `session_info.txt` for the full session information.

## License

[Add license]

## Citation

[Add citation when published]
