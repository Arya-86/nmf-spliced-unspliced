# nmf-spliced-unspliced
Code for NMF cross-validation of spliced/unspliced scRNA-seq data

# Joint Modeling of Spliced and Unspliced scRNA-seq Data using NMF

This repository contains analysis scripts and figures supporting the manuscript  
**"Jointly Modeling Spliced and Unspliced Single-Cell Transcriptomes Increases Revealable Rank Relative to Standard Methods."**

All analyses were performed in **R** using the [`singlet`](https://github.com/zdebruine/singlet) package and follow a reproducible workflow across multiple 10x Genomics datasets (PBMC, Brain, Neuron, Heart).

---

## Repository Structure


```plaintext
nmf-spliced-unspliced/
├── code/                  # R scripts and notebooks
│   └── final_code-4.Rmd   # Full reproducible workflow (CV, PCA, GSEA, violin plots)
├── figures/               # Output figures (PNG, PDF)
├── data/                  # Example data or rank files (optional)
├── Table1_mean_sd.csv     # Summary of revealable NMF ranks (mean ± SD)
└── README.md
```
## Analysis Overview

### **1. Cross-validation (CV) of NMF ranks**

Script: `final_code-4.Rmd` (section: *CV figure by seeds*)

This section estimates **test reconstruction error across ranks (k)** and **replicates (seeds)** for each modality:
- Spliced
- Unspliced
- Added (s + u)
- Concatenated (s ⊕ u)

### 2. PCA Comparison

This section uses PCA to compute variance explained per component across modalities and produce scree plots.

Script snippet:

get_pca_scree()  # Helper to compute PCA variance explained
ggplot(...)      # Faceted scree plot by modality

### 3. Gene Set Enrichment Analysis (GSEA)

Performs GSEA using MSigDB C5 collections on modality-specific NMF loadings.

Key steps:

Runs fgsea on spliced/unspliced weights for each factor

Identifies enriched biological processes and contrasts spliced vs. unspliced pathways

Generates comparison plots for Factor 9 and Factor 28

### 4. Violin Plots

Shows distribution of spliced and unspliced weights per NMF factor.

ggplot(pbmc_scores, aes(...)) + geom_violin() + geom_jitter()

### 5. Heatmap and Additional Figures

Additional figures such as pbmc_cv_figure_with_concat.png visualize comparative rank performance and modality integration.

## Reproducibility
To reproduce all figures and outputs:

rmarkdown::render("code/final_code-4.Rmd")

### Required R packages:
install.packages(c("dplyr", "ggplot2", "cowplot", "Matrix", "Seurat",
                   "fgsea", "msigdbr", "irlba", "reticulate", "dplyr", 
                   "viridis", "RccpML"))
remotes::install_github("zdebruine/singlet")

## Data Availability
All datasets used in this analysis are publicly available.
Paired spliced and unspliced count matrices were obtained from:
Gorin, G., & Chari, T. (2022). 10 scRNA-seq datasets processed using 3 unspliced counting pipelines (1.0) [Data set]. CaltechDATA. https://doi.org/10.22002/D1.20030

These datasets provide processed 10x Genomics single-cell RNA-seq matrices with spliced and unspliced read annotations, used to construct modality-specific inputs for NMF analysis.


