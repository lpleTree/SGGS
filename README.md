# SGGS: Systematic GWAS-Assisted Genomic Selection Framework

[![R](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-Academic%20Use%20Only-lightgrey)]()

## Overview

**SGGS** (Systematic GWAS-assisted Genomic Selection) is an open-source R framework that enhances complex trait prediction by systematically integrating GWAS-derived priors into genomic selection (GS) models. Unlike conventional single-algorithm approaches, SGGS synthesizes **cross-algorithm consensus signals** from seven GWAS models and deploys a **data-driven P-value gradient scan** to adaptively select optimal prior information—eliminating subjective threshold decisions and intrinsic data leakage.

Validated across four diverse plant species (Catalpa, Norway spruce, wheat, and maize), SGGS provides a modular, reproducible pipeline from raw genotype/phenotype data to publication-ready comparative evaluation.

## Key Features

- **7-GWAS-model consensus** — GLM, MLM, CMLM, MLMM, SUPER, FarmCPU, and BLINK run in parallel via GAPIT3
- **Data-driven prior selection** — Gradient P-value scanning (10 thresholds from 0.01 to 1.0) replaces arbitrary cutoffs
- **Training-set-only design** — GWAS and parameter tuning are strictly confined to the training fold; test data is fully masked, eliminating leakage
- **9 GS models** — BLUP-based (rrBLUP, GBLUP, cBLUP, sBLUP) and Bayesian (BRR, BayesA, BayesB, BayesCπ, Bayesian Lasso)
- **8 GWAS-assisted strategies** — Baseline, fixed-effect integration, and multi-model union approaches
- **Cross-species validation** — Tested on four species spanning annual crops and perennial forest trees
- **Reproducible** — 5× repeated 10-fold cross-validation with fixed random seed

## Repository Structure

```
SGGS/
├── README.md
├── .gitignore
├── .gitattributes                          # Git LFS tracking config
├── data/
│   ├── catalpa/
│   │   ├── Catalpa_geno.txt                # 180 genotypes × 44,474 SNPs (012 coding)
│   │   ├── Catalpa_pheno.txt               # 3 traits: TH, DBH, DR
│   │   └── Catalpa_SNPinfor.txt            # SNP map: SNP, Chromosome, Position
│   ├── maize/
│   │   ├── Maize_geno.txt                  # CUBIC population: 1,404 lines × 20,000 SNPs
│   │   ├── Maize_pheno.txt                 # 3 traits: DTA, EW, EL
│   │   └── Maize_SNPinfor.txt              # SNP map
│   └── wheat/
│       ├── Wheat_geno.txt                  # CIMMYT: 599 lines × 1,279 DArT SNPs
│       ├── Wheat_pheno.txt                 # Grain yield across 4 environments
│       └── Wheat_SNPinfor.txt              # SNP map
└── scripts/
    ├── Differ_SNP_data_sets_11_7LPL.R      # Main analysis pipeline
    └── gapit_functionsw.txt                # Modified GAPIT R package functions
```

## Data Description

### Genotype Data Format

All genotype files use **012 numeric coding** (tab-separated), where:
- `0` = homozygous reference allele
- `1` = heterozygous
- `2` = homozygous alternative allele

| Species   | Population          | Individuals | Markers   | Genotype File Size | LFS Tracked |
|-----------|---------------------|-------------|-----------|---------------------|-------------|
| Catalpa   | *C. bungei* germplasm  | 180         | 44,474    | ~16 MB              | No          |
| Maize     | CUBIC inbred lines  | 1,404       | 20,000    | ~54 MB              | Yes         |
| Wheat     | CIMMYT inbred lines | 599         | 1,279     | ~1.5 MB             | No          |

### Phenotype Data

| Species | Traits | Description |
|---------|--------|-------------|
| Catalpa | TH, DBH, DR | Tree height (cm), diameter at breast height (mm), disease resistance (score 0–4) |
| Maize   | DTA, EW, EL | Days to anthesis, ear weight (g), ear length (cm) |
| Wheat   | GY (W1–W4)  | Grain yield across four environments (CIMMYT global program) |

### SNP Map

Each `*_SNPinfor.txt` is tab-separated with three columns:
- `SNP` — Marker identifier
- `Chromosome` — Chromosome number
- `Position` — Physical position (bp)

## Methods Workflow

SGGS implements a modular five-stage pipeline:

### Stage 1: Quality Control & Phenotypic Profiling

Trait distributions, correlations, and missingness are characterized. Genomic QC includes MAF filtering, call rate thresholds, and LD-based pruning. All QC outputs include publication-ready visualizations.

### Stage 2: Multi-Model GWAS

Seven GAPIT3 models are executed in parallel on the full dataset:

| Model   | Type               | Key Characteristics                                          |
|---------|--------------------|--------------------------------------------------------------|
| **GLM** | General Linear     | PCA covariates only; fast baseline                           |
| **MLM** | Mixed Linear       | PCA + kinship (K); standard approach                         |
| **CMLM**| Compressed MLM     | Clustered kinship for computational efficiency               |
| **MLMM**| Multi-Locus MLM    | Stepwise forward-backward regression                         |
| **SUPER**| Settlement of MLM | P-value threshold optimization with LD-based binning         |
| **FarmCPU** | Fixed & Random | Iterative fixed-effect testing + random-effect prediction    |
| **BLINK** | Bayesian           | Bayesian-information LD Iteratively Nested Keyway            |

The model detecting the most significant associations is automatically identified as the **best model** for downstream prior construction.

### Stage 3: P-Value Gradient Scanning

From the best GWAS model, SNP subsets are constructed at 10 P-value thresholds: **0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0**. This gradient scan replaces arbitrary threshold selection with a data-driven exploration of the full prior-strength spectrum.

### Stage 4: GWAS-Assisted Genomic Selection

Within each cross-validation fold, GWAS is re-run strictly on the training set. Eight strategies are evaluated across three GS model families:

**GS Models:**

| Model    | Full Name        | Description                                      |
|----------|------------------|--------------------------------------------------|
| **cBLUP**| Compressed BLUP  | BLUP with compressed relationship matrix         |
| **gBLUP**| Genomic BLUP     | Standard genomic BLUP with G matrix              |
| **sBLUP**| Subset BLUP      | BLUP using GWAS-filtered marker subsets          |

**Integration Strategies (8 total):**

| Strategy | Category     | Description                                                     |
|----------|--------------|-----------------------------------------------------------------|
| 1        | Baseline     | Full SNP set without prior information                          |
| 2        | Baseline     | Gradient P-value SNP subsets from best model only               |
| 3        | Fixed-effect | Peak SNP (most significant) as fixed-effect covariate           |
| 4        | Fixed-effect | Peak SNP + gradient P-value subsets                             |
| 5        | Fixed-effect | Multi-model union SNPs as fixed effects                         |
| 6        | Fixed-effect | Multi-model union + gradient subsets                            |
| 7        | Fixed-effect | Consensus SNPs (≥2 models) as fixed effects                     |
| 8        | Fixed-effect | Consensus SNPs + gradient subsets                               |

### Stage 5: Comparative Evaluation

Predictive accuracy is measured as the **Pearson correlation (r)** between predicted and observed phenotypes under 5× repeated 10-fold cross-validation. Results are systematically compared across models, strategies, P-value thresholds, and species to identify optimal prior integration approaches.

## Dependencies

### R Environment

- **R** ≥ 4.0
- Recommended: R ≥ 4.3

### Required R Packages

```r
library(openxlsx)      # Excel I/O for result export
library(tidyverse)     # Data manipulation and visualization
library(data.table)    # High-performance data reading
library(ggplot2)       # Publication-quality plotting
library(caret)         # Cross-validation fold generation
library(dplyr)         # Data manipulation
library(Matrix)        # Sparse matrix operations
```

### GAPIT (Genome Association and Prediction Integrated Tool)

The pipeline uses a **customized version** of GAPIT3 functions loaded from `gapit_functionsw.txt`. Key modifications include adjusted output formatting, modified model parameter defaults, and customized GWAS result parsing. Standard GAPIT is available at [zzlab.net/GAPIT](http://zzlab.net/GAPIT/).

## Quick Start

### 1. Clone the Repository

```bash
git clone git@github.com:lpleTree/SGGS.git
cd SGGS
```

### 2. Install Git LFS (for Maize genotype file)

```bash
git lfs install
git lfs pull
```

### 3. Prepare Your Data

Place your genotype (012 coding), phenotype, and SNP map files in the appropriate `data/<species>/` directory. The expected format:

- **Genotype**: Tab-separated, first column `Taxa` (sample ID), remaining columns are SNP markers with values 0/1/2
- **Phenotype**: Tab-separated, `Taxa` + trait columns
- **SNP Map**: Tab-separated, `SNP` + `Chromosome` + `Position`

### 4. Configure the Script

Open `scripts/Differ_SNP_data_sets_11_7LPL.R` in R or RStudio and modify the data paths:

```r
GD  <- read.table("data/wheat/Wheat_geno.txt", header = TRUE)
myY <- read.table("data/wheat/Wheat_pheno.txt", header = TRUE)[, c(1,2)]  # trait column
GM  <- read.table("data/wheat/Wheat_SNPinfor.txt", header = TRUE)
```

### 5. Run the Analysis

```bash
Rscript scripts/Differ_SNP_data_sets_11_7LPL.R
```

Or source line-by-line in RStudio for interactive exploration.

### Key Parameters

At the bottom of the script, the main function call controls the analysis scope:

```r
result <- GWAS_assisted_GS(GD, GM, myY,
    models    = c("sBLUP", "gBLUP", "cBLUP"),
    p.levels  = c(0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0),
    n_repeats = 5,
    n_folds   = 10)
```

| Parameter  | Description                                              |
|------------|----------------------------------------------------------|
| `models`   | GS models to evaluate (`cBLUP`, `gBLUP`, `sBLUP`)       |
| `p.levels` | P-value thresholds for SNP subsetting                   |
| `n_repeats`| Number of cross-validation repeats (default: 5)          |
| `n_folds`  | Number of cross-validation folds (default: 10)           |

### Output

- **Phase 1**: `GAPIT.Association.Filter_GWAS_results.csv` — Significant SNPs from all seven GWAS models with per-model annotations
- **Phase 2**: `differ-stra/` directory containing prediction accuracy (Pearson *r*) for each strategy × model × P-value combination, exported as CSV

## Analysis Workflow Diagram

```
┌─────────────┐    ┌──────────────────┐    ┌───────────────────────┐
│  Genotype   │    │  7-Model GWAS    │    │  Best Model Selection │
│  Phenotype  │───▶│  (GAPIT3)        │───▶│  (most sig. signals)  │
│  SNP Map    │    │  GLM/MLM/CMLM/   │    └───────────┬───────────┘
└─────────────┘    │  MLMM/SUPER/     │                │
                   │  FarmCPU/BLINK   │                ▼
                   └──────────────────┘    ┌───────────────────────┐
                                           │  P-Value Gradient     │
                                           │  Scan (0.01 → 1.0)   │
                                           └───────────┬───────────┘
                                                       │
┌──────────────────────────────────────────────────────┘
│
▼
┌──────────────────────────────────────────────────────────────────┐
│              5× Repeated 10-Fold Cross-Validation                 │
│                                                                  │
│  ┌──────────┐   ┌──────────────┐   ┌───────────────────────────┐ │
│  │ Training │   │ GWAS re-run  │   │ 8 Strategies × 3 Models   │ │
│  │   Set    │──▶│ (fold-level) │──▶│ cBLUP / gBLUP / sBLUP     │ │
│  │  GWAS    │   │              │   │                           │ │
│  └──────────┘   └──────────────┘   └───────────┬───────────────┘ │
│                                                │                 │
│  ┌──────────┐                    ┌─────────────▼───────────────┐ │
│  │   Test   │                    │  Predict & Evaluate (r)     │ │
│  │   Set    │◀───────────────────│  per strategy/model/P-value │ │
│  │ (Masked) │                    └─────────────────────────────┘ │
│  └──────────┘                                                    │
└──────────────────────────────────────────────────────────────────┘
```

## Notes

- The pipeline processes **one trait at a time**; modify the column index in `myY` to switch traits
- GWAS is re-executed within each training fold to **prevent data leakage**; do not use genome-wide significant SNPs from the full dataset in prediction
- The `differ-stra/` output directory is created automatically; ensure write permissions in the working directory
- Random seed is fixed (`set.seed(123)`) for reproducibility
- The maize genotype file exceeds 50 MB and is tracked with **Git LFS** — run `git lfs pull` after cloning

## Citation

If you use SGGS in your research, please cite:

> Li P, Yang X, Yin T, Chen S, Lv S, Shi C, Ma W, Xie G, Du J, Zhai W, Zhang M\*, Wang J\*. SGGS: A Systematic GWAS-Assisted Genomic Selection Framework for Enhancing Prediction Accuracy Across Diverse Plant Species. *In preparation*.

And the GAPIT software:

> Wang J, Zhang Z. (2021). GAPIT Version 3: Boosting Power and Accuracy for Genomic Association and Prediction. *Genomics, Proteomics & Bioinformatics*, 19(4), 629–640.

## License

This project is shared for academic and research purposes. All rights reserved. Please contact the corresponding authors for permissions regarding commercial use or redistribution.

---

**Contact**: wangjh@caf.ac.cn | State Key Laboratory of Tree Genetics and Breeding, Chinese Academy of Forestry, Beijing, China
