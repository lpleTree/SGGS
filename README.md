# SGGS: Systematic GWAS-Assisted Genomic Selection Framework

**SGGS** is an open-source R package that optimizes complex trait prediction by integrating GWAS-derived priors into genomic selection (GS) models. It synthesizes consensus signals from seven GWAS algorithms, deploys a data-driven P-value gradient scan for adaptive prior selection, and evaluates diverse prior-assisted strategies across multiple GS models.

## Framework Overview

SGGS implements a modular end-to-end pipeline:

1. **Quality Control** – Phenotypic profiling and genomic QC with publication-quality visualizations
2. **Multi-Model GWAS** – Parallel association mapping across seven GAPIT3 algorithms (GLM, MLM, CMLM, MLMM, SUPER, FarmCPU, BLINK)
3. **Prior Stratification** – Gradient P-value-based SNP stratification with consensus filtering
4. **Genomic Prediction** – Cross-validated GS modeling with 9 models
5. **Comparative Evaluation** – Benchmarking GWAS-assisted vs. baseline strategies

## Repository Structure

`
SGGS/
├── README.md
├── .gitignore
├── data/
│   ├── catalpa/       # Catalpa bungei (楸树)
│   ├── maize/         # Zea mays CUBIC population (玉米)
│   └── wheat/         # Triticum aestivum CIMMYT lines (小麦)
└── scripts/
    ├── Differ_SNP_data_sets_11_7LPL.R  # Differential SNP analysis
    └── gapit_functionsw.txt            # GAPIT utility functions
`

## Data Description

| Species     | Population              | Samples | Traits                          | Genotype Format |
|-------------|-------------------------|---------|---------------------------------|-----------------|
| Catalpa     | C. bungei germplasm     | 180     | TH, DBH, DR                     | Numeric (-1/0/1)|
| Maize       | CUBIC inbred lines      | 1,404   | DTA, EW, EL                     | Numeric (-1/0/1)|
| Wheat       | CIMMYT inbred lines     | 599     | GY (4 environments)             | Numeric (-1/0/1)|

## GS Models

Nine genomic selection models are supported:
- **BLUP-based**: rrBLUP, GBLUP, cBLUP, sBLUP
- **Bayesian**: BRR, BayesA, BayesB, BayesCπ, Bayesian Lasso

## GWAS-Assisted Strategies

Eight strategies for incorporating GWAS priors:
- Baseline: full-SNP GBLUP/cBLUP, P-value gradient scanning
- Fixed-effect: peak-SNP integration, multi-model union, consensus filtering

All strategies use a strict training-set-only design to prevent data leakage.

## Software Dependencies

- **R** ≥ 4.0
- **GAPIT3** (for GWAS)
- **BGLR**, **rrBLUP** (for GS)
- Supporting packages: igmemory, EMMREML, multtest, gplots, scatterplot3d, genetics, pe, lme4

## Citation

Li P, Yang X, Yin T, et al. SGGS: A Systematic GWAS-Assisted Genomic Selection Framework for Enhancing Prediction Accuracy Across Diverse Plant Species. *In preparation*.

## License

All rights reserved. For academic use only.
