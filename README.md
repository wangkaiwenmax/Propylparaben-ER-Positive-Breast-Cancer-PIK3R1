# Propylparaben ER-Positive Breast Cancer PIK3R1

## Overview

This repository contains R scripts supporting the study:

**Propylparaben Increases the Risk of Estrogen Receptor–Positive Breast Cancer via PIK3R1: A Comprehensive Integrative Analysis**

The project integrates network toxicology, TCGA transcriptomic data, machine learning, Mendelian Randomization (MR), molecular docking, and immune infiltration analysis to investigate the carcinogenic effects of propylparaben (PP) in estrogen receptor–positive breast cancer.

---

## Workflow

* Identification of PP-related targets (network toxicology)
* TCGA differential expression analysis
* PPI network construction
* Machine learning–based core gene screening
* ROC analysis for diagnostic performance
* Mendelian Randomization for causal inference
* Molecular docking between PP and PIK3R1
* Immune infiltration analysis using ssGSEA

---

## Requirements

### Environment
* **R (≥ 4.2.0)**

### Required packages
* `limma`
* `GSVA`
* `GSEABase`
* `data.table`
* `clusterProfiler`
* `org.Hs.eg.db`
* `ggplot2`
* `pROC`
* `TwoSampleMR`
* `STRINGdb`
* `reshape2`

> **Note:** Additional packages may be required depending on your specific environment.

---

## Usage

### 1. Clone this repository
```bash
git clone https://github.com/yourname/propylparaben-erpos-breast-cancer-pik3r1.git

```

### 2. Open R and run:

```r
source("main.R")

```

*Please modify file paths inside the script before execution.*

---

## Input Data

* TCGA breast cancer expression matrix
* PP-related target genes
* eQTL summary statistics
* GWAS summary data for ER-positive breast cancer

> **Notice:** Due to privacy restrictions, raw data are not included in this repository.

---

## Output

* Differentially expressed genes
* Core genes identified by machine learning
* ROC curves
* MR results
* Docking scores
* Immune infiltration scores

*All results will be saved locally after execution.*

---

## Code Availability

All scripts used in this study are publicly available in this repository to ensure reproducibility.

## Citation

If you use this code, please cite:

> **Propylparaben Increases the Risk of Estrogen Receptor–Positive Breast Cancer via PIK3R1: A Comprehensive Integrative Analysis**

## Contact

For questions, please contact the corresponding author.
