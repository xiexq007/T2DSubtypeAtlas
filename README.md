# Single-cell profiling of pancreatic islets reveals subtype-specific molecular perturbation landscape in type 2 diabetes

![Main Figure](./main_figure.png)

## Overview

This repository contains scripts to analysis data used to generate manuscript figures presented in our manuscript. Our study provides a comprehensive single-cell molecular landscape of pancreatic islets across 
various type 2 diabetes (T2D) subtypes, identifying specific molecular perturbations in severe insulin-deficient diabetes (SIDD), severe insulin-resistant diabetes (SIRD), mild obesity-related diabetes (MOD), and mild age-related 
diabetes (MARD) groups.

---

## Repository Structure

### `scripts/`
The analysis pipeline is organized by the corresponding figures in the manuscript:
* **`Preprocessing_for_scdata.r`**: Initial quality control (including SoupX), normalization, and data integration.
* **`Figure1.r` / `Figure1.py`**: Clustering, cell type identification, and clinical feature correlations.
* **`Figure2.r`**: Compositional analysis and cell proportion shifts across subtypes.
* **`Figure3.r`**: Differential expression analysis (DEGs) and molecular disturbance profiles.
* **`Figure4.r`**: Pathway enrichment (GSEA/GSVA), radar plots, and PPI network preparations.
* **`Figure5.r`**: Advanced analysis including Slingshot trajectory, Monocle 2 lineages, and functional module scoring (Maturity/ER stress).

### `data/`
This folder contains the clinical metadata required for analysis:
* **Subtype Information**: Clinical classification for the T2D subgroups used in initial clustering.
* **Donor Metadata**: Comprehensive clinical details for all donors (ND and all four T2D subtypes).

---

## Data Availability
The raw sequencing data used in this study were obtained from the **Human Pancreas Analysis Program (HPAP)** database: [https://hpap.pmacs.upenn.edu/](https://hpap.pmacs.upenn.edu/).

### Donor Cohort Summary
Our study analyzed a large-scale cohort from HPAP. The table below summarizes the donor distribution:

| Subgroup | Total Donors (Metadata) | Donors with scRNA-seq Data |
| :--- | :---: | :---: |
| **ND** (Normal Donor) | 18 | 9 |
| **MARD** (Mild Age-Related Diabetes) | 17 | 11 |
| **MOD** (Mild Obesity-Related Diabetes) | 13 | 8 |
| **SIDD** (Severe Insulin-Deficient Diabetes) | 10 | 3 |
| **SIRD** (Severe Insulin-Resistant Diabetes) | 3 | 1 |

---

## System Requirements
* **R Version**: >= 4.1.0
* **Key R Packages**: `Seurat` (v5.0), `ComplexHeatmap`, `clusterProfiler`, `Monocle2`, `Slingshot`, `ggplot2`, `dplyr`.
* **Python Version**: >= 3.8 (for specific network/clustering scripts).

---

## Citation
If you use the code or data in this repository, please cite our paper:
> *Single-cell profiling of pancreatic islets reveals subtype-specific molecular perturbation landscape in type 2 diabetes. (Journal Name, Year)*

---

## Contact
For any questions regarding the analysis or data, please open an issue in this repository or contact the corresponding author.
