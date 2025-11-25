# Transcriptomic-Profiling-GSE65194-DGE-Limma-Pathway-Visualization-Breast-Cancer-Subtypes
Integrated Transcriptomic Profiling of GSE65194 Differential Gene Expression Limma Pathway Enrichment and Interactive Visualization of HER2 vs TNBC Breast Cancer Subtypes
# Integrated Transcriptomic Analysis Pipeline for GSE65194

## 🔬 Project Title
**Integrated Transcriptomic Profiling of GSE65194: Differential Gene Expression (Limma), Functional Enrichment, and Interactive Visualization of HER2 vs. TNBC Breast Cancer**

---

## 📋 Overview
This repository contains a complete, robust R pipeline for analyzing the **GSE65194** Affymetrix microarray dataset. The core analysis focuses on identifying key molecular differences and perturbed pathways between **HER2-positive** and **Triple-Negative Breast Cancer (TNBC)** tumor samples.

### Key Features
* **Methodology:** Uses the gold-standard **Limma** package for accurate differential expression analysis on microarray data.
* **Filtering:** Implements **smart filtering logic** to correctly isolate tumor samples (HER2 and TNBC) while excluding technical replicates, normal tissues, and contaminating cell lines.
* **Annotation:** Utilizes **`hgu133plus2.db`** (GPL570 platform) for reliable Probe ID to Entrez ID mapping.
* **Visualization:** Generates a complete set of files for both static reporting and interactive exploration.

---

## 🚀 Execution & Usage

### 1. Requirements
This script requires **R** and the following packages (automatically installed if missing):
* **Bioconductor:** `GEOquery`, `limma`, `clusterProfiler`, `ReactomePA`, `hgu133plus2.db`
* **CRAN:** `ggplot2`, `plotly`, `DT`, `heatmaply`

### 2. Run the Analysis
1.  Ensure your output path (default: `D:/DOWNLOADS/GSE65194_Results`) is accessible.
2.  Run the entire R script file.

The script will automatically download the necessary data from GEO, perform normalization, run the Limma model, and execute all visualization and pathway analysis steps.

---

## 📂 Results Output

All generated files will be automatically saved into the specified directory:
`D:/DOWNLOADS/GSE65194_Results`

The output is a hybrid collection of static and interactive files, all prefixed with **`GSE65194_`**.

### Interactive HTML Files (Open in Web Browser)
| Filename Prefix | Content | Feature |
| :--- | :--- | :--- |
| `GSE65194_1_Interactive_Table.html` | Full DGE Results Table | **Searchable & Filterable** DGE results. |
| `GSE65194_2_Interactive_3D_PCA.html` | PCA Plot | **Rotatable 3D view** of sample clustering. |
| `GSE65194_3_Interactive_Heatmap.html` | Top 50 Variable Genes | **Zoomable, hoverable** heatmap of expression data. |
| `GSE65194_4_Interactive_Volcano.html` | Volcano Plot | **Hover** to view Gene Symbol, LogFC, and FDR for each probe. |
| `GSE65194_6_Interactive_GO.html` | GO Dot Plot | **Hover** to see count and p-value for Biological Process terms. |

### Static & Data Files (PNG, CSV)
| Filename Prefix | Content | Purpose |
| :--- | :--- | :--- |
| `GSE65194_Full_DGE_Results.csv` | Master table with **LogFC, FDR, Symbol, and Entrez ID**. | Input for other tools (e.g., GSEA) and record keeping. |
| `GSE65194_Significant_DEGs.csv` | List of genes passing the threshold (FDR < 0.05, \|LogFC\| > 1). | Focused list of key results. |
| `GSE65194_2_PCA_Plot.png` | Static image of the PCA plot. | Reporting and presentations. |
| `GSE65194_7_KEGG_Results.csv` | Enriched KEGG pathways. | Functional validation of gene sets. |
| `GSE65194_6_GO_Dotplot.png` | Static image of the GO enrichment results. | Publication-quality figure. |

---

## 🧬 Methodology Summary

1.  **Data Curation:** Smart filtering isolates samples based on `title`, `source_name`, and `characteristics` columns to ensure only HER2 and TNBC **tumor tissue** samples are included (Section 2.2).
2.  **Normalization:** Quantile normalization is applied using `limma::normalizeBetweenArrays`.
3.  **DGE:** Differential expression is calculated using `limma` (Empirical Bayes method).
4.  **Annotation:** Probe IDs are translated to Gene Symbols and Entrez IDs using `hgu133plus2.db` (Section 3.3).
5.  **Enrichment:** Over-Representation Analysis (ORA) is performed on significant genes using:
    * **GO (BP):** `enrichGO`
    * **KEGG:** `enrichKEGG`
    * **Reactome:** `enrichPathway` (Section 5).
