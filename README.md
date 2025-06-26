# scRNAseq Analysis Code — Benamar et al.

This repository contains all the code used in the analysis for the manuscript:

**"Notch3+ Regulatory T cells Drive Autoimmune Neuroinflammation in Multiple Sclerosis"**  
_Benamar et al._

- **Manuscript:** [Link] (DOI: [..])
- **Contact:** Dr. Mehdi Benamar — Mehdi.Benamar@childrens.harvard.edu

---

## Overview

This repository collects the scripts for processing and analyzing single-cell RNA sequencing (scRNAseq) data as described in the manuscript.  
Analyses were performed using [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) (10x Genomics) and [Seurat](https://satijalab.org/seurat/) v4, with additional supporting R packages.

---

## Analysis Workflow

1. **Pre-processing with Cellranger**
    - `cellranger count` — Quantification of gene expression for each individual sample.
    - `cellranger aggr` — Aggregation of multiple samples.

2. **Quality Control and Processing (Seurat)**
    - See: `1_QC and processing.R`
    - Import filtered matrices, quality control, filtering, normalization, integration, and clustering.

3. **Differential Gene Expression**
    - See: `2_DE genes with pseudobulk.R`
    - Identification of differentially expressed genes using pseudobulk analysis.

4. **Visualization and Figure Generation**
    - See: `3_Plots.R`
    - Generation of all figures and plots included in the manuscript.

---

## Data Availability

- **Raw sequencing data:** Available at GEO/SRA: [GSE253602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253602)
- **Processed Seurat object:** Available on [Zenodo](https://zenodo.org/) (DOI: ...)

---

## Requirements

All R package dependencies and the R version used for this analysis are specified in the `renv.lock` file included in this repository.

To recreate the exact software environment, install the [renv](https://rstudio.github.io/renv/) package in R and run:

```r
install.packages("renv")
renv::restore()
```

This will install the correct R packages and versions as used for the original analysis.

> **Note:**  
> The analysis was performed using R version 4.3.2  
> For details about all package versions, please refer to the `renv.lock` file.

---

## Usage

1. Clone the repository:
    ```bash
    git clone https://github.com/tuo-utente/nome-repo.git
    cd nome-repo
    ```

2. Follow the scripts in the order described above.  
   Please refer to comments within each script for details on requirements and parameters.

**Note:**  
For Cellranger usage, refer to the [10x Genomics documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest).

---

## Issues & Contact

If you encounter any problems or have questions, please open an issue in this repository or contact Dr. Mehdi Benamar at [Mehdi.Benamar@childrens.harvard.edu](mailto:Mehdi.Benamar@childrens.harvard.edu).

---

