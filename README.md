<a href="https://doi.org/10.5281/zenodo.14211887"><img src="https://zenodo.org/badge/795955983.svg" alt="DOI"></a>
# SCCAF-D: Single-Cell Clustering Assessment Framework optimised reference for Deconvolution
**SCCAF-D** is a novel framework for **cell type deconvolution** using single-cell RNA-seq data as reference. By integrating multiple single-cell datasets and employing advanced machine learning techniques, SCCAF-D produces optimized reference data to ensure reliable and accurate deconvolution results. This approach provides valuable insights into the cell composition and heterogeneity of biological samples.

## Metadata Requirements

To perform deconvolution, the metadata of the single-cell datasets should contain the following essential columns:

- **cellID**: Unique identifiers for individual cells in the dataset.
- **cellType**: Annotations or labels indicating the type of each cell.
- **sampleID**: Identifiers for biological samples to facilitate comparative analysis across different conditions.

These three columns are critical for accurate cell type deconvolution using SCCAF-D.

## Parameters

SCCAF-D requires the following parameters to perform deconvolution:
   - **Bulk**: The bulk RNA-seq data used to be deconvolved.
   - **Reference**: The reference data used for deconvolution.
   - **Transformation** (string): The method for transforming both the bulk and reference data. Options include `'none'` (default), `'log'`, `'sqrt'`, and `'vst'`.
   - **Deconv_type** (string): Deconvolution methods are categorized into bulk (`'bulk'`) and single-cell (`'sc'`) approaches based on the reference source. Bulk methods, such as `CIBERSORT` and `FARDEEP`, use a reference signature from sorted cell types or a marker gene list. In contrast, single-cell methods, such as `MuSiC` and `DWLS`, use single-cell datasets as reference data.
   - **Normalization_C** (string): Normalization for reference data. Eighteen normalization methods are supported, including `'column'`, `'row'`, `'mean'`, `'column z-score'`, `'global z-score'`, `'column min-max'`, `'global min-max'`, `'LogNormalize'`, `'none'`, `'QN'`, `'TMM'`, `'UQ'`, `'median ratios'`, `'TPM'`, `'SCTransform'`, `'scran'`, `'scater'`, and `'Linnorm'`.

   - **Normalization_T** (string): Normalization for bulk data. The same normalization methods as listed for reference data. **Suggestion**: If using the single-cell (`sc`) method for deconvolution, the choice for this parameter should match the reference data normalization.

   - **Method** (string): Twenty-five deconvolution algorithms are available, including `'DWLS'`, `'FARDEEP'`, `'MuSiC'`, `'nnls'`, `'RLR'`, `'EpiDISH'`, `'OLS'`, `'EPIC'`, `'elasticNet'`, `'lasso'`, `'proportionsInAdmixture'`, `'ridge'`, `'CIBERSORT'`, `'SCDC'`, `'BisqueRNA'`, `'CDSeq'`, `'CPM'`, `'DCQ'`, `'DSA'`, `'DeconRNASeq'`, `'TIMER'`, `'deconf'`, `'dtangle'`, `'ssFrobenius'`, and `'ssKL'`.
   - **Number_cells** (integer): The number of cells to select when preparing simulated 'pseudobulk' from single-cell data.
   - **To_remove** (string): Specify any cell types to exclude from the reference data.
   - **Num_cores** (integer): The number of cores to use for parallelization during deconvolution.
   - **NormTrans** (logical): Whether to perform data normalization or transformation first. `TRUE` indicates normalization first, while `FALSE` indicates transformation first.
   - **Return_expr** (logical): Whether to return the estimated expression matrix of the bulk data. Default is `FALSE`.
   - **Batch_key** (string): The parameter used to calculate highly variable genes in SCANPY.
   - **Span** (numeric): The fraction of cells used when estimating variance in the loess model fit in SCANPY (when `flavor = 'seurat_v3'`).
   - **Python_home** (string): The path to the Python executable.

----

## Installation

To install the required [SCCAF](https://github.com/SCCAF/sccaf) package, use one of the following commands:

```shell
conda install sccaf
```

or

```shell
pip install sccaf
```

The installation of SCCAF takes about a few minutes.

Additionally, several R and Bioconductor packages are needed:

```R
# CRAN packages
packages <- c("devtools", "BiocManager", "data.table", "ggplot2", "tidyverse", 
              "reticulate", "pheatmap", "Matrix", "matrixStats", "gtools",
              "foreach", "doMC", "doSNOW", "Seurat", "sctransform", "nnls", 
              "MASS", "glmnet","reshape2","quadprog","reshape","e1071","Seurat","ROCR",
              "varhandle")
install.packages(packages)

# Bioconductor packages
bioc_packages <- c('limma', 'edgeR', 'DESeq2', 'pcaMethods', 'BiocParallel', 
                   'preprocessCore', 'scater', 'SingleCellExperiment', 'Linnorm',
                   'DeconRNASeq', 'multtest', 'GSEABase', 'annotate', 'genefilter', 
                   'graph', 'MAST', 'Biobase', 'sparseMatrixStats')
BiocManager::install(bioc_packages)

# github packages
devtools::install_github("cellgeni/sceasy")
```

**Note**: Installing the required packages may take some time, depending on the user's R environment and network conditions. The installation process could take several minutes or even longer, potentially up to an hour.

**Note**: if user wants to use other 24 deconvolution methods (except for DWLS) in SCCAF-D framework, please install packages mentioned in <a href="https://github.com/rnacentre/SCCAF-D/blob/main/install.md"> install.md </a> file. 
## Example Usage in R

```R
# Set working directory to the location of your data
setwd('/data/')

# Load the SCCAF-D R function
source('/SCCAF-D/SCCAF_D.R')

# Specify the path to the Python environment
python_home <- '/home/miniconda3/envs/sccaf/bin/python'

# Run SCCAF-D deconvolution
results <- SCCAF_D(bulk = 'pseudobulk_Baron_T.rds',reference = 'integrated_baron.rds',python_home='/home/miniconda3/envs/sccaf/bin/python')
```
We provide a detailed example of **SCCAF-D** usage in a Jupyter Notebook, which can be found [here](https://github.com/rnacentre/SCCAF-D/blob/main/SCCAF-D%20example.ipynb). This example demonstrates the deconvolution process using integrated single-cell data, including **three** single-cell datasets, consisting of a total of **32523** cells, as the reference, and **five** bulk RNA-seq samples for deconvolution. Example data is available <a href="https://drive.google.com/drive/folders/1nMCtnaRN-5c-Tv5KxXl4QbXqOW1faxVz?usp=drive_link"> here</a>. Processing this demo will take approximately **30** minutes on a server with 32 CPUs and 256 GB of RAM. For larger datasets, increased memory may be required to ensure efficient computation.

## Scalability

To assess the scalability of the SCCAF-D method, we analyzed 20 distinct datasets, each comprising 1000 synthetically simulated bulk samples. The scalability of each method was evaluated by quantifying both the running time in seconds and the memory usage in gigabytes. We calculated the maximum memory usage and the total running time across all datasets, with detailed data provided in Table 1(/data/Table 1). Notably, the SCCAF-D method showed a maximum memory usage of 33 GB, with reference data derived from the integrated dataset of Dominguez et al., Madissoon et al., and He et al., which comprised 156,334 cells, and bulk data from the Tabula dataset. We observed the longest total running time for SCCAF-D on the dataset from Habermann et al., which contained 81,372 cells, reaching a maximum of 37603.691 seconds (10.45 hours)(**Fig. 1**). This evaluation was conducted on a normal workstation with a setup comprising 376 GB of RAM and 20 CPUs.
![memory](https://github.com/rnacentre/SCCAF-D/blob/main/data/github_time_memory.png)
**Fig. 1** Memory usage and runtime of selected deconvolution methods. (**a**) Maximum memory usage of different deconvolution methods. Each point represents a dataset, with colours indicating the deconvolution method. (**b**) The runtime of different simulated bulk dataset using SCCAF-D. The X-axis denotes the names of the simulated bulk datasets.



## Citation
To cite SCCAF-D, please refer to the following:

> Shuo Feng, Liangfeng Huang, Anna Vathrakokili Pournara, Ziliang Huang, Xinlu Yang, Yongjian Zhang, Alvis Brazma, Ming Shi, Irene Papatheodorou, Zhichao Miao "Alleviating batch effects in cell type deconvolution." (<a href='https://www.nature.com/articles/s41467-024-55213-x'> Article link</a>)
