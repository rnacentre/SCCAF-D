SCCAF-D
=====

**SCCAF-D** is novel approach in single-cell data analysis, particularly in **cell type deconvolution**. This method leverages the integration of multiple single-cell datasets and employs sophisticated machine learning techniques to generate an optimised reference. By doing so, SCCAF-D ensures the production of reliable and accurate deconvolution results.

In SCCAF-D, the metadata associated with single-cell data must contain three essential columns:

**cellID**: This column contains unique identifiers for individual cells in the dataset. Each cell is assigned a specific ID that distinguishes it from others in the dataset.

**cellType**: The cellType column provides information about the type or identity of each cell. It assigns a specific label or annotation to indicate the cell type it represents. These cell type annotations are crucial for downstream analysis and interpretation.

**sampleID**: This column specifies the sample. It helps organize the data and enables researchers to analyze cell type composition across different biological samples.

By including these three key metadata columns in the single-cell datasets, SCCAF-D can effectively perform cell type deconvolution and provide valuable insights into the cellular composition and heterogeneity within the analyzed samples.

**The concrete parameters include a list of 11 parameters**:

**[1] Bulk name**: The name of bulk data used.
**[2] Reference dataset name**: The name of reference data used.
**[3] Transformation**: The four methods available for transforming both the bulk and reference data are none (default), log, sqrt, and vst.
**[4] Deconvolution type**: Deconvolution methods are categorised into bulk and single-cell (sc) approaches based on the reference source. Bulk methods, such as CIBERSORT and FARDEEP, use a reference signature from sorted cell types or a marker gene list, while sc methods, like MuSiC and DWLS, rely on single-cell datasets.
**[5] Normalization for C, normalization**: C denotes single-cell reference data. 18 normalisation methods are supported, including column, row, mean, column z-score, global z-score, column min-max, global min-max, LogNormalize, none, QN, TMM, UQ, median ratios, TPM, SCTransform, scran, scater, and Linnorm.
**[6] Normalization for T, marker strategy**: T refers to bulk data, with normalization methods for bulk data as outlined in [5]. If using the sc method for deconvolution, this selection should be the same as in reference data.
**[7] Deconvolution method**: Twenty-five deconvolution algorithms are available for selection, including DWLS, FARDEEP, MuSiC, nnls, RLR, EpiDISH, OLS, EPIC, elasticNet, lasso, proportionsInAdmixture, ridge, CIBERSORT, SCDC, BisqueRNA, CDSeq, CPM, DCQ, DSA, DeconRNASeq, TIMER, deconf, dtangle, ssFrobenius, and ssKL.

**[8] Number of cells used**: The number of cells selected during the preparation of simulated 'pseudobulk' samples from single-cell transcriptomic data.

**[9] Remove cell type or not**: Whether to remove any cell types from the reference data.

**[10] Number of cores used**: Select the number of cores to utilize for conducting deconvolution.

**[11] Normalization first (T) or Transformation first (F)**: Whether to perform data normalization or transformation first.


----

**To use SCCAF-D ensure the SCCAF package is installed:**
```
conda install sccaf
```
Or

```
pip install sccaf
```

**You will also need to install reticulate package:**
```
install.packages('reticulate')
```
Example
Usage within R environment
--
```
#Set a working directory and place the files used for analysis in this directory
setwd('***')

#load SCCAF_D.R function
source('./SCCAF_D.R')

#selection the parameters to conduct deconvolution
param=c("bulk.rds",'single-reference.rds',"none","sc","TMM","TMM",'DWLS',10000,"none",1,'T')

#python_home:Specify the python to use
results <- SCCAF_D(param,python_home='/home/feng_shuo/miniconda3/envs/sccaf/bin/python')
```

How to cite
-----
For citation please refer to this:

Feng *et al*. Alleviating batch effects in cell type deconvolution.

The current version will be soon available at research square.



