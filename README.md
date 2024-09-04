SCCAF-D
=====

**SCCAF-D** is novel approach in single-cell data analysis, particularly in **cell type deconvolution**. This method leverages the integration of multiple single-cell datasets and employs sophisticated machine learning techniques to generate an optimised reference. By doing so, SCCAF-D ensures the production of reliable and accurate deconvolution results.

In SCCAF-D, the metadata associated with single-cell data must contain three essential columns:

**cellID**: This column contains unique identifiers for individual cells in the dataset. Each cell is assigned a specific ID that distinguishes it from others in the dataset.

**cellType**: The cellType column provides information about the type or identity of each cell. It assigns a specific label or annotation to indicate the cell type it represents. These cell type annotations are crucial for downstream analysis and interpretation.

**sampleID**: This column specifies the sample. It helps organize the data and enables researchers to analyze cell type composition across different biological samples.

By including these three key metadata columns in the single-cell datasets, SCCAF-D can effectively perform cell type deconvolution and provide valuable insights into the cellular composition and heterogeneity within the analyzed samples.

**The concrete parameters include a list of 11 parameters**:

**Param[1]**: string. The name of bulk data used to be deconvolved.
**Param[2]**: string. The name of reference data used for deconvolution.
**Param[3]**: string. The four methods available for transforming both the bulk and reference data are none (default), log, sqrt, and vst.
**Param[4]**: string. Deconvolution methods are categorised into bulk and single-cell (sc) approaches based on the reference source. Bulk methods, such as CIBERSORT and FARDEEP, use a reference signature from sorted cell types or a marker gene list, while ‘sc’ methods, like MuSiC and DWLS, use single-cell datasets as reference data.
**Param[5]**: string. Normalization for reference data. 18 normalisation methods are supported, including column, row, mean, column z-score, global z-score, column min-max, global min-max, LogNormalize, none, QN, TMM, UQ, median ratios, TPM, SCTransform, scran, scater, and Linnorm.
**Param[6]**: string. Normalization for bulk data. With normalization methods for bulk data as outlined in [5]. If using the sc method for deconvolution, this selection should be the same as in reference data. If using the bulk methods, this selection should be ‘all’.
**Param[7]**: string. Twenty-five deconvolution algorithms are available for selection, including DWLS, FARDEEP, MuSiC, nnls, RLR, EpiDISH, OLS, EPIC, elasticNet, lasso, proportionsInAdmixture, ridge, CIBERSORT, SCDC, BisqueRNA, CDSeq, CPM, DCQ, DSA, DeconRNASeq, TIMER, deconf, dtangle, ssFrobenius, and ssKL.

**Param[8]**: string. The number of cells selected during the preparation of simulated 'pseudobulk' from single-cell data.

**Param[9]**: string. Whether to remove any cell types from the reference data.

**Param[10]**: string. Select the number of cores to be used for deconvolution.
**Param[11]**: string. Whether to perform data normalization or transformation first. T means the normalization first. F means the transformation first.



----

**To use SCCAF-D ensure the SCCAF package is installed:**
```
conda install sccaf
```
Or

```
pip install sccaf
```

**You will also need to install another packages:**
```
packages <- c("devtools", "BiocManager","data.table","ggplot2","tidyverse","reticulate","pheatmap",
			  "Matrix","matrixStats",
			  "gtools",
			  "foreach","doMC","doSNOW", 
			  "Seurat","sctransform", 
			  "nnls","MASS","glmnet") 

for (i in packages){ install.packages(i, character.only = TRUE)}

packages1 = c('limma','edgeR','DESeq2','pcaMethods','BiocParallel','preprocessCore','scater','SingleCellExperiment','Linnorm','DeconRNASeq','multtest','GSEABase','annotate','genefilter','preprocessCore','graph','MAST','Biobase','sparseMatrixStats')
for (i in packages1){ BiocManager::install(i, character.only = TRUE)}

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



