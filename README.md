SCCAF-D
=====

**SCCAF-D** is novel approach in single-cell data analysis, particularly in **cell type deconvolution**. This method leverages the integration of multiple single-cell datasets and employs sophisticated machine learning techniques to generate an optimised reference. By doing so, SCCAF-D ensures the production of reliable and accurate deconvolution results.

In SCCAF-D, the metadata associated with single-cell data must contain three essential columns:

**cellID**: This column contains unique identifiers for individual cells in the dataset. Each cell is assigned a specific ID that distinguishes it from others in the dataset.

**cellType**: The cellType column provides information about the type or identity of each cell. It assigns a specific label or annotation to indicate the cell type it represents. These cell type annotations are crucial for downstream analysis and interpretation.

**sampleID**: This column specifies the sample. It helps organize the data and enables researchers to analyze cell type composition across different biological samples.

By including these three key metadata columns in the single-cell datasets, SCCAF-D can effectively perform cell type deconvolution and provide valuable insights into the cellular composition and heterogeneity within the analyzed samples.

**The concrete parameters include a list of 11 parameters**:

1: **Bulk name**: The bulk dataset used to conduct deconvolution.

2: **Reference dataset name**: The single-cell dataset used to conduct deconvolution.

3: **Transformation**: none (defalt), log, sqrt, vst. The choice whether to transform the dataset both bulk and single-cell dataset.

4: **Deconvolution type**: bulk, sc. The different deconvolution types. We divide the methods into bulk and sc methods refer to the reference source, 
like bulk methods includes the CIBERSORT, FARDEEP (reference is either a signature of sorted cell types or a marker gene list), 
or the sc methods includes MuSiC, DWLS using the single-cell dataset.

5: **Normalization for C, normalization**. Normalize the single-cell dataset.

6: **Normalization for T, marker strategy**. Normalize the bulk dataset, if deconvolution type choose the bulk, this parameter should set 'all',
if choose the 'sc',the parameter should set the same normalization as the bulk dataset.

7: **Deconvolution method**. The different deconvolution algorithms.

8: **Number of cells used**. How many cells used to produce the pseudobulk. this parameter is used in pseudobulk process.

9: **Remove cell type or not (none: default)**. Remove a cell type in bulk and singl-cell dataset before deconvolution.

10: **Number of cores used**. 

11: **Normalize first (T) or Transform first(F)**.



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



