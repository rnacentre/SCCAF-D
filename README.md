SCCAF-D
=====

**SCCAF-D** is novel approach in single-cell data analysis, particularly in **cell type deconvolution**. This method leverages the integration of multiple single-cell datasets and employs sophisticated machine learning techniques to generate an optimised reference. By doing so, SCCAF-D ensures the production of reliable and accurate deconvolution results.

In SCCAF-D, the metadata associated with single-cell data must contain three essential columns:

**cellID**: This column contains unique identifiers for individual cells in the dataset. Each cell is assigned a specific ID that distinguishes it from others in the dataset.

**cellType**: The cellType column provides information about the type or identity of each cell. It assigns a specific label or annotation to indicate the cell type it represents. These cell type annotations are crucial for downstream analysis and interpretation.

**sampleID**: This column specifies the sample. It helps organize the data and enables researchers to analyze cell type composition across different biological samples.

By including these three key metadata columns in the single-cell datasets, SCCAF-D can effectively perform cell type deconvolution and provide valuable insights into the cellular composition and heterogeneity within the analyzed samples.


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

Feng,S *et al*. Alleviating batch effects in cell type deconvolution.

The current version is not acceptable.



