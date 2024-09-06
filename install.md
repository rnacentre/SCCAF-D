This is the installation instruction for all the deconvolution methods used in the benchmarking.

# Notice
We suggest to install all the programs in a [conda](https://anaconda.org/) environment. 
The R programs can be classifed into THREE categories: from source codes, from [conda](https://anaconda.org/), from [bioconductor](bioconductor.org) and from [github](https://github.com/).

# Source code packages
```
# DSA
wget https://github.com/zhandong/DSA/raw/master/Package/version_1.0/DSA_1.0.tar.gz
R CMD INSTALL DSA_1.0.tar.gz
```

# conda dependencies
```
# glmnet
conda install -c r r-glmnet

# deconrnaseq
conda install -c bioconda bioconductor-deconrnaseq

# ADAPTS dependent
conda install -c conda-forge r-foreign
conda install -c conda-forge r-rmisc
conda install -c conda-forge libcurl

# CAMTHC  # this package is not used in the pipeline, instead debCAM is used
conda install -c r r-rjava
conda install -c conda-forge r-pcapp
conda install -c conda-forge r-dmwr

# CDSeq
conda install -c conda-forge r-slam

# Seurat
conda install -c conda-forge r-seuratdisk
conda install -c conda-forge r-seurat






```

# Bioconductor packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")
BiocManager::install("GO.db")
BiocManager::install("Rhtslib")

packages <- c("devtools", "BiocManager","data.table","ggplot2","tidyverse",
			  "Matrix","matrixStats",
			  "gtools",'stringi','future','rlang','multcomp','Rtsne','irlba','ROCR',
			  "foreach","doMC","doSNOW", #for parallelism
			  "Seurat","sctransform", #sc-specific normalization
			  "nnls","FARDEEP","MASS","glmnet","ComICS","dtangle") #bulk deconvolution methods

for (i in packages){ install.packages(i, character.only = TRUE)}

packages3 = c('debCAM','limma','edgeR','DESeq2','pcaMethods','BiocParallel','preprocessCore','scater','SingleCellExperiment','Linnorm','DeconRNASeq','multtest','GSEABase','annotate','genefilter','preprocessCore','graph','MAST','Biobase') #last two are required by DWLS and MuSiC, respectively.
for (i in packages3){ BiocManager::install(i, character.only = TRUE)}

# Dependencies for CellMix: 'NMF', 'csSAM', 'GSEABase', 'annotate', 'genefilter', 'preprocessCore', 'limSolve', 'corpcor', 'graph', 'BiocInstaller'
packages2 = c('NMF','csSAM','limSolve','corpcor','pkgmaker', 'NMF', 'csSAM', 'registry', 'stringr', 'GSEABase', 'Biobase', 'BiocGenerics', 'annotate', 'matrixStats', 'genefilter', 'AnnotationDbi', 'RSQLite', 'DBI', 'preprocessCore', 'limSolve', 'quadprog', 'corpcor', 'xtable', 'gtools', 'beeswarm', 'graph', 'BiocInstaller', 'bibtex', 'digest', 'ggplot2', 'plyr')
for (i in packages2){ install.packages(i, character.only = TRUE)}


install.packages('remotes')
```

# GitHub packages


```
devtools::install_github("LTLA/BiocNeighbors")  
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE) #requires knitr
devtools::install_github("xuranw/MuSiC") 
devtools::install_github("meichendong/SCDC")
devtools::install_github("rosedu1/deconvSeq")
devtools::install_github("cozygene/bisque")
devtools::install_github("dviraran/SingleR@v1.0")
devtools::install_github("jingshuw/descend")
install.packages("scBio")
devtools::install_github("amitfrish/scBio")
BiocManager::install("EpiDISH")
devtools::install_bitbucket("yuanlab/dwls")  # this package is not used in the pipeline, instead the /DWLS.R is used
devtools::install_github("Danko-Lab/TED/TED")
devtools::install_github("chichaumiau/CellMix")
devtools::install_github("Lululuella/CAMTHC")  # this package is not used in the pipeline
devtools::install_github("kkang7/CDSeq_R_Package")
devtools::install_github('shenorrlab/bseqsc')

```
# More
the functions for the methods below can be found in the irrespective files in this repository: 
TIMER : /TIMER.R  more info about the package at http://cistrome.org/TIMER/
DWLS : /DWLS.R   more info about the package at https://github.com/dtsoucas/DWLS

CIBERSORT 
The source code for CIBERSORT can be asked from the authors at https://cibersort.stanford.edu/download.php
#bseqsc method is also dependent on CIBERSORT



# SessionInfo()
```
R version 3.6.1 (2019-07-05)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.4 (Maipo)
Matrix products: default
BLAS/LAPACK: /nfs/leia/research/ma/chichau/bin/anaconda3/lib/R/lib/libRblas.so
locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
attached base packages:
 [1] grid      stats4    compiler  parallel  stats     graphics  grDevices
 [8] utils     datasets  methods   base     
other attached packages:
 [1] BisqueRNA_1.0.5      ADAPTS_1.0.6         scBio_0.1.6         
 [4] descend_1.0.0        EpiDISH_2.2.2        TED_1.1             
 [7] SCDC_0.0.0.9000      MuSiC_0.1.1          deconvSeq_0.3.0     
[10] DeconRNASeq_1.28.0   ggplot2_3.3.3        pcaMethods_1.78.0   
[13] limSolve_1.5.6       dtangle_2.0.9        MASS_7.3-53.1       
[16] FARDEEP_1.0.1        EPIC_1.1.5           nnls_1.4            
[19] CellMix_1.6.2        GSEABase_1.48.0      graph_1.64.0        
[22] annotate_1.64.0      XML_3.99-0.3         AnnotationDbi_1.48.0
[25] IRanges_2.20.2       S4Vectors_0.24.4     stringr_1.4.0       
[28] csSAM_1.2.4          NMF_0.23.0           Biobase_2.46.0      
[31] BiocGenerics_0.32.0  cluster_2.1.1        rngtools_1.5        
[34] pkgmaker_0.32.2      registry_0.5-1       glmnet_4.1-1        
[37] Matrix_1.3-2         Seurat_3.1.2        
loaded via a namespace (and not attached):
  [1] rsvd_1.0.3                  Hmisc_4.5-0                
  [3] ica_1.0-2                   corpcor_1.6.9              
  [5] class_7.3-18                Rsamtools_2.2.3            
  [7] foreach_1.5.1               lmtest_0.9-38              
  [9] crayon_1.4.1                rbibutils_2.0              
 [11] nlme_3.1-152                backports_1.2.1            
 [13] rlang_0.4.10                XVector_0.26.0             
 [15] ROCR_1.0-11                 irlba_2.3.3                
 [17] SparseM_1.81                nloptr_1.2.2.2             
 [19] limma_3.42.2                scater_1.14.6              
 [21] L1pack_0.38.196             BiocParallel_1.20.1        
 [23] bit64_4.0.5                 glue_1.4.2                 
 [25] pheatmap_1.0.12             sctransform_0.3.2          
 [27] vipor_0.4.5                 mcmc_0.9-7                 
 [29] tidyselect_1.1.0            SummarizedExperiment_1.16.1
 [31] fitdistrplus_1.1-3          tidyr_1.1.3                
 [33] zoo_1.8-9                   GenomicAlignments_1.22.1   
 [35] xtable_1.8-4                MatrixModels_0.5-0         
 [37] magrittr_2.0.1              evaluate_0.14              
 [39] bibtex_0.4.2.3              Rdpack_2.1.1               
 [41] zlibbioc_1.32.0             sn_1.6-2                   
 [43] rstudioapi_0.13             fastmatrix_0.3-81          
 [45] rpart_4.1-15                mathjaxr_1.4-0             
 [47] locfdr_1.1-8                BiocSingular_1.2.2         
 [49] xfun_0.22                   askpass_1.1                
 [51] multtest_2.42.0             caTools_1.18.1             
 [53] pbdZMQ_0.3-5                methylKit_1.12.0           
 [55] tibble_3.1.0                lpSolve_5.6.15             
 [57] quantreg_5.85               ggrepel_0.9.1              
 [59] ape_5.4-1                   listenv_0.8.0              
 [61] Biostrings_2.54.0           png_0.1-7                  
 [63] future_1.21.0               withr_2.4.1                
 [65] bitops_1.0-6                plyr_1.8.6                 
 [67] e1071_1.7-4                 dqrng_0.2.1                
 [69] coda_0.19-4                 pillar_1.5.1               
 [71] gplots_3.1.1                cachem_1.0.4               
 [73] multcomp_1.4-16             DelayedMatrixStats_1.8.0   
 [75] vctrs_0.3.6                 ellipsis_0.3.1             
 [77] generics_0.1.0              metap_1.4                  
 [79] tools_3.6.1                 foreign_0.8-71             
 [81] beeswarm_0.3.1              munsell_0.5.0              
 [83] DelayedArray_0.12.3         fastmap_1.1.0              
 [85] rtracklayer_1.46.0          plotly_4.9.3               
 [87] GenomeInfoDbData_1.2.2      gridExtra_2.3              
 [89] MCMCpack_1.5-0              edgeR_3.28.1               
 [91] lattice_0.20-41             mutoss_0.1-12              
 [93] utf8_1.2.1                  dplyr_1.0.5                
 [95] BiocFileCache_1.10.2        jsonlite_1.7.2             
 [97] scales_1.1.1                pbapply_1.4-3              
 [99] genefilter_1.68.0           lazyeval_0.2.2             
[101] doParallel_1.0.16           latticeExtra_0.6-29        
[103] R.utils_2.10.1              reticulate_1.18            
[105] checkmate_2.0.0             sandwich_3.0-0             
[107] cowplot_1.1.1               statmod_1.4.35             
[109] Rtsne_0.15                  uwot_0.1.10                
[111] igraph_1.2.6                survival_3.2-7             
[113] numDeriv_2016.8-1.1         plotrix_3.8-1              
[115] htmltools_0.5.1.1           memoise_2.0.0              
[117] locfit_1.5-9.4              quadprog_1.5-8             
[119] viridisLite_0.3.0           digest_0.6.27              
[121] assertthat_0.2.1            rappdirs_0.3.3             
[123] repr_1.1.3                  emdbook_1.3.12             
[125] RSQLite_2.2.4               future.apply_1.7.0         
[127] Rmisc_1.5                   data.table_1.14.0          
[129] blob_1.2.1                  R.oo_1.24.0                
[131] preprocessCore_1.48.0       splines_3.6.1              
[133] Formula_1.2-4               RCurl_1.98-1.2             
[135] hms_1.0.0                   colorspace_2.0-0           
[137] base64enc_0.1-3             BiocManager_1.30.10        
[139] mnormt_2.0.2                ggbeeswarm_0.6.0           
[141] SDMTools_1.1-221            GenomicRanges_1.38.0       
[143] shape_1.4.5                 tmvnsim_1.0-2              
[145] nnet_7.3-15                 Rcpp_1.0.6                 
[147] mclust_5.4.7                RANN_2.6.1                 
[149] mvtnorm_1.1-1               fansi_0.4.2                
[151] conquer_1.0.2               truncnorm_1.0-8            
[153] parallelly_1.23.0           IRdisplay_1.0              
[155] R6_2.5.0                    ggridges_0.5.3             
[157] lifecycle_1.0.0             TFisher_0.2.0              
[159] curl_4.3                    minqa_1.2.4                
[161] xbioc_0.1.19                leiden_0.3.7               
[163] qvalue_2.18.0               RcppAnnoy_0.0.18           
[165] TH.data_1.0-10              RColorBrewer_1.1-2         
[167] iterators_1.0.13            fastseg_1.32.0             
[169] htmlwidgets_1.5.3           biomaRt_2.42.1             
[171] purrr_0.3.4                 globals_0.14.0             
[173] openssl_1.4.3               htmlTable_2.1.0            
[175] bdsmatrix_1.3-4             codetools_0.2-18           
[177] matrixStats_0.58.0          gtools_3.8.2               
[179] prettyunits_1.1.1           psych_2.0.12               
[181] SingleCellExperiment_1.8.0  dbplyr_2.1.0               
[183] gridBase_0.4-7              R.methodsS3_1.8.1          
[185] GenomeInfoDb_1.22.1         gtable_0.3.0               
[187] tsne_0.1-3                  DBI_1.1.1                  
[189] httr_1.4.2                  KernSmooth_2.23-18         
[191] stringi_1.5.3               progress_1.2.2             
[193] reshape2_1.4.4              uuid_0.1-4                 
[195] viridis_0.5.1               bbmle_1.0.23.1             
[197] boot_1.3-27                 IRkernel_1.1.1             
[199] BiocNeighbors_1.9.4         lme4_1.1-26                
[201] geneplotter_1.64.0          DESeq2_1.26.0              
[203] scran_1.14.6                bit_4.0.4                  
[205] jpeg_0.1-8.1                pkgconfig_2.0.3            
[207] Rsolnp_1.16                 knitr_1.31    
```

# footnote
```
ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libcrypto.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib

ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libcurl.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib

ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libz.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib

ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libssl.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib

ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libgfortran.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib

ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libgomp.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib/

ln -s /nfs/leia/research/ma/chichau/bin/anaconda3/lib/libquadmath.so \
/nfs/leia/research/ma/chichau/bin/anaconda3/x86_64-conda_cos6-linux-gnu/lib/

cd /nfs/leia/research/ma/chichau/bin/anaconda3/lib/R/etc
cp ../../../pkgs/r-base-3.6.1-h9bb98a2_1/lib/R/etc/Makeconf .

```
