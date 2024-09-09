#' @description
#' Use a single-cell references to deconvolve a bulk data. 
#' @details
#' 
#' @param param: a list of 11 parameters:
#' 1: bulk name: The bulk dataset used to conduct deconvolution.
#' 2: reference name: The single-cell dataset used to conduct deconvolution.
#' 3: Transformation: none (defalt), log, sqrt, vst. The choice whether to transform the dataset both bulk and single-cell dataset.
#' 4: deconvolution type: bulk, sc. The different deconvolution types. We divide the methods into bulk and sc methods refer to the reference source, 
#' like bulk methods includes the CIBERSORT, FARDEEP (reference is either a signature of sorted cell types or a marker gene list), 
#' or the sc methods includes MuSiC, DWLS using the single-cell dataset.

#' 5: Normalization for C, normalization. Normalize the single-cell dataset.
#' 6: Normalization for T, marker strategy. Normalize the bulk dataset, if deconvolution type choose the bulk, this parameter should set 'all',
#' if choose the 'sc',the parameter should set the same normalization as the bulk dataset.
#' 7: Deconvolution method. The different deconvolution algorithms.
#' 8: number of cells used. How many cells used to produce the pseudobulk. this parameter is used in pseudobulk process.
#' 9: remove cell type or not (none: default). Remove a cell type in bulk and singl-cell dataset before deconvolution.
#' 10: number of cores used. 
#' 11: Normalize first (T) or Transform first(F).
#' 12: Return_expr: logical. Whether to return the estimated expression matrix of bulk data. Default (FALSE).
#' 13: Batch_key: string. The parameter used to calculate the highly variable genes in SCANPY.
#' 14: Span: numeric. The parameter used to calculate the highly variable genes in SCANPY. The fraction of the data (cells) used when estimating the variance in the loess model fit if flavor = 'seurat_v3'.
#' 15: Python_home: string. The path where python is located.

SCCAF_D=function(bulk=bulk,reference=reference,transformation='none',deconv_type='sc',normalization_C='TMM',normalization_T='TMM',marker_strategy='all',method='DWLS',number_cells=10000,to_remove='none',num_cores=1,NormTrans=TRUE,return_expr=FALSE,batch_key='sampleID',span=0.3,python_home = Sys.which("python")) {
    reticulate::use_python(python_home)
    source("benchmark1.R")
    source("deconvolution1.R")
    source("Frame.R")
    source("DWLS.R")
    source("expr.R")
    #####
    if (!reticulate::py_module_available("SCCAF")) {
        stop("python module SCCAF does not seem to be installed; - try running 'pip install SCCAF'")
    }
    reticulate::source_python("/home/feng_shuo/deconvolution/SCCAF-D/scanpy_workflow.py")
    #####
    pathway=getwd()
    ####
    reference1=strsplit(reference,'\\.')[[1]][1]
    ####
    X2_1 = readRDS(reference)
    sceasy::convertFormat(X2_1, from = "seurat", to = "anndata",main_layer = "counts",outFile = paste(reference1,'.h5',sep=''))
    ###
    ad=paste(pathway,'/',reference1,'.h5',sep='')
    #####selection cells--scanpy
    scanpy_workflow(ad,pathway,batch_key=batch_key,span=span)
    ###selection top100 cells for each cell type
    reference = selection_cells(pathway,reference1,reference)
    ###
    print(table(reference$cellType))
    ####
    ####read data
    X1 = read_bulk(bulk)
    X2 = read_data(paste(pathway,'/',reference1,'_',"sccaf-reference",".rds",sep=""))
    ####
    flag = FALSE
    if (FALSE) {
        X2 <- QC(X2)
    }
    to_keep = intersect(rownames(X1$data), rownames(X2$data))
    print(paste0("number of intersect features: ", length(to_keep)))
    ####
    X1$data = X1$data[to_keep, ]
    X2$data = X2$data[to_keep, ]
    Xtrain = prepare_train(X2$data, X2$original_cell_names)
    pDataC = X2$pData
    P <- as.data.frame(matrix(0.1, nrow = length(unique(colnames(X2$data))), 
        ncol = dim(X1$data)[2]))
    rownames(P) <- unique(colnames(X2$data))
    colnames(P) <- colnames(X1$data)
    Xtest <- list(T = X1$data, P = P)
    #####
    cell_proportion=Framework(deconv_type, NormTrans, Xtest, Xtrain, normalization, 
        normalization_T, normalization_C, transformation, 
        marker_strategy, to_remove, Xtest$P, method, pDataC)
    ###expression matrix
    ref=X2$data
    ref = t(as.matrix(ref))
    subBulk=X1$data
    subBulk = t(subBulk)
    reshaped_data <- dcast(cell_proportion, formula = tissue ~ CT, value.var = "observed_values")
    # set CT rowname
    rownames(reshaped_data) <- reshaped_data$tissue
    reshaped_data$tissue <- NULL  #remove CT
    expr=deconvExp(reference = ref, cell.state.labels =pDataC$cellType, bluk = subBulk, cell.state.fraction = reshaped_data)
    
    if (!return_expr){
        return(cell_proportion)
    }else {
        return(list('cell_pro'=cell_proportion,'expr'=expr))
    }
    
}
