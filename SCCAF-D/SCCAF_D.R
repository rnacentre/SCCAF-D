#' @description
#' Use a single-cell references to deconvolve a bulk data. 
#' @details
#' 
#' @param param: a list of 11 parameters:
#' 1: bulk name: The bulk dataset used to conduct deconvolution.
#' 2: reference dataset name: The single-cell dataset used to conduct deconvolution.
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
#' 
#' @return
#' The deconvolution results as a table.
#' @example
#' @batch_key scanpy: scanpy.pp.highly_variable_genes
#' @span: scanpy: scanpy.pp.highly_variable_genes
#' @ python_home: Specify the python to use
SCCAF_D = function(param, batch_key='sampleID',span=0.3,python_home = Sys.which("python")) {
    reticulate::use_python(python_home)
    source("./benchmark1.R")
    source("./deconvolution1.R")
    source("./Frame.R")
    source("./DWLS.R")
    #####
    if (!reticulate::py_module_available("SCCAF")) {
        stop("python module SCCAF does not seem to be installed; - try running 'pip install SCCAF'")
    }
    reticulate::source_python("./scanpy_workflow.py")
    #####
    if (length(param) != 11) {
        print("Please check that all required parameters are indicated or are correct")
        print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
        print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
        stop()
    }
    flag = FALSE
    bulk = param[1]
    dataset = param[2]
    transformation = param[3]
    deconv_type = param[4]
    if (deconv_type == "bulk") {
        normalization = param[5]
        marker_strategy = param[6]
    }
    else if (deconv_type == "sc") {
        normalization_scC = param[5]
        normalization_scT = param[6]
    }
    else {
        print("Please enter a valid deconvolution framework")
        stop()
    }
    method = param[7]
    number_cells = round(as.numeric(param[8]), digits = -2)
    to_remove = param[9]
    num_cores = min(as.numeric(param[10]), parallel::detectCores() - 
        1)
    if (param[11] == "T") {
        NormTrans = TRUE
    }
    else {
        NormTrans = FALSE
    }
    ###
    pathway=getwd()
    ####
    dataset1=strsplit(dataset,'\\.')[[1]][1]
    ####
    X2_1 = readRDS(dataset)
    sceasy::convertFormat(X2_1, from = "seurat", to = "anndata",main_layer = "counts",outFile = paste(dataset1,'.h5',sep=''))
    ###
    ad=paste(pathway,'/',dataset1,'.h5',sep='')
    #####selection cells--scanpy
    scanpy_workflow(ad,pathway,batch_key=batch_key,span=span)
    ###selection top100 cells for each cell type
    reference = selection_cells(pathway,dataset1,dataset)
    ###
    print(table(reference$cellType))
    ####
    ####read data
    X1 = read_bulk(bulk)
    X2 = read_data(paste(pathway,'/',dataset1,'_',"sccaf-reference",".rds",sep=""))
    ####
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
    return(Framework(deconv_type, NormTrans, Xtest, Xtrain, normalization, 
        normalization_scT, normalization_scC, transformation, 
        marker_strategy, to_remove, Xtest$P, method, pDataC))
}
