# This file includes functions in the benchmarking process. 
# TODO: "LinSeed"
# TODO: document the functions

require(Biobase)

#################################################
##########    DECONVOLUTION METHODS    ##########
# T is pseudo-bulk
Deconvolution <- function(T, C, method, phenoDataC, P = NULL, elem = NULL, STRING = NULL, marker_distrib, refProfiles.var){ 

    bulk_methods = c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elasticNet","lasso","ridge","EPIC",
					 "DSA","ssKL","ssFrobenius","dtangle", "deconf", "proportionsInAdmixture", "svmdecon", "EpiDISH","CAMmarker","CDSeq" )
    sc_methods = c("MuSiC","BisqueRNA","DWLS","deconvSeq","SCDC","bseqsc","CPM","TIMER")
    
    keep = intersect(rownames(C),rownames(T)) 
    C = C[keep,]
    T = T[keep,]

    ########## Using marker information for bulk_methods
    if(method %in% bulk_methods){
#         C = C[rownames(C) %in% marker_distrib$gene,]
#         T = T[rownames(T) %in% marker_distrib$gene,]
        refProfiles.var = refProfiles.var[rownames(refProfiles.var) %in% marker_distrib$gene,]

    } else { 
        sample_column = grep("sampleID",colnames(phenoDataC)) 
        colnames(phenoDataC)[sample_column] = "SubjectName"
        rownames(phenoDataC) = phenoDataC$cellID

        phenoDataC<-phenoDataC[colnames(C),]
        
        require(SingleCellExperiment) 
        T.eset=as.matrix(T)
        C.eset <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(C)),
                                    colData = phenoDataC)
    }
    
    fun<-get(sprintf('run_%s', method))
    RESULTS = fun(T = T, 
				  C = C, 
				  T.eset = T.eset, 
				  C.eset = C.eset, 
				  phenoDataC = phenoDataC, 
				  marker_distrib = marker_distrib,
				  refProfiles.var = refProfiles.var, 
				  STRING = STRING,
				  elem = elem)
#     saveRDS(RESULTS, "RESULTS0.rds")
    RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)),]                  
    RESULTS = data.table::melt(RESULTS)                   
	colnames(RESULTS) <-c("CT","tissue","observed_values")
#     print(RESULTS)
#     saveRDS(RESULTS, "RESULTS.rds")
#     saveRDS(P, "P.rds")
	if(!is.null(P)){
		P = P[gtools::mixedsort(rownames(P)),]
		P$CT = rownames(P)
		P = data.table::melt(P, id.vars="CT")
		colnames(P) <-c("CT","tissue","expected_values")
		RESULTS = merge(RESULTS,P)
		RESULTS$expected_values <-round(RESULTS$expected_values,3)
		RESULTS$observed_values <-round(RESULTS$observed_values,3)
	}
#     print(RESULTS)
    return(RESULTS) 

}

########################################################################
# bulk Methods
run_EpiDISH<-function(T, C, ...){
	#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default

	require(EpiDISH)
	RESULTS = t(EpiDISH::epidish(beta.m = T, ref.m = C, method = "RPC")$estF)
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}

run_FARDEEP<-function(T, C, ...){
	require(FARDEEP)
	RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}

run_RLR<-function(T, C, ...){
	# RLR = robust linear regression
	require(MASS)
	RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
	rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}

run_nnls<-function(T, C, ...){
	require(nnls)
#     saveRDS(T,"Tnnls.rds")
#     saveRDS(C,"Cnnls.rds")
    print(dim(T))
    print(dim(C))
	RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
	rownames(RESULTS) <- colnames(C)
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}

run_MuSiC<-function(T.eset, C.eset, ...){
	
	#DEBUG# MuSiC requires ExpressionSet format input, which requires dense matrix. Thus, may cause memory issue. 
    
	require(MuSiC)
	RESULTS = t(MuSiC::music_prop(bulk.mtx = T.eset, sc.sce = C.eset, clusters = 'cellType',
										markers = NULL, normalize = FALSE, samples = 'SubjectName', 
										verbose = F)$Est.prop.weighted)
	return(RESULTS)
}

run_DWLS<-function(T, C, phenoDataC, elem, ...){
# 	require(DWLS)
	# source('/home/feng_shuo/deconvolution/Real_bulk/DWLS.R')
	path=paste(getwd(),"/DWLS_",sep="")

	if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created
		dir.create(path)
	} 
	Signature <- buildSignatureMatrixMAST(scdata = C, id = phenoDataC[,"cellType"], path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)
	
	RESULTS <- apply(T,2, function(x){
		b = setNames(x, rownames(T))
		tr <- trimData(Signature, b)
		RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
	})

	rownames(RESULTS) <- as.character(unique(phenoDataC$cellType))
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}

