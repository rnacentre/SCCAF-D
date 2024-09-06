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

    RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)),]                  
    RESULTS = data.table::melt(RESULTS)                   
	colnames(RESULTS) <-c("CT","tissue","observed_values")

	if(!is.null(P)){
		P = P[gtools::mixedsort(rownames(P)),]
		P$CT = rownames(P)
		P = data.table::melt(P, id.vars="CT")
		colnames(P) <-c("CT","tissue","expected_values")
		RESULTS = merge(RESULTS,P)
		RESULTS$expected_values <-round(RESULTS$expected_values,3)
		RESULTS$observed_values <-round(RESULTS$observed_values,3)
	}

    return(RESULTS) 

}

########################################################################
# bulk Methods
run_CAMmarker<-function(T, C, marker_distrib, ...){
	library(debCAM)
	 #Full version, irrespective of C
	ML = CellMix::MarkerList()
	ML@.Data <- tapply(as.character(marker_distrib$gene),as.character(marker_distrib$CT),list)
	RESULTS = t(AfromMarkers(T, ML))
	colnames(RESULTS) <- colnames(T)
	rownames(RESULTS) <- names(ML)
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}
run_CDSeq<-function(T, C, marker_distrib, ...){

	
	RESULTS<-CDSeq::CDSeq(bulk_data = T, reference_gep = C, cell_type_number = length(unique(marker_distrib$CT)), mcmc_iterations = 1000,
                                  cpu_number=10,block_number=6,gene_subset_size=15)$estProp
	return(RESULTS)
}
run_CIBERSORT<-function(T, C, ...){
	source('CIBERSORT.R')
	write.table(C, file = xf <- 'reference.tsv', sep = "\t", row.names = TRUE, col.names = NA)
	write.table(T, file = yf <- 'mixture.tsv', sep = "\t", row.names = TRUE, col.names = NA)
	message("* Running CIBERSORT ... ", appendLF = FALSE)
	RESULTS = CIBERSORT(sig_matrix = xf, mixture_file = yf, QN = FALSE) 
	RESULTS = t(RESULTS[,1:(ncol(RESULTS)-3)])
	return(RESULTS)
}

run_DCQ<-function(T, C, ...){
	#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default
	require(ComICS)
	RESULTS = t(ComICS::dcq(reference_data = C, mix_data = T, marker_set = as.data.frame(row.names(C)) , alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average)
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
	return(RESULTS)
}

run_DSA<-function(T, C, marker_distrib, ...){
	#DSA algorithm assumes that the input mixed data are in linear scale; If log = FALSE the data is left unchanged

	require(CellMix)

	ML = CellMix::MarkerList()
	ML@.Data <- tapply(as.character(marker_distrib$gene),as.character(marker_distrib$CT),list)
	RESULTS = CellMix::ged(as.matrix(T), ML, method = "DSA", log = FALSE)@fit@H
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x))
	return(RESULTS)
}

run_DeconRNASeq<-function(T, C, ...){
	#nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
	#datasets and reference matrix: signatures, need to be non-negative. 
	#"use.scale": whether the data should be centered or scaled, default = TRUE
	unloadNamespace("Seurat") #needed for PCA step
	library(pcaMethods) #needed for DeconRNASeq to work
	RESULTS = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(T), 
										 signatures = as.data.frame(C), 
										 proportions = NULL, 
										 checksig = FALSE, 
										 known.prop = FALSE, 
										 use.scale = FALSE, 
										 fig = FALSE)$out.all)
	colnames(RESULTS) = colnames(T)
# 	require(Seurat)
	return(RESULTS)
}

run_EPIC<-function(T, C, marker_distrib, refProfiles.var, ...){
	

	require(EPIC)
	marker_distrib = marker_distrib[marker_distrib$gene %in% rownames(C),]
	markers = as.character(marker_distrib$gene)
	C_EPIC <- list()

	common_CTs <- intersect(colnames(C),colnames(refProfiles.var))

	C_EPIC[["sigGenes"]] <- rownames(C[markers,common_CTs])
	C_EPIC[["refProfiles"]] <- as.matrix(C[markers,common_CTs])
	C_EPIC[["refProfiles.var"]] <- refProfiles.var[markers,common_CTs]

	RESULTS <- t(EPIC::EPIC(bulk=as.matrix(T), reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) 

	RESULTS = RESULTS[!rownames(RESULTS) %in% "otherCells",]
	RESULTS[is.na(RESULTS)] <- 0
	return(RESULTS)
}

run_EpiDISH<-function(T, C, ...){
	#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default

	require(EpiDISH)
	RESULTS = t(EpiDISH::epidish(beta.m = T, ref.m = C, method = "RPC")$estF)
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
	return(RESULTS)
}

run_FARDEEP<-function(T, C, ...){
	require(FARDEEP)
	RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
	return(RESULTS)
}

run_OLS<-function(T, C, ...){
	# OLS = ordinary least squares
	
	RESULTS = apply(T,2,function(x) lm(x ~ as.matrix(C))$coefficients[-1])
	rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
	RESULTS[is.na(RESULTS)] <- 0
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x))
    return(RESULTS)
}

run_RLR<-function(T, C, ...){
	# RLR = robust linear regression
	require(MASS)
	RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
	rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x))
	return(RESULTS)
}

run_deconf<-function(T, C, marker_distrib, ...){

	##source("deconf.R")

	require(CellMix)
	#Full version, irrespective of C
	ML = CellMix::MarkerList()
	ML@.Data <- tapply(as.character(marker_distrib$gene),as.character(marker_distrib$CT),list)

	res <- CellMix::ged(T, x=length(unique(as.character(marker_distrib$CT))), method = "deconf", maxIter = 500, verbose= TRUE) 
	RESULTS<- CellMix::match.nmf(res, ML)@fit@H 
    
	return(RESULTS)
}

run_dtangle<-function(T, C, marker_distrib, ...){
	#Only works if T & C are log-transformed

	require(dtangle)
	mixture_samples = t(T)
	reference_samples = t(C)

	marker_distrib <- tidyr::separate_rows(marker_distrib,"CT",sep="\\|")    
	marker_distrib = marker_distrib[marker_distrib$gene %in% rownames(C),] 
	MD <- tapply(marker_distrib$gene,marker_distrib$CT,list)
	MD <- lapply(MD,function(x) sapply(x, function(y) which(y==rownames(C))))

	RESULTS = t(dtangle::dtangle(Y=mixture_samples, reference=reference_samples, markers = MD)$estimates)
	return(RESULTS)
}

run_elasticNet<-function(T, C, ...){
	#standardize = TRUE by default. lambda=NULL by default 

	require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
	RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0.2, standardize = TRUE, 
														  lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
	return(RESULTS)
}

run_lasso<-function(T, C, ...){ 
	#alpha=1; shrinking some coefficients to 0. 

	require(glmnet)
	RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 1, standardize = TRUE, 
														  lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x))
	RESULTS[is.na(RESULTS)] <- 0 
	return(RESULTS)
}

run_ridge<-function(T, C, ...){
	#alpha=0

	require(glmnet)
	RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0, standardize = TRUE, 
														  lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
	return(RESULTS)
}

run_nnls<-function(T, C, ...){
	require(nnls)
    print(dim(T))
    print(dim(C))
	RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
	rownames(RESULTS) <- colnames(C)
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
	return(RESULTS)
}

run_proportionsInAdmixture<-function(T, C, ...){
	#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default

	require(ADAPTS)
	RESULTS = ADAPTS::estCellPercent(refExpr = C, geneExpr = T, method="proportionsInAdmixture")
	RESULTS[is.na(RESULTS)] <- 0  ####Anna## convert NAs to zeros so you can apply sum to one constraint
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
	RESULTS <- RESULTS[-nrow(RESULTS),]  
	return(RESULTS)
}


run_ssFrobenius<-function(T, C, marker_distrib, ...){

	require(CellMix) # require NMF 0.23.0, but NMF 0.30.1 does not work.
	md = marker_distrib #Full version, irrespective of C
	ML = CellMix::MarkerList()
	ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
	RESULTS <- CellMix::ged(as.matrix(T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE)@fit@H #equivalent to coef(CellMix::ged(T,...)
	return(RESULTS)
}

run_ssKL<-function(T, C, marker_distrib, ...){ 

	require(CellMix)
	md = marker_distrib #Full version, irrespective of C
	ML = CellMix::MarkerList()
	ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
	
	res <- CellMix::ged(as.matrix(T), x=length(unique(as.character(marker_distrib$CT))),
														method = "ssKL",markers= "semi", 
														sscale = FALSE, maxIter=500, 
														log = FALSE, 
														verbose= TRUE)  
    RESULTS<- CellMix::match.nmf(res, ML)@fit@H  
	return(RESULTS)
}




########################################################################
# Single-cell Methods
run_BisqueRNA<-function(T.eset, C.eset, ...){
	#By default, BisqueRNA uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter
	# Anna multisubject method, needs more that one individual in the reference C
	#DEBUG# BisqueRNA requires ExpressionSet format input, which requires dense matrix. Thus, may cause memory issue. 
	
	require(BisqueRNA)
	RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, 
													  C.eset, 
													  markers=NULL, 
													  use.overlap=FALSE)$bulk.props 
	return(RESULTS)
}

run_CPM<-function(T, C, phenoDataC, ...){
	#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default

	require(scBio)
	TReduced = T - rowMeans(T)
	p <- prcomp(t(C), center = TRUE,scale. = TRUE)$x[,1:2]
	celltypes.sc = as.character(phenoDataC$cellType)
	
	RESULTS = CPM(C, celltypes.sc, TReduced, p, quantifyTypes = TRUE, no_cores = 6)$cellTypePredictions
	RESULTS = t(RESULTS)
	RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) 
	RESULTS = apply(RESULTS,2,function(x) x/sum(x)) 
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
run_bseqsc<-function(T, C, T.eset, C.eset, ...){

	require(bseqsc)
	bseqsc_config('CIBERSORT.R')
	cells = unique(C.eset@phenoData@data$cellType)
	keep = intersect(rownames(C),rownames(T)) 
	markers<-list()
	for(cell in cells){
		markers[cell]<-list(keep)
	}
	B.eset = bseqsc_basis(C.eset, markers, clusters = 'cellType', samples = 'SubjectName', ct.scale = TRUE)
	fit <- bseqsc_proportions(T.eset, B.eset, verbose = TRUE)
	RESULTS = coef(fit)
	return(RESULTS)
}
run_DWLS<-function(T, C, phenoDataC, elem, ...){
    path=paste(getwd(),"/sccaf-d/",sep="")
	if(!dir.exists(path)){ 
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

run_deconvSeq<-function(T, C, T.eset, C.eset, phenoDataC, ...){
      
	singlecelldata = C.eset 
	celltypes.sc = as.character(phenoDataC$cellType) #To avoid "Design matrix not of full rank" when removing 1 CT 
	tissuedata = T.eset 

	design.singlecell = model.matrix(~ -1 + as.factor(celltypes.sc))
	colnames(design.singlecell) = levels(as.factor(celltypes.sc))
	rownames(design.singlecell) = colnames(singlecelldata)

	dge.singlecell = deconvSeq::getdge(singlecelldata,design.singlecell, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
	b0.singlecell = deconvSeq::getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min =1, nsamp.min = 4)
	dge.tissue = deconvSeq::getdge(tissuedata, NULL, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
	
	RESULTS = t(deconvSeq::getx1.rnaseq(NB0 = "all",b0.singlecell, dge.tissue)$x1) #all genes
	return(RESULTS)
}
run_SCDC<-function(T, C, T.eset, C.eset, phenoDataC, ...){
	##Proportion estimation with traditional deconvolution + >1 subject
	require(SCDC)
	RESULTS <- t(SCDC::SCDC_prop(bulk.eset = T.eset, 
								 sc.eset = C.eset, 
								 ct.varname = "cellType", 
								 sample = "SubjectName", 
								 ct.sub = unique(as.character(phenoDataC$cellType)), 
								 iter.max = 200)$prop.est.mvw)
	return(RESULTS)
}
run_TIMER<-function(T, C, phenoDataC, ...){
	source('TIMER.R')
	ref_anno <- phenoDataC$cellID
	names(ref_anno)<- phenoDataC$cellType
	RESULTS = t(TIMER_deconv(mix = T, ref = C, curated.cell.types = ref_anno, sig = rownames(T)))
	return(RESULTS)
}
