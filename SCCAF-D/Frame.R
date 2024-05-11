#' @title Framework
#' 
#' @description
#' This function contains the whole framework after data preprocessing and before deconvolution.
#' @details
#' This function is wrapped, because self_reference, cross_reference and bulk_reference can all reuse it. 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' 
#' @return
#' The final deconvolution results. 
#' @example
#' 
Framework<-function(deconv_type,
					NormTrans,
					Xtest,
					Xtrain,
					normalization,
					normalization_scT,
					normalization_scC,
					transformation,
					marker_strategy,
					to_remove,
					P,
					method,
					pDataC){
	#-------------------------------------------------------
	### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:
	if(deconv_type == "bulk"){
		if(NormTrans){
			T = Scaling(Xtest$T, normalization)
			C = Scaling(Xtrain$C, normalization)
			
			T = Transformation(T, transformation)
			C = Transformation(C, transformation)
		}else{
			T = Transformation(Xtest$T, transformation)
			C = Transformation(Xtrain$C, transformation)

			T = Scaling(T, normalization)
			C = Scaling(C, normalization)
		}

		# marker selection (on training data) 
		marker_distrib = marker_strategies(Xtrain$markers, marker_strategy, C)

		#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
		if(to_remove != "none"){

			T <- T[,P[to_remove,] != 0]
			C <- C[, colnames(C) %in% rownames(P) & (!colnames(C) %in% to_remove)]
			P <- P[!rownames(P) %in% to_remove, colnames(T)]
			Xtrain$ref = Xtrain$ref[,colnames(Xtrain$ref) %in% rownames(P) & (!colnames(Xtrain$ref) %in% to_remove)]
			marker_distrib <- marker_distrib[marker_distrib$CT %in% rownames(P) & (marker_distrib$CT != to_remove),]
		}
	} else if (deconv_type == "sc"){
		if(NormTrans){
			T = Scaling(Xtest$T, normalization_scT)
			C = Scaling(Xtrain$train_cellID, normalization_scC)
			
			T = Transformation(T, transformation)
			C = Transformation(C, transformation)
		}else{
			T = Transformation(Xtest$T, transformation)
			C = Transformation(Xtrain$train_cellID, transformation)

			T = Scaling(T, normalization_scT)
			C = Scaling(C, normalization_scC)
		}
		#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
		if(to_remove != "none"){

			T <- T[,P[to_remove,] != 0]
			C <- C[,pDataC$cellType != to_remove]
			P <- P[!rownames(P) %in% to_remove, colnames(T)]
			pDataC <- pDataC[pDataC$cellType != to_remove,]

		}
		marker_distrib = NULL
	}
	## 
	# T is the pseudo-bulk exprs
	# C is the reference exprs (average for each cell type)
	# P is the preassigned proportion for the pseudo-bulk
	# pDataC is the meta-data for the reference
#     saveRDS(T, 'T.rds')
#     saveRDS(C, 'C.rds')
#     saveRDS(pDataC, 'pDataC.rds')
	RESULTS = Deconvolution(T = T, 
							C = C, 
							method = method, 
							P = P, 
							elem = to_remove, 
							refProfiles.var = Xtrain$ref, 
							STRING = as.character(sample(1:10000, 1)), 
# 							STRING = "text", 
							marker_distrib = marker_distrib, 
							phenoDataC = pDataC) 
    
	return(RESULTS)
}
###
selection_cells=function(pathway,dataset1,dataset){
    data <- read.csv(paste(pathway,'/sccaf-result.csv',sep=''))
    test <- read.csv(paste(pathway,'/y_test-result.csv',sep=''))
    data <- data[,-1]
    sccaf <- cbind(test,data)
    rownames(sccaf) <- sccaf[,1]
    sccaf <- sccaf[,-1]
    index <- sccaf$cellType==sccaf$predict
    sccaf_index <- sccaf[index,]
    library(dplyr)
    mydata <- data.frame(matrix(ncol=2,nrow=0))
    for(i in 3:ncol(sccaf_index)){
        name <- arrange(sccaf_index,desc(sccaf_index[,i]))
        name <- name[c(1:100),]
        name <- name[,c(1,i)]
        colnames(name) <- c('celltype','score')
        mydata <- rbind(mydata,name)
    }  
    dataset1=dataset1
    save(mydata,file=paste(pathway,'/',dataset1,'_',"sccaf-cells",".RData",sep=""))
    dataset=dataset
    X2=readRDS(paste(pathway,'/',dataset,sep=''))
    X2=X2[,rownames(mydata)]
    saveRDS(X2,file = paste(pathway,'/',dataset1,'_',"sccaf-reference",".rds",sep=""))
    return(mydata)
}
