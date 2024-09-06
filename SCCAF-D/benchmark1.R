# This file includes functions in the data preprocessing of the benchmarking process. 
# TODO: document the functions

options(stringsAsFactors = FALSE)
set.seed(24)
require(limma); require(dplyr); require(pheatmap); require(Matrix); require(sparseMatrixStats);require(reshape2)

# Added parameters: sampleCT = FALSE, propsample = TRUE, pct.var=30
# Pseudo bulk generator
Generator <- function(sce, phenoData, sampleCT = FALSE, propsample = TRUE, pct.var=30, Num.mixtures = 1000, 
					  pool.size = 100, min.percentage = 1, max.percentage = 99, seed = 24, mode = 2){ 

	CT = unique(phenoData$cellType)
	?stopifnot(length(CT) >= 2)

	set.seed(seed)
	require(dplyr)
	require(gtools)
	cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE) 
	colnames(cell.distribution) = c("CT","max.n")

	Tissues = list()
	Proportions = list()

	if(mode == 1){
		for(y in 1:Num.mixtures){

			#Only allow feasible mixtures based on cell distribution
			while(!exists("P")){
				if(sampleCT){
					num.CT.mixture = sample(x = 3:length(CT),1) # more than 3 cell types to fit a curve.
					selected.CT = sample(CT, num.CT.mixture, replace = FALSE)
				}else{
					num.CT.mixture = length(CT)
					selected.CT = CT
				}
				if(propsample){
					x = round(runif(num.CT.mixture, 100-pct.var, 100+pct.var))/100.0
					P = cell.distribution[selected.CT,]$max.n*x
				}else{
					P = runif(num.CT.mixture, min.percentage, max.percentage)
				}
				P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
				P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)

				missing.CT = CT[!CT %in% selected.CT]
				missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)
                
				P = rbind.data.frame(P, missing.CT)
				potential.mix = merge(P, cell.distribution)
				potential.mix$size = potential.mix$expected * pool.size

				if(sum(P$expected) != 1){
					rm(list="P") 
				}
			  
			}

			# Using info in P to build T simultaneously
			chosen_cells <- sapply(which(P$expected != 0), function(x){
				n.cells = P$expected[x] * pool.size
				chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]],
							  n.cells, replace = TRUE)
				chosen
			}) %>% unlist()
			T <- Matrix::rowSums(sce[,as.vector(chosen_cells)]) %>% as.data.frame()
			colnames(T) = paste("mix",y,sep="")

			P = P[,c("CT","expected")]
			P$mix = paste("mix",y,sep="")

			Tissues[[y]] <- T
			Proportions[[y]] <- P

			rm(list=c("T","P","chosen_cells","missing.CT"))

		}
	}
	else if(mode ==2){
		for(y in 1:Num.mixtures){
			chosen_cells = sample(phenoData$cellID, pool.size, replace = TRUE)
			T <- Matrix::rowSums(sce[,chosen_cells]) %>% as.data.frame()
			P1 = phenoData$cellType
			names(P1)<-phenoData$cellID
			P = data.frame(table(P1[chosen_cells]),stringsAsFactors = FALSE) 
			colnames(P) = c("CT","expected")
            colnames(T) = paste("mix",y,sep="")
			P$mix = paste("mix",y,sep="")
			Tissues[[y]] <- T
			Proportions[[y]] <- P
		}
	}

	P = do.call(rbind.data.frame, Proportions)
	T = do.call(cbind.data.frame, Tissues)
	P = data.table::dcast(P, CT ~ mix, 
						value.var = "expected",
						fun.aggregate = sum) %>% data.frame(.,row.names = 1) 

	P = P[,gtools::mixedsort(colnames(P))]
    
	return(list(T = T, P = P))
  
} 

Transformation <- function(matrix, option){
    
    #############################################################
    ##########    DATA TRANSFORMATION (on full data)   ##########
    if(option=="none"){
		
        matrix = matrix
        
    }else if(option=="log"){
		
        matrix = log1p(matrix)
        
    }else if(option=="sqrt"){
		
        matrix = sqrt(matrix)

    }else if(option=="vst"){
		
		#DEBUG# This function does not work for sparse matrix and may cause memory issue. 
        
        matrix = DESeq2::varianceStabilizingTransformation(as.matrix(matrix))

    }

    return(matrix)

}

Scaling <- function(matrix, option, phenoDataC=NULL){

    ##########    Remove rows & columns full of zeroes   ##########
    #Error: NMF::nmf - Input matrix x contains at least one null or NA-filled row.
    matrix = matrix[rowSums(matrix)!=0,]
    #OLS with error if all elements within a row are equal (e.g. all 0, or all a common value after log/sqrt/vst transformation)
    matrix = matrix[!apply(matrix, 1, function(x) var(x) == 0),]
    matrix = matrix[,colSums(matrix)!=0]

    if(option=="column"){
        
        matrix = apply(matrix,2,function(x) x/sum(x)) 

    } else if(option=="row"){ 
        
        matrix = t(apply(matrix,1,function(x) x/sum(x))) 

    } else if(option=="mean"){ 
        
        matrix = apply(matrix,2,function(x) x - mean(x)) 

    } else if(option=="column_z-score"){ 
        
        matrix = scale(matrix, center = TRUE, scale = TRUE)

    } else if(option=="global_z-score"){
		
		#DEBUG# global_z-score desifies sparseMatrix and may cause memory issue.
        
        matrix = (matrix - mean(matrix)) / sd(matrix)

    } else if(option=="column_min-max"){
        
        matrix = apply(matrix, 2, function(x) (x - min(x))/(max(x) - min(x)))

    } else if(option=="global_min-max"){
        
        matrix = (matrix - min(matrix))/(max(matrix) - min(matrix))

    } else if (option=="LogNormalize"){
		
		#DEBUG# Seurat::LogNormalize desifies sparseMatrix and may cause memory issue.
        
        matrix = as.matrix(expm1(Seurat::LogNormalize(matrix, verbose = FALSE))) #for v3

    } else if (option=="none"){
      
      #DEBUG# Seurat::LogNormalize desifies sparseMatrix and may cause memory issue.
      
      matrix = as.matrix(matrix) #for v3
      
    } else if (option=="QN"){
		
		#DEBUG# QN requires matrix format input and may cause memory issue.

        matrix_rownames <- rownames(matrix); matrix_colnames <- colnames(matrix)

        matrix = preprocessCore::normalize.quantiles(as.matrix(matrix))

        rownames(matrix) <- matrix_rownames; colnames(matrix) <- matrix_colnames

    } else if (option=="TMM"){# CPM counts coming from TMM-normalized library sizes; https://support.bioconductor.org/p/114798/

        if(!is.null(phenoDataC)){#use CT info for scRNA-seq

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }
        matrix <- edgeR::DGEList(counts=matrix, group=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(data))
        if(length(CtrlGenes)>1){
            spikes <- data[CtrlGenes,]
            spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
            matrix$samples$norm.factors <- spikes$samples$norm.factors
        
        } else {
      
            matrix <- edgeR::calcNormFactors(matrix, method = "TMM")  
        
        }
        matrix <- edgeR::cpm(matrix)

    } else if (option=="UQ"){

        if(!is.null(phenoDataC)){#use CT info for scRNA-seq

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }

        matrix <- edgeR::DGEList(counts=matrix, group=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(data))
  
        if(length(CtrlGenes)>1){
            
            spikes <- data[CtrlGenes,]
            spikes <- edgeR::calcNormFactors(spikes, method = "upperquartile") 
            matrix$samples$norm.factors <- spikes$samples$norm.factors
        
        } else {
        
            matrix <- edgeR::calcNormFactors(matrix, method = "upperquartile")  
        
        }

        matrix <- edgeR::cpm(matrix)

    } else if (option=="median_ratios"){#requires integer values
        
        if(!is.null(phenoDataC)){#use CT info for scRNA-seq

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }

        metadata <- data.frame(Celltype=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(matrix))
        matrix = DESeq2::DESeqDataSetFromMatrix(matrix, colData = metadata, design = ~ Celltype)

        if(length(CtrlGenes)>1 & sum(rowSums(DESeq2::counts(matrix[CtrlGenes,]) != 0) >= 0.5*(ncol(matrix))) >= 2){

            dds <- DESeq2::estimateSizeFactors(matrix, type = "ratio", controlGenes = CtrlGenes)

        } else {

            dds <- DESeq2::estimateSizeFactors(matrix, type = "ratio")
            
        }

        matrix <- DESeq2::counts(dds, normalized=TRUE)
    
    } else if (option=="TPM"){

        require(SingleR)
        data(human_lengths)

        # Doesn't work with Ensembl IDs:
        if(length(grep("ENSG000",rownames(matrix))) > 100){
            
            suppressMessages(library("AnnotationDbi"))
            suppressMessages(library("org.Hs.eg.db"))
            temp = mapIds(org.Hs.eg.db,keys=names(human_lengths),column="ENSEMBL",keytype="SYMBOL",multiVals="first")
            names(human_lengths) = as.character(temp)

        }
        
        matrix = SingleR::TPM(counts = matrix, lengths = human_lengths)
        rownames(matrix) = toupper(rownames(matrix))
        detach(package:SingleR, unload=TRUE)
    
    ####################################################################################
    ## scRNA-seq specific  

    } else if (option=="SCTransform"){# SCTransform = RNBR

        #Model formula is y ~ log_umi
        ##following line needed to solve "Wrong R type for mapped matrix"
        #https://github.com/ChristophH/sctransform/issues/24
        matrix = as(matrix, "dgCMatrix")
        matrix = sctransform::vst(matrix, return_corrected_umi=TRUE, verbosity = FALSE)$umi_corrected

    } else if (option=="scran"){
        sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts=matrix))
        sce <- scran::computeSumFactors(sce, clusters=NULL) 
		sce <- scater::logNormCounts(sce, log=FALSE)
		matrix = sce@assays@data$normcounts
        
    } else if (option=="scater"){  

        size_factors = scater::librarySizeFactors(matrix)
        matrix <- scater::normalizeCounts(matrix, size_factors = size_factors)

    } else if (option=="Linnorm"){
		matrix = expm1(Linnorm::Linnorm(matrix))
    }

    return(matrix)

}

########################################################################
read_data<-function(dataset){
    
    ad = readRDS(dataset)
    if (class(ad)== "Seurat") {
        data = ad@assays$RNA@counts 
        pData = ad@meta.data
    } else if (class(ad) =="ExpressionSet"){
        data = ad@assayData$exprs
        pData = ad@phenoData@data
               
    } else {
		print("the format of the input data is neither Seurat nor ExpressionSet")
    }
    rownames(pData)<-NULL
    
    colnames(data)<-pData$cellID
	original_cell_names = colnames(data)
    colnames(data) <- as.character(pData$cellType)
	cell_counts = table(colnames(data))
	return(list('data'=data, 'pData'=pData, 'original_cell_names'=original_cell_names, 'cell_counts'=cell_counts))
}

# Read bulk and pseudobulk data
read_bulk<-function(dataset){
    if(endsWith(dataset, ".rds")){
        data = readRDS(dataset)
    }else{
        data = read.csv(dataset, row.names=1)
    }
	
	if (class(data)== "list") {
		return(data)
	} else{
		return(list('data'=data))
	}
}
                       
QC<-function(X){
	# filter features
	keep <- which(Matrix::rowSums(data > 0) >= 3)
	X$data = X$data[keep,]
	# Keep CTs with >= 4 cells after QC
	to_keep = names(X$cell_counts)[X$cell_counts >= 4]
	X$pData <- X$pData[X$pData$cellType %in% to_keep,]
	to_keep = which(colnames(X$data) %in% to_keep)   
	X$data <- X$data[,to_keep]
	X$original_cell_names <- X$original_cell_names[to_keep]
	return(X)
}

ref_mat<-function(train){
	cellType <- colnames(train)
	group = list()
	for(i in unique(cellType)){ 
		group[[i]] <- which(cellType %in% i)
	}
	#C should be made with the mean (not sum) to agree with the way markers were selected
	C = lapply(group,function(x) Matrix::rowMeans(train[,x]))
	C = round(do.call(cbind.data.frame, C))

	refProfiles.var = lapply(group,function(x) train[,x])
	refProfiles.var = lapply(refProfiles.var, function(x) sparseMatrixStats::rowSds(x))
	refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
	rownames(refProfiles.var) <- rownames(train)
	return(list('C'=C, 'ref'=refProfiles.var))
}

# This function is used in the `marker_selection` function. 
NormalizationTMM <- function(data){

    data <- edgeR::DGEList(data)
    CtrlGenes <- grep("ERCC-",rownames(data))

    if(length(CtrlGenes)>1){

        spikes <- data[CtrlGenes,]
        spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
        data$samples$norm.factors <- spikes$samples$norm.factors

    } else {

        data$samples$lib.size[data$samples$lib.size<1] = 100
        data <- edgeR::calcNormFactors(data, method = "TMM")  

    }

    return(data)

}

marker.fc <- function(fit2, cont.matrix, log2.threshold = 1, output_name = "markers"){
  
	topTable_RESULTS = limma::topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2.threshold)
    
    if(dim(topTable_RESULTS)[1]<1){
        topTable_RESULTS = limma::topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", 
                                       p.value = 0.05, lfc = 0.1)
        log2.threshold = 0.1
    }
	AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
	topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]
	if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }

	markers <- apply(topTable_RESULTS,1,function(x){
	temp = sort(x)
	((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)

	})
	topTable_RESULTS = topTable_RESULTS[markers,]
	markers <- cbind.data.frame(rownames(topTable_RESULTS),
	                                   t(apply(topTable_RESULTS, 1, function(x){
	                                     temp = max(x)
	                                     if(temp < log2.threshold){
	                                       temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
	                                     } else {
	                                       temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
	                                     } 
	                                     temp
	                                   })))
	colnames(markers) <- c("gene","log2FC","CT")
	markers$log2FC = as.numeric(as.character(markers$log2FC))
	markers <- markers %>% dplyr::arrange(CT,desc(log2FC)) 
	markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
	markers$gene <- as.character(markers$gene)
	markers$CT <- as.character(markers$CT)
	return(markers)
  
}

marker_strategies <- function(marker_distrib, marker_strategy, C){

    set.seed(4)

    if(marker_strategy == "all"){
        
        #using all markers that were found
        markers = marker_distrib

    } else if (marker_strategy == "pos_fc"){

        # using only markers with positive FC (=over-expressed in cell type of interest)
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% as.data.frame()

    } else if (marker_strategy == "top_50p_logFC"){

        # top 50% of markers (per CT) based on logFC
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "bottom_50p_logFC"){

        # bottom 50% of markers based on logFC
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(floor(n()*-0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "top_50p_AveExpr"){

        # top 50% of markers based on average gene expression (baseline expression)
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(AveExpr)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "bottom_50p_AveExpr"){

        # low 50% based on average gene expression.
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(AveExpr)) %>% group_by(CT) %>% dplyr::top_n(floor(n()*-0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "top_n2"){

        # using the top 2 genes/CT with highest log2FC
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(2, wt=log2FC) %>% as.data.frame()
    
    } else if (marker_strategy == "random5"){

        # using 5 random markers for each different cell types
        markers = marker_distrib[1:(ncol(C)*5),]
        markers$CT = rep(colnames(C),5) #labelling purposes: important for semi-supervised
        markers$gene = sample(rownames(C), nrow(markers), replace = FALSE)

    }

    return(markers) 

}

marker_selection<-function(train){
	#-------------------------------------------------------
	#Normalization of "train" followed by marker selection 
	#for marker selection, keep genes where at least 30% of cells within a cell type have a read/UMI count different from 0
	cellType = colnames(train) 
	keep <- sapply(unique(cellType), function(x) {
		CT_hits = which(cellType %in% x)
		size = ceiling(0.3*length(CT_hits))
		Matrix::rowSums(train[,CT_hits,drop=FALSE] != 0) >= size
	})
#     print("====================================1")
	train = train[Matrix::rowSums(keep) > 0,]
	train2 = NormalizationTMM(train)
	# INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account 
	#[compare one group with average expression of all other groups]
	annotation = factor(colnames(train2))
    if(length(unique(annotation))<2){
        print("Error in cell type annotation!!!")
    }
    
	design <- model.matrix(~0+annotation)
	colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
	cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
	colnames(cont.matrix) <- colnames(design)
	diag(cont.matrix) <- (ncol(design)-1)/ncol(design)
	v <- limma::voom(train2, design=design, plot=FALSE) 
	fit <- limma::lmFit(v, design)
	fit2 <- limma::contrasts.fit(fit, cont.matrix)
	fit2 <- limma::eBayes(fit2, trend=TRUE)
	markers = marker.fc(fit2, cont.matrix=cont.matrix, log2.threshold = log2(2))
	return(markers)
}

prepare_train<-function(train, train_cell_names){
	# Generate phenodata for reference matrix C
	train_cellID = train
	colnames(train_cellID) = train_cell_names
	# reference matrix (C) + refProfiles.var from TRAINING dataset
	cellType <- colnames(train)
	group = list()
	for(i in unique(cellType)){ 
		group[[i]] <- which(cellType %in% i)
	}
	#C should be made with the mean (not sum) to agree with the way markers were selected
	C = lapply(group, function(x) if(class(train[,x])=='numeric') train[,x] else Matrix::rowMeans(train[,x]) )
	C = round(do.call(cbind.data.frame, C))
	refProfiles.var = lapply(group,function(x) train[,x])
    refProfiles.var = lapply(refProfiles.var, function(x) if(class(x)=='numeric') x else sparseMatrixStats::rowSds(x)) 
	refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
	rownames(refProfiles.var) <- rownames(train)
	markers = marker_selection(train)
	return(list('train_cellID'=train_cellID, 'C'=C, 'ref' = refProfiles.var, 'markers'=markers))
}

QC_feature<-function(data){
	
	keep <- which(Matrix::rowSums(data > 0) >= 3)
	data = data[keep,]
}
