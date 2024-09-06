#' An S4 class to represent the input data and parameters for Gibbs sampling in deconvolution analysis.
#' This class encapsulates the reference matrix, bulk matrix, cell type fractions, and control parameters
#' for the Gibbs sampling process.
#'
#' @slot reference A matrix of dimension KxG to denote the reference data (single cell data).
#' @slot X A matrix of dimension NxG to denote the bulk data (observed data).
#' @slot X.fraction A matrix of dimension NxK to denote the observed values in bulk data.
#' @slot gibbs.control A list containing parameters for controlling the Gibbs sampler.
#' @docType class
#' @name gibbsSampler
setClass("gibbsSampler",
         slots = c(
           reference = "matrix",
           X = "matrix",
		   X.fraction = "data.frame",
           gibbs.control = "list"
         ),
         prototype = list(
           reference = matrix(),
           X = matrix(),
		   X.fraction = data.frame(),
           gibbs.control = list(chain.length = NULL,
           					    burn.in = NULL,
           					    thinning = NULL,
           					    n.cores = NULL,
           					    seed = NULL,
           					    alpha = NULL)
         ),
         validity = function(object){
         	errors <- character()
         	if (length(errors) == 0) TRUE else errors
         }
)

#' function to return the index set of MCMC chain
#' @param gibbs.control a list containing parameters of the Gibbs sampler 
#' @return a numeric vector of the index of samples to be retained from MCMC chain 
get.gibbs.idx <- function(gibbs.control){
	chain.length <- gibbs.control$chain.length
	burn.in <- gibbs.control$burn.in
	thinning <- gibbs.control$thinning
	all.idx <- 1: chain.length
	burned.idx <-  all.idx[-(1: burn.in)]
	thinned.idx <- burned.idx[seq(1,length(burned.idx), by= thinning)]
	thinned.idx
}

#' function that generates one sample from Dirichlet distribution
#' param alpha a numeric vector to represent the parameter of dirichlet distribution
#' return a numeric vector (one sample from Dirichlet distribution)
rdirichlet <- function(alpha){
    l <- length(alpha)
    x <- rgamma(length(alpha), alpha)
    return(x/sum(x))
}

#' function to run Gibbs sampling for Z and theta on each bulk
#' @param X_n a numeric vector of reads count of the nth bulk sample
#' @param phi an array of dimension K*G to denote reference matrix
#' @param alpha a numeric value to denote the symmetric Dirichlet prior
#' @param cs_fraction a numeric vector of the cell type fraction of the nth bulk sample
#' @param gibbs.idx a numeric vector to denote the index of samples to be retained from MCMC chain
#' return a list containing the posterior mean of Z_n and theta_n
sample.Z <- function(X_n,
                     phi,
                     alpha,
					 cs_fraction,
                     gibbs.idx) {
    G <- ncol(phi)
    K <- nrow(phi)

	if (!is.null(cs_fraction)) {
        theta_n.i <- as.numeric(cs_fraction)
    } else {
        theta_n.i <- rep(1 / K, K)
    }

    Z_n.i <- array(NA, c(G, K))

    Z_n.sum <- array(0, c(G, K))

    for (i in 1:max(gibbs.idx)) {
        # sample Z for patient n
        prob.mat <- phi * theta_n.i # multiply by each column

        for (g in 1:G) {
            Z_n.i[g, ] <- rmultinom(
                n = 1,
                size = X_n[g],
                prob = prob.mat[, g]
            )
        }
        # sample theta for patient n
        Z_nk.i <- colSums(Z_n.i) # total count for each cell type
        theta_n.i <- rdirichlet(alpha = Z_nk.i + alpha)

        if (i %in% gibbs.idx) {
            # collect sample and compute posterior sum
            Z_n.sum <- Z_n.sum + Z_n.i
        }

        if ((i %% 50) == 0) gc()
    }

    samples.size <- length(gibbs.idx)
    # gibbs.constant <- multinom.coef + z.logtheta + theta.dirichlet.alpha

    Z_n <- Z_n.sum / samples.size

    return(list(Z_n = Z_n))
}

#' function to run initial Gibbs sampling if reference is of the class refPhi
#'
#' @param gibbsSampler.obj, a gibbsSampler object
#' @import snowfall
#' @return Z
run.gibbs.refPhi.ini <- function(gibbsSampler.obj) {
    phi <- gibbsSampler.obj@reference
    X <- gibbsSampler.obj@X
	fraction = gibbsSampler.obj@X.fraction
    gibbs.control <- gibbsSampler.obj@gibbs.control
    alpha <- gibbs.control$alpha

    gibbs.idx <- get.gibbs.idx(gibbs.control)
    seed <- gibbs.control$seed

    cat("Start run... \n")

    if (gibbs.control$n.cores > 1) {
        # parallel using snowfall
        sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK")

        cpu.fun <- function(n) {
            if (!is.null(seed)) set.seed(seed)
            # load nth mixture from disk
            file.name.X_n <- paste(tmp.dir, "/mixture_", n, ".rdata", sep = "")
            load(file.name.X_n)

            sample.Z(
                X_n = X_n,
                phi = phi,
                alpha = alpha,
				cs_fraction = fraction[n, ],
                gibbs.idx = gibbs.idx
            )
        }
        tmp.dir <- tempdir(check = TRUE)
        sfExport(
            "phi", "alpha", "gibbs.idx", "seed",
            "sample.Z", "tmp.dir", "fraction"
        )
        environment(cpu.fun) <- globalenv()
        gibbs.list <- sfLapply(1:nrow(X), cpu.fun)
        sfStop()
    } else {
        # single thread
        cpu.fun <- function(n) {
            if (!is.null(seed)) set.seed(seed)
            cat(n, " ")
            sample.Z(
                X_n = X[n, ], phi = phi, alpha = alpha,
				cs_fraction = fraction[n, ],
                gibbs.idx = gibbs.idx
            )
        }
        gibbs.list <- lapply(1:nrow(X), cpu.fun)
        cat("\n")
    }
    bulkID <- rownames(X)
    geneID <- colnames(X)
    cellType <- rownames(phi)

    N <- length(bulkID)
    G <- length(geneID)
    K <- length(cellType)

    stopifnot(length(gibbs.list) == N)

    Z <- array(NA,
        dim = c(N, G, K),
        dimnames = list(bulkID, geneID, cellType)
    )

    for (n in 1:N) {
        Z[n, , ] <- gibbs.list[[n]]$Z_n
    }
    return(Z)
}

#' function to run Gibbs sampling
#'
#' @param gibbsSampler.obj, a gibbsSampler object
#' @import snowfall
#' @return a gibbsSampler object with Z and theta entries
run.gibbs <- function(gibbsSampler.obj) {
    cat("Run Gibbs sampling... \n")
    return(run.gibbs.refPhi.ini(gibbsSampler.obj = gibbsSampler.obj))
}

#' function to validate opt.control
#' @param control a named list of parameters required to control optimization
valid.opt.control <- function(control) {
    ctrl <- list(
        maxit = 100000, maximize = FALSE, trace = 0, eps = 1e-07,
        dowarn = TRUE, tol = 0, maxNA = 500, n.cores = 1, optimizer = "MAP", sigma = 2
    )
    namc <- names(control)

    if (!all(namc %in% names(ctrl))) {
        stop("unknown names in opt.control: ", namc[!(namc %in% names(ctrl))])
    }
    ctrl[namc] <- control

    if (!ctrl$optimizer %in% c("MAP", "MLE")) {
        stop("unknown names of optimizer: ", ctrl$optimizer)
    }

    stopifnot(length(ctrl$optimizer) == 1)

    if (ctrl$optimizer == "MAP") {
        if (!is.numeric(ctrl$sigma)) {
            stop("sigma needs to be a numeric variable")
        } else {
            if (ctrl$sigma < 0) {
                stop("sigma needs to be positive")
            }
        }
    }

    return(ctrl)
}


#' function to validate gibbs.control
#' @param control a named list of parameters required to control optimization
valid.gibbs.control <- function(control) {
    ctrl <- list(chain.length = 1000, burn.in = 500, thinning = 2, n.cores = 1, seed = 123, alpha = 1)
    namc <- names(control)

    if (!all(namc %in% names(ctrl))) {
        stop("unknown names in gibbs.control: ", namc[!(namc %in% names(ctrl))])
    }
    ctrl[namc] <- control

    if (ctrl$alpha < 0) stop("alpha needs to be positive")

    return(ctrl)
}

#' function to summing up reads for each level in labels
#' @param ref a matrix of N*G
#' @param labels a character vector of length N
#' return a matrix of K*G, with K=number of unique levels in labels
collapse <- function(ref, labels) {
    stopifnot(nrow(ref) == length(labels))

    # remove NA in labels
    non.na.idx <- !is.na(labels)
    if (sum(!non.na.idx) > 0) print("Warning: NA found in the cell type/state labels. These cells will be excluded!")
    labels <- labels[non.na.idx]
    ref <- ref[non.na.idx, ]

    labels.uniq <- unique(labels)

    ref.collapsed <- do.call(
        rbind,
        lapply(
            labels.uniq,
            function(label.i) {
                colSums(ref[labels == label.i, , drop = F])
            }
        )
    )

    rownames(ref.collapsed) <- labels.uniq

    return(ref.collapsed)
}

#' function to normalize expression matrix, s.t. it sum up to one for each row, with the zero entries = pseudo.min
#' if no zero entries, return direct normalization (no pseudo.min returned)
#' @param ref a unnormalized matrix of dimension K*G (with rownames and colnames supplied)
#' @param pseudo.min the desired min values to replace zero after normalization
#' return a normalized matrix of the same dimension
norm.to.one <- function(ref,
                        pseudo.min) {
    G <- ncol(ref)

    phi <- ref / rowSums(ref) * (1 - pseudo.min * G) + pseudo.min

    # if the minimum value is greater than zero. simply normalize by total depth
    min.value <- apply(ref, 1, min)
    which.row <- min.value > 0
    if (any(which.row)) {
        # cat("One or more cell types have all genes with non-zero expression. pseudo.min is not applied to these cell types. \n")
        phi[which.row, ] <- ref[which.row, , drop = F] / rowSums(ref[which.row, , drop = F])
    }

    return(phi)
}

#' This function performs make single cell reference and bulk expression input.
#' @param reference A matrix containing the expression data of single cells.
#' @param mixture A matrix containing the expression data of bulk samples to be deconvoluted.
#' @param cell.state.labels A character vector where each element corresponds to the cell state label for each cell in the reference data.
#' @param cell.state.fraction the bulk-celltype matrix in the format of sample IDs (rows) by cell type fraction (columns).
#' @param outlier.cut A numeric value used to define the threshold for identifying outlier genes.
#'                     Genes with expression values greater than 'outlier.cut' in more than
#'                     'outlier.fraction' of the bulk samples are considered outliers.
#' @param outlier.fraction A numeric value representing the fraction of bulk samples
#'                        in which a gene's expression must exceed 'outlier.cut'
#'                        to be considered an outlier.
#' @param pseudo.min A small pseudocount added to the expression values to improve numerical stability.
#' @return A list containing the processed data matrices as input.
#' @export
make_input <- function(reference, mixture, cell.state.labels, cell.state.fraction, outlier.cut = 0.01, outlier.fraction = 0.1, pseudo.min = 1E-8) {
    # 确保 reference 和 mixture 是矩阵
    if (is.data.frame(reference)) reference <- as.matrix(reference)
    if (is.data.frame(mixture)) mixture <- as.matrix(mixture)

    # 处理 cell.state.labels
    cell.state.labels <- as.character(cell.state.labels)

    # 计算每个 cell state 的细胞数量
    cell_counts <- table(cell.state.labels)
    if (min(cell_counts) < 20) {
        cat("Recommend to have sufficient number of cells in each cell state\n")
        insufficient_states <- names(cell_counts[cell_counts < 20])
        cat("The following cell states have fewer than 20 cells:", paste(insufficient_states, collapse = ", "), "\n")
    }

    # 处理 mixture 数据
    if (is.null(dim(mixture))) mixture <- matrix(mixture, nrow = 1, dimnames = list("mixture-1", names(mixture)))

    # 处理 mixture 行名
    if (is.null(rownames(mixture))) rownames(mixture) <- paste("mixture", 1:nrow(mixture), sep = "-")

    # 过滤 bulk 中的异常值
    # mixture <- filter.bulk.outlier(mixture=mixture, outlier.cut=outlier.cut, outlier.fraction=outlier.fraction)

    # 移除 reference 中的零基因
    # reference <- reference[, colSums(reference) > 0]

    # 检查 reference 和 mixture 之间的基因注释是否匹配
    gene.shared <- intersect(colnames(reference), colnames(mixture))
    if (length(gene.shared) == 0) {
        stop("Error: gene names of reference and mixture do not match!")
    }
    if (length(gene.shared) < 100) {
        warning("Warning: very few gene from reference and mixture match! Please double check your gene names.")
    }

    # 折叠 reference
    ref.cs <- collapse(ref = reference, labels = cell.state.labels)

    # 对齐 reference 和 mixture
    cat("Aligning reference and mixture... \n")
    ref.cs <- ref.cs[, gene.shared, drop = F]
    mixture <- mixture[, gene.shared, drop = F]

	# 对应细胞比例，保证行对应X sample，列对应sc cell type
	fraction = cell.state.fraction[rownames(mixture), rownames(ref.cs)]

    # 归一化
    cat("Normalizing reference... \n")
    ref.cs <- norm.to.one(ref = ref.cs, pseudo.min = pseudo.min)

    return(list(ref.cs = ref.cs, mixture = mixture, fraction = fraction))
}

#' main function to run deconvolution and get cell type expressions.
#' @param reference expression matrix of single cell reference, should be raw counts. Must in the format of cell IDs (rows) by genes (columns).
#' @param cell.state.labels a character vector to denote cell state of each cell.
#' @param mixture the bulk RNA-seq matrix in the format of sample IDs (rows) by genes (columns).
#' 		colnames(mixture) need to match the annotation of colnames(reference) and will intersect colnames(reference) and colnames(mixture) and perform deconvolution using the intersected gene set.
#' @param cell.state.fraction the bulk-celltype matrix in the format of sample IDs (rows) by cell type fraction (columns).
#' @param n.cores number of cores to use. Default =1.
#' 		If needs to set different number of cores for Gibbs sampling and optimization,
#' 		supply n.cores in gibbs.control and/or opt.control
#' @param gibbs.control a list containing parameters of the Gibbs sampler
#' 		chain.length: length of MCMC chain. Default=1000;
#' 		burn.in: length of burn in period. Default=500;
#' 		thinning: retain every # of MCMC samples after the burn in period to reduce auto-correlation. Default=2;
#' 		n.cores: number of cores to use. Default uses n.cores in the main argument.
#' 		seed: seed number to use for repoducibility. Default = 123. set to NULL if use pseudo random.
#' 		alpha: a numeric vector to represent the parameter of dirichlet distribution. Default=1. (1E-8 may yield theta=0 due to underflow. causes issue when psuedo.min=0)
#' @param opt.control a list containing parameters for the optimization step:
#' 		maxit: maximum number of cycles to interate. Default=100000
#' 		sigma: hyper-parameter of the prior if optimizer="MAP". Default=2.
#' 		optimizer a character string to denote which algorithm to use
#' 			"MAP": the one used by the BayesPrism paper, with cell type-specific gamma under a log-normal prior
#' 			"MLE": the new algorithm that models a single gamma across cell types. Useful when some cell types are low in Z_k, e.g. spatial data
#' 	 		default to "MAP"
#' return a three-dimensional array Z, structured as samples * genes * cell types.
#' @export
deconvExp <- function(reference,
                      cell.state.labels,
                      bluk,
					  cell.state.fraction,
                      pseudo.min = 1E-8,
                      n.cores = 1,
                      gibbs.control = list(),
                      opt.control = list()) {
    if (!"n.cores" %in% names(gibbs.control)) gibbs.control$n.cores <- n.cores
    if (!"n.cores" %in% names(opt.control)) opt.control$n.cores <- n.cores
    stopifnot(is.numeric(n.cores) & length(n.cores) == 1)

    opt.control <- valid.opt.control(opt.control)
    gibbs.control <- valid.gibbs.control(gibbs.control)

    if (pseudo.min == 0) {
        gibbs.control$alpha <- max(1, gibbs.control$alpha)
    }

    # make input
    input <- make_input(reference, bluk, cell.state.labels, cell.state.fraction)

    bluk <- input$mixture
    ref.cs <- input$ref.cs

    # write bluk to disk and load each X_n to the corresponding node to save memory
    tmp.dir <- tempdir(check = TRUE)
    for (n in 1:nrow(bluk)) {
        X_n <- bluk[n, ]
        file.name <- paste(tmp.dir, "/bluk_", n, ".rdata", sep = "")
        save(X_n, file = file.name)
    }

    # sampling cell states (cs)
    gibbsSampler.ini.cs <- new("gibbsSampler",
        reference = ref.cs,
        X = bluk,
		X.fraction = input$fraction,
        gibbs.control = gibbs.control
    )

    Z <- run.gibbs(gibbsSampler.ini.cs)

    unlink(tmp.dir, recursive = TRUE)

    return(Z)
}
