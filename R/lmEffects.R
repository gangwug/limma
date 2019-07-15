.lmEffects <- function(y,design=NULL,contrast=ncol(design),array.weights=NULL,gene.weights=NULL,weights=NULL,block=NULL,correlation)
#	Compute matrix of effects from genewise linear models
#	Gordon Smyth
#	Created 11 Apr 2016.  Last modified 4 Feb 2018.
{
#	Extract components from y
	y <- getEAWP(y)
	ngenes <- nrow(y$exprs)
	n <- ncol(y$exprs)

#	Check y
	if(anyNA(y$exprs)) stop("All y values must be finite and non-NA")

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		stop("design matrix not specified")
	} else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	p <- ncol(design)

	if(n <= p) stop("No residual degrees of freedom")

#	Check contrast
	if(is.character(contrast)) {
		if(length(contrast)>1L) {
			warning("using only first entry for contrast")
			contrast <- contrast[1]
		}
		contrast <- which(contrast==colnames(design))
		if(length(contrast)==0L) stop("coef ",contrast," not found")
	}
	if(all(contrast == 0)) stop("contrast all zero")

#	Reform design matrix so that contrast is last coefficient
	if(length(contrast) == 1L) {
		contrast <- as.integer(contrast)
		if(contrast < p)
			X <- cbind(design[,-contrast,drop=FALSE],design[,contrast,drop=FALSE])
		else
			X <- design
	} else {
		if(length(contrast) != p) stop("length of contrast must match column dimension of design")
		X <- contrastAsCoef(design, contrast, first=FALSE)$design
	}
	CoefName <- colnames(X)[p]
	if(is.null(CoefName)) CoefName <- "Contrast"

#	Allow array.weights to be alternatively passed via 'weights', as per lmFit documentation
	if(is.null(array.weights) && length(weights)==n) {
		array.weights <- weights
		weights <- NULL
	}

#	Check array.weights
	if(!is.null(array.weights)) {
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of arrays")
		AnyNeg <- any(array.weights <= 0)
		if(anyNA(AnyNeg) || AnyNeg) stop("array.weights must be positive")
	}

#	Allow gene.weights to be alternatively passed via 'weights', as per lmFit documentation
	if(is.null(gene.weights) && length(weights)==ngenes) {
		gene.weights <- weights
		weights <- NULL
	}

#	Check gene.weights
	if(!is.null(gene.weights)) {
		if(length(gene.weights) != ngenes) stop("Length of gene.weights doesn't match number of genes")
		AnyNeg <- any(gene.weights <= 0)
		if(anyNA(AnyNeg) || AnyNeg) stop("gene.weights must be positive")
	}

#	Check observation weights
	if(is.null(weights)) weights <- y$weights
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		dimw <- dim(weights)
		if(dimw[1] != ngenes || dimw[2] != n) stop("weights must have same dimensions as y")
		AnyNeg <- any(weights <= 0)
		if(anyNA(AnyNeg) || AnyNeg) stop("weights must be positive")
	}

#	Reduce to numeric expression matrix
	y <- y$exprs
	geneid <- rownames(y)
	if(is.null(geneid)) geneid <- 1:ngenes

#	Divide out array weights
	if(!is.null(array.weights)) {
		ws <- sqrt(array.weights)
		X <- X*ws
		y <- .matvec(y,ws)
		array.weights <- NULL
	}

#	Correlation matrix
	if(!is.null(block)) {
		if(missing(correlation) || is.null(correlation)) stop("correlation must be set")
		block <- as.vector(block)
		if (length(block) != n) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,n,nblocks) == matrix(ub,n,nblocks,byrow = TRUE)
		cormatrix <- Z %*% (correlation * t(Z))
		diag(cormatrix) <- 1
		R <- chol(cormatrix)

#		Divide out correlations if possible
		if(is.null(weights)) {	
			y <- t(backsolve(R, t(y), transpose = TRUE))
			X <- backsolve(R, X, transpose = TRUE)
		}
 	}

	qrX <- qr(X)
	if(qrX$rank < p) stop("design must be full column rank")

#	Compute effects for contrasts and residuals
	if(is.null(weights)) {
		Effects <- t(qr.qty(qrX,t(y)))
		signc <- sign(qrX$qr[p,p])
#		Remove model effects other than the required contrast
		if(p>1) Effects <- Effects[,p:n,drop=FALSE]
#		Preserve sign of estimated effect
		if(signc<0) Effects[,1] <- signc*Effects[,1]
	} else {
		Effects <- matrix(0,ngenes,n)
		signc <- rep.int(0,ngenes)
		ws <- sqrt(weights)
		for (g in 1:ngenes) {			
			wX <- X*ws[g,]
			wy <- y[g,]*ws[g,]
			if(!is.null(block)) {
				wy <- backsolve(R,wy,transpose=TRUE)
				wX <- backsolve(R,wX,transpose=TRUE)
			}
			qrX <- qr(wX)
			signc[g] <- sign(qrX$qr[p,p])
			Effects[g,] <- qr.qty(qrX,wy)
		}
#		Remove model effects other than the required contrast
		if(p>1) Effects <- Effects[,p:n,drop=FALSE]
#		Preserve sign of estimated effect
		Effects[,1] <- signc*Effects[,1]
	}

#	Apply gene weights
	if(!is.null(gene.weights)) Effects <- sqrt(gene.weights) * Effects

#	Dimension names
	EffectNames <- p:n
	EffectNames[1] <- CoefName
	dimnames(Effects) <- list(geneid,c(EffectNames))

	Effects
}
