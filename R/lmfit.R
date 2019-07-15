#  LINEAR MODELS

lmFit <- function(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation,weights=NULL,method="ls",...)
#	Fit genewise linear models
#	Gordon Smyth
#	30 June 2003.  Last modified 6 Oct 2015.
{
#	Extract components from y
	y <- getEAWP(object)

#	Check design matrix
	if(is.null(design)) design <- y$design
	if(is.null(design))
		design <- matrix(1,ncol(y$exprs),1)
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
		if(nrow(design) != ncol(y$exprs)) stop("row dimension of design doesn't match column dimension of data object")
	}
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")

#	Weights and spacing arguments can be specified in call or stored in y
#	Precedence for these arguments is
#	1. Specified in function call
#	2. Stored in object
#	3. Default values
	if(missing(ndups) && !is.null(y$printer$ndups)) ndups <- y$printer$ndups
	if(missing(spacing) && !is.null(y$printer$spacing)) spacing <- y$printer$spacing
	if(missing(weights) && !is.null(y$weights)) weights <- y$weights

#	Check method
	method <- match.arg(method,c("ls","robust"))

#	If duplicates are present, reduce probe-annotation and Amean to correct length
	if(ndups>1) {
		if(!is.null(y$probes)) y$probes <- uniquegenelist(y$probes,ndups=ndups,spacing=spacing)
		if(!is.null(y$Amean)) y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean),ndups=ndups,spacing=spacing),na.rm=TRUE)
	}

#	Dispatch fitting algorithms
	if(method=="robust")
		fit <- mrlm(y$exprs,design=design,ndups=ndups,spacing=spacing,weights=weights,...)
	else
		if(ndups < 2 && is.null(block))
			fit <- lm.series(y$exprs,design=design,ndups=ndups,spacing=spacing,weights=weights)
		else {
			if(missing(correlation)) stop("the correlation must be set, see duplicateCorrelation")
			fit <- gls.series(y$exprs,design=design,ndups=ndups,spacing=spacing,block=block,correlation=correlation,weights=weights,...)
		}

#	Possible warning on missing coefs
	if(NCOL(fit$coef)>1) {
		n <- rowSums(is.na(fit$coef))
		n <- sum(n>0 & n<NCOL(fit$coef))
		if(n>0) warning("Partial NA coefficients for ",n," probe(s)",call.=FALSE) 
	}

#	Output
	fit$genes <- y$probes
	fit$Amean <- y$Amean
	fit$method <- method
	fit$design <- design
	new("MArrayLM",fit)
}

lm.series <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL)
#	Fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	18 Apr 2002. Revised 26 June 2015.
{
#	Check expression matrix
	M <- as.matrix(M)
	narrays <- ncol(M)

#	Check design
	if(is.null(design))
		design <- matrix(1,narrays,1)
	else
		design <- as.matrix(design)
	nbeta <- ncol(design)
	coef.names <- colnames(design)
	if(is.null(coef.names)) coef.names <- paste("x",1:nbeta,sep="")

#	Check weights
	if(!is.null(weights)) {
		weights <- asMatrixWeights(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}

#	Reform duplicated rows into columns
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}

#	Initialize standard errors
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))

#	Check whether QR-decomposition is constant for all genes
#	If so, fit all genes in one sweep
	NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights,"arrayweights")))
	if(NoProbeWts) {
		if(is.null(weights))
			fit <- lm.fit(design, t(M))
		else {
			fit <- lm.wfit(design, t(M), weights[1,])
			fit$weights <- NULL
		}
		if(fit$df.residual>0) {
			if(is.matrix(fit$effects))
				fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 1):narrays,,drop=FALSE]^2))
			else
				fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 1):narrays]^2))
		} else
			fit$sigma <- rep(NA,ngenes)
		fit$fitted.values <- fit$residuals <- fit$effects <- NULL
		fit$coefficients <- t(fit$coefficients)
		fit$cov.coefficients <- chol2inv(fit$qr$qr,size=fit$qr$rank)
		est <- fit$qr$pivot[1:fit$qr$rank]
		dimnames(fit$cov.coefficients) <- list(coef.names[est],coef.names[est])
		stdev.unscaled[,est] <- matrix(sqrt(diag(fit$cov.coefficients)),ngenes,fit$qr$rank,byrow = TRUE)
		fit$stdev.unscaled <- stdev.unscaled
		fit$df.residual <- rep.int(fit$df.residual,ngenes)
		dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
		fit$pivot <- fit$qr$pivot
		return(fit)
	}

#	Genewise QR-decompositions are required, so iterate through genes
	beta <- stdev.unscaled
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		if(sum(obs) > 0) {
			X <- design[obs,,drop=FALSE]
			y <- y[obs]
			if(is.null(weights))
				out <- lm.fit(X,y)
			else {
				w <- as.vector(weights[i,obs])
				out <- lm.wfit(X,y,w)
			}
			est <- !is.na(out$coef)
			beta[i,] <- out$coef
			stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
			df.residual[i] <- out$df.residual
			if(df.residual[i] > 0) sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
		}
	}

#	Correlation matrix of coefficients
	QR <- qr(design)
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])

	list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,df.residual=df.residual,cov.coefficients=cov.coef,pivot=QR$pivot,rank=QR$rank)
}

mrlm <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL,...)
#	Robustly fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	20 Mar 2002.  Last revised 26 June 2015.
{
	if(!requireNamespace("MASS",quietly=TRUE)) stop("MASS package required but is not available")
	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	coef.names <- colnames(design)
	nbeta <- ncol(design)
	if(!is.null(weights)) {
		weights <- asMatrixWeights(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		X <- design[obs,,drop=FALSE]
		y <- y[obs]
		if(is.null(weights))
			w <- rep(1,length(y))
		else
			w <- as.vector(weights[i,obs])
		if(length(y) > nbeta) {
			out <- MASS::rlm(x=X,y=y,weights=w,...)
			beta[i,] <- coef(out)
			stdev.unscaled[i,] <- sqrt(diag(chol2inv(out$qr$qr)))
			df.residual[i] <- length(y) - out$rank
			if(df.residual[i] > 0) sigma[i] <- out$s
		}
	}
	QR <- qr(design)
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])
	list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,df.residual=df.residual,cov.coefficients=cov.coef,pivot=QR$pivot,rank=QR$rank)
}

gls.series <- function(M,design=NULL,ndups=2,spacing=1,block=NULL,correlation=NULL,weights=NULL,...)
#	Fit linear model for each gene to a series of microarrays.
#	Fit is by generalized least squares allowing for correlation between duplicate spots.
#	Gordon Smyth
#	11 May 2002.  Last revised 29 Dec 2015.
{
#	Check M
	M <- as.matrix(M)
	ngenes <- nrow(M)
	narrays <- ncol(M)

#	Check design
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	if(nrow(design) != narrays) stop("Number of rows of design matrix does not match number of arrays")
	nbeta <- ncol(design)
	coef.names <- colnames(design)

#	Check correlation
	if(is.null(correlation)) correlation <- duplicateCorrelation(M,design=design,ndups=ndups,spacing=spacing,block=block,weights=weights,...)$consensus.correlation
	if(abs(correlation) >= 1) stop("correlation is 1 or -1, so the model is degenerate")

#	Check weights
	if(!is.null(weights)) {
		weights[is.na(weights)] <- 0
		weights <- asMatrixWeights(weights,dim(M))
		M[weights < 1e-15 ] <- NA
		weights[weights < 1e-15] <- NA
	}

#	Unwrap duplicates (if necessary) and make correlation matrix
	if(is.null(block)) {
#		Correlated within-array probes
		if(ndups<2) {
			warning("No duplicates (ndups<2)")
			ndups <- 1
			correlation <- 0
		}
		cormatrix <- diag(rep(correlation,len=narrays),nrow=narrays,ncol=narrays) %x% array(1,c(ndups,ndups))
		if(is.null(spacing)) spacing <- 1
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		colnames(design) <- coef.names
		ngenes <- nrow(M)
		narrays <- ncol(M)
	} else {
#		Correlated samples
		ndups <- spacing <- 1
		block <- as.vector(block)
		if(length(block)!=narrays) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,narrays,nblocks)==matrix(ub,narrays,nblocks,byrow=TRUE)
		cormatrix <- Z%*%(correlation*t(Z))
	}
	diag(cormatrix) <- 1

	stdev.unscaled <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))

#	If weights and missing values are absent, use a fast computation
	NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights,"arrayweights")))
	if(NoProbeWts) {
		V <- cormatrix
		if(!is.null(weights)) {
			wrs <- 1/sqrt(weights[1,])
			V <- wrs * t(wrs * t(V))
		}
		cholV <- chol(V)
		y <- backsolve(cholV,t(M),transpose=TRUE)
		dimnames(y) <- rev(dimnames(M))
		X <- backsolve(cholV,design,transpose=TRUE)
		dimnames(X) <- dimnames(design)
		fit <- lm.fit(X,y)
		if(fit$df.residual>0) {
			if(is.matrix(fit$effects))
				fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank),,drop=FALSE]^2))
			else
				fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
		} else
			fit$sigma <- rep(NA,ngenes)
		fit$fitted.values <- fit$residuals <- fit$effects <- NULL
		fit$coefficients <- t(fit$coefficients)
		fit$cov.coefficients <- chol2inv(fit$qr$qr,size=fit$qr$rank)
		est <- fit$qr$pivot[1:fit$qr$rank]
		dimnames(fit$cov.coefficients) <- list(coef.names[est],coef.names[est])
		stdev.unscaled[,est] <- matrix(sqrt(diag(fit$cov.coefficients)),ngenes,fit$qr$rank,byrow = TRUE)
		fit$stdev.unscaled <- stdev.unscaled
		fit$df.residual <- rep.int(fit$df.residual,ngenes)
		dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
		fit$pivot <- fit$qr$pivot
		fit$ndups <- ndups
		fit$spacing <- spacing
		fit$block <- block
		fit$correlation <- correlation
		return(fit)
	}

#	Weights or missing values are present, to have to iterate over probes
	beta <- stdev.unscaled
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		y <- y[o]
		n <- length(y)
		if(n > 0) {
			X <- design[o,,drop=FALSE]
			V <- cormatrix[o,o]
			if(!is.null(weights)) {
				wrs <- 1/sqrt(drop(weights[i,o]))
				V <- wrs * t(wrs * t(V))
			}
			cholV <- chol(V)
			y <- backsolve(cholV,y,transpose=TRUE)
			if(all(X==0)) {
				df.residual[i] <- n
				sigma[i] <- sqrt( array(1/n,c(1,n)) %*% y^2 )
			} else {
				X <- backsolve(cholV,X,transpose=TRUE)
				out <- lm.fit(X,y)
				est <- !is.na(out$coefficients)
				beta[i,] <- out$coefficients
				stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
				df.residual[i] <- out$df.residual
				if(df.residual[i] > 0)
					sigma[i] <- sqrt( array(1/out$df.residual,c(1,n)) %*% out$residuals^2 )
			}
		}
	}
	cholV <- chol(cormatrix)
	QR <- qr(backsolve(cholV,design,transpose=TRUE))
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])
	list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,df.residual=df.residual,ndups=ndups,spacing=spacing,block=block,correlation=correlation,cov.coefficients=cov.coef,pivot=QR$pivot,rank=QR$rank)
}

is.fullrank <- function(x)
#	Check whether a numeric matrix has full column rank
#	Gordon Smyth
#	18 August 2003.  Last modified 9 March 2004.
{
	x <- as.matrix(x)
	e <- eigen(crossprod(x),symmetric=TRUE,only.values=TRUE)$values
	e[1] > 0 && abs(e[length(e)]/e[1]) > 1e-13
}

nonEstimable <- function(x)
#	Check whether a numeric matrix has full column rank
#	If not, return names of redundant columns
#	Gordon Smyth
#	10 August 2004
{
	x <- as.matrix(x)
	p <- ncol(x)
	QR <- qr(x)
	if(QR$rank < p) {
		n <- colnames(x)
		if(is.null(n)) n <- as.character(1:p)
		notest <- n[QR$pivot[(QR$rank+1):p]]
		blank <- notest==""
		if(any(blank)) notest[blank] <- as.character(((QR$rank+1):p)[blank])
		return(notest)
	} else {
		return(NULL)
	}
}

fitted.MArrayLM <- function(object,...)
#	Fitted values from MArray linear model fit
#	Gordon Smyth
#	29 November 2005.  Last modified 12 Feb 2019.
{
	if(!is.null(object$contrasts)) stop("Object contains contrasts rather than coefficients, so fitted values cannot be computed.")
	object$coefficients %*% t(object$design)
}

residuals.MArrayLM <- function(object,y,...)
#	Residuals from MArray linear model fit
#	Gordon Smyth
#	29 November 2005
{
	as.matrix(y) - fitted(object)
}

getEAWP <- function(object)
#	Given any microarray data object, extract basic information needed
#	for linear modelling.
#	Gordon Smyth
#	9 March 2008. Last modified 18 July 2016.
{
	y <- list()
	
	if(is(object,"list")) {
#		Method for MAList (classed or unclassed) or EList objects
		if(is(object,"EList")) {
			y$exprs <- as.matrix(object$E)
			y$Amean <- rowMeans(y$exprs,na.rm=TRUE)
		} else {
			if(is(object,"EListRaw")) stop("EListRaw object: please normalize first")
			if(is(object,"RGList")) stop("RGList object: please normalize first")
			y$printer <- object$printer
			if(is.null(object$M)) stop("data object isn't of a recognized data class")
			y$exprs <- as.matrix(object$M)
			if(!is.null(object$A)) y$Amean <- rowMeans(as.matrix(object$A),na.rm=TRUE)
		}
		y$weights <- object$weights
		y$probes <- object$genes
		y$design <- object$design
	} else {
	if(is(object,"ExpressionSet")) {
		if(!requireNamespace("Biobase",quietly=TRUE)) stop("Biobase package required but is not available")
		y$exprs <- Biobase::exprs(object)
		if(length(object@featureData@data)) y$probes <- object@featureData@data
		y$Amean <- rowMeans(y$exprs,na.rm=TRUE)
		if("weights" %in% Biobase::assayDataElementNames(object)) y$weights <- Biobase::assayDataElement(object,"weights")
	} else {
	if(is(object,"PLMset")) {
		y$exprs <- object@chip.coefs
		if(length(y$exprs)==0) stop("chip.coefs has length zero")
		if(length(object@se.chip.coefs)) y$weights <- 1/pmax(object@se.chip.coefs,1e-5)^2
		y$Amean <- rowMeans(y$exprs,na.rm=TRUE)
	} else {
	if(is(object,"marrayNorm")) {
		y$exprs <- object@maM
		if(length(object@maW)) y$weights <- object@maW
		if(length(object@maGnames@maInfo)) {
			y$probes <- object@maGnames@maInfo
			attr(y$probes, "Notes") <- object@maGnames@maNotes
		}
		if(length(object@maA)) y$Amean <- rowMeans(object@maA,na.rm=TRUE)
	} else {
	if(is(object,"eSet")) {
		if(!requireNamespace("Biobase",quietly=TRUE)) stop("Biobase package required but is not available")
		y$exprs <- Biobase::assayDataElement(object,"exprs")
		if(length(object@featureData@data)) y$probes <- object@featureData@data
		y$Amean <- rowMeans(y$exprs,na.rm=TRUE)
		if("weights" %in% Biobase::assayDataElementNames(object)) y$weights <- Biobase::assayDataElement(object,"weights")
	} else {
#		Default method for matrices, data.frames, vsn objects etc.
		if(is.vector(object))
			y$exprs <- matrix(object,nrow=1)
		else
			y$exprs <- as.matrix(object)
		y$Amean <- rowMeans(y$exprs,na.rm=TRUE)
	}}}}}

#	Check expression values are numeric
	if(mode(y$exprs) != "numeric") stop("Data object doesn't contain numeric expression values")

#	Get rownames from probes?
	if(is.null(rownames(y$exprs)) && !is.null(row.names(y$probes))) rownames(y$exprs) <- row.names(y$probes)

#	Check rownames are unique
#	rn <- rownames(y$exprs)
#	if(is.null(rn))
#		rownames(y$exprs) <- 1:nrow(y$exprs)
#	else
#		if(anyDuplicated(rn)>0) {
#			rownames(y$exprs) <- 1:nrow(y$exprs)
#			if(is.null(y$probes))
#				y$probes <- data.frame(ID=rn,stringsAsFactors=FALSE)
#			else
#				if("ID" %in% names(y$probes))
#					y$probes$ID0 <- rn
#				else
#					y$probes$ID <- rn
#		}

	y
}

