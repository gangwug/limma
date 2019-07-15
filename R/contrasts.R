#  CONTRASTS

contrasts.fit <- function(fit,contrasts=NULL,coefficients=NULL)
#	Convert coefficients and std deviations in fit object to reflect contrasts of interest
#	Note: does not completely take probe-wise weights into account
#	because this would require refitting the linear model for each probe
#	Gordon Smyth
#	Created 13 Oct 2002.  Last modified 14 Feb 2019.
{
#	Check number of arguments
	if(is.null(contrasts) == is.null(coefficients)) stop("Must specify exactly one of contrasts or coefficients")

#	Remove test statistics in case eBayes() has previously been run on the fit object
	fit$t <- NULL
	fit$p.value <- NULL
	fit$lods <- NULL
	fit$F <- NULL
	fit$F.p.value <- NULL

#	Number of coefficients in fit
	ncoef <- NCOL(fit$coefficients)

#	Check contrasts. If coefficients are specified, convert into contrast matrix.
	if(!is.null(contrasts)) {
		contrasts <- as.matrix(contrasts)
		rn <- rownames(contrasts)
		cn <- colnames(fit$coefficients)
		if(!is.null(rn) && !is.null(cn) && any(rn != cn)) warning("row names of contrasts don't match col names of coefficients")
	} else {
		ncont <- length(coefficients)
		contrasts <- diag(ncoef)
		rownames(contrasts) <- colnames(contrasts) <- colnames(fit$coefficients)
		contrasts <- contrasts[,coefficients,drop=FALSE]
	}
	if(nrow(contrasts)!=ncoef) stop("Number of rows of contrast matrix must match number of coefficients in fit")
	if(anyNA(contrasts)) stop("NAs not allowed in contrasts")

	fit$contrasts <- contrasts
	if(is.null(fit$cov.coefficients)) {
		warning("no coef correlation matrix found in fit - assuming orthogonal")
		cormatrix <- diag(ncoef)
	} else
		cormatrix <- cov2cor(fit$cov.coefficients)

#	If design matrix was singular, reduce to estimable coefficients
	r <- nrow(cormatrix)
	if(r < ncoef) {
		if(is.null(fit$pivot)) stop("cor.coef not full rank but pivot column not found in fit")
		est <- fit$pivot[1:r]
		if(any(contrasts[-est,]!=0)) stop("trying to take contrast of non-estimable coefficient")
		contrasts <- contrasts[est,,drop=FALSE]
		fit$coefficients <- fit$coefficients[,est,drop=FALSE]
		fit$stdev.unscaled <- fit$stdev.unscaled[,est,drop=FALSE]
		ncoef <- r
	}
#	fit$coefficients <- fit$coefficients %*% contrasts
	fit$coefficients <- .zeroDominantMatrixMult(fit$coefficients,contrasts)

#	Test whether design was orthogonal
	if(length(cormatrix) < 2) {
		orthog <- TRUE
	} else {
		orthog <- all(abs(cormatrix[lower.tri(cormatrix)]) < 1e-14)
	}

#	Contrast correlation matrix
	R <- chol(fit$cov.coefficients)
	fit$cov.coefficients <- crossprod(R %*% contrasts)
	fit$pivot <- NULL

	if(orthog)
		fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
	else {
		R <- chol(cormatrix)
		ngenes <- NROW(fit$stdev.unscaled)
		ncont <- NCOL(contrasts)
		U <- matrix(1,ngenes,ncont,dimnames=list(rownames(fit$stdev.unscaled),colnames(contrasts)))
		o <- array(1,c(1,ncoef))
		for (i in 1:ngenes) {
			RUC <- R %*% .vecmat(fit$stdev.unscaled[i,],contrasts)
			U[i,] <- sqrt(o %*% RUC^2)
		}
		fit$stdev.unscaled <- U
	}
	fit
}

.zeroDominantMatrixMult <- function(A,B)
#	Computes A %*% B, except that a zero in B will always produce
#	zero even when multiplied by an NA in A
#	Gordon Smyth
#	Created 16 Feb 2018.
{
	HasZero <- (rowSums(B==0) > 0L)
	if(any(HasZero)) {
		if(mean(HasZero) > 0.4) {
#			If the matrix is big, it's much quicker to check the whole matrix than to subset it
			HasNA <- anyNA(A)
		} else {
			HasNA <- anyNA(A[,HasZero])
		}
	} else {
		HasNA <- FALSE
	}
	if(HasZero && HasNA) {
		D <- matrix(0,nrow(A),ncol(B))
		for (j in 1:ncol(B)) {
			z <- B[,j]==0
			if(any(z))
				D[,j] <- A[,!z,drop=FALSE] %*% B[!z,j,drop=FALSE]
			else
				D[,j] <- A %*% B[,j]
		}
		dimnames(D) <- list(rownames(A),colnames(B))
	} else {
		D <- A %*% B
	}
	D
}
