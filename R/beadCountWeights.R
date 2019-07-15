beadCountWeights <- function(y, x, design=NULL, bead.stdev=NULL, bead.stderr=NULL, nbeads=NULL, array.cv=TRUE, scale=FALSE)
#	Compute weights for BeadChips based on bead-level counts and standard deviations
#	Charity Law and Gordon Smyth
#	Created 4 August 2010. Last modified 4 Oct 2017.
{
#	Matrix of normalized expression values
	E <- as.matrix(y)
	P <- nrow(E)
	A <- ncol(E)

#	Matrix of non-normalized expression values on unlog scale
	E.raw <- as.matrix(x)
	if(nrow(E.raw) != P) stop("dimensions don't match")
	if(ncol(E.raw) != A) stop("dimensions don't match")

#	Check design matrix
	if(is.null(design)) design <- y$design
	if(is.null(design))
		design <- matrix(1,A,1)
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
		if(nrow(design) != A) stop("row dimension of design doesn't match column dimension of data object")
	}

#	Number of beads per probe
	if(is.null(nbeads)) {
		nbeads <- y$other$Avg_NBEADS
		if(is.null(nbeads)) stop("nbeads not provided and Avg_NBEADS not found in data object")
	}	
	if(nrow(nbeads) != P) stop("dimensions don't match")
	if(ncol(nbeads) != A) stop("dimensions don't match")

#	Standard deviation of beads per probe
	if(is.null(bead.stdev)) {
		if(is.null(bead.stderr)) {
			if(is.null(y$other$BEAD_STDEV)) {
				if(is.null(y$other$BEAD_STDERR)) stop("Neither bead.stdev nor bead.stderr were provided (and neither could be found in the data object y). At least one is required.")
				bead.stdev <- y$other$BEAD_STDERR*sqrt(nbeads)
			} else {
				bead.stdev <- y$other$BEAD_STDEV
			}
		} else {
			bead.stdev <- bead.stderr*sqrt(nbeads)
		}
	}
	if(nrow(bead.stdev) != P) stop("dimensions don't match")
	if(ncol(bead.stdev) != A) stop("dimensions don't match")

#	Coefficient of variation of bead-level observations
#	Array-specific or constant. 
	cv <- bead.stdev/E.raw
	cv.array <- apply(sqrt(cv), 2, function(k) mean(k, trim=0.125, na.rm=TRUE))^2
	cv.constant <- mean(sqrt(cv), trim=0.125, na.rm=TRUE)^2
	if(array.cv) {cv <- rep(cv.array, each=nrow(y))} 
	else cv <- cv.constant
	
#	Predicted variance of normalized probe-level values
	tv <- log(cv^2/nbeads+1)/log(2)^2

# 	Squared-residuals to calculate biological variance
	qrX <- qr(design)
	res <- qr.resid(qrX, t(E))
	r2 <- (t(res))^2

# 	Leverages
	Q <- qr.Q(qrX) 
	h <- rowSums(Q^2)
	h <- matrix(rep(h, each=P), nrow=P, ncol=A)

	bv <- .ilmn.biological.variance(var.from.beads=tv, sq.residuals=r2, leverage=h)

	out <- list()
	out$var.biological <- bv
	out$var.technical <- tv
	out$cv.array <- cv.array
	out$cv.constant <- cv.constant
	out$weights <- 1/(tv + bv)
	out$weights <- out$weights/rowMeans(out$weights)
	if(scale) out$weights <- out$weights*rowMeans(1/tv)
	out
}


.ilmn.biological.variance <- function(var.from.beads, sq.residuals, leverage)
#	Find the component of the between-array variance
#	not explained by bead-level variability
#	Gordon Smyth and Charity Law
#	28 July 2010. Modified 16 July 2012.
{
	if(any(var.from.beads < 0)) stop("negative variances not allowed")
	if(any(sq.residuals < 0)) stop("negative variances not allowed")
	tv <- as.matrix(var.from.beads)
	r2 <- as.matrix(sq.residuals)
	h <- as.matrix(leverage)
	P <- nrow(tv)
	A <- ncol(tv)
	if(nrow(r2) != P) stop("dimensions don't match")
	if(ncol(r2) != A) stop("dimensions don't match")
	if(nrow(h) != P) stop("dimensions don't match")
	if(ncol(h) != A) stop("dimensions don't match")
	
	# Newton's iteration
	F <- rowMeans(r2/(2*tv^2)-(1-h)/(2*tv))	
	bv <- rep(0,length=P)
	i <- (F > 0)
	iter <- 0
	while(any(i)) {
		iter <- iter+1
		if(iter > 200) {
			warning("More than 200 iterations required")
			return(bv)
		}
		Fdash <- - rowMeans(r2[i,,drop=FALSE]/(bv[i]+tv[i,,drop=FALSE])^3 - (1-h[i,,drop=FALSE])/(2*(bv[i]+tv[i,,drop=FALSE])^2))
		step <- - F[i]/Fdash
		bv[i] <- bv[i] + step
		F[i] <- rowMeans(r2[i,,drop=FALSE]/(2*(bv[i]+tv[i,,drop=FALSE])^2) - (1-h[i,,drop=FALSE])/(2*(bv[i]+tv[i,,drop=FALSE])))
		i[i] <- (step > 1e-5)
		print(sum(i))
	}
	bv
}


