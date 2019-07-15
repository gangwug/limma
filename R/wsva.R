wsva <- function(y, design, n.sv=1L, weight.by.sd=FALSE, plot=FALSE, ...)
#	Weighted surrogate variable analysis
#	Yifang Hu and Gordon Smyth
#	Created 26 Nov 2015.  Last modified 2 Mar 2017.
{
	y <- as.matrix(y)
	ngenes <- nrow(y)
	narrays <- ncol(y)
	p <- ncol(design)
	d <- narrays-p
	n.sv <- max(n.sv,1L)
	n.sv <- min(n.sv, d)
	if(n.sv <= 0L) stop("No residual df")

	if(weight.by.sd) {
		if(plot) message("Plot not available with weight.by.sd=TRUE")
		for(i in 1L:n.sv) {
			Effects <- .lmEffects(y, design, ...)[,-1L]
			s <- sqrt(rowMeans(Effects^2))
			Effects <- s * Effects
			u <- drop(svd(Effects,nu=1L,nv=0L)$u)
			u <- u*s
			sv <- colSums(u*y)
			design <- cbind(design, sv)
		}
		SV <- t(design[,-(1:p),drop=FALSE])
	} else {
		Effects <- .lmEffects(y, design, ...)[,-1L]
		SVD <- svd(Effects,nu=n.sv,nv=0L)
		SV <- crossprod(SVD$u,y)
		if(plot) {
			lambda <- SVD$d^2
			lambda <- lambda/sum(lambda)
			plot(lambda,xlab="Surrogate variable number",ylab="Proportion variance explained")
		}
	}

	A <- rowMeans(SV^2)
	SV <- t( SV / sqrt(A) )
	colnames(SV) <- paste0("SV",1L:n.sv)
	SV
}
