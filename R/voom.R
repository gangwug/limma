voom <- function(counts,design=NULL,lib.size=NULL,normalize.method="none",span=0.5,plot=FALSE,save.plot=FALSE,...) 
#	Linear modelling of count data with mean-variance modelling at the observational level.
#	Creates an EList object for entry to lmFit() etc in the limma pipeline.
#	Gordon Smyth and Charity Law
#	Created 22 June 2011.  Last modified 19 April 2019.
{
	out <- list()

#	Check counts
	if(is(counts,"DGEList")) {
		out$genes <- counts$genes
		out$targets <- counts$samples
		if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
		if(is.null(lib.size)) lib.size <- counts$samples$lib.size*counts$samples$norm.factors
		counts <- counts$counts
	} else {
		isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
		if(isExpressionSet) {
			if(length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
			if(length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
			counts <- Biobase::exprs(counts)
		} else {
			counts <- as.matrix(counts)
		}
	}

	n <- nrow(counts)
	if(n < 2L) stop("Need at least two genes to fit a mean-variance trend")

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(counts),1)
		rownames(design) <- colnames(counts)
		colnames(design) <- "GrandMean"
	}

#	Check lib.size
	if(is.null(lib.size)) lib.size <- colSums(counts)

#	Fit linear model to log2-counts-per-million
	y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
	y <- normalizeBetweenArrays(y,method=normalize.method)
	fit <- lmFit(y,design,...)
	if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)

#	If no replication found, set all weight to 1
	NWithReps <- sum(fit$df.residual > 0L)
	if(NWithReps < 2L) {
		if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
		if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
		out$E <- y
		out$weights <- y
		out$weights[] <- 1
		out$design <- design
		if(is.null(out$targets))
			out$targets <- data.frame(lib.size=lib.size)
		else
			out$targets$lib.size <- lib.size
		return(new("EList",out))
	}

#	Fit lowess trend to sqrt-standard-deviations by log-count-size
	sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)
	sy <- sqrt(fit$sigma)
	allzero <- rowSums(counts)==0
	if(any(allzero)) {
		sx <- sx[!allzero]
		sy <- sy[!allzero]
	}
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lines(l,col="red")
	}

#	Make interpolating rule
#	Special treatment of zero counts is now removed;
#	instead zero counts get same variance as smallest gene average.
#	l$x <- c(0.5^0.25, l$x)
#	l$x <- c(log2(0.5), l$x)
#	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
#	var0 <- max(var0,1e-6)
#	l$y <- c(var0, l$y)
	f <- approxfun(l, rule=2, ties=list("ordered",mean))

#	Find individual quarter-root fitted counts
	if(fit$rank < ncol(design)) {
		j <- fit$pivot[1:fit$rank]
		fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
	} else {
		fitted.values <- fit$coef %*% t(fit$design)
	}
	fitted.cpm <- 2^fitted.values
	fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
	fitted.logcount <- log2(fitted.count)

#	Apply trend to individual observations
	w <- 1/f(fitted.logcount)^4
	dim(w) <- dim(fitted.logcount)

#	Output
	out$E <- y
	out$weights <- w
	out$design <- design
	if(is.null(out$targets))
		out$targets <- data.frame(lib.size=lib.size)
	else
		out$targets$lib.size <- lib.size
	if(save.plot) {
		out$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )")
		out$voom.line <- l
	}

	new("EList",out)
}

