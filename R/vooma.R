vooma <- function(y,design=NULL,correlation,block=NULL,plot=FALSE,span=NULL)
# Linear modelling of microarray data with mean-variance modelling at the observational level.
# Creates an EList object for entry to lmFit() etc in the limma pipeline.
# Gordon Smyth and Charity Law
# Created 31 July 2012.  Last modified 17 May 2019.
{
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	narrays <- ncol(y)
	ngenes <- nrow(y)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		design <- matrix(1,narrays,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "GrandMean"
	}

#	Fit linear model
	if(is.null(block)) {
		fit <- lm.fit(design,t(y$E))
		mu <- fit$fitted.values
	} else {
		block <- as.vector(block)
		if(length(block)!=narrays) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,narrays,nblocks)==matrix(ub,narrays,nblocks,byrow=TRUE)
		cormatrix <- Z%*%(correlation*t(Z))
		diag(cormatrix) <- 1
		cholV <- chol(cormatrix)
		z <- backsolve(cholV,t(y$E),transpose=TRUE)
		X <- backsolve(cholV,design,transpose=TRUE)
		fit <- lm.fit(X,z)
		mu <- crossprod(cholV,fit$fitted.values)
	}
	s2 <- colMeans(fit$effects[-(1:fit$rank),,drop=FALSE]^2)

#	Fit lowess trend to sqrt-standard-deviations by ave log intensity
	sx <- rowMeans(y$E)
	sy <- sqrt(sqrt(s2))
	if(is.null(span)) if(ngenes<=10) span <- 1 else span <- 0.3+0.7*(10/ngenes)^0.5
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="Average log2 expression",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lines(l,col="red")
	}

#	Make interpolating rule
	f <- approxfun(l, rule=2, ties=list("ordered",mean))

#	Apply trend to individual observations
	w <- 1/f(t(mu))^4
	dim(w) <- dim(y)
	colnames(w) <- colnames(y)
	rownames(w) <- rownames(y)

#	Output
	y$meanvar.trend <- list(x=sx,y=sy)
	y$weights <- w
	y$design <- design
	y$span <- span
	y
}

voomaByGroup <- function(y,group,design=NULL,correlation,block=NULL,plot=FALSE,span=NULL,col=NULL,lwd=1,alpha=0.5,pch=16,cex=0.3,legend="topright")
#	Vooma by group
#	Linear modelling of microarray data with mean-variance modelling at the observational level by fitting group-specific trends.
#	Creates an EList object for entry to lmFit() etc in the limma pipeline.
#	Charity Law and Gordon Smyth
#	Created 13 Feb 2013.  Modified 8 Sept 2014.
{
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	ngenes <- nrow(y)
	narrays <- ncol(y)

#	Check group
	group <- as.factor(group)
	intgroup <- as.integer(group)
	levgroup <- levels(group)
	ngroups <- length(levgroup)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) design <- model.matrix(~group)

#	Check color
	if(is.null(col))
		if(ngroups==2L)
			col <- c("red","blue")
		else
			col <- 1L+1L:ngroups
	col <- rep_len(col,ngroups)

	w <- y$E
	sx <- sy <- matrix(0, nrow=nrow(y), ncol=nlevels(group))
	colnames(sx) <- levgroup
	rownames(sx) <- rownames(y)
	for (lev in 1L:ngroups) {
		i <- intgroup==lev
		yi <- y[,i]
		designi <- design[i,,drop=FALSE]
		voomi <- vooma(y=yi,design=designi,correlation=correlation, block=block[i], plot=FALSE, span=span)
		w[,i] <- voomi$weights
		sx[,lev] <- voomi$meanvar.trend$x
		sy[,lev] <- voomi$meanvar.trend$y	
	}
	span <- voomi$span

# 	Voom plot	
	if(plot) {
		RGB <- col2rgb(col)/255
		plot(sx,sy,xlab="Average log2 expression",ylab="Sqrt( standard deviation )",main="voom: Mean-variance trend",type="n")
		for (lev in 1:nlevels(group)) {
			coli.transparent <- rgb(RGB[1,lev],RGB[2,lev],RGB[3,lev],alpha=alpha)
			points(sx[,lev],sy[,lev],pch=pch,cex=cex,col=coli.transparent)
			l <- lowess(sx[,lev],sy[,lev],f=span)
			lines(l,col=col[lev],lwd=lwd)
		}
		if(is.character(legend)) legend(legend, levels(group), col=col, lty=1, lwd=lwd)
	}

# 	Output	
	y$meanvar.trend <- list(x=sx,y=sy)
	y$weights <- w
	y$design <- design
	y$span <- span
	y
}
