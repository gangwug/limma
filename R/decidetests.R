#  DECIDETESTS.R

setClass("TestResults",representation("matrix"))

summary.TestResults <- function(object,...)
#	Gordon Smyth
#	Created 26 Feb 2004.  Last modified 6 Jan 2017.
{
	Levels <- attr(object,"levels")
	if(is.null(Levels)) Levels <- c(-1L,0L,1L)
	nlevels <- length(Levels)
	tab <- matrix(0L,nlevels,ncol(object))
	Labels <- attr(object,"labels")
	if(is.null(Labels)) Labels <- as.character(Levels)
	dimnames(tab) <- list(Labels,colnames(object))
	for (i in 1:nlevels) tab[i,] <- colSums(object==Levels[i],na.rm=TRUE)
	class(tab) <- "table"
	tab
}

setMethod("show","TestResults",function(object) {
	cat("TestResults matrix\n")
	printHead(object@.Data)
})

levels.TestResults <- function(x) attr(x,"levels")

labels.TestResults <- function(object,...) attr(object,"labels")


decideTests <- function(object,...) UseMethod("decideTests")

decideTests.default <- function(object,method="separate",adjust.method="BH",p.value=0.05,lfc=0,coefficients=NULL,cor.matrix=NULL,tstat=NULL,df=Inf,genewise.p.value=NULL,...)
#	Accept or reject hypothesis tests across genes and contrasts
#	from a matrix of p-values
#	Gordon Smyth
#	17 Aug 2004. Last modified 13 Dec 2017.
{
	method <- match.arg(method,c("separate","global","hierarchical","nestedF"))
	if(method=="nestedF") stop("nestedF adjust method requires an MArrayLM object",call.=FALSE)

	adjust.method <- match.arg(adjust.method,c("none","bonferroni","holm","BH","fdr","BY"))
	if(adjust.method=="fdr") adjust.method <- "BH"

	p <- as.matrix(object)
	if(any(p>1) || any(p<0)) stop("object doesn't appear to be a matrix of p-values")

	switch(method,

	  separate={
		for (i in 1:ncol(p)) p[,i] <- p.adjust(p[,i],method=adjust.method)

	},global={
		p[] <- p.adjust(p[],method=adjust.method)

	},hierarchical={
		if(is.null(genewise.p.value)) {
#			Apply Simes' method by rows to get genewise p-values
			genewise.p.value <- rep_len(1,nrow(p))
			ngenes <- nrow(p)
			ncontrasts <- ncol(p)
			Simes.multiplier <- ncontrasts/(1:ncontrasts)
			for (g in 1:ngenes) {
				op <- sort(p[g,],na.last=TRUE)
				genewise.p.value[g] <- min(op*Simes.multiplier,na.rm=TRUE)
			}
		}
#		Adjust genewise p-values
		DEgene <- p.adjust(genewise.p.value,method=adjust.method) <= p.value
#		Adjust row-wise p-values
		p[!DEgene,] <- 1
		gDE <- which(DEgene)
		for (g in gDE) p[g,] <- p.adjust(p[g,],method=adjust.method)
#		Adjust p-value cutoff for number of DE genes
		nDE <- length(gDE)
		a <- switch(adjust.method,
			none=1,
			bonferroni=1/ngenes,
			holm=1/(ngenes-nDE+1),
			BH=nDE/ngenes,
			BY=nDE/ngenes/sum(1/(1:ngenes))
		)
		p.value <- a*p.value
	},nestedF={
		stop("nestedF adjust method requires an MArrayLM object",call.=FALSE)
	})

	isDE <- array(0L,dim(p),dimnames=dimnames(p))
	isDE[p <= p.value] <- 1L
	if(is.null(coefficients)) coefficients <- tstat
	if(is.null(coefficients)) {
		attr(isDE,"levels") <- c(0L,1L)
		attr(isDE,"labels") <- c("NotSig","Sig")
	} else {
		attr(isDE,"levels") <- c(-1L,0L,1L)
		attr(isDE,"labels") <- c("Down","NotSig","Up")
		coefficients <- as.matrix(coefficients)
		if( !all(dim(coefficients)==dim(p)) ) stop("dim(object) disagrees with dim(coefficients)")
		i <- coefficients<0
		isDE[i] <- -isDE[i]
		if(lfc>0) isDE[ abs(coefficients)<lfc ] <- 0L
	}

	new("TestResults",isDE)
}

decideTests.MArrayLM <- function(object,method="separate",adjust.method="BH",p.value=0.05,lfc=0,...)
#	Accept or reject hypothesis tests across genes and contrasts
#	Gordon Smyth
#	17 Aug 2004. Last modified 13 Dec 2017.
{
	if(is.null(object$p.value)) object <- eBayes(object)
	method <- match.arg(method,c("separate","global","hierarchical","nestedF"))
	adjust.method <- match.arg(adjust.method,c("none","bonferroni","holm","BH","fdr","BY"))
	if(adjust.method=="fdr") adjust.method <- "BH"
	switch(method,separate={
		p <- as.matrix(object$p.value)
		for (j in 1:ncol(p)) {
			o <- !is.na(p[,j])
			p[o,j] <- p.adjust(p[o,j],method=adjust.method)
		}
		s <- sign(as.matrix(object$coefficients))
		results <- new("TestResults",s*(p<p.value))
	},global={
		p <- as.matrix(object$p.value)
		o <- !is.na(p)
		p[o] <- p.adjust(p[o],method=adjust.method)
		s <- sign(as.matrix(object$coefficients))
		results <- new("TestResults",s*(p<p.value))
	},hierarchical={
		if(any(is.na(object$F.p.value))) stop("Can't handle NA p-values yet")
		sel <- p.adjust(object$F.p.value,method=adjust.method) < p.value
		i <- sum(sel,na.rm=TRUE)
		n <- sum(!is.na(sel))
		a <- switch(adjust.method,
			none=1,
			bonferroni=1/n,
			holm=1/(n-i+1),
			BH=i/n,
			BY=i/n/sum(1/(1:n))
		)
		results <- new("TestResults",array(0,dim(object$t)))
		dimnames(results) <- dimnames(object$coefficients)
		if(any(sel)) results[sel,] <- .classifyTestsP(object[sel,],p.value=p.value*a,method=adjust.method)
	},nestedF={
		if(any(is.na(object$F.p.value))) stop("nestedF method can't handle NA p-values",call.=FALSE)
		sel <- p.adjust(object$F.p.value,method=adjust.method) < p.value
		i <- sum(sel,na.rm=TRUE)
		n <- sum(!is.na(sel))
		a <- switch(adjust.method,
			none=1,
			bonferroni=1/n,
			holm=1/(n-i+1),
			BH=i/n,
			BY=i/n/sum(1/(1:n))
		)
		results <- new("TestResults",array(0,dim(object$t)))
		dimnames(results) <- dimnames(object$coefficients)
		if(any(sel)) results[sel,] <- classifyTestsF(object[sel,],p.value=p.value*a)
	})
	if(lfc>0) {
		if(is.null(object$coefficients))
			warning("lfc ignored because coefficients not found")
		else
			results@.Data <- results@.Data * (abs(object$coefficients)>lfc)
	}
	attr(results,"levels") <- c(-1L,0L,1L)
	attr(results,"labels") <- c("Down","NotSig","Up")
	results
}

classifyTestsF <- function(object,cor.matrix=NULL,df=Inf,p.value=0.01,fstat.only=FALSE) {
#	Use F-tests to classify vectors of t-test statistics into outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 6 June 2009.

#	Method intended for MArrayLM objects but accept unclassed lists as well
	if(is.list(object)) {
		if(is.null(object$t)) stop("tstat cannot be extracted from object")
		if(is.null(cor.matrix) && !is.null(object$cov.coefficients)) cor.matrix <- cov2cor(object$cov.coefficients)
		if(missing(df) && !is.null(object$df.prior) && !is.null(object$df.residual)) df <- object$df.prior+object$df.residual
		tstat <- as.matrix(object$t)
	} else {
		tstat <- as.matrix(object)
	}
	ngenes <- nrow(tstat)
	ntests <- ncol(tstat)
	if(ntests == 1) {
		if(fstat.only) {
			fstat <- tstat^2
			attr(fstat,"df1") <- 1
			attr(fstat,"df2") <- df
			return(fstat)
		} else {
			p <- 2 * pt(abs(tstat), df, lower.tail=FALSE)
			return(new("TestResults", sign(tstat) * (p < p.value) ))
		}
	}

#	cor.matrix is estimated correlation matrix of the coefficients
#	and also the estimated covariance matrix of the t-statistics
	if(is.null(cor.matrix)) {
		r <- ntests
		Q <- diag(r)/sqrt(r)
	} else {
		E <- eigen(cor.matrix,symmetric=TRUE)
		r <- sum(E$values/E$values[1] > 1e-8)
		Q <- .matvec( E$vectors[,1:r], 1/sqrt(E$values[1:r]))/sqrt(r)
	}

#	Return overall moderated F-statistic only
	if(fstat.only) {
		fstat <- drop( (tstat %*% Q)^2 %*% array(1,c(r,1)) )
		attr(fstat,"df1") <- r
		attr(fstat,"df2") <- df
		return(fstat)
	}

#	Return TestResults matrix
	qF <- qf(p.value, r, df, lower.tail=FALSE)
	if(length(qF)==1) qF <- rep(qF,ngenes) 
	result <- matrix(0,ngenes,ntests,dimnames=dimnames(tstat))
	for (i in 1:ngenes) {
		x <- tstat[i,]
		if(any(is.na(x)))
			result[i,] <- NA
		else
			if( crossprod(crossprod(Q,x)) > qF[i] ) {
				ord <- order(abs(x),decreasing=TRUE)
				result[i,ord[1]] <- sign(x[ord[1]])
				for (j in 2:ntests) {
					bigger <- ord[1:(j-1)]
					x[bigger] <- sign(x[bigger]) * abs(x[ord[j]])
					if( crossprod(crossprod(Q,x)) > qF[i] )
						result[i,ord[j]] <- sign(x[ord[j]])
					else
						break
				}
			}
	}
	new("TestResults",result)
}

#FStat <- function(object,cor.matrix=NULL)
##	Compute overall F-tests given a matrix of t-statistics
##	Gordon Smyth
##	24 February 2004.  Last modified 21 July 2004.
#{
#	m <- as.list(match.call())
#	m[[1]] <- as.name("classifyTestsF")
#	m$fstat.only <- TRUE
#	eval(as.call(m))
#}

.classifyTestsP <- function(object,df=Inf,p.value=0.05,method="holm") {
#	TestResults by rows for a matrix t-statistics using adjusted p-values
#	Gordon Smyth
#	12 July 2003.  Last modified 23 March 2004.

#	Method intended for MArrayLM objects but accept unclassed lists as well
	if(is.list(object)) {
		if(is.null(object$t)) stop("tstat cannot be extracted from object")
		tstat <- object$t
		if(!is.null(object$df.residual)) df <- object$df.residual
		if(!is.null(object$df.prior)) df <- df+object$df.prior
	} else {
		tstat <- object
	}
	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	ngenes <- nrow(tstat)
	P <- 2*pt(-abs(tstat),df=df)
	result <- tstat
	for (i in 1:ngenes) {
		P[i,] <- p.adjust(P[i,],method=method)
		result[i,] <- sign(tstat[i,])*(P[i,]<p.value)
	}
	new("TestResults",result)
}
