##  CAMERA.R

interGeneCorrelation <- function(y, design)
#	Estimate variance-inflation factor for means of correlated genes
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 11 Feb 2012
{
	m <- nrow(y)
	qrdesign <- qr(design)
	y <- qr.qty(qrdesign, t(y))[-(1:qrdesign$rank),]
#	Gives same result as the following
#	ny <- t(y) / sqrt(colSums(y^2))
#	cormatrix <- tcrossprod(ny)
#	correlation <- mean(cormatrix[lower.tri(cormatrix)])
#	1+correlation*(n-1)
	y <- t(y) / sqrt(colMeans(y^2))
	vif <- m * mean(colMeans(y)^2)
	correlation <- (vif-1)/(m-1)
	list(vif=vif,correlation=correlation)
}

camera <- function(y,...) UseMethod("camera")

camera.default <- function(y,index,design=NULL,contrast=ncol(design),weights=NULL,use.ranks=FALSE,allow.neg.cor=FALSE,inter.gene.cor=0.01,trend.var=FALSE,sort=TRUE,...)
#	Competitive gene set test allowing for correlation between genes
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 30 July 2017.
{
#	Issue warning if extra arguments found
	dots <- names(list(...))
	if(length(dots)) warning("Extra arguments disregarded: ",sQuote(dots))

#	Extract components from y
	y <- getEAWP(y)
	G <- nrow(y$exprs)
	n <- ncol(y$exprs)
	ID <- rownames(y$exprs)
	if(G<3) stop("Two few genes in dataset: need at least 3")

#	Check index
	if(!is.list(index)) index <- list(set1=index)
	nsets <- length(index)
	if(nsets==0L) stop("index is empty")

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design))
		stop("design matrix not specified")
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	p <- ncol(design)
	df.residual <- n-p
	if(df.residual < 1) stop("No residual df: cannot compute t-tests")

#	Check weights
	if(is.null(weights)) weights <- y$weights

#	Check inter.gene.cor
	fixed.cor <- !(is.na(inter.gene.cor) || is.null(inter.gene.cor))

#	Set df for camera tests
	if(fixed.cor) {
		if(use.ranks)
			df.camera <- Inf
		else
			df.camera <- G-2
	} else {
		df.camera <- min(df.residual,G-2)
	}

#	Reduce to numeric expression matrix
	y <- y$exprs

#	Check weights
	if(!is.null(weights)) {
		if(any(weights<=0)) stop("weights must be positive")
		if(length(weights)==n) {
			sw <- sqrt(weights)
			y <- t(t(y)*sw)
			design <- design*sw
			weights <- NULL
		}
	}
	if(!is.null(weights)) {
		if(length(weights)==G) weights <- matrix(weights,G,n)
		weights <- as.matrix(weights)
		if(any( dim(weights) != dim(y) )) stop("weights not conformal with y")
	}

#	Reform design matrix so that contrast of interest is last column
	if(is.character(contrast)) {
		contrast <- which(contrast==colnames(design))
		if(length(contrast)==0) stop("coef ",contrast," not found")
	}
	if(length(contrast)==1) {
		j <- c((1:p)[-contrast], contrast)
		if(contrast<p) design <- design[,j]
	} else {
		QR <- qr(contrast)
		design <- t(qr.qty(QR,t(design)))
		if(sign(QR$qr[1,1]<0)) design[,1] <- -design[,1]
		design <- design[,c(2:p,1)]
	}

#	Compute effects matrix
	if(is.null(weights)) {
		QR <- qr(design)
		if(QR$rank<p) stop("design matrix is not of full rank")
		effects <- qr.qty(QR,t(y))
		unscaledt <- effects[p,]
		if(QR$qr[p,p]<0) unscaledt <- -unscaledt
	} else {
		effects <- matrix(0,n,G)
		colnames(effects) <- ID
		unscaledt <- rep.int(0,G)
		names(unscaledt) <- ID
		sw <- sqrt(weights)
		yw <- y*sw
		for (g in 1:G) {
			xw <- design*sw[g,]
			QR <- qr(xw)
			if(QR$rank<p) stop("weighted design matrix not of full rank for gene ",g)
			effects[,g] <- qr.qty(QR,yw[g,])
			unscaledt[g] <- effects[p,g]
			if(QR$qr[p,p]<0) unscaledt[g] <- -unscaledt[g]
		}
	}

#	Standardized residuals
	U <- effects[-(1:p),,drop=FALSE]
	sigma2 <- colMeans(U^2)
	U <- t(U) / sqrt(pmax(sigma2,1e-8))

#	Moderated t
	if(trend.var) A <- rowMeans(y) else A <- NULL
	sv <- squeezeVar(sigma2,df=df.residual,covariate=A)
	modt <- unscaledt / sqrt(sv$var.post)
	if(use.ranks)
		Stat <- modt
	else {
		df.total <- min(df.residual+sv$df.prior, G*df.residual)
		Stat <- zscoreT(modt, df=df.total, approx=TRUE)
	}

#	Global statistics
	meanStat <- mean(Stat)
	varStat <- var(Stat)

	tab <- matrix(0,nsets,5)
	rownames(tab) <- names(index)
	colnames(tab) <- c("NGenes","Correlation","Down","Up","TwoSided")
	for (i in 1:nsets) {
		iset <- index[[i]]
		if(is.character(iset)) iset <- which(ID %in% iset)
		StatInSet <- Stat[iset]
		m <- length(StatInSet)
		m2 <- G-m
		if(fixed.cor) {
			correlation <- inter.gene.cor
			vif <- 1+(m-1)*correlation
		} else {
			if(m>1) {
				Uset <- U[iset,,drop=FALSE]
				vif <- m * mean(colMeans(Uset)^2)
				correlation <- (vif-1)/(m-1)
			} else {
				vif <- 1
				correlation <- NA
			}
		}

		tab[i,1] <- m
		tab[i,2] <- correlation
		if(use.ranks) {
			if(!allow.neg.cor) correlation <- max(0,correlation)
			tab[i,3:4] <- rankSumTestWithCorrelation(iset,statistics=Stat,correlation=correlation,df=df.camera)
		} else {	
			if(!allow.neg.cor) vif <- max(1,vif)
			meanStatInSet <- mean(StatInSet)
			delta <- G/m2*(meanStatInSet-meanStat)
			varStatPooled <- ( (G-1)*varStat - delta^2*m*m2/G ) / (G-2)
			two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
			tab[i,3] <- pt(two.sample.t,df=df.camera)
			tab[i,4] <- pt(two.sample.t,df=df.camera,lower.tail=FALSE)
		}
	}
	tab[,5] <- 2*pmin(tab[,3],tab[,4])

#	New column names (Jan 2013)
	tab <- data.frame(tab,stringsAsFactors=FALSE)
	Direction <- rep.int("Up",nsets)
	Direction[tab$Down < tab$Up] <- "Down"
	tab$Direction <- Direction
	tab$PValue <- tab$TwoSided
	tab$Down <- tab$Up <- tab$TwoSided <- NULL

#	Remove correlation column if it was not estimated
	if(fixed.cor) tab$Correlation <- NULL

#	Add FDR
	if(nsets>1) tab$FDR <- p.adjust(tab$PValue,method="BH")

#	Sort by p-value
	if(sort && nsets>1) {
		o <- order(tab$PValue)
		tab <- tab[o,]
	}

	tab
}
