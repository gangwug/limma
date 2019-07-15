#  TOPTABLE.R

topTable <- function(fit,coef=NULL,number=10,genelist=fit$genes,adjust.method="BH",sort.by="B",resort.by=NULL,p.value=1,lfc=0,confint=FALSE)
#	Summary table of top genes, object-orientated version
#	Gordon Smyth
#	4 August 2003.  Last modified 1 Feb 2018.
{
#	Check fit
	if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM object")
	if(is.null(fit$t) && is.null(fit$F)) stop("Need to run eBayes or treat first")
	if(is.null(fit$coefficients)) stop("coefficients not found in fit object")
	if(confint && is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")

	if(is.null(coef)) {
		if(is.null(fit$treat.lfc)) {
			coef <- 1:ncol(fit)
			cn <- colnames(fit)
			if(!is.null(cn)) {
				i <- which(cn == "(Intercept)")
				if(length(i)) {
					coef <- coef[-i]
					message("Removing intercept from test coefficients")
				}
			}
		} else
			coef <- ncol(fit)
	}
	if(length(coef)>1) {
		if(!is.null(fit$treat.lfc)) stop("Treat p-values can only be displayed for single coefficients")
		coef <- unique(coef)
		if(length(fit$coef[1,coef]) < ncol(fit)) fit <- fit[,coef]
		if(sort.by=="B") sort.by <- "F"
		return(topTableF(fit,number=number,genelist=genelist,adjust.method=adjust.method,sort.by=sort.by,p.value=p.value,lfc=lfc))
	}
	fit <- unclass(fit)
	ebcols <- c("t","p.value","lods")
	if(confint) ebcols <- c("s2.post","df.total",ebcols)
	.topTableT(fit=fit[c("coefficients","stdev.unscaled")],
		coef=coef,
		number=number,
		genelist=genelist,
		A=fit$Amean,
		eb=fit[ebcols],
		adjust.method=adjust.method,
		sort.by=sort.by,
		resort.by=resort.by,
		p.value=p.value,
		lfc=lfc,
		confint=confint)
}

topTableF <- function(fit,number=10,genelist=fit$genes,adjust.method="BH",sort.by="F",p.value=1,lfc=0)
#	Summary table of top genes by F-statistic
#	Gordon Smyth
#	27 August 2006. Last modified 24 June 2014.
{
#	Check fit
	if(is.null(fit$coefficients)) stop("Coefficients not found in fit")
	M <- as.matrix(fit$coefficients)
	rn <- rownames(M)
	if(is.null(colnames(M))) colnames(M) <- paste("Coef",1:ncol(M),sep="")
	Amean <- fit$Amean
	Fstat <- fit$F
	Fp <- fit$F.p.value
	if(is.null(Fstat)) stop("F-statistics not found in fit")

#	Ensure genelist is a data.frame
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ProbeID=genelist,stringsAsFactors=FALSE)

#	Check rownames
	if(is.null(rn))
		rn <- 1:nrow(M)
	else
		if(anyDuplicated(rn)) {
			if(is.null(genelist))
				genelist <- data.frame(ID=rn,stringsAsFactors=FALSE)
			else
				if("ID" %in% names(genelist))
					genelist$ID0 <- rn
				else
					genelist$ID <- rn
			rn <- 1:nrow(M)
		}

#	Check sort.by
	sort.by <- match.arg(sort.by,c("F","none"))

#	Apply multiple testing adjustment
	adj.P.Value <- p.adjust(Fp,method=adjust.method)

#	Thin out fit by lfc and p.value thresholds
	if(lfc > 0 || p.value < 1) {
		if(lfc>0)
			big <- rowSums(abs(M)>lfc,na.rm=TRUE)>0
		else
			big <- TRUE
		if(p.value<1) {
			sig <- adj.P.Value <= p.value
			sig[is.na(sig)] <- FALSE
		} else
			sig <- TRUE
		keep <- big & sig
		if(!all(keep)) {
			M <- M[keep,,drop=FALSE]
			rn <- rn[keep]
			Amean <- Amean[keep]
			Fstat <- Fstat[keep]
			Fp <- Fp[keep]
			genelist <- genelist[keep,,drop=FALSE]
			adj.P.Value <- adj.P.Value[keep]
		}
	}

#	Enough rows left?
	if(nrow(M) < number) number <- nrow(M)
	if(number < 1) return(data.frame())

#	Find rows of top genes
	if(sort.by=="F")
		o <- order(Fp,decreasing=FALSE)[1:number]
	else
		o <- 1:number

#	Assemble data.frame
	if(is.null(genelist))
		tab <- data.frame(M[o,,drop=FALSE])
	else
		tab <- data.frame(genelist[o,,drop=FALSE],M[o,,drop=FALSE])
	tab$AveExpr <- Amean[o]
	tab <- data.frame(tab,F=Fstat[o],P.Value=Fp[o],adj.P.Val=adj.P.Value[o])
	rownames(tab) <- rn[o]
	tab
}

.topTableT <- function(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="BH",sort.by="B",resort.by=NULL,p.value=1,lfc=0,confint=FALSE,...)
#	Summary table of top genes for a single coefficient
#	Gordon Smyth
#	21 Nov 2002. Forked from toptable() 1 Feb 2018. Last revised 1 Feb 2018.
{
#	Check fit
	fit$coefficients <- as.matrix(fit$coefficients)
	rn <- rownames(fit$coefficients)

#	Check coef is length 1
	if(length(coef)>1) {
		coef <- coef[1]
		warning("Treat is for single coefficients: only first value of coef being used")
	}

#	Ensure genelist is a data.frame
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ID=genelist,stringsAsFactors=FALSE)

#	Check rownames
	if(is.null(rn))
		rn <- 1:nrow(fit$coefficients)
	else
		if(anyDuplicated(rn)) {
			if(is.null(genelist))
				genelist <- data.frame(ID=rn,stringsAsFactors=FALSE)
			else
				if("ID" %in% names(genelist))
					genelist$ID0 <- rn
				else
					genelist$ID <- rn
			rn <- 1:nrow(fit$coefficients)
		}

#	Check sort.by
	sort.by <- match.arg(sort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B","none"))
	if(sort.by=="M") sort.by="logFC"
	if(sort.by=="A" || sort.by=="Amean") sort.by <- "AveExpr"
	if(sort.by=="T") sort.by <- "t"
	if(sort.by=="p") sort.by <- "P"

#	Check resort.by
	if(!is.null(resort.by)) {
		resort.by <- match.arg(resort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B"))
		if(resort.by=="M") resort.by <- "logFC"
		if(resort.by=="A" || resort.by=="Amean") resort.by <- "AveExpr"
		if(resort.by=="p") resort.by <- "P"
		if(resort.by=="T") resort.by <- "t"
	}

#	Check A
	if(is.null(A)) {
		if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
	} else {
		if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
	}

#	Check for lods component
	if(is.null(eb$lods)) {
		if(sort.by=="B") stop("Trying to sort.by B, but B-statistic (lods) not found in MArrayLM object",call.=FALSE)
		if(!is.null(resort.by)) if(resort.by=="B") stop("Trying to resort.by B, but B-statistic (lods) not found in MArrayLM object",call.=FALSE)
		include.B <- FALSE
	} else {
		include.B <- TRUE
	}

#	Extract statistics for table
	M <- fit$coefficients[,coef]
	tstat <- as.matrix(eb$t)[,coef]
	P.Value <- as.matrix(eb$p.value)[,coef]
	if(include.B) B <- as.matrix(eb$lods)[,coef]

#	Apply multiple testing adjustment
	adj.P.Value <- p.adjust(P.Value,method=adjust.method)

#	Thin out fit by p.value and lfc thresholds	
	if(p.value < 1 | lfc > 0) {
		sig <- (adj.P.Value <= p.value) & (abs(M) >= lfc)
		if(any(is.na(sig))) sig[is.na(sig)] <- FALSE
		if(all(!sig)) return(data.frame())
		genelist <- genelist[sig,,drop=FALSE]
		M <- M[sig]
		A <- A[sig]
		tstat <- tstat[sig]
		P.Value <- P.Value[sig]
		adj.P.Value <- adj.P.Value[sig]
		if(include.B) B <- B[sig]
		rn <- rn[sig]
	}

#	Are enough rows left?
	if(length(M) < number) number <- length(M)
	if(number < 1) return(data.frame())

#	Select top rows
	ord <- switch(sort.by,
		logFC=order(abs(M),decreasing=TRUE),
		AveExpr=order(A,decreasing=TRUE),
		P=order(P.Value,decreasing=FALSE),
		t=order(abs(tstat),decreasing=TRUE),
		B=order(B,decreasing=TRUE),
		none=1:length(M)
	)
	top <- ord[1:number]

#	Assemble output data.frame
	if(is.null(genelist))
		tab <- data.frame(logFC=M[top])
	else {
		tab <- data.frame(genelist[top,,drop=FALSE],logFC=M[top],stringsAsFactors=FALSE)
	}
	if(confint) {
		if(is.numeric(confint)) alpha <- (1+confint[1])/2 else alpha <- 0.975
		margin.error <- sqrt(eb$s2.post[top])*fit$stdev.unscaled[top,coef]*qt(alpha,df=eb$df.total[top])
		tab$CI.L <- M[top]-margin.error
		tab$CI.R <- M[top]+margin.error
	}
	if(!is.null(A)) tab$AveExpr <- A[top]
	tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],adj.P.Val=adj.P.Value[top])
	if(include.B) tab$B <- B[top]
	rownames(tab) <- rn[top]

#	Resort table
	if(!is.null(resort.by)) {
		ord <- switch(resort.by,
			logFC=order(tab$logFC,decreasing=TRUE),
			AveExpr=order(tab$AveExpr,decreasing=TRUE),
			P=order(tab$P.Value,decreasing=FALSE),
			t=order(tab$t,decreasing=TRUE),
			B=order(tab$B,decreasing=TRUE)
		)
		tab <- tab[ord,]
	}

	tab
}

toptable <- function(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="BH",sort.by="B",resort.by=NULL,p.value=1,lfc=0,confint=FALSE,...)
#	Summary table of top genes for a single coefficient
#	Old non-object-orientated version, marked as deprecated from 1 Feb 2018.
#	Gordon Smyth
#	21 Nov 2002. Last revised 1 Feb 2018.
{
	if(is(fit,"MArrayLM")) {
		message("toptable() is deprecated and will be removed in the future version of limma. Please use topTable() instead.")
	}

#	Check fit
	fit <- unclass(fit)
	fit$coefficients <- as.matrix(fit$coefficients)
	rn <- rownames(fit$coefficients)

#	Check coef is length 1
	if(length(coef)>1) {
		coef <- coef[1]
		warning("Treat is for single coefficients: only first value of coef being used")
	}

#	Ensure genelist is a data.frame
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ID=genelist,stringsAsFactors=FALSE)

#	Check rownames
	if(is.null(rn))
		rn <- 1:nrow(fit$coefficients)
	else
		if(anyDuplicated(rn)) {
			if(is.null(genelist))
				genelist <- data.frame(ID=rn,stringsAsFactors=FALSE)
			else
				if("ID" %in% names(genelist))
					genelist$ID0 <- rn
				else
					genelist$ID <- rn
			rn <- 1:nrow(fit$coefficients)
		}

#	Check sort.by
	sort.by <- match.arg(sort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B","none"))
	if(sort.by=="M") sort.by="logFC"
	if(sort.by=="A" || sort.by=="Amean") sort.by <- "AveExpr"
	if(sort.by=="T") sort.by <- "t"
	if(sort.by=="p") sort.by <- "P"

#	Check resort.by
	if(!is.null(resort.by)) {
		resort.by <- match.arg(resort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B"))
		if(resort.by=="M") resort.by <- "logFC"
		if(resort.by=="A" || resort.by=="Amean") resort.by <- "AveExpr"
		if(resort.by=="p") resort.by <- "P"
		if(resort.by=="T") resort.by <- "t"
	}

#	Check A
	if(is.null(A)) {
		if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
	} else {
		if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
	}

#	If eb is NULL, look for components in fit
	if(is.null(eb)) {
		ebcols <- c("t","p.value","lods")
		if(confint) ebcols <- c("s2.post","df.total",ebcols)
		if(all(ebcols %in% names(fit))) eb <- fit[ebcols]
#		If still NULL, compute using ebayes just for the one column required
		if(is.null(eb)) {
			message("Argument 'eb' not provided, so running ebayes() to get moderated t-statistics and p-values.")
			fit$coefficients <- fit$coefficients[,coef,drop=FALSE]
			fit$stdev.unscaled <- as.matrix(fit$stdev.unscaled)[,coef,drop=FALSE]
			eb <- .ebayes(fit,...)
			coef <- 1
		}
	}

#	Check for lods component
	if(is.null(eb$lods)) {
		if(sort.by=="B") stop("Trying to sort.by B, but B-statistic (lods) not found in MArrayLM object")
		if(!is.null(resort.by)) if(resort.by=="B") stop("Trying to resort.by B, but B-statistic (lods) not found in MArrayLM object",call.=FALSE)
		include.B <- FALSE
	} else {
		include.B <- TRUE
	}

#	Extract statistics for table
	M <- fit$coefficients[,coef]
	tstat <- as.matrix(eb$t)[,coef]
	P.Value <- as.matrix(eb$p.value)[,coef]
	if(include.B) B <- as.matrix(eb$lods)[,coef]

#	Apply multiple testing adjustment
	adj.P.Value <- p.adjust(P.Value,method=adjust.method)

#	Thin out fit by p.value and lfc thresholds	
	if(p.value < 1 | lfc > 0) {
		sig <- (adj.P.Value <= p.value) & (abs(M) >= lfc)
		if(any(is.na(sig))) sig[is.na(sig)] <- FALSE
		if(all(!sig)) return(data.frame())
		genelist <- genelist[sig,,drop=FALSE]
		M <- M[sig]
		A <- A[sig]
		tstat <- tstat[sig]
		P.Value <- P.Value[sig]
		adj.P.Value <- adj.P.Value[sig]
		if(include.B) B <- B[sig]
		rn <- rn[sig]
	}

#	Are enough rows left?
	if(length(M) < number) number <- length(M)
	if(number < 1) return(data.frame())

#	Select top rows
	ord <- switch(sort.by,
		logFC=order(abs(M),decreasing=TRUE),
		AveExpr=order(A,decreasing=TRUE),
		P=order(P.Value,decreasing=FALSE),
		t=order(abs(tstat),decreasing=TRUE),
		B=order(B,decreasing=TRUE),
		none=1:length(M)
	)
	top <- ord[1:number]

#	Assemble output data.frame
	if(is.null(genelist))
		tab <- data.frame(logFC=M[top])
	else {
		tab <- data.frame(genelist[top,,drop=FALSE],logFC=M[top],stringsAsFactors=FALSE)
	}
	if(confint) {
		if(is.numeric(confint)) alpha <- (1+confint[1])/2 else alpha <- 0.975
		margin.error <- sqrt(eb$s2.post[top])*fit$stdev.unscaled[top,coef]*qt(alpha,df=eb$df.total[top])
		tab$CI.L <- M[top]-margin.error
		tab$CI.R <- M[top]+margin.error
	}
	if(!is.null(A)) tab$AveExpr <- A[top]
	tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],adj.P.Val=adj.P.Value[top])
	if(include.B) tab$B <- B[top]
	rownames(tab) <- rn[top]

#	Resort table
	if(!is.null(resort.by)) {
		ord <- switch(resort.by,
			logFC=order(tab$logFC,decreasing=TRUE),
			AveExpr=order(tab$AveExpr,decreasing=TRUE),
			P=order(tab$P.Value,decreasing=FALSE),
			t=order(tab$t,decreasing=TRUE),
			B=order(tab$B,decreasing=TRUE)
		)
		tab <- tab[ord,]
	}

	tab
}

