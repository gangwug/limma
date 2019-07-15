fry <- function(y,...) UseMethod("fry")

fry.default <- function(y,index=NULL,design=NULL,contrast=ncol(design),geneid=NULL,standardize="posterior.sd",sort="directional",...)
#	Quick version of roast gene set test assuming equal variances between genes
#	The up and down p-values are equivalent to those from roast with nrot=Inf
#	in the special case of prior.df=Inf.
#	Gordon Smyth and Goknur Giner
#	Created 30 January 2015.  Last modified 11 May 2016
{
#	Partial matching of extra arguments
#	array.weights and trend.var included for (undocumented) backward compatibility with code before 9 May 2016.
	Dots <- list(...)
	PossibleArgs <- c("array.weights","weights","block","correlation","trend.var","robust","winsor.tail.p")
	if(!is.null(names(Dots))) {
		i <- pmatch(names(Dots),PossibleArgs)
		names(Dots) <- PossibleArgs[i]
	}

#	Defaults for extra arguments
	if(is.null(Dots$trend.var)) {
		Dots$trend <- FALSE
	} else {
		Dots$trend <- Dots$trend.var
		Dots$trend.var <- NULL
	}
	if(is.null(Dots$robust)) Dots$robust <- FALSE
	if(Dots$robust & is.null(Dots$winsor.tail.p)) Dots$winsor.tail.p <- c(0.05,0.1)
	if(!is.null(Dots$block) & is.null(Dots$correlation)) stop("Intra-block correlation must be specified")

#	Covariate for trended eBayes
	covariate <- NULL
	if(Dots$trend) covariate <- rowMeans(as.matrix(y))

#	Compute effects matrix (with df.residual+1 columns)
	Effects <- .lmEffects(y=y,design=design,contrast=contrast,
		array.weights=Dots$array.weights,
		weights=Dots$weights,
		block=Dots$block,
		correlation=Dots$correlation)

#	Divide out genewise standard deviations
	standardize <- match.arg(standardize, c("none","residual.sd","posterior.sd","p2"))
	if(!standardize=="none") {

#		Estimate genewise sds robustly
		OK <- requireNamespace("statmod",quietly=TRUE)
		if(!OK) stop("statmod package required but isn't installed (or can't be loaded)")
		gq <- statmod::gauss.quad.prob(128,"uniform")
		df.residual <- ncol(Effects)-1
		Eu2max <- sum( (df.residual+1)*gq$nodes^df.residual*qchisq(gq$nodes,df=1)*gq$weights )
		u2max <- apply(Effects^2,1,max)
		s2.robust <- (rowSums(Effects^2)-u2max) / (df.residual+1-Eu2max)

#		Empirical Bayes moderation method I: estimate hyperparameters from robust variances
		if(standardize=="p2") {
			sv <- squeezeVar(s2.robust, df=0.92*df.residual,
				covariate=covariate,
				robust=Dots$robust,
				winsor.tail.p=Dots$winsor.tail.p)
			s2.robust <- sv$var.post
		}

#		Empirical Bayes moderation method II: estimate hyperparameters from residual variances but apply squeezing to robust variances
		if(standardize=="posterior.sd") {
			s2 <- rowMeans(Effects[,-1,drop=FALSE]^2)
			if(Dots$robust) {
				fit <- fitFDistRobustly(s2, df1=df.residual, covariate=covariate, winsor.tail.p=Dots$winsor.tail.p)
				df.prior <- fit$df2.shrunk
			} else {
				fit <- fitFDist(s2, df1=df.residual, covariate=covariate)
				df.prior <- fit$df2
			}
			s2.robust <- .squeezeVar(s2.robust,df=0.92*df.residual,var.prior=fit$scale,df.prior=df.prior)
		}

		Effects <- Effects/sqrt(s2.robust)
	}

#	Check geneid (can be a vector of gene IDs or an annotation column name)
	if(is.null(geneid)) {
		geneid <- rownames(Effects)
	} else {
		geneid <- as.character(geneid)
		if(length(geneid)==1) 
			geneid <- as.character(y$genes[,geneid])
		else
			if(length(geneid) != nrow(y)) stop("geneid vector should be of length nrow(y)")
	}

	.fryEffects(effects=Effects,index=index,geneid=geneid,sort=sort)
}


.fryEffects <- function(effects,index=NULL,geneid=rownames(effects),sort=TRUE)
#	fry given the effects matrix
#	Gordon Smyth and Goknur Giner
#	Created 30 January 2015.  Last modified 28 June 2016
{
	G <- nrow(effects)
	neffects <- ncol(effects)
	df.residual <- neffects-1L

#	Check index
	if(is.null(index)) index <- list(set1=1L:G)
	if(is.data.frame(index) || !is.list(index)) index <- list(set1=index)
	nsets <- length(index)
	if(nsets==0L) stop("index is empty")
	if(is.null(names(index)))
		names(index) <- paste0("set",formatC(1L:nsets,width=1L+floor(log10(nsets)),flag="0"))
	else
		if(anyDuplicated(names(index))) stop("Gene sets don't have unique names",call. =FALSE)

#	Global statistics
	NGenes <- rep.int(0L,nsets)
	PValue.Mixed <- t.stat <- rep.int(0,nsets)
	for (i in 1:nsets) {
		iset <- index[[i]]
		if(is.data.frame(iset)) {
			if(ncol(iset)>1 && is.numeric(iset[,2])) {
				iw <- iset[,2]
				iset <- iset[,1]
				if(is.factor(iset)) iset <- as.character(iset)
				if(is.character(iset)) {
					if(anyDuplicated(iset)) stop("Duplicate gene ids in set ",i)
					j <- match(geneid,iset)
					inset <- which(!is.na(j))
					EffectsSet <- effects[inset,,drop=FALSE]
					iw <- iw[j[inset]]
				} else {
					EffectsSet <- effects[iset,,drop=FALSE]
				}
				MeanEffectsSet <- iw %*% EffectsSet
			} else {
				stop("index ",i," is a data.frame but doesn't contain gene weights")
			}
		} else {
			if(is.factor(iset)) iset <- as.character(iset)
			if(is.character(iset)) iset <- which(geneid %in% iset)
			EffectsSet <- effects[iset,,drop=FALSE]
			MeanEffectsSet <- colMeans(EffectsSet)
		}
		t.stat[i] <- MeanEffectsSet[1] / sqrt(mean(MeanEffectsSet[-1]^2))
		NGenes[i] <- nrow(EffectsSet)

		if(NGenes[i]>1) {
			SVD <- svd(EffectsSet,nu=0)
			A <- SVD$d^2
			d1 <- length(A)
			d <- d1-1L
			beta.mean <- 1/d1
			beta.var <- d/d1/d1/(d1/2+1)
			Fobs <- (sum(EffectsSet[,1]^2)-A[d1]) / (A[1]-A[d1])
			Frb.mean <- (sum(A) * beta.mean - A[d1]) / (A[1]-A[d1])
			COV <- matrix(-beta.var/d,d1,d1)
			diag(COV) <- beta.var
			Frb.var <- (A %*% COV %*% A ) / (A[1]-A[d1])^2
			alphaplusbeta <- Frb.mean*(1-Frb.mean)/Frb.var-1
			alpha <- alphaplusbeta*Frb.mean
			beta <- alphaplusbeta-alpha
			PValue.Mixed[i] <- pbeta(Fobs,shape1=alpha,shape2=beta,lower.tail=FALSE)
		}
	}

	Direction <- rep.int("Up",nsets)
	Direction[t.stat<0] <- "Down"
	PValue <- 2*pt(-abs(t.stat),df=df.residual)
	PValue.Mixed[NGenes==1] <- PValue[NGenes==1]

#	Add FDR
	if(nsets>1) {
		FDR <- p.adjust(PValue,method="BH")
		FDR.Mixed <- p.adjust(PValue.Mixed,method="BH")
		tab <- data.frame(NGenes=NGenes,Direction=Direction,PValue=PValue,FDR=FDR,PValue.Mixed=PValue.Mixed,FDR.Mixed=FDR.Mixed)
	} else {
		tab <- data.frame(NGenes=NGenes,Direction=Direction,PValue=PValue,PValue.Mixed=PValue.Mixed)
	}
	rownames(tab) <- names(index)

#	Sort results
	if(is.logical(sort)) if(sort) sort <- "directional" else sort <- "none"
	sort <- match.arg(sort,c("directional","mixed","none"))
	if(sort=="none") return(tab)
	if(sort=="directional") {
		o <- order(tab$PValue,-tab$NGenes,tab$PValue.Mixed)
	} else {
		o <- order(tab$PValue.Mixed,-tab$NGenes,tab$PValue)
	}
	tab[o,,drop=FALSE]
}
