cameraPR <- function(statistic,...) UseMethod("cameraPR")

cameraPR.default <- function(statistic,index,use.ranks=FALSE,inter.gene.cor=0.01,sort=TRUE,...)
#	Competitive gene set test allowing for correlation between genes: pre-ranked statistic.
#	Gordon Smyth
#	Created 18 April 2017.
{
#	Issue warning if extra arguments found
	dots <- names(list(...))
	if(length(dots)) warning("Extra arguments disregarded: ",sQuote(dots))

#	Check statistic
	if(is.list(statistic)) stop("statistic should be a numeric vector")
	storage.mode(statistic) <- "numeric"
	if(anyNA(statistic)) stop("NA values for statistic not allowed")
	G <- length(statistic)
	ID <- names(statistic)
	if(G<3) stop("Two few genes in dataset: need at least 3")

#	Check index
	if(!is.list(index)) index <- list(set1=index)
	nsets <- length(index)

#	Check inter.gene.cor
	if(anyNA(inter.gene.cor)) stop("NA inter.gene.cor not allowed")
	if(any(abs(inter.gene.cor) >= 1)) stop("inter.gene.cor too large or small")
	if(length(inter.gene.cor) > 1L) {
		if(length(inter.gene.cor) != nsets) stop("Length of inter.gene.cor doesn't match number of sets")
		fixed.cor <- FALSE
	} else {
		fixed.cor <- TRUE
		inter.gene.cor <- rep_len(inter.gene.cor,nsets)
	}

#	Set df
	if(use.ranks)
		df.camera <- Inf
	else
		df.camera <- G-2L

#	Global statistics
	meanStat <- mean(statistic)
	varStat <- var(statistic)

	NGenes <- Down <- Up <- rep_len(0,nsets)
	for (i in 1:nsets) {
		iset <- index[[i]]
		if(is.character(iset)) iset <- which(ID %in% iset)
		StatInSet <- statistic[iset]
		m <- length(StatInSet)
		NGenes[i] <- m
		if(use.ranks) {
			p.value <- rankSumTestWithCorrelation(iset,statistics=statistic,correlation=inter.gene.cor[i],df=df.camera)
			Down[i] <- p.value[1]
			Up[i] <- p.value[2]
		} else {	
			vif <- 1+(m-1)*inter.gene.cor[i]
			m2 <- G-m
			meanStatInSet <- mean(StatInSet)
			delta <- G/m2*(meanStatInSet-meanStat)
			varStatPooled <- ( (G-1L)*varStat - delta^2*m*m2/G ) / (G-2L)
			two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
			Down[i] <- pt(two.sample.t,df=df.camera)
			Up[i] <- pt(two.sample.t,df=df.camera,lower.tail=FALSE)
		}
	}
	TwoSided <- 2*pmin(Down,Up)

#	Assemble into data.frame
	D <- (Down < Up)
	Direction <- rep_len("Up",nsets)
	Direction[D] <- "Down"
	if(fixed.cor)
		tab <- data.frame(NGenes=NGenes,Direction=Direction,PValue=TwoSided,stringsAsFactors=FALSE)
	else
		tab <- data.frame(NGenes=NGenes,Correlation=inter.gene.cor,Direction=Direction,PValue=TwoSided,stringsAsFactors=FALSE)
	rownames(tab) <- names(index)

#	Add FDR
	if(nsets>1L) tab$FDR <- p.adjust(tab$PValue,method="BH")

#	Sort by p-value
	if(sort && nsets>1L) {
		o <- order(tab$PValue)
		tab <- tab[o,]
	}

	tab
}
