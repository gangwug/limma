##  GENESET.R

geneSetTest <- function(index,statistics,alternative="mixed",type="auto",ranks.only=TRUE,nsim=9999)
#	Competitive gene set test using either rank sum test or simulation.
#	Gordon Smyth
#	3 September 2004. Last modified 11 Feb 2012.
{
	alternative <- match.arg(alternative,c("mixed","either","down","up","less","greater","two.sided"))
	if(alternative=="two.sided") alternative <- "either"
	if(alternative=="less") alternative <- "down"
	if(alternative=="greater") alternative <- "up"
	type <- match.arg(tolower(type),c("auto","t","f"))
	allsamesign <- all(statistics >= 0) || all(statistics <= 0)
	if(type=="auto") {
		if(allsamesign)
			type <- "f"
		else
			type <- "t"
	}
	if(type=="f" & alternative!="mixed") stop("Only alternative=\"mixed\" is possible with F-like statistics.")
	if(alternative=="mixed") statistics <- abs(statistics)
	if(alternative=="down") {
		statistics <- -statistics
		alternative <- "up"
	}
	if(ranks.only) {
#		The test statistic is the mean rank of the selected genewise statistics
#		and the p-value is obtained explicitly from the rank sum test
		pvalues <- rankSumTestWithCorrelation(index=index,statistics=statistics,df=Inf)
		p.value <- switch(alternative,
			"down" = pvalues["less"],
			"up" = pvalues["greater"],
			"either" = 2*min(pvalues),
			"mixed" = pvalues["greater"])
	} else {
#		The global test statistic is the mean of the selected genewise statistics
#		and the p-value is obtained by random permutation
		ssel <- statistics[index]
		ssel <- ssel[!is.na(ssel)]
		nsel <- length(ssel)
		if(nsel==0) return(1)
		stat <- statistics[!is.na(statistics)]
		msel <- mean(ssel)
		if(alternative=="either")
			posstat <- abs
		else
			posstat <- function(x) x
		msel <- posstat(msel)
		ntail <- 1
		for (i in 1:nsim) if(posstat(mean(sample(stat,nsel))) >= msel) ntail <- ntail+1
		p.value <- ntail/(nsim+1)
	}
	as.vector(p.value)
}

wilcoxGST <- function(index,statistics,...)
#	Mean-rank gene set test using Wilcox test.
#	Gordon Smyth
#	27 July 2009.  Last modified 3 Sep 2011.
geneSetTest(index=index,statistics=statistics,...,ranks.only=TRUE)

