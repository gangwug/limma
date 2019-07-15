rankSumTestWithCorrelation <- function(index,statistics,correlation=0,df=Inf)
#	Rank sum test as for two-sample Wilcoxon-Mann-Whitney test,
#	but allowing for correlation between members of test set.
#	Gordon Smyth and Di Wu
#	Created 2007.  Last modified 24 Feb 2012.
{
	n <- length(statistics)
	r <- rank(statistics)
	r1 <- r[index]
	n1 <- length(r1)
	n2 <- n-n1
	U <- n1*n2 + n1*(n1+1)/2 - sum(r1)
	mu <- n1*n2/2

	if(correlation==0 || n1==1) {
		sigma2 <- n1*n2*(n+1)/12
	} else {
		sigma2 <- asin(1)*n1*n2 + asin(0.5)*n1*n2*(n2-1) + asin(correlation/2)*n1*(n1-1)*n2*(n2-1) + asin((correlation+1)/2)*n1*(n1-1)*n2
		sigma2 <- sigma2/2/pi
	}

	TIES <- (length(r) != length(unique(r)))
	if(TIES) {
		NTIES <- table(r)
		adjustment <- sum(NTIES*(NTIES+1)*(NTIES-1)) / (n*(n+1)*(n-1))
		sigma2 <- sigma2 * (1 - adjustment)
	}
	zlowertail <- (U+0.5-mu)/sqrt(sigma2)
	zuppertail <- (U-0.5-mu)/sqrt(sigma2)

#	Lower and upper tails are reversed on output
#	because R's ranks are the reverse of Mann-Whitney's ranks
	pvalues <- c(less=pt(zuppertail,df=df,lower.tail=FALSE), greater=pt(zlowertail,df=df))
	pvalues	
}

