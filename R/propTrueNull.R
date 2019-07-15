propTrueNull <- function(p, method="lfdr", nbins=20, ...)
#	Estimate proportion of null p-values
#	Belinda Phipson and Gordon Smyth
#	Created 23 May 2012. Last revised 2 Dec 2012.
{
	method <- match.arg(method, c("lfdr","mean","hist","convest"))
	switch(method,
		lfdr = .propTrueNullByLocalFDR(p),
		mean = .propTrueNullByMeanP(p),
		hist = .propTrueNullFromHistogram(p, nbins=nbins),
		convest = convest(p, ...)
	)
}

.propTrueNullByLocalFDR <- function(p)
#	Estimate proportion of null p-values
#	by average local FDR
#	Belinda Phipson and Gordon Smyth
#	23 May 2012. Last revised 30 July 2012.
{
	n <- length(p)
	i <- n:1L
	p <- sort(p, decreasing = TRUE)
	q <- pmin(n/i * p, 1)
	n1 <- n+1L
	sum(i*q) / n/n1*2
}

.propTrueNullByMeanP <- function(p)
#	Estimate proportion of null p-values
#	by mean p-value
#	Gordon Smyth
#	12 Feb 2013.
{
	n <- length(p)
	i <- 1:n
	p <- sort(p)
	q <- pmin(p, (i-0.5)/n)
	2*mean(q)
}

.propTrueNullFromHistogram <- function(p, nbins=20) 
#	Estimate proportion of null p-values
#	by histogram method

#	Adapted by Gordon Smyth from the function estimate.m0 from
#	http://www.public.iastate.edu/~dnett/microarray/multtest.txt
#	Accessed March 2012

#	Created 2 Dec 2012. Last revised 9 May 2016.
{
	bin <- c(-0.1, (1:nbins)/nbins)
	bin.counts <- rep_len(0L, length.out=nbins)
	tab <- tabulate(cut(p,bin))
	bin.counts[1:length(tab)] <- tab
	tail.means <- rev(cumsum(rev(bin.counts))/(1:nbins))
	index <- which(tail.means >= bin.counts)[1]
	tail.means[index]/tail.means[1]
}


