#  SCORE.R

zscore <- function(q, distribution=NULL, ...) 
#  Z-score equivalents for deviates from specified distribution
#  Gordon Smyth
#  13 June 2012
{
	z <- q
	n <- length(q)
	pdist <- get(paste("p",as.character(distribution),sep=""))
	pupper <- pdist(q,...,lower.tail=FALSE,log.p=TRUE)
	plower <- pdist(q,...,lower.tail=TRUE,log.p=TRUE)
	up <- pupper<plower
	if(any(up)) z[up] <- qnorm(pupper[up],lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(plower[!up],lower.tail=TRUE,log.p=TRUE)
	z
}

zscoreGamma <- function(q, shape, rate = 1, scale = 1/rate) 
#  Z-score equivalents for gamma deviates
#  Gordon Smyth
#  1 October 2003
{
	z <- q
	n <- length(q)
	shape <- rep(shape,length.out=n)
	scale <- rep(scale,length.out=n)
	up <- (q > shape*scale)
	if(any(up)) z[up] <- qnorm(pgamma(q[up],shape=shape[up],scale=scale[up],lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(pgamma(q[!up],shape=shape[!up],scale=scale[!up],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	z
}

zscoreT <- function(x, df, approx=FALSE)
#  Z-score equivalents for t distribution deviates
#  Gordon Smyth
#  24 August 2003. Last modified 3 June 2014.
{
	if(approx) {
		if(all(df > 10000)) return(x)

#		Approximation from Hill (1970)
		A <- df-0.5
		B <- 48*A*A
		W <- A*log1p(x*x/df)
		z <- (((((-0.4*W-3.3)*W-24)*W-85.5)/(0.8*W*W+100+B)+W+3)/B+1)*sqrt(W)
		z[x<0] <- -z[x<0]
	} else {
		z <- x
		df <- rep(df,length.out=length(x))
		pos <- x>0
		if(any(pos)) z[pos] <- qnorm(pt(x[pos],df=df[pos],lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE) 
		if(any(!pos)) z[!pos] <- qnorm(pt(x[!pos],df=df[!pos],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	}
	z
}

tZscore <- function(x, df)
#  t-statistics equivalents for z-scores deviates
#  Gordon Smyth
#  1 June 2004
{
	z <- x
	df <- rep(df,length.out=length(x))
	pos <- x>0
	if(any(pos)) z[pos] <- qt(pnorm(x[pos],lower.tail=FALSE,log.p=TRUE),df=df[pos],lower.tail=FALSE,log.p=TRUE) 
	if(any(!pos)) z[!pos] <- qt(pnorm(x[!pos],lower.tail=TRUE,log.p=TRUE),df=df[!pos],lower.tail=TRUE,log.p=TRUE)
	z
}
