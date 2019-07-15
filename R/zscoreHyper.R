zscoreHyper <- function(q,m,n,k) 
#  Z-score equivalents for deviates from hypergeometric distribution
#  Gordon Smyth
#  12 Aug 2012
{
	z <- q
	suppressWarnings(d <- dhyper(q,m,n,k,log=TRUE)-log(2))
	suppressWarnings(pupper <- phyper(q,m,n,k,lower.tail=FALSE,log.p=TRUE))
	suppressWarnings(plower <- phyper(q-1,m,n,k,lower.tail=TRUE,log.p=TRUE))
	d[is.na(d)] <- -Inf
	pupper[is.na(pupper)] <- -Inf
	plower[is.na(plower)] <- -Inf

#	Add half point probability to upper tail probability preserving log-accuracy
	a <- pupper
	b <- d-pupper
	a[b>0] <- d[b>0]
	b <- -abs(b)
	pmidupper <- a+log1p(exp(b))
	pmidupper[is.infinite(a)] <- a[is.infinite(a)]

#	Similarly for lower tail probability preserving log-accuracy
	a <- plower
	b <- d-plower
	a[b>0] <- d[b>0]
	b <- -abs(b)
	pmidlower <- a+log1p(exp(b))
	pmidlower[is.infinite(a)] <- a[is.infinite(a)]

	up <- pmidupper<pmidlower
	if(any(up)) z[up] <- qnorm(pmidupper[up],lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(pmidlower[!up],lower.tail=TRUE,log.p=TRUE)
	z
}

