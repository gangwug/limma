logcosh <- function(x)
#	Compute log(cosh(x)) without floating over or underflow
#	Gordon Smyth  8 April 2007
{
	y <- abs(x)-log(2)
	i <- abs(x) < 1e-4
	y[i] <- 0.5*x[i]^2
	i <- !i & (abs(x) < 17)
	y[i] <- log(cosh(x[i]))
	y
}

logsumexp <- function(x,y)
#	Compute log(exp(x) + exp(y)) without floating overflow or underflow
#	Gordon Smyth
#	1 May 2018
{
#	pmin and pmax
	mi <- pmin(x,y)
	ma <- pmax(x,y)

#	Special cases
	IsNA <- is.na(ma)
	Infma <- ma == Inf
	Infmi <- mi == -Inf
	Special <- IsNA | Infma | Infmi
	if(any(Special)) {
		s <- x
		Infma[IsNA] <- FALSE
		Infmi[IsNA] <- TRUE
		s[Infma] <- Inf
		s[Infmi] <- ma[Infmi]
		i <- which(!Special)
		s[i] <- Recall(x[i],y[i])
		return(s)
	}

#	From here all values are finite
	m <- (x+y)/2
	m + logcosh(ma-m) + log(2)
}
