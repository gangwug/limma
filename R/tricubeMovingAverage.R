tricubeMovingAverage <- function(x,span=0.5,power=3)
#	Moving average filter with tricube weights for a time series
#	Gordon Smyth and Yifang Hu
#	Created 5 Feb 2014.  Last modified 6 Nov 2015.
{
#	Check span
	if(length(span)>1) {
		warning("only first value of span used")
		span <- span[1]
	}
	if(span>1) span <- 1
	if(span<=0) return(x)

#	Check power
	if(length(power)>1) {
		warning("only first value of power used")
		span <- span[1]
	}
	if(power<0) power <- 0

#	Convert span to width of moving window
	n <- length(x)
	width <- span*n

#	Round width of smoothing window to nearest odd number
	hwidth <- as.integer(width %/% 2L)
	width <- 2L * hwidth + 1L

#	Make sure window width can't be greater than n
	if(width>n) {
		width <- width-2L
		hwidth <- hwidth-1L
	}

#	If width is 1 or smaller, return original series
	if(hwidth <= 0L) return(x)

#	Tricube weights with all positive values
	u <- seq(from=-1,to=1,length=width) * width / (width+1)
	tricube.weights <- (1-abs(u)^3)^power
	tricube.weights <- tricube.weights/sum(tricube.weights)

#	Extend x with zeros and compute weighted moving averages
	z <- numeric(hwidth)
	x <- as.vector(filter(c(z,x,z),tricube.weights),mode="numeric")[(hwidth+1):(n+hwidth)]

#	Rescale boundary values to remove influence of the outside values of zero
	cw <- cumsum(tricube.weights)
	x[1:hwidth] <- x[1:hwidth] / cw[(width-hwidth):(width-1)]
	x[(n-hwidth+1):n] <- x[(n-hwidth+1):n] / cw[(width-1):(width-hwidth)]

	x
}
