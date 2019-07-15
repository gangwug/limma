#  ROC.R

auROC <- function(truth, stat=NULL)
#	Area under Receiver Operating Curve for empirical data
#	Gordon Smyth
#	21 Dec 2003. Last modified 28 April 2015.
{
#	Don't allow any NA
	if(any(is.na(truth))) return(NA)

#	Make logical and integer vectors
	ntests <- length(truth)
	truthl <- as.logical(truth)
	truthi <- as.integer(truthl)

#	Return NA if truth is always the same
	npos <- sum(truthi)
	if(npos==0 || npos==ntests) return(NA)

	if(is.null(stat)) {
		sensitivity <- cumsum(truthi) / npos
		return( mean(sensitivity[!truthl]) )
	}
	
#	From here, stat is not NULL

#	Check stat
	stat <- as.vector(stat)
	if(length(stat) != ntests) stop("lengths differ")
	if(any(is.na(stat))) return(NA)

#	From here, stat is not NA

#	Sort truth by stat
	o <- order(stat,decreasing=TRUE)
	truthi <- truthi[o]
	truthl <- truthl[o]
	stat <- stat[o]
	sensitivity <- cumsum(truthi) / npos

#	Replace sensitivities with averages for tied stat
	i <- 1:(ntests-1)
	iseq2prev <- stat[i]==stat[i+1]
	if(any(iseq2prev)) {
		iseq2prev <- c(FALSE,iseq2prev)
		tied.first <- which(!iseq2prev)
		tied.last <- c(tied.first[-1]-1,ntests)
		sensitivity.last <- sensitivity[tied.last]
		sensitivity.previous <- c(0,sensitivity.last[-length(sensitivity.last)])
		sensitivity.average <- (sensitivity.last+sensitivity.previous)/2
		sensitivity <- rep(sensitivity.average, tied.last-tied.first+1)
	}

	mean(sensitivity[!truthl])
}
