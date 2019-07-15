detectionPValues <- function(x,...) UseMethod("detectionPValues")

detectionPValues.EListRaw <- function(x,status=NULL,...)
#	Detection p-values from negative controls for EListRaw
#	Gordon Smyth
#	Created 13 June 2016
{
	if(is.null(status)) status <- x$genes$Status
	detectionPValues(x=x$E, status=status, ...)
}

detectionPValues.default <- function(x, status, negctrl="negative",...)
#	Detection p-values from negative controls
#	Gordon Smyth
#	Created 13 June 2016.  Last modified 14 June 2016.
{
	x <- as.matrix(x)

#	Identify negative control rows
	if(missing(status) || is.null(status)) stop("need to specify status vector")
	isneg <- as.integer(status==negctrl)
	notneg <- 1L-isneg
	nneg <- sum(isneg)

#	Count the number of negative controls greater than each value
#	cs1 counts number of negative controls >= the observed value
#	cs2 counts number of negative controls > the observed value
#	By taking the average of cs1 and cs2, we count ties as 1/2
	DetectionPValue <- x
	cs1 <- cs2 <- rep_len(0,length.out=nrow(x))
	for (j in 1:ncol(x)) {
		o1 <- order(x[,j],isneg,decreasing=TRUE)
		o2 <- order(x[,j],notneg,decreasing=TRUE)
		cs1[o1] <- cumsum(isneg[o1])
		cs2[o2] <- cumsum(isneg[o2])
		DetectionPValue[,j] <- cs1+cs2
	}

#	Convert to p-values
	DetectionPValue/(2*nneg)
}
