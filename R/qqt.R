qqt <- function(y,df=Inf,ylim=range(y),main="Student's t Q-Q Plot",xlab="Theoretical Quantiles",ylab="Sample Quantiles",plot.it=TRUE,...)
#	Student's t probability plot
#	Gordon Smyth
#	3 Oct 2002
{
	y <- y[!is.na(y)]
	if(0 == (n <- length(y))) stop("y is empty")
	x <- qt(ppoints(n),df=df)[order(order(y))]
	if (plot.it) plot(x,y,main=main,xlab=xlab,ylab=ylab,ylim=ylim,...)
	invisible(list(x=x,y=y))
}

qqf <- function(y, df1, df2, ylim=range(y), main="F Distribution Q-Q Plot", xlab="Theoretical Quantiles", ylab="Sample Quantiles", plot.it=TRUE,...)
#	F probability plot
#	Belinda Phipson
#	14 May 2012
{
	y <- y[!is.na(y)]
	if(0 == (n <- length(y))) stop("y is empty")
	x <- qf(ppoints(n),df1=df1,df2=df2)[order(order(y))]
	if (plot.it) plot(x,y,main=main,xlab=xlab,ylab=ylab,ylim=ylim,...)
	invisible(list(x=x,y=y))
}
