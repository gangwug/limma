convest <- function(p,niter=100,plot=FALSE,report=FALSE,file="",tol=1e-06)
# Estimates pi0 using a convex decreasing density estimate
# Input: p=observed p-values,k=number of iterations,	
# plot=TRUE plots and report=TRUE writes results of each iteration.	
# Returns: An estimate of pi0 
# The methology underlying the function can be found in the
# preprint: "Estimating the proportion of true null hypotheses,
# with application to DNA microarray data." that can be downloaded
# from http://www.math.ntnu.no/~mettela/
#	Written by Egil Ferkingstad.
#	Received from Mette Langaas 26 Jun 2004.
#	Modified for limma by Gordon Smyth, 29 Oct 2004, 28 May 2005.
#	Report modified by Marcus Davy, 24 June 2007. Implemented and edited by Gordon Smyth, 9 Sep 2012.
{
	if(!length(p)) return(NA)
	if(any(is.na(p))) stop("Missing values in p not allowed")
	if(any(p<0 | p>1)) stop("All p-values must be between 0 and 1")
	k <- niter
#	accuracy of the bisectional search for finding a new
#	convex combination of the current iterate and the mixing density
	ny <- tol
	p <- sort(p)
	m <- length(p)
	p.c <- ceiling(100*p)/100
	p.f <- floor(100*p)/100
	t.grid <- (1:100)/100
	x.grid <- (0:100)/100
	t.grid.mat <- matrix(t.grid,ncol=1)
	f.hat <- rep(1,101) #f.hat at the x-grid
	f.hat.p <- rep(1,m) #f.hat at the p-values
	theta.hat <- 0.01*which.max(apply(t.grid.mat,1,function(theta) sum((2*(theta-p)*(p<theta)/theta^2))))
	f.theta.hat <- 2*(theta.hat-x.grid)*(x.grid<theta.hat)/theta.hat^2 # f.theta.hat at the x-grid
	f.theta.hat.p <- 2*(theta.hat-p)*(p<theta.hat)/theta.hat^2 # f.theta.hat at the p-values
	i<-1
	j<-0
	thetas <- numeric()
	if(report) cat("j\tpi0\ttheta.hat\t\tepsilon\tD\n", file=file, append=FALSE)
	for (j in 1:k) {
		f.hat.diff <- (f.hat.p - f.theta.hat.p)
		if (sum(f.hat.diff/f.hat.p) > 0)
			eps <- 0
		else {
			l <- 0
			u <- 1
			while (abs(u-l)>ny) {
				eps <- (l+u)/2
				if(sum((f.hat.diff/((1-eps)*f.hat.p+eps*f.theta.hat.p))[f.hat.p>0])<0)
					l <- eps
				else
					u <- eps
			}
		}
		f.hat <- (1-eps)*f.hat + eps*f.theta.hat
		pi.0.hat <- f.hat[101]
		d <- sum(f.hat.diff/f.hat.p)
		if(report) cat(j, "\t",pi.0.hat, "\t",theta.hat,"\t",eps, "\t",d, "\n", file=file, append=TRUE)
		f.hat.p <- 100*(f.hat[100*p.f+1]-f.hat[100*p.c+1])*(p.c-p)+f.hat[100*p.c+1]
		theta.hat <- 0.01*which.max(apply(t.grid.mat,1,function(theta) sum((2*(theta-p)*(p<theta)/theta^2)/f.hat.p)))
		f.theta.hat <- 2*(theta.hat-x.grid)*(x.grid<theta.hat)/theta.hat^2
		f.theta.hat.p <- 2*(theta.hat-p)*(p<theta.hat)/theta.hat^2
		if (sum(f.theta.hat.p/f.hat.p)<sum(1/f.hat.p)) { # check if the Unif[0,1]-density is the new "f.theta.hat"
			theta.hat <- 0
			f.theta.hat <- rep(1,101)
			f.theta.hat.p <- rep(1,m)
		}
		if (sum(thetas==theta.hat)==0) {
			thetas[i] <- theta.hat
			thetas <- sort(thetas)
			i <- i + 1
		}
#		pi.0.hat <- f.hat[101]
		if (plot) {
			plot(x.grid,f.hat,type="l",main=paste(format(round(pi.0.hat,5),digits=5)),ylim=c(0,1.2))
			points(thetas,f.hat[100*thetas+1],pch=20,col = "blue")
		}
	}
	pi.0.hat
}

