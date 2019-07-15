#  BACKGROUND-NORMEXP.R

normexp.signal <- function(par,x)
#	Expected value of signal given foreground in normal + exponential model
#	Gordon Smyth
#	24 Aug 2002. Last modified 24 February 2012.
{
	mu <- par[1]
	sigma <- exp(par[2])
	sigma2 <- sigma*sigma
	alpha <- exp(par[3])
#	cat(c(mu,sigma,alpha),"\n")
	if(alpha <= 0) stop("alpha must be positive")
	if(sigma <= 0) stop("sigma must be positive")
	mu.sf <- x-mu-sigma2/alpha
	signal <- mu.sf + sigma2 * exp(dnorm(0,mean=mu.sf,sd=sigma,log=TRUE) - pnorm(0,mean=mu.sf,sd=sigma,lower.tail=FALSE,log.p=TRUE))
	o <- !is.na(signal)
	if(any(signal[o]<0)) {
		warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
		signal[o] <- pmax(signal[o],1e-6)
	}
	signal
}

normexp.fit <- function(x, method="saddle", n.pts=NULL, trace=FALSE)
#	Estimate parameters of normal+exponential convolution model
#	Pure R Version Gordon Smyth 24 Aug 2002.
#	Version with C by Jeremy Silver 29 Oct 2007.
#	Last modified 14 January 2015.
{
	isna <- is.na(x)
	if(any(isna)) x <- x[!isna]
	if(length(x)<4) stop("Not enough data: need at least 4 non-missing corrected intensities")
	if(trace) cat("trace not currently implemented\n")

	method <- match.arg(method,c("mle","saddle","rma","rma75","mcgee","nlminb","nlminblog"))
#	Backward compatility with old names
	if(method=="mcgee") method <- "rma75"
	if(method=="nlminb") method <- "mle"
	if(method=="nlminblog") method <- "mle"

	if(method=="rma") {
		if(!requireNamespace("affy",quietly=TRUE)) stop("affy package required but is not available")
		out <- affy::bg.parameters(x)
		return(list(par=c(out$mu,log(out$sigma),-log(out$alpha))))
	}

	if(method=="rma75") {
		out <- .bg.parameters.rma75(x)
		return(list(par=c(out$mu,log(out$sigma),-log(out$alpha))))
	}

#	Starting values for parameters mu, alpha and sigma
	q <- quantile(x, c(0,0.05,0.1,1), na.rm = TRUE, names = FALSE)
	if(q[1]==q[4]) return(list(par=c(q[1],-Inf,-Inf),m2loglik=NA,convergence=0))
	if(q[2] > q[1]) {
		mu <- q[2]
	} else {
		if(q[3] > q[1]) {
			mu <- q[3]
		} else {
			mu <- q[1] + 0.05*(q[4]-q[1])
		}
	}
	sigma2 <- mean((x[x<mu]-mu)^2, na.rm = TRUE)
	alpha <- mean(x,na.rm = TRUE) - mu
	if(alpha <= 0) alpha <- 1e-6
	Par0 <- c(mu,log(sigma2)/2,log(alpha))

#	Use a maximum of n.pts points for the fit
	if(!is.null(n.pts)) if(n.pts >= 4 & n.pts < length(x)) {
		a <- 0.5
		x <- quantile(x,((1:n.pts)-a)/n.pts,type=5)
	}

#	Maximize saddlepoint approximation to likelihood
	out1 <- .C("fit_saddle_nelder_mead",
		par = as.double(Par0), 
		X = as.double(x), 
		N = as.integer(length(x)), 
		convergence = as.integer(0), 
		fncount = as.integer(0), 
		m2loglik = as.double(0),
		PACKAGE="limma")
	out1$X <- out1$N <- NULL

	if(method=="saddle") return(out1)
	
	Par1 <- out1$par
#	Convert from log-sd to log-var parametrization
	Par1[2] <- 2*Par1[2]
	LL1 <- .normexp.m2loglik(Par1, f = x)
	out2 <- nlminb(start = Par1,
		objective = .normexp.m2loglik,
		gradient = .normexp.gm2loglik,
		hessian = .normexp.hm2loglik,
		f = x,
		scale = median(abs(Par1))/abs(Par1))
#	Convert back to log-sd parametrization
	out2$par[2] <- out2$par[2]/2
	out2$m2loglik <- out2$objective
	out2$objective <- NULL

#	Check whether nlminb helped
	if(out2$m2loglik >= LL1) return(out1)
	out2
}

.normexp.m2loglik <- function(theta,f)
#	normexp minus-twice log-likelihood
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 25 Sept 2008.
{
  mu <- theta[1]
  s2 <- exp(theta[2])
  al <- exp(theta[3])
  
  .C("normexp_m2loglik",
    mu = as.double(mu), 
    s2 = as.double(s2), 
    al = as.double(al), 
    n = as.integer(length(f)), 
    f = as.double(f), 
    m2LL = double(1),
    PACKAGE="limma"
  )$m2LL
}

.normexp.gm2loglik <- function(theta,f)
#	Gradient of normexp m2loglik
#	with respect to mu, log(sigma^2) and log(alpha)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 25 Sept 2008.
{
  mu <- theta[1]
  s2 <- exp(theta[2])
  al <- exp(theta[3])

  .C("normexp_gm2loglik",
    mu = as.double(mu), 
    s2 = as.double(s2), 
    al = as.double(al), 
    n = as.integer(length(f)), 
    f = as.double(f), 
    dm2LL = double(3),
    PACKAGE = "limma"
  )$dm2LL

}

.normexp.hm2loglik <- function(theta,f)
#	Hessian of normexp m2loglik
#	with respect to mu, log(sigma^2) and log(alpha)
#  Jeremy Silver
#  29 Oct 2007.
#	Last modified 25 Sept 2008.
{
  mu <- theta[1]
  s2 <- exp(theta[2])
  al <- exp(theta[3])

  matrix(.C("normexp_hm2loglik",
    mu = as.double(mu), 
    s2 = as.double(s2), 
    al = as.double(al), 
    n = as.integer(length(f)), 
    f = as.double(f), 
    d2m2LL = double(9),
    PACKAGE="limma"
  )$d2m2LL,3,3)
}


.bg.parameters.rma75 <- function(pm,n.pts = 2^14)
#	Estimate normexp parameters
#	This code is extracted without alteration from the RMA-75 function of
#	McGee, M. and Chen, Z. (2006). Parameter estimation for the
#	exponential-normal convolution model for background correction
#	of Affymetrix GeneChip data.
#	Stat Appl Genet Mol Biol, 5(1), Article 24.
{
##	mu-correction function
	mu.est.correct <- function(m,s,a) { 
		f <- function(x) (dnorm(x-s*a)-s*a*(pnorm(x-s*a)+pnorm(m/s+s*a)-1))
		t <- uniroot(f, c(-5, 10), tol = 1e-12)$root
		t <- m-s*t
		return(t)
	}

##	getting mode function
	max.density <- function(x, n.pts) {
		aux <- density(x, kernel = "epanechnikov", n = n.pts, na.rm = TRUE)
		aux$x[order(-aux$y)[1]]
	}

	pmbg <- max.density(pm, n.pts)
	bg.data <- pm[pm < pmbg]		   
	pmbg <- max.density(bg.data, n.pts) 
	mubg <- pmbg   ## the mode
	bg.data <- pm[pm < pmbg]
	bg.data <- bg.data - pmbg
	bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2) ## estimate sigma
	sig.data<-pm[pm > pmbg]
	sig.data <- sig.data - pmbg   
	q75 <- 0.75
	alpha3 <- -(quantile(pm,q75)-pmbg)/log(1-q75) ## 75th quantile estimation

##	mode-correction
	mu3 <- mu.est.correct(m=mubg,s=bgsd,a=1/alpha3)
	mu3 <- (mu3+mubg)/2  ## take ave
	bg.data3<- pm[pm < mu3] 
	bg.data3 <- bg.data3 - mu3
	bgsd3 <- sqrt(sum(bg.data3^2)/(length(bg.data3) - 1)) * sqrt(2)
	sig.data3 <- pm[pm > mu3]
	alpha3<- -(quantile(pm,q75)-mu3)/log(1-q75)
	list(alpha = 1/alpha3, mu = mu3, sigma = bgsd3)
}

