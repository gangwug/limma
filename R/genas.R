##  GENAS.R

genas <- function(fit,coef=c(1,2),subset="all",plot=FALSE,alpha=0.4)
#	Genuine association of gene expression profiles
#	Belinda Phipson and Gordon Smyth
#	21 September 2009. Last modified 11 April 2015.
{
	if(is.null(fit$cov.coefficients))
		if(is.null(fit$qr)) {
			stop("The fit object appears to have been fitted with observation weights or missing values. genas is not yet implemented for these situations.")
		} else {
			fit$cov.coefficients <- chol2inv(fit$qr$qr[,1:fit$qr$rank,drop=FALSE])
		}

	out <- list(
		technical.correlation=NA,
		covariance.matrix=matrix(NA,2,2),
		biological.correlation=NA,
		deviance=0,
		p.value=1,
		n=0
	)

#	Check fit
	if(nrow(fit)<1) {
		message("fit object has zero rows")
		return(out)
	}
	if(is.null(fit$s2.post)) fit <- eBayes(fit)
	trend <- (length(fit$s2.prior) > 1)
	robust <- (length(fit$df.prior) > 1)

#	Check coef
	if(length(coef) != 2) stop("coef should contain two column numbers")

#	Check subset
	if(subset=="n") subset <- "all"  # for backward compatibility
	subset <- match.arg(subset, c("all","Fpval","p.union","p.int","logFC","predFC"))

#	Keep only the two fit contrasts to be correlated
	if(ncol(fit)>2) fit <- fit[,coef]
	fit.plot <- fit

	x1 <- fitGammaIntercept(fit$coeff[,1]^2/fit$s2.post,offset=fit$cov.coeff[1,1])
	x2 <- fitGammaIntercept(fit$coeff[,2]^2/fit$s2.post,offset=fit$cov.coeff[2,2])
	if(x1 > 0 && x2 > 0) {
		v0null <- matrix(c(x1,0,0,x2),2,2)
		C <- chol(v0null)
		x <- log(diag(C))
	} else
		x <- c(0,0)
	m <- 2

#	Subset to genes that show some differential expression
	if(subset != "all") {
		fit <- .whichGenes(fit,subset)
		if(nrow(fit)<1) {
			message("No genes met criterion for inclusion in analysis")
			return(out)
		}
		fit <- eBayes(fit,trend=trend,robust=robust)
	}

	Q2 <- optim(x, .multTLogLikNull, fit = fit, m = m)
	Q1 <- optim(c(Q2$par[1], Q2$par[2], 0),.multTLogLik,fit=fit,m=m)
		
	L <- matrix(c(1,Q1$par[3],0,1),2,2)
	D <- matrix(c(exp(Q1$par[1]),0,0,exp(Q1$par[2])),2,2)
	V0 <- L%*%D%*%t(L)
	rhobiol <- V0[2,1]/sqrt(V0[1,1]*V0[2,2])

	V <- fit$cov.coefficients
	rhotech <- V[2,1]/sqrt(V[1,1]*V[2,2])

	if(plot) {
		if(!requireNamespace("ellipse",quietly=TRUE)) stop("ellipse package required but is not available")
		lim <- mean(c(sd(fit.plot$coeff[,1]),sd(fit.plot$coeff[,2])))
		if(nrow(fit)<500) lim <- 1.5*lim else lim <- 2*lim
		max <- max(c(fit.plot$coeff[,1],fit.plot$coeff[,2]))
		min <- min(c(fit.plot$coeff[,1],fit.plot$coeff[,2]))
		max <- sign(max)*max(abs(min),abs(max))
		min <- sign(min)*max(abs(min),abs(max))
		if(abs(rhobiol)>abs(rhotech)){
			plot(fit.plot$coeff[,1],fit.plot$coeff[,2],pch=16,cex=0.4,ylim=c(min,max),xlim=c(min,max),xlab=colnames(fit.plot$coeff)[1],ylab=colnames(fit.plot$coeff)[2])
			polygon(ellipse::ellipse(rhotech,scale=c(lim,lim)),col=rgb(0,0,1,alpha=alpha),border=rgb(0,0,1,alpha=alpha))
			polygon(ellipse::ellipse(rhobiol,scale=c(lim,lim)),col=rgb(0,1,0,alpha=alpha),border=rgb(0,1,0,alpha=alpha))
			points(fit.plot$coeff[,1],fit.plot$coeff[,2],pch=16,cex=0.4)
			abline(h=0,v=0)
			legend(min,max,legend=bquote(rho[biol]==.(round(rhobiol,3))),col=rgb(0,1,0,alpha=alpha),pch=16,bty="n",cex=0.8)
			legend(min,max-lim/2,legend=bquote(rho[tech]==.(round(rhotech,3))),col=rgb(0,0,1,alpha=alpha),pch=16,bty="n",cex=0.8)
		} else {
			plot(fit.plot$coeff[,1],fit.plot$coeff[,2],pch=16,cex=0.4,ylim=c(min,max),xlim=c(min,max),xlab=colnames(fit.plot$coeff)[1],ylab=colnames(fit.plot$coeff)[2])
			polygon(ellipse::ellipse(rhobiol,scale=c(lim,lim)),col=rgb(0,1,0,alpha=alpha),border=rgb(0,1,0,alpha=alpha))
			polygon(ellipse::ellipse(rhotech,scale=c(lim,lim)),col=rgb(0,0,1,alpha=alpha),border=rgb(0,0,1,alpha=alpha))
			points(fit.plot$coeff[,1],fit.plot$coeff[,2],pch=16,cex=0.4)
			abline(h=0,v=0)
			legend(min,max,legend=bquote(rho[biol]==.(round(rhobiol,3))),col=rgb(0,1,0,alpha=alpha),pch=16,bty="n",cex=0.8)
			legend(min,max-lim/2,legend=bquote(rho[tech]==.(round(rhotech,3))),col=rgb(0,0,1,alpha=alpha),pch=16,bty="n",cex=0.8)
		}
	}

	D <- abs(2*(Q2$value-Q1$value))
	p.val <- pchisq(D,df=1,lower.tail=FALSE)

	list(technical.correlation=rhotech,covariance.matrix=V0,biological.correlation=rhobiol,deviance=D,p.value=p.val,n=nrow(fit))
}

.multTLogLikNull <- function(x,fit,m) 
#	Calculate the log-likelihood under the null hypothesis of no biological correlation
#	Belinda Phipson and Gordon Smyth
#	21 September 2009. Last modified 2 December 2013.
{
	df.total <- fit$df.total
	s <- fit$s2.post
	B <- fit$coefficients
	V <- fit$cov.coefficients
	a1 <- x[1]
	a2 <- x[2]
	chol <- matrix(c(exp(a1),0,0,exp(a2)),2,2)
	V0 <- t(chol) %*% chol
	
	R <- chol(V0+V)
	Second <- sum(log(diag(R)))

	W <- backsolve(R,t(B),transpose=TRUE)
	Q <- colSums(W^2)

	Third <- 0.5*(m+df.total)*log(1+Q/s/df.total)
	
	sum(Second+Third)
}


.multTLogLik <- function(x,fit,m) 
#	Calculate the log-likelihood with biological correlation
#	Belinda Phipson and Gordon Smyth
#	21 September 2009. Last modified 2 December 2013.
{	
	df.total <- fit$df.total
	s <- fit$s2.post
	B <- fit$coefficients
	V <- fit$cov.coefficients
	a1 <- x[1]
	a2 <- x[2]
	b <- x[3]
	L <- matrix(c(1,b,0,1),2,2)
	D <- matrix(c(exp(a1),0,0,exp(a2)),2,2)
	V0 <- L %*% D %*% t(L)

	R <- chol(V0+V)
	Second <- sum(log(diag(R)))

	W <- backsolve(R,t(B),transpose=TRUE)
	Q <- colSums(W^2)

	Third <- 0.5*(m+df.total)*log(1+Q/s/df.total)

	sum(Second+Third)
}


.whichGenes <- function(fit,subset)
{
	if(subset=="all") return(genes)

	if(subset=="Fpval") {
		p <- 1-propTrueNull(fit$F.p.value)
		R <- rank(fit$F.p.value)
		cut <- p*nrow(fit)
		genes <- (R <= cut)
	}

	if(subset=="p.union"){
		p1 <- 1-propTrueNull(fit$p.value[,1])
		p2 <- 1-propTrueNull(fit$p.value[,2])
		cut1 <- p1*nrow(fit)
		cut2 <- p2*nrow(fit)
		if(p1==0 & p2==0){ 
			genes <- FALSE
		} else {
			R1 <- rank(fit$p.value[,1])
			R2 <- rank(fit$p.value[,2])
			genes <- R1 <= cut1 | R2 <= cut2
		}
	}

	if(subset=="p.int"){
		p1 <- 1-propTrueNull(fit$p.value[,1])
		p2 <- 1-propTrueNull(fit$p.value[,2])
		R1 <- rank(fit$p.value[,1])
		R2 <- rank(fit$p.value[,2])
		cut1 <- p1*nrow(fit)
		cut2 <- p2*nrow(fit)
		genes <- R1 <= cut1 & R2 <= cut2
	}
 
 	if(subset=="logFC") {
		q1 <- quantile(abs(fit$coeff[,1]),probs=0.9)
		q2 <- quantile(abs(fit$coeff[,2]),probs=0.9)
		genes <- abs(fit$coeff[,1]) > q1 | abs(fit$coeff[,2]) > q2
		fit$coeff[,1] <- sign(fit$coeff[,1]) * (abs(fit$coeff[,1])-q1)
		fit$coeff[,2] <- sign(fit$coeff[,2]) * (abs(fit$coeff[,2])-q2)
	}
 
 	if(subset=="predFC") {
		pfc1 <- predFCm(fit, coef=1)
		pfc2 <- predFCm(fit, coef=2)
		q1 <- quantile(abs(pfc1),probs=0.9)
		q2 <- quantile(abs(pfc2),probs=0.9)
		genes <- abs(pfc1) > q1 | abs(pfc2) > q2
		fit$coeff[,1] <- pfc1
		fit$coeff[,2] <- pfc2
		fit$coeff[,1] <- sign(fit$coeff[,1]) * (abs(fit$coeff[,1])-q1)
		fit$coeff[,2] <- sign(fit$coeff[,2]) * (abs(fit$coeff[,2])-q2)
	}

	fit[genes,]
}

