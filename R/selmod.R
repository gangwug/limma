##  selmod.R

selectModel <- function(y, designlist, criterion="aic", df.prior=0, s2.prior=NULL, s2.true=NULL, ...)
# y is a data matrix to be fitted, with rows as genes and columns as arrays.
# designlist is a list of design matrices to be compared.
# The function returns AIC or BIC values for each design
# and the prefered model for each gene with minimim criterion value.
# Written 17/7/08 Alicia Oshlack
# Last revised 2 Oct 2008, Gordon Smyth
{
	ym <- as.matrix(y)
	if(any(is.na(ym))) stop("NAs not allowed")
	narrays <- ncol(ym)
	rm(ym)

	nmodels <- length(designlist)
	models <- names(designlist)
	if(is.null(models)) models <- as.character(1:nmodels)
	if(df.prior>0 & is.null(s2.prior)) stop("s2.prior must be set")
	if(df.prior==0) s2.prior <- 0
	criterion <- match.arg(criterion,c("aic","bic","mallowscp"))

	if(criterion=="mallowscp") {
		if(is.null(s2.true)) stop("Need s2.true values")
		for(i in 1:nmodels) {
			fit <- lmFit(y, designlist[[i]], ...)
			npar <- narrays-fit$df.residual[1]
			if(i==1) {
				IC <- matrix(nrow=nrow(fit),ncol=nmodels,dimnames=list(Probes=rownames(fit),Models=models))
				if(length(s2.true)!=nrow(fit) && length(s2.true)!=1) stop("s2.true wrong length")
			}
			IC[,i] <- fit$df.residual*fit$sigma^2/s2.true+npar*2-narrays
		}
	} else {
		ntotal <- df.prior+narrays
		penalty <- switch(criterion,bic=log(narrays),aic=2)
		for(i in 1:nmodels) {
			fit <- lmFit(y, designlist[[i]], ...)
			npar <- narrays-fit$df.residual[1]+1
			s2.post <- (df.prior*s2.prior+fit$df.residual*fit$sigma^2)/ntotal
			if(i==1) IC <- matrix(nrow=nrow(fit),ncol=nmodels,dimnames=list(Probes=rownames(fit),Models=models))
			IC[,i] <- ntotal*log(s2.post)+npar*penalty
		}
	}

	pref <- factor(apply(IC,1,which.min),levels=1:nmodels,labels=models)
	list(IC=IC,pref=pref,criterion=criterion)
}

