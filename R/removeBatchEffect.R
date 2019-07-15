#  removeBatchEffect.R

#  A refinement would be to empirical Bayes shrink
#  the batch effects before subtracting them.

removeBatchEffect <- function(x,batch=NULL,batch2=NULL,covariates=NULL,design=matrix(1,ncol(x),1),...)
#  Remove batch effects from matrix of expression data
#  Gordon Smyth and Carolyn de Graaf
#  Created 1 Aug 2008. Last revised 1 June 2014.
{
	if(is.null(batch) && is.null(batch2) && is.null(covariates)) return(as.matrix(x))
	if(!is.null(batch)) {
		batch <- as.factor(batch)
		contrasts(batch) <- contr.sum(levels(batch))
		batch <- model.matrix(~batch)[,-1,drop=FALSE]
	}
	if(!is.null(batch2)) {
		batch2 <- as.factor(batch2)
		contrasts(batch2) <- contr.sum(levels(batch2))
		batch2 <- model.matrix(~batch2)[,-1,drop=FALSE]
	}
	if(!is.null(covariates)) covariates <- as.matrix(covariates)
	X.batch <- cbind(batch,batch2,covariates)
	fit <- lmFit(x,cbind(design,X.batch),...)
	beta <- fit$coefficients[,-(1:ncol(design)),drop=FALSE]
	beta[is.na(beta)] <- 0
	as.matrix(x) - beta %*% t(X.batch)
}
