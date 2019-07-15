voomWithQualityWeights <- function(counts, design=NULL, lib.size=NULL, normalize.method="none", plot=FALSE, span=0.5, var.design=NULL, method="genebygene", maxiter=50, tol=1e-10, trace=FALSE, col=NULL, ...)
#	Combine voom weights with sample-specific weights estimated by arrayWeights() function for RNA-seq data
#	Matt Ritchie, Cynthia Liu, Gordon Smyth
#	Created 22 Sept 2014. Last modified 21 Dec 2017.
{
	if(plot) {
		oldpar <- par(mfrow=c(1,2))
		on.exit(par(oldpar))
	}
	v <- voom(counts, design=design, lib.size=lib.size, normalize.method=normalize.method, plot=FALSE, span=span, ...)
	aw <- arrayWeights(v, design=design, method=method, maxiter=maxiter, tol=tol, var.design=var.design)
	v <- voom(counts, design=design, weights=aw, lib.size=lib.size, normalize.method=normalize.method, plot=plot, span=span, ...)
	aw <- arrayWeights(v, design=design, method=method, maxiter=maxiter, tol=tol, trace=trace, var.design=var.design)

	v$weights <- t(aw * t(v$weights))
	v$targets$sample.weights <- aw

	if(plot) {
		barplot(aw, names=1:length(aw), main="Sample-specific weights", ylab="Weight", xlab="Sample", col=col)
		abline(h=1, col=2, lty=2)
	}

	v		
}

