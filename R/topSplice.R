topSplice <- function(fit, coef=ncol(fit), test="simes", number=10L, FDR=1, sort.by="p")
#	Collate diffSplice results into data.frame, ordered from most significant at top
#	Gordon Smyth and Yunshun Chen
#	Created 18 Dec 2013.  Last modified 29 Sep 2017.
{
#	Check fit is as produced by diffSplice
	if(is.null(fit$gene.genes$NExons)) stop("fit should be fit object produced by diffSplice")

#	Can only specify one coefficient
	coef <- coef[1]

	test <- match.arg(test,choices=c("simes","F","f","t"))
	if(test=="f") test <- "F"

	sort.by <- match.arg(sort.by,choices=c("p","none","logFC","NExons"))
	if(sort.by=="logFC" & test!="t") stop("Sorting by logFC only available with Simes test")
	if(sort.by=="NExons" & test=="t") stop("Sorting by NExons only available with gene-level tests")

#	Assemble data.frame of results for this coef
	switch(test,
	t = {
		out <- fit$genes
		out$logFC <- as.matrix(fit$coefficients)[,coef]
		out$t <- as.matrix(fit$t)[,coef]
		out$P.Value <- as.matrix(fit$p.value)[,coef]
	},

	F = {
		out <- fit$gene.genes
		out$F <- as.matrix(fit$gene.F)[,coef]
		out$P.Value <- as.matrix(fit$gene.F.p.value)[,coef]
	},

	simes = {
		out <- fit$gene.genes
		out$P.Value <- as.matrix(fit$gene.simes.p.value)[,coef]
	}
	)
	out$FDR <- p.adjust(out$P.Value, method="BH")

#	Reduce to significant genes
	if(FDR<1) out <- out[out$FDR <= FDR,]

#	Is the number of rows requested more than number available?
	number <- min(number, nrow(out))
	if(number <= 1L) return(out)

#	Sort rows
	o <- switch(sort.by,
			p = order(out$P.Value, decreasing=FALSE),
			logFC = order(abs(out$logFC), decreasing=TRUE),
			NExons = order(out$NExons, -out$P.Value, decreasing=TRUE),
			none = 1:nrow(out)
		)
	o <- o[1:number]
	out[o,]
}
