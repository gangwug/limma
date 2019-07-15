plotSplice <- function(fit, coef=ncol(fit), geneid=NULL, genecolname=NULL, rank=1L, FDR=0.05)
#	Plot exons of chosen gene
#	fit is output from diffSplice
#	Gordon Smyth, Yifang Hu and Yunshun Chen
#	Created 3 Jan 2014.  Last modified 28 Sep 2017.
{
	if(is.null(genecolname)) 
		genecolname <- fit$genecolname
	else
		genecolname <- as.character(genecolname)

	if(is.null(geneid)) {
#		Find gene from specified rank 
		if(rank==1L)
			i <- which.min(fit$gene.F.p.value[,coef])
		else
			i <- order(fit$gene.F.p.value[,coef])[rank]
		geneid <- paste(fit$gene.genes[i,genecolname], collapse=".")
	} else {
#		Find gene from specified name
		geneid <- as.character(geneid)
		i <- which(fit$gene.genes[,genecolname]==geneid)[1]
		if(!length(i)) stop(paste("geneid",geneid,"not found"))
	}

#	Row numbers containing exons
	j <- fit$gene.firstexon[i]:fit$gene.lastexon[i]

	exoncolname <- fit$exoncolname

#	Get strand if possible		
	strcol <- grepl("strand", colnames(fit$gene.genes), ignore.case=TRUE)
	if(any(strcol)) geneid <- paste0(geneid, " (", as.character(fit$gene.genes[i, strcol])[1], ")")

	if(is.null(exoncolname)) {
		plot(fit$coefficients[j, coef], xlab="Exon", ylab="logFC (this exon vs rest)", main=geneid, type="b")
	} else {
		exon.id <- fit$genes[j, exoncolname]
		xlab <- paste("Exon", exoncolname, sep=" ")
		plot(fit$coefficients[j, coef], xlab="", ylab="logFC (this exon vs rest)", main=geneid, type="b", xaxt="n")
		axis(1, at=1:length(j), labels=exon.id, las=2, cex.axis=0.5)
		mtext(xlab, side=1, padj=5.2)
	}

#	Mark the topSpliced exons
	top <- topSplice(fit, coef=coef, number=Inf, test="t", sort.by="none")
	m <- which(top[, genecolname] %in% fit$gene.genes[i, genecolname])
	fdr <- top$FDR[m]
	sig <- fdr < FDR
	if(any(sig)){
		fdr.sig <- fdr[sig]
		if(length(unique(fdr.sig))==1)
			cex <- 1.5
		else {
			abs.fdr <- abs(log10(fdr.sig))
			from <- range(abs.fdr)
			to <- c(1,2)
			cex <- (abs.fdr - from[1])/diff(from) * diff(to) + to[1]
		}
		points((1:length(j))[sig], fit$coefficients[j[sig], coef], col="red", pch=16, cex=cex)
	}

	abline(h=0,lty=2)
	invisible()
}
