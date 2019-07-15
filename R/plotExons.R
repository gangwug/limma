plotExons <- function(fit, coef = ncol(fit), geneid = NULL, genecolname = "GeneID", exoncolname = NULL, rank = 1L, FDR = 0.05)
#	Plot log2 fold changes of the exons and mark the differentially expressed exons
#	Yifang Hu and Gordon Smyth.
#	Created 1 May 2014. Last modified 7 October 2014.
{

	# Check input
	if(!is(fit, "MArrayLM")) stop("fit is not MArrayLM object")
	if(is.null(fit$p.value)) stop("fit object should contain p value from eBayes")
	if(is.null(fit$method)) stop("fit object should have fitting method from least square or robust regression")

	genecolname <- as.character(genecolname)[1]
	if(! (genecolname %in% colnames(fit$genes))) stop(paste("genecolname", genecolname, "not found"))

	if( ! is.null(exoncolname)) {

		exoncolname <- as.character(exoncolname)[1]
		if( ! (exoncolname %in% colnames(fit$genes))) stop(paste("exoncolname", exoncolname, "not found"))

	}

	# Gene to plot using rank
	if (is.null(geneid)) {

		if (rank == 1L) igene <- which.min(fit$p.value[, coef])

		else {

			o.p <- order(fit$p.value[, coef])
			fit.o <- fit[o.p,]
			fit.uniq <- fit.o[!duplicated(fit.o$genes[, genecolname]),]
			igene <- match(fit.uniq$genes[rank, genecolname], fit$genes[, genecolname])

		}

	# Gene to plot using geneid
	} else {

		geneid <- as.character(geneid)
		igene <- match(geneid[1], as.character(fit$genes[, genecolname]))
		if (is.na(igene)) stop(paste("geneid", geneid, "not found"))

	}

	# Gene annotation
	ilab <- grep(paste0(genecolname, "|geneid|symbol"), colnames(fit$genes), ignore.case = TRUE)
	genecollab <- colnames(fit$genes)[ilab]
	genelab <- paste(sapply(fit$genes[igene,genecollab], as.character), collapse = ", ")

	# Get strand if possible
	strcol <- grepl("strand", colnames(fit$genes), ignore.case = TRUE)
	if(any(strcol)) genelab <- paste0(genelab, " (", as.character(fit$genes[igene, strcol])[1], ")")

	# Exon level DE results
	fdr <- p.adjust(fit$p.value[, coef], method = "BH")
	de <- data.frame(fit$genes, logFC = fit$coefficients[, coef], fdr)
	m <- which(de[, genecolname] %in% fit$genes[igene, genecolname])
	exon.de <- de[m,]

	# Add exon identifier
	if(is.null(exoncolname)) {

		exoncolname <- "ID"
		exon.de$ID <- 1:nrow(exon.de)

	}

	# Order exon level results by exon identifier
	exon.de <- exon.de[order(exon.de[,exoncolname]),]

	# Plot exons
	plot(exon.de$logFC, xlab = "", ylab = "Log2 Fold Change", main = genelab, type = "b", xaxt = "n")
	axis(1, at = 1:nrow(exon.de), labels = exon.de[, exoncolname], las = 2, cex.axis = 0.5)
	mtext(text = paste("Exon", exoncolname), side = 1, padj = 5.2)

	# Mark DE exons
	mark <- exon.de$fdr < FDR

	if (sum(mark) > 0) {

		col <- rep(NA, nrow(exon.de))
		col[mark & (exon.de$logFC > 0)] <- "red"
		col[mark & (exon.de$logFC < 0)] <- "blue"

		if (sum(mark) == 1) cex <- 1.5 else {

			abs.fdr <- abs(log10(exon.de$fdr[mark]))
			from <- range(abs.fdr)
			to <- c(1, 2)
			cex <- (abs.fdr - from[1])/diff(from) * diff(to) + to[1]

		}

		points((1:nrow(exon.de))[mark], exon.de$logFC[mark], col = col[mark], pch = 16, cex = cex)

	}

	abline(h = mean(exon.de$logFC), lty = 2)

}
