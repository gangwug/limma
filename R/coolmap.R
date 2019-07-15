coolmap <- function(x, cluster.by="de pattern", col=NULL, linkage.row="complete", linkage.col="complete", show.dendrogram="both", ...)
#	Interface to heatmap.2 with useful defaults for log2-expression data
#	Gordon Smyth
#	Created 23 Sep 2016. Last modified 18 Sep 2018.
{
#	Check arguments
	x <- as.matrix(x)
	nsamples <- ncol(x)
	if(nsamples < 2L) stop("Need at least 2 rows and 2 columns")
	cluster.by <- match.arg(cluster.by, c("de pattern","expression level"))
	if(is.null(col)) {
		if(cluster.by=="de pattern") col <- "redblue" else col <- "yellowblue"
	} else {
		if(length(col)==1L) col <- match.arg(col, c("redblue","redgreen","yellowblue","whitered"))
	}
	show.dendrogram <- match.arg(show.dendrogram, c("both","row","column","none"))

#	Require gplots package
	suppressPackageStartupMessages(OK <- requireNamespace("gplots",quietly=TRUE))
	if(!OK) stop("gplots package required but not installed (or can't be loaded)")

#	Linkage method for rows (genes)
	linkage.row <- as.character(linkage.row)
	if(linkage.row %in% c("w","wa","war","ward")) linkage.row <- "ward.D2"
	METHODS <- c("none", "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
	linkage.row <- match.arg(linkage.row, METHODS)

#	Linkage method for columns (samples)
	linkage.col <- as.character(linkage.col)
	if(linkage.col %in% c("w","wa","war","ward")) linkage.col <- "ward.D2"
	linkage.col <- match.arg(linkage.col, METHODS)

#	Dendrogram for columns (samples) uses Euclidean distance without scaling
	if(linkage.col=="none") {
		hc <- FALSE
		if(show.dendrogram=="both") show.dendrogram <- "row"
		if(show.dendrogram=="column") show.dendrogram <- "none"
	} else {
		hc <- as.dendrogram(hclust(dist(t(x), method="euclidean"), method=linkage.col))
	}

#	Standardize rows if clustering by patterns
	if(cluster.by=="de pattern") {
		M <- rowMeans(x, na.rm=TRUE)
		DF <- nsamples - 1L
		IsNA <- is.na(x)
		if(any(IsNA)) {
			mode(IsNA) <- "integer"
			DF <-  DF - rowSums(IsNA)
			DF[DF==0L] <- 1L
		}
		x <- x-M
		V <- rowSums(x^2L, na.rm=TRUE) / DF
		x <- x / sqrt(V+0.01)
		sym <- TRUE
		key.xlab <- "Z-Score"
	} else {
		sym <- FALSE
		key.xlab <- "log2(expression)"
	}

#	Dendrogram for rows (genes) uses Euclidean distance
#	If rows are scaled, then this is equivalent to clustering using distance = 1-correlation.
	if(linkage.row=="none") {
		hr <- FALSE
		if(show.dendrogram=="both") show.dendrogram <- "column"
		if(show.dendrogram=="row") show.dendrogram <- "none"
	} else {
		hr <- as.dendrogram(hclust(dist(x, method="euclidean"), method=linkage.row))
	}

#	Set color pallate
	if(length(col)==1L) col <- switch(col,
		"redblue"=gplots::colorpanel(256,"blue2","white","red2"),
		"redgreen"=gplots::colorpanel(256,"green","black","red"),
		"yellowblue"=gplots::colorpanel(256,"blue2","white","yellow2"),
		"whitered"=gplots::colorpanel(256,low="white",high="red2")
	)

#	Plot heatmap
	gplots::heatmap.2(x, Rowv=hr, Colv=hc, scale="none", density.info="none", trace="none", col=col, symbreaks=sym, symkey=sym, dendrogram=show.dendrogram, key.title=NA, key.xlab=key.xlab, ...)
}