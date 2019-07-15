kegga <- function(de,...) UseMethod("kegga")

kegga.MArrayLM <- function(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, ...)
#	KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway analysis of DE genes from linear model fit
#	Gordon Smyth and Yifang Hu
#	Created 15 May 2015. Last modified 3 June 2015.
{
#	Avoid argument collision with default method
	dots <- names(list(...))
	if("universe" %in% dots) stop("kegga.MArrayLM defines its own universe",call.=FALSE)
	if((!is.logical(trend) || trend) && "covariate" %in% dots) stop("kegga.MArrayLM defines it own covariate",call.=FALSE)

#	Check fit
	if(is.null(de$coefficients)) stop("de does not appear to be a valid MArrayLM fit object.")
	if(is.null(de$p.value)) stop("p.value not found in fit object, perhaps need to run eBayes first.")	
	if(length(coef) != 1) stop("Only one coef can be specified.")
	ngenes <- nrow(de)

#	Check geneid
#	Can be either a vector of gene IDs or an annotation column name
	geneid <- as.character(geneid)
	if(length(geneid) == ngenes) {
		universe <- geneid
	} else {
		if(length(geneid) == 1L) {
			universe <- de$genes[[geneid]]
			if(is.null(universe)) stop("Column ",geneid," not found in de$genes")
		} else
			stop("geneid of incorrect length")
	}

#	Check trend
#	Can be logical, or a numeric vector of covariate values, or the name of the column containing the covariate values
	if(is.logical(trend)) {
		if(trend) {
			covariate <- de$Amean
			if(is.null(covariate)) stop("Amean not found in fit")
		}
	} else {
		if(is.numeric(trend)) {
			if(length(trend) != ngenes) stop("If trend is numeric, then length must equal nrow(de)")
			covariate <- trend
			trend <- TRUE
		} else {
			if(is.character(trend)) {
				if(length(trend) != 1L) stop("If trend is character, then length must be 1")
				covariate <- de$genes[[trend]]
				if(is.null(covariate)) stop("Column ",trend," not found in de$genes")
				trend <- TRUE
			} else
				stop("trend is neither logical, numeric nor character")
		}
	}

#	Check FDR
	if(!is.numeric(FDR) | length(FDR) != 1) stop("FDR must be numeric and of length 1.")
	if(FDR < 0 | FDR > 1) stop("FDR should be between 0 and 1.")

#	Get up and down DE genes
	fdr.coef <- p.adjust(de$p.value[,coef], method = "BH")
	EG.DE.UP <- universe[fdr.coef < FDR & de$coef[,coef] > 0]
	EG.DE.DN <- universe[fdr.coef < FDR & de$coef[,coef] < 0]
	DEGenes <- list(Up=EG.DE.UP, Down=EG.DE.DN)

#	If no DE genes, return data.frame with 0 rows
	if(length(EG.DE.UP)==0 && length(EG.DE.DN)==0) {
		message("No DE genes")
		return(data.frame())
	}

	if(trend)
		kegga(de=DEGenes, universe = universe, covariate=covariate, ...)
	else
		kegga(de=DEGenes, universe = universe, ...)
}

kegga.default <- function(de, universe=NULL,  restrict.universe=FALSE, species="Hs", species.KEGG=NULL, convert=FALSE, gene.pathway=NULL, pathway.names = NULL,prior.prob=NULL, covariate=NULL, plot=FALSE, ...)
#	KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway analysis of DE genes
#	Gordon Smyth and Yifang Hu
#	Created 18 May 2015.  Modified 19 May 2019.
{
#	Ensure de is a list
	if(is.list(de)) {
		if(is.data.frame(de)) stop("de should be a list of character vectors. It should not be a data.frame.")
	} else {
		de <- list(DE = de)
	}
	nsets <- length(de)

#	Stop if components of de are not vectors
	if(!all(vapply(de,is.vector,TRUE))) stop("components of de should be vectors")

#	Ensure de gene IDs are unique and of character mode
	for (s in 1:nsets) de[[s]] <- unique(as.character(de[[s]]))

#	Ensure de components have unique names
	names(de) <- trimWhiteSpace(names(de))
	NAME <- names(de)
	i <- which(NAME == "" | is.na(NAME))
	NAME[i] <- paste0("DE",i)
	names(de) <- makeUnique(NAME)

#	Select species
	if(is.null(species.KEGG)) {
		species <- match.arg(species, c("Ag", "At", "Bt", "Ce", "Dm", "Dr", "EcK12", "EcSakai", "Gg", "Hs", "Mm", "Mmu", "Pf", "Pt", "Rn", "Ss", "Xl"))
#		Convert from Bioconductor to KEGG species codes
		species.KEGG <- switch(species, "Ag"="aga", "At"="ath", "Bt"="bta", "Ce"="cel", "Cf"="cfa", "Dm"="dme", "Dr"="dre", "EcK12"="eco", "EcSakai"="ecs", "Gg"="gga", "Hs"="hsa", "Mm"="mmu", "Mmu"="mcc", "Pf"="pfa", "Pt"="ptr", "Rn"="rno", "Ss"="ssc", "Xl"="xla")
	}

#	Get pathway annotation
	if(is.null(gene.pathway))
		GeneID.PathID <- getGeneKEGGLinks(species.KEGG, convert=convert)
	else {
		GeneID.PathID <- gene.pathway
		d <- dim(GeneID.PathID)
		if(is.null(d)) stop("gene.pathway must be data.frame or matrix")
		if(d[2] < 2) stop("gene.pathway must have at least 2 columns")
		isna <- rowSums(is.na(GeneID.PathID[,1:2])) > 0.5
		GeneID.PathID <- GeneID.PathID[!isna,]

#		Remove repeated rows
		ID.ID <- paste(GeneID.PathID[,1],GeneID.PathID[,2],sep=".")
		if(anyDuplicated(ID.ID)) {
			d <- duplicated(ID.ID)
			GeneID.PathID <- GeneID.PathID[!d,]
		}
	}

#	Get pathway names
	if(is.null(pathway.names))
		PathID.PathName <- getKEGGPathwayNames(species.KEGG, remove.qualifier=TRUE)
	else {
		PathID.PathName <- pathway.names
		d <- dim(PathID.PathName)
		if(is.null(d)) stop("pathway.names must be data.frame or matrix")
		if(d[2] < 2) stop("pathway.names must have at least 2 columns")
		isna <- rowSums(is.na(PathID.PathName[,1:2])) > 0.5
		PathID.PathName <- PathID.PathName[!isna,]
	}

#	Universe defaults to all annotated genes
#	prior.prob and covariate must have same length as universe
#	Ensure universe unique
	if(is.null(universe)) {
		universe <- unique(GeneID.PathID[,1])
		prior.prob <- covariate <- NULL
	} else {
		universe <- as.character(universe)
		lu <- length(universe)
		if(!lu) stop("No genes in universe")
		if(!is.null(prior.prob) && length(prior.prob)!=lu) stop("universe and prior.prob must have same length")
		if(!is.null(covariate) && length(covariate)!=lu) stop("universe and covariate must have same length")
		if(restrict.universe) {
			i <- universe %in% GeneID.PathID[,1]
			universe <- universe[i]
			if(!is.null(prior.prob)) prior.prob <- prior.prob[i]
			if(!is.null(covariate)) covariate <- covariate[i]
		}
	}

#	Consolidate any replicated entries in universe, averaging corresponding prior.probs
	if(anyDuplicated(universe)) {
		d <- duplicated(universe)
		if(!is.null(covariate)) {
			covariate <- rowsum(covariate,group=universe,reorder=FALSE)
			n <- rowsum(rep_len(1L,length(universe)),group=universe,reorder=FALSE)
			covariate <- covariate/n
		}
		if(!is.null(prior.prob)) {
			prior.prob <- rowsum(prior.prob,group=universe,reorder=FALSE)
			n <- rowsum(rep_len(1L,length(universe)),group=universe,reorder=FALSE)
			prior.prob <- prior.prob/n
		}
		universe <- universe[!d]
	}
	NGenes <- length(universe)
	if(NGenes<1L) stop("No annotated genes found in universe")

#	Restrict DE genes to universe
	for (s in 1:nsets) de[[s]] <- de[[s]][de[[s]] %in% universe]

#	Restrict pathways to universe
	i <- GeneID.PathID[,1] %in% universe
	if(sum(i)==0L) stop("Pathways do not overlap with universe")
	GeneID.PathID <- GeneID.PathID[i,]

#	Get prior.prob trend in DE with respect to the covariate, combining all de lists
	if(!is.null(covariate)) {
		if(!is.null(prior.prob)) message("prior.prob being recomputed from covariate")
		covariate <- as.numeric(covariate)
		isDE <- (universe %in% unlist(de))
		o <- order(covariate)
		prior.prob <- covariate
		span <- approx(x=c(20,200),y=c(1,0.5),xout=sum(isDE),rule=2,ties=list("ordered",mean))$y
		prior.prob[o] <- tricubeMovingAverage(isDE[o],span=span)
		if(plot) barcodeplot(covariate, index=isDE, worm=TRUE, span.worm=span, main="DE status vs covariate")
	}

#	Overlap pathways with DE genes
#	Create incidence matrix (X) of gene.pathway by DE sets
	if(is.null(prior.prob)) {
		X <- matrix(1,nrow(GeneID.PathID),nsets+1)
		colnames(X) <- c("N",names(de))
	} else {
		names(prior.prob) <- universe
		X <- matrix(1,nrow(GeneID.PathID),nsets+2)
		X[,nsets+2] <- prior.prob[GeneID.PathID[,1]]
		colnames(X) <- c("N",names(de),"PP")
	}
	for (s in 1:nsets) X[,s+1] <- (GeneID.PathID[,1] %in% de[[s]])

#	Count #genes and #DE genes and sum prior.prob for each pathway
	S <- rowsum(X, group=GeneID.PathID[,2], reorder=FALSE)

#	Overlap tests
	PValue <- matrix(0,nrow=nrow(S),ncol=nsets)
	colnames(PValue) <- paste("P", names(de), sep=".")
	nde <- lengths(de, use.names=FALSE)
	if(!is.null(prior.prob)) {

#		Probability ratio for each pathway vs rest of universe
		SumPP <- sum(prior.prob)
		M2 <- NGenes-S[,"N"]
		Odds <- S[,"PP"] / (SumPP-S[,"PP"]) * M2 / S[,"N"]

#		Wallenius' noncentral hypergeometric test
#		Note that dWNCHypergeo() is more accurate than pWNCHypergeo(), hence use of 2-terms
		if(!requireNamespace("BiasedUrn",quietly=TRUE)) stop("BiasedUrn package required but is not installed (or can't be loaded)")
		for(j in seq_len(nsets)) for(i in seq_len(nrow(S)))
			PValue[i,j] <- BiasedUrn::pWNCHypergeo(S[i,1L+j], S[i,"N"], M2[i], nde[j], Odds[i], lower.tail=FALSE) +
			               BiasedUrn::dWNCHypergeo(S[i,1L+j], S[i,"N"], M2[i], nde[j], Odds[i])

#		Remove sum of prob column, not needed for output
		S <- S[,-ncol(S)]

	} else {

#		Fisher's exact test
		for(j in seq_len(nsets))
			PValue[,j] <- phyper(S[,1L+j]-0.5, nde[j], NGenes-nde[j], S[,"N"], lower.tail=FALSE)

	}

#	Assemble output
	g <- rownames(S)
	m <- match(g, PathID.PathName[,1])
	Results <- data.frame(Pathway = PathID.PathName[m,2], S, PValue, stringsAsFactors=FALSE)
	rownames(Results) <- g

	Results
}

getGeneKEGGLinks <- function(species.KEGG="hsa", convert=FALSE)
#	Read pathway-gene mapping from KEGG website
#	Gordon Smyth
#	Created 7 Jan 2016.  Last modified 11 Feb 2016.
{
	URL <- paste("http://rest.kegg.jp/link/pathway",species.KEGG,sep="/")
	Path.Gene <- read.table(URL,sep="\t",quote="\"",fill=TRUE,comment.char="",stringsAsFactors=FALSE)
	colnames(Path.Gene) <- c("GeneID","PathwayID")

#	Human, mouse, rat and chimp IDs are already Entrez, so don't need to convert
	if(convert && species.KEGG %in% c("hsa","mmu","rno","ptr")) convert <- FALSE

	if(convert) {
#		Convert KEGG IDs to Entrez Gene IDs
		URL <- paste("http://rest.kegg.jp/conv",species.KEGG,"ncbi-geneid",sep="/")
		EntrezID.KEGGGeneID <- read.table(URL,sep="\t",quote="\"",fill=TRUE,comment.char="",stringsAsFactors=FALSE)
		m <- match(Path.Gene[,1],EntrezID.KEGGGeneID[,2])
		Path.Gene[,1] <- EntrezID.KEGGGeneID[m,1]
		Path.Gene[,1] <- sub("^ncbi-geneid:", "", Path.Gene[,1])
	} else {
		Path.Gene[,1] <- sub(paste0("^",species.KEGG,":"), "", Path.Gene[,1])
	}
	Path.Gene
}

getKEGGPathwayNames <- function(species.KEGG=NULL, remove.qualifier=FALSE)
#	Read pathways from KEGG website
#	Gordon Smyth
#	Created 7 Jan 2016.  Last modified 8 Jan 2016.
{
	URL <- "http://rest.kegg.jp/list/pathway"
	if(!is.null(species.KEGG)) URL <- paste(URL,species.KEGG,sep="/")
	PathID.PathName <- read.table(URL,sep="\t",quote="\"",fill=TRUE,comment.char="",stringsAsFactors=FALSE)
	colnames(PathID.PathName) <- c("PathwayID","Description")
	if(remove.qualifier) PathID.PathName[,2] <- removeExt(PathID.PathName[,2], sep=" - ")
	PathID.PathName
}

topKEGG <- function(results, sort = NULL, number = 20L, truncate.path=NULL)
#	Extract top KEGG pathways from kegga output 
#	Gordon Smyth and Yifang Hu
#	Created 15 May 2015. Modified 25 Jan 2016.
{
#	Check results
	if(!is.data.frame(results)) stop("results should be a data.frame.")
	dimres <- dim(results)

#	Check number
	if(!is.numeric(number)) stop("number should be a positive integer")
	if(number > dimres[1L]) number <- dimres[1]
	if(number < 1L) return(results[integer(0),])

#	Number of gene lists for which results are reported
#	Lists are either called "Up" and "Down" or have user-supplied names
	nsets <- (dimres[2L]-2L) %/% 2L
	if(nsets < 1L) stop("results has wrong number of columns")
	setnames <- colnames(results)[3L:(2L+nsets)]

#	Check sort. Defaults to all gene lists.
	if(is.null(sort)) {
		isort <- 1L:nsets
	} else {
		sort <- as.character(sort)
		isort <- which(tolower(setnames) %in% tolower(sort))
		if(!length(isort)) stop("sort column not found in results")
	}

#	Sort by minimum p-value for specified gene lists
	P.col <- 2L+nsets+isort
	if(nsets==1L)
		P <- results[,P.col]
	else
		P <- do.call(pmin,results[,P.col,drop=FALSE])
	o <- order(P,results$N,results$Pathway)
	tab <- results[o[1L:number],,drop=FALSE]

#	Truncate Pathway column for readability
	if(!is.null(truncate.path)) {
		truncate.path <- as.integer(truncate.path[1])
		truncate.path <- max(truncate.path,5L)
		truncate.path <- min(truncate.path,1000L)
		tm2 <- truncate.path-3L
		i <- (nchar(tab$Pathway) > tm2)
		tab$Pathway[i] <- paste0(substring(tab$Pathway[i],1L,tm2),"...")
	}

	tab
}
