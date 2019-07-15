goana <- function(de,...) UseMethod("goana")

goana.MArrayLM <- function(de, coef = ncol(de), geneid = rownames(de), FDR = 0.05, trend = FALSE, ...)
#	Gene ontology analysis of DE genes from linear model fit
#	Gordon Smyth and Yifang Hu
#	Created 20 June 2014.  Last modified 1 May 2015.
{
#	Avoid argument collision with default method
	dots <- names(list(...))
	if("universe" %in% dots) stop("goana.MArrayLM defines its own universe",call.=FALSE)
	if((!is.logical(trend) || trend) && "covariate" %in% dots) stop("goana.MArrayLM defines it own covariate",call.=FALSE)

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
		goana(de=DEGenes, universe = universe, covariate=covariate, ...)
	else
		goana(de=DEGenes, universe = universe, ...)
}

goana.default <- function(de, universe = NULL, species = "Hs", prior.prob = NULL, covariate=NULL, plot=FALSE, ...)
#	Gene ontology analysis of DE genes
#	Gordon Smyth and Yifang Hu
#	Created 20 June 2014.  Last modified 19 May 2019.
{
#	Get access to package of GO terms
	suppressPackageStartupMessages(OK <- requireNamespace("GO.db",quietly=TRUE))
	if(!OK) stop("GO.db package required but not installed (or can't be loaded)")

#	Get access to required annotation functions
	suppressPackageStartupMessages(OK <- requireNamespace("AnnotationDbi",quietly=TRUE))
	if(!OK) stop("AnnotationDbi package required but not installed (or can't be loaded)")

#	Load appropriate organism package
	orgPkg <- paste0("org.",species,".eg.db")
	suppressPackageStartupMessages(OK <- requireNamespace(orgPkg,quietly=TRUE))
	if(!OK) stop(orgPkg," package required but not not installed (or can't be loaded)")

#	Get GO to Entrez Gene mappings
	obj <- paste0("org.",species,".egGO2ALLEGS")
	egGO2ALLEGS <- tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE)
	if(is.logical(egGO2ALLEGS)) stop("Can't find gene ontology mappings in package ",orgPkg)

#	Convert gene-GOterm mappings to data.frame
#	Remove duplicate entries from both mappings and universe
	if(is.null(universe)) {
		GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id","go_id","Ontology")]
		i <- !duplicated(GeneID.PathID[,c("gene_id", "go_id")])
		GeneID.PathID <- GeneID.PathID[i, ]
		universe <- unique(GeneID.PathID[,1])
		prior.prob <- covariate <- NULL
	} else {
		universe <- as.character(universe)
		lu <- length(universe)
		if(!is.null(prior.prob) && length(prior.prob)!=lu) stop("universe and prior.prob must have same length")
		if(!is.null(covariate) && length(covariate)!=lu) stop("universe and covariate must have same length")
		if(anyDuplicated(universe)) {
			i <- !duplicated(universe)
			if(!is.null(covariate)) covariate <- covariate[i]
			if(!is.null(prior.prob)) prior.prob <- prior.prob[i]
			universe <- universe[i]
		}
#		Make universe and set of all annotated genes agree
		i <- (universe %in% AnnotationDbi::Lkeys(egGO2ALLEGS))
		universe <- universe[i]
		if(!is.null(covariate)) covariate <- covariate[i]
		if(!is.null(prior.prob)) prior.prob <- prior.prob[i]
		AnnotationDbi::Lkeys(egGO2ALLEGS) <- universe
#		Convert GO annotation to data.frame
		GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id","go_id","Ontology")]
		d <- duplicated(GeneID.PathID[,c("gene_id", "go_id")])
		GeneID.PathID <- GeneID.PathID[!d, ]
	}

#	From here, code is mostly the same as kegga.default

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

#	Check universe isn't empty
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
	GOID <- rownames(S)
	TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,keys=GOID,columns="TERM"))
	m <- match(GOID,GeneID.PathID[,2])
	Ont <- GeneID.PathID[m,3]
	Results <- data.frame(Term=TERM[,2], Ont=Ont, S, PValue, stringsAsFactors=FALSE)

	Results
}

topGO <- function(results, ontology = c("BP", "CC", "MF"), sort = NULL, number = 20L, truncate.term=NULL)
#	Extract top GO terms from goana output 
#	Gordon Smyth and Yifang Hu
#	Created 20 June 2014. Last modified 23 June 2016.
{
#	Check results
	if(!is.data.frame(results)) stop("results should be a data.frame.")

#	Check ontology
	ontology <- match.arg(unique(ontology), c("BP", "CC", "MF"), several.ok = TRUE)

#	Limit results to specified ontologies
	if(length(ontology) < 3L) {
		sel <- results$Ont %in% ontology
		results <- results[sel,]
	}
	dimres <- dim(results)

#	Check number
	if(!is.numeric(number)) stop("number should be a positive integer")
	if(number > dimres[1L]) number <- dimres[1]
	if(number < 1L) return(results[integer(0),])

#	Number of gene lists for which results are reported
#	Lists are either called "Up" and "Down" or have user-supplied names
	nsets <- (dimres[2L]-3L) %/% 2L
	if(nsets < 1L) stop("results has wrong number of columns")
	setnames <- colnames(results)[4L:(3L+nsets)]

#	Check sort. Defaults to all gene lists.
	if(is.null(sort)) {
		isort <- 1L:nsets
	} else {
		sort <- as.character(sort)
		isort <- which(tolower(setnames) %in% tolower(sort))
		if(!length(isort)) stop("sort column not found in results")
	}

#	Sort by minimum p-value for specified gene lists
	P.col <- 3L+nsets+isort
	if(length(P.col)==1L) {
		P <- results[,P.col]
	} else {
		P <- do.call("pmin",as.data.frame(results[,P.col,drop=FALSE]))
	}
	o <- order(P,results$N,results$Term)
	tab <- results[o[1L:number],,drop=FALSE]

#	Truncate Term column for readability
	if(!is.null(truncate.term)) {
		truncate.term <- as.integer(truncate.term[1])
		truncate.term <- max(truncate.term,5L)
		truncate.term <- min(truncate.term,1000L)
		tm2 <- truncate.term-3L
		i <- (nchar(tab$Term) > tm2)
		tab$Term[i] <- paste0(substring(tab$Term[i],1L,tm2),"...")
	}

	tab
}
