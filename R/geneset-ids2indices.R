ids2indices <- function(gene.sets, identifiers, remove.empty=TRUE)
# Make a list of gene identifiers into a list of indices for gene sets
# Gordon Smyth and Yifang Hu
# 25 March 2009.  Last modified 6 May 2015.
{
	if(!is.list(gene.sets)) gene.sets <- list(Set1=gene.sets)
	identifiers <- as.character(identifiers)
	index <- lapply(gene.sets, function(x) which(identifiers %in% as.character(x)))
	if(remove.empty) {
		Empty <- which(lengths(index)==0L)
		if(length(Empty)) for (i in length(Empty):1) index[[Empty[i]]] <- NULL
	}
	index
}
