cumOverlap <- function(ol1, ol2)
#	Cumulative overlap analysis.
#	Test whether two ordered lists of IDs are significantly overlapped. 
#	Di Wu and Gordon Smyth
#	Createdin 2007. Last modified 20 Sep 2018
{
#	Check for duplicates
	if(anyDuplicated(ol1)) stop("Duplicate IDs found in ol1")
	if(anyDuplicated(ol2)) stop("Duplicate IDs found in ol2")

#	Reduce to IDs found in both lists
	m <- match(ol1,ol2)
	redo <- FALSE
	if(anyNA(m)) {
		ol1 <- ol1[!is.na(m)]
		redo <- TRUE
	}
	m2 <- match(ol2,ol1)
	if(anyNA(m2)) {
		ol2 <- ol2[!is.na(m2)]
		redo <- TRUE
	}

#	Match ol1 to ol2
	if(redo) m <- match(ol1,ol2)

#	Count overlaps
	ngenes <- length(ol1)
	if(ngenes == 0L) return(list(n.total=0L))
	i <- noverlap <- 1:ngenes
	for (j in i) noverlap[j] <- sum(m[1:j] <= j)

#	Hypergeometric p-valules
	p <- phyper(noverlap-0.5,m=i,n=ngenes-i,k=i,lower.tail=FALSE)

#	Directed Bonferroni adjustment, starting from top of list
	p.b <- p*i
	nmin <- which.min(p.b)
	p.b <- pmin(p.b,1)

#	Which ids contribute to overlap?
	idoverlap <- ol1[which(m[1:nmin] <= nmin)]

	list(n.total=ngenes,n.min=nmin,p.min=p.b[nmin],n.overlap=noverlap,id.overlap=idoverlap,p.value=p,adj.p.value=p.b)
}
