plotExonJunc <- function(fit, coef=ncol(fit), geneid, genecolname=NULL, FDR=0.05, annotation=NULL)
#	Assuming the data for diffSplice analysis contains both exons and junctions,
#	'fit' is a MArrayLM object produced by diffSplice().
#	To distinguish between exons and junctions, 'fit$genes$Length' are set to 1 for all the junctions.
#	Since the diffSplice analysis is usually performed after filtering,
#	the full annotation (e.g. the inbuilt annotation used by featureCounts)
#	is required for producing the plot.
#	Yunshun Chen and Gordon Smyth.
#	Created 9 March 2018. Last modified 29 Oct 2018.
{
	if(is.null(genecolname)) 
		genecolname <- fit$genecolname
	else
		genecolname <- as.character(genecolname)

	geneid <- as.character(geneid)
	i <- fit$genes[, genecolname]==geneid
	if(!any(i)) stop(paste0(geneid, " not found."))

#	Subsetting
	fdr <- p.adjust(fit$p.value[,coef], method="BH")
	fdr <- fdr[i]
	genes <- fit$genes[i,]
	strand <- genes$Strand[1]
	p.value <- fit$p.value[i,coef]
	coefficients <- fit$coefficients[i,coef]

#	Sorting	exons and junctions by their start positions
	o <- order(genes$Start, genes$End)
	genes <- genes[o,]
	p.value <- p.value[o]
	coefficients <- coefficients[o]
	
#	Which ones are exons? (Junctions are assigned length of 1 prior to the diffSplice analysis)
	IsExon <- genes$Length > 1L
	genes.e <- genes[IsExon, ]

#	Check the format of the annotation file.
	if(!is.null(annotation)){
		if(is.null(annotation$GeneID)) stop("Annotation file must contain Entrez gene ids.")
		if(is.null(annotation$Start) | is.null(annotation$End)) stop("Annotation file must contain start-end positions of exons.")
		if(is.null(annotation$Length)) annotation$Length <- abs(annotation$End - annotation$Start) + 1L

#		Retrieve annotation information for the exons that have been filtered prior to the diffSplice analysis
		sel <- annotation$GeneID == genes$GeneID[1]
		genes.e2 <- annotation[sel, ]
		genes.e2 <- genes.e2[order(genes.e2$Start), ]
		m <- match(genes.e$Start, genes.e2$Start)
		genes.e <- genes.e2
	}

#	Get the start-end positions and the length for all the introns
	Start.i <- genes.e$End[-nrow(genes.e)] + 1L
	End.i <- genes.e$Start[-1] - 1L
	genes.i <- genes.e[-1,]
	genes.i$Start <- Start.i
	genes.i$End <- End.i
	genes.i$Length <- genes.i$End - genes.i$Start + 1L

#	Get the start-end positions for all the junctions
	if(any(!IsExon)){
		genes.j <- genes[!IsExon, ]
		
#		Extend the plotting range for the gene in case there are junctions outside of the gene range.
		if(min(genes.j$Start) < min(genes.e$Start)){
			intron <- genes.i[1,,drop=FALSE]
			intron$Start <- min(genes.j$Start)
			intron$End <- min(genes.e$Start) - 1L
			intron$Length <- intron$End - intron$Start + 1L
			genes.i <- rbind(intron, genes.i)
		}
		if(max(genes.j$End) > max(genes.e$End)){
			intron <- genes.i[1,,drop=FALSE]
			intron$Start <- max(genes.e$End) + 1L
			intron$End <- max(genes.j$End)
			intron$Length <- intron$End - intron$Start + 1L
			genes.i <- rbind(genes.i, intron)
		}
	}

#	Combine introns and exons for plotting
	genes.ie <- rbind(cbind(genes.e, Flag="Exon"), cbind(genes.i, Flag="Intron"))
	genes.ie <- genes.ie[order(genes.ie$Start), ]

#	Scale the length of intron/exon segments for better visualization
	pseudo.length <- (genes.ie$Length)^.5
	pseudo.pos <- cumsum((genes.ie$Length)^.5)
	pseudo.start <- c(0, pseudo.pos[-nrow(genes.ie)])
	pseudo.end <- pseudo.pos
	genes.ie <- cbind(genes.ie, pseudo.start=pseudo.start, pseudo.end=pseudo.end, pseudo.length=pseudo.length)

#	Update start-end postions for junctions on the pseudo scale
	if(any(!IsExon)){
		genes.j <- cbind(genes.j, pseudo.start=0, pseudo.end=0)
		for(j in 1:nrow(genes.j)){
			k <- which(genes.j$Start[j] <= genes.ie$End)[1]
			genes.j$pseudo.start[j] <- genes.ie$pseudo.end[k] - (genes.ie$End[k] - genes.j$Start[j]) / genes.ie$Length[k] * genes.ie$pseudo.length[k]
			k <- which(genes.j$End[j] <= genes.ie$End)[1]
			genes.j$pseudo.end[j] <- genes.ie$pseudo.end[k] - (genes.ie$End[k] - genes.j$End[j]) / genes.ie$Length[k] * genes.ie$pseudo.length[k]
		}
	}

#	Setup the plot
	GeneStart <- min(genes.ie$pseudo.start)
	GeneEnd <- max(genes.ie$pseudo.end)
	gene.length <- GeneEnd - GeneStart
	plot.new()
	plot.window(xlim=c(GeneStart, GeneEnd), ylim=c(-0.7, 0.7))
	title(main=paste0(geneid, " (", genes$Strand[1], ")"))

#	Plot gene range
	rect(xleft=GeneStart, xright=GeneEnd, ybottom=-0.02, ytop=0.02, col="gray", border="gray")
	if(strand=="+"){
		tx.left <- "5'"
		tx.right <- "3'"
	} else {
		tx.left <- "3'"
		tx.right <- "5'"
	}
	text(x=-0.02*gene.length, y=0.1, labels=tx.left)
	text(x=1.02*gene.length, y=0.1, labels=tx.right)
	
#	Direction and significance of the diffSplice results
	up <- coefficients > 0
	down <- coefficients < 0
	IsSig <- fdr < FDR
	#IsSig <- p.adjust(p.value, method="holm") < FDR
	IsSig.j <- IsSig[!IsExon]
	down.j <- down[!IsExon]
	
#	Colouring
	col <- rep("black", sum(i))
	col[up & IsSig] <- "red"
	col[down & IsSig] <- "dodgerblue"
	col.e <- col[IsExon]
	col.j <- col[!IsExon]
	
#	Filtered exons are coloured in grey
	if(!is.null(annotation)){
		col.e2 <- rep("grey", nrow(genes.e))
		col.e2[m] <- col.e
		col.e <- col.e2
	}

#	Plot exons
	ex <- genes.ie$Flag=="Exon"
	rect(xleft=genes.ie$pseudo.start[ex], xright=genes.ie$pseudo.end[ex], ybottom=-0.1,ytop=0.1, col=col.e, border=col.e)

#	Plot junctions
	if(any(!IsExon)){
		MidPoint <- (genes.j$pseudo.start + genes.j$pseudo.end)/2
		y0 <- rep(0.11, sum(!IsExon))
		y1 <- rep(0.4, sum(!IsExon))
		y1[IsSig.j] <- 0.6
		y0[down.j] <- -y0[down.j]
		y1[down.j] <- -y1[down.j]
		segments(x0=genes.j$pseudo.start, x1=MidPoint, y0=y0, y1=y1, col=col.j, lwd=2)
		segments(x0=MidPoint, x1=genes.j$pseudo.end, y0=y1, y1=y0, col=col.j, lwd=2)
	}

#	Label axis
	if(genes$Strand[1]=="+")
		labels <- paste0("Exon.", 1:length(col.e))
	else
		labels <- paste0("Exon.", length(col.e):1)
	axis(side=1, at=(genes.ie$pseudo.start[ex]+genes.ie$pseudo.end[ex])/2, las=2, labels=labels)

	invisible()
}




