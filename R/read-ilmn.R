#  READ-ILMN.R

read.ilmn <- function(files=NULL, ctrlfiles=NULL, path=NULL, ctrlpath=NULL, probeid="Probe", annotation=c("TargetID", "SYMBOL"), expr="AVG_Signal", other.columns="Detection",sep="\t", quote="\"", verbose=TRUE, ...)
#	Read one or more files of Illumina BeadStudio output
#	Wei Shi and Gordon Smyth.
#	Created 15 July 2009. Last modified 21 July 2014.
{
	if(!is.null(files)){
		f <- unique(files)
		if(!is.null(path)) f <- file.path(path, f)
		n <- length(f)
		for(i in 1:n){
			if(verbose) cat("Reading file", f[i], "... ...\n")
			elist1 <- .read.oneilmnfile(f[i], probeid, annotation, expr, other.columns, sep, quote, verbose, ...)
			if(i==1)
				elist <- elist1
			else
				elist <- cbind(elist, elist1)
		}
		if(!is.null(ctrlfiles))
			elist$genes$Status <- "regular"
	}
		
	if(!is.null(ctrlfiles)){
		cf <- unique(ctrlfiles)
		if(is.null(ctrlpath)) ctrlpath <- path
		if(!is.null(ctrlpath)) cf <- file.path(ctrlpath, cf)
		n <- length(cf)
		for(i in 1:n){
			if(verbose) cat("Reading file", cf[i], "... ...\n")
			elist.ctrl1 <- .read.oneilmnfile(cf[i], probeid, annotation, expr, other.columns, sep, quote, verbose, ...)
			if(i==1)
				elist.ctrl <- elist.ctrl1
			else
				elist.ctrl <- cbind(elist.ctrl, elist.ctrl1)
		}
		
		if(is.null(elist.ctrl$genes)) elist.ctrl$genes <- data.frame(Status = rep("negative", nrow(elist.ctrl)))
		else {
			ctrl.status.negative <- apply(elist.ctrl$genes, 2, function(x) sum(tolower(x) %in% "negative"))
			STATUS.col <- which.max(ctrl.status.negative)
			if(ctrl.status.negative[STATUS.col] > 0) elist.ctrl$genes$Status <- elist.ctrl$genes[[STATUS.col]]
			else elist.ctrl$genes$Status <- "negative"
		}
	}
	
	if(!is.null(files))
		if(!is.null(ctrlfiles)){
			REG.col <- setdiff(colnames(elist$genes), colnames(elist.ctrl$genes))
			if(length(REG.col)) for(i in REG.col) elist.ctrl$genes[[i]] <- NA
			CTRL.col <- setdiff(colnames(elist.ctrl$genes), colnames(elist$genes))
			if(length(CTRL.col)) for(i in CTRL.col) elist$genes[[i]] <- NA
			return(rbind(elist, elist.ctrl))
		}
		else
			return(elist)
	else
		if(!is.null(ctrlfiles))
			return(elist.ctrl)
	}
	
read.ilmn.targets <- function(targets, ...)
#  Read Illumina BeadStudio output using targets frame
#  Wei Shi
#  15 July 2009
{
	if(is.null(targets$files) && is.null(targets$ctrlfiles))
		stop("Can not find column \"files\" or \"ctrlfiles\" in targets\n") 
	x <- read.ilmn(targets$files, targets$ctrlfiles, ...) 
	x$targets <- targets
	x
}

.read.oneilmnfile <- function(fname, probeid, annotation, expr, other.columns, sep, quote, verbose, ...)
#	Read a single file of Illumina BeadStudio output
#	Wei Shi and Gordon Smyth
#	Created 15 July 2009. Last modified 21 July 2014.
{
	h <- readGenericHeader(fname,columns=expr,sep=sep)
	skip <- h$NHeaderRecords
	header <- h$ColumnNames

	elist <- new("EListRaw")
	elist$source <- "illumina"
	reqcol <- header[grep(tolower(paste(c(probeid, annotation, expr, other.columns), collapse="|")), tolower(header))]
	reqcol <- trimWhiteSpace(reqcol)

	x <- read.columns(file=fname, required.col=reqcol, skip=skip, sep=sep, quote=quote, stringsAsFactors=FALSE,	...)
	nprobes <- nrow(x)

#	Match column names to find column numbers
	cn <- tolower(colnames(x))
	idcol <- grep(tolower(probeid), cn)
	anncol <- grep(tolower(paste(annotation,collapse="|")), cn)
	exprcol <- grep(tolower(expr), cn)

#	Probe IDs
	pids <- x[,idcol]

#	Sample names
	snames <- colnames(x)[exprcol]
	snames <- unlist(strsplit(snames, paste("[.]*", expr, "-*", sep="")))
	snames <- snames[snames != ""]
	nsamples <- length(snames)

#	Expression matrix	
	elist$E <- data.matrix(x[,exprcol])
	colnames(elist$E) <- snames
	rownames(elist$E) <- pids

#	Add probe annotation	
	if(length(anncol)) {
		elist$genes <- x[,anncol,drop=FALSE]
		if(length(pids) & !any(duplicated(pids))) row.names(elist$genes) <- pids
	}

#	elist$targets <- data.frame(SampleNames=snames, stringsAsFactors=FALSE)

#	Add other columns if required	
	if(!is.null(other.columns)){
		elist$other <- vector("list", length(other.columns))
		for(i in 1:length(other.columns)){
			if(length(co <- grep(tolower(other.columns[i]), tolower(colnames(x)))))
				elist$other[[i]] <- as.matrix(x[, co])
			else
				elist$other[[i]] <- matrix(NA, nprobes, nsamples)
			colnames(elist$other[[i]]) <- snames
			rownames(elist$other[[i]]) <- pids
		}
		names(elist$other) <- other.columns
	}

	elist
}
