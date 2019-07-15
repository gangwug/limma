#  PRESENTATION PLOTS FROM FITTED MODEL OBJECTS

volcanoplot <- function(fit,coef=1L,style="p-value",highlight=0L,names=fit$genes$ID,hl.col="blue",xlab="Log2 Fold Change",ylab=NULL,pch=16,cex=0.35, ...)
#	Volcano plot of log-fold-change vs significance (p-value or B-statistic)
#	Gordon Smyth
#	Created 27 Oct 2006. Last modified 14 Oct 2018.
{
	if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM")
	x <- as.matrix(fit$coef)[,coef]
	style <- match.arg(tolower(style), c("p-value","b-statistic"))
	if(style=="p-value") {
		if(is.null(fit$p.value)) stop("No p-values found in linear model fit object")
		y <- as.matrix(fit$p.value)[,coef]
		y <- -log10(y)
		if(is.null(ylab)) ylab="-log10(P-value)"
	} else {
		if(is.null(fit$lods)) stop("No B-statistics found in linear model fit object")
		y <- as.matrix(fit$lods)[,coef]
		if(is.null(ylab)) ylab="Log Odds of Differential Expression"
	}
	plot(x,y,xlab=xlab,ylab=ylab,pch=pch,cex=cex,...)
	if(highlight>0) {
		if(is.null(names)) names <- 1:length(x)
		names <- as.character(names)
		o <- order(y,decreasing=TRUE)
		i <- o[1:highlight]
		text(x[i],y[i],labels=substring(names[i],1,8),cex=0.8,col=hl.col)
	}
	invisible()
}

heatDiagram <- function(results,coef,primary=1,names=NULL,treatments=colnames(coef),limit=NULL,orientation="landscape",low="green",high="red",cex=1,mar=NULL,ncolors=123,...) {
#	Heat diagram to display fold changes of genes under different conditions
#	Gordon Smyth
#	Created 27 Oct 2002. Last revised 11 Oct 2004.

#	Check input
	results <- as.matrix(results)
	results[is.na(results)] <- 0
	coef <- as.matrix(coef)
	if(!identical(dim(results),dim(coef))) stop("results and coef must be the same size")
	nt <- ncol(results)
	if(is.null(names)) names <- as.character(1:nrow(coef))
	names <- substring(names,1,15)
	if(is.null(treatments)) treatments <- as.character(1:nt)
	orientation <- match.arg(orientation,c("landscape","portrait"))

#	Sort coefficients
	DE <- (abs(results[,primary]) > 0.5)
	ng <- sum(DE)
	if(ng == 0) {
		warning("Nothing significant to plot")
		return(invisible())
	}
	results <- results[DE,,drop=FALSE]
	coef <- coef[DE,,drop=FALSE]
	coef[abs(results) < 0.5] <- NA
	names <- names[DE]
	ord <- order(coef[,primary],decreasing=TRUE)

#	Truncate coefficients if limit is preset
	if(!is.null(limit))
		if(limit > 0) {
			coef[coef < -limit] <- -limit
			coef[coef > limit] <- limit
		} else
			warning("limit ignored because not positive")

#	Check colours
	if(is.character(low)) low <- col2rgb(low)/255
	if(is.character(high)) high <- col2rgb(high)/255
	r <- range(coef,na.rm=TRUE)
	r <- r/max(abs(r))
	low2 <- low + (high-low)*(1+r[1])/2
	high2 <- high + (low-high)*(1-r[2])/2
	col <- rgb( seq(low2[1],high2[1],len=ncolors), seq(low2[2],high2[2],len=ncolors), seq(low2[3],high2[3],len=ncolors) )

#	Output dataframe
	coef <- coef[ord,,drop=FALSE]
	names <- names[ord]
	out <- coef
	rownames(out) <- names

#	Insert white space between up and down
	nup <- sum(coef[,primary]>=0)
	if(nup>0 && nup<ng) {
		coef <- rbind(coef[1:nup,,drop=FALSE],matrix(NA,1,nt),coef[(nup+1):ng,,drop=FALSE])
		names <- c(names[1:nup],"",names[(nup+1):ng])
		ng <- ng+1
	}
	if(orientation=="portrait") {
		coef <- t(coef)
		coef <- coef[,ng:1,drop=FALSE]
	}

#	Heat plot
	on.exit(par(old.par))
	if(orientation=="portrait") {
		if(is.null(mar)) mar <- cex*c(1,1,4,3)
		old.par <- par(mar=mar)
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,8/nt)
		mtext(paste(" ",treatments,sep=""),side=3,las=2,at=(cext-1)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,40/ng)
		mtext(paste(" ",names,sep=""),side=4,las=2,at=(1-cex)*0.005+((ng-1):0)/(ng-1),cex=cex)
	} else {
		if(is.null(mar)) mar <- cex*c(5,6,1,1)
		old.par <- par(mar=mar)
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,12/nt)
		mtext(paste(treatments," ",sep=""),side=2,las=1,at=(1-cext)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,60/ng)
		mtext(paste(names," ",sep=""),side=1,las=2,at=(cex-1)*0.005+(0:(ng-1))/(ng-1),cex=cex)
	}
	invisible(out)
}

heatdiagram <- function(stat,coef,primary=1,names=NULL,treatments=colnames(stat),critical.primary=4,critical.other=3,limit=NULL,orientation="landscape",low="green",high="red",cex=1,mar=NULL,ncolors=123,...) {
#	Heat diagram to display fold changes of genes under different conditions
#	Similar to heatDiagram with classifyTestsT(stat,t1=critical.primary,t2=critical.other)
#	except that heatdiagram() requires primary column to be significant at the first step-down level
#	Gordon Smyth
#	27 Oct 2002. Last revised 25 Feb 2003.  mar added 11 Oct 2004

#	Check input
	stat <- as.matrix(stat)
	coef <- as.matrix(coef)
	if(any(dim(stat) != dim(coef))) stop("stat and coef must be the same size")
	nt <- ncol(stat)
	if(is.null(names)) names <- as.character(1:nrow(stat))
	names <- substring(names,1,15)
	if(is.null(treatments)) treatments <- as.character(1:nt)

#  Sort coefficients
	DE <- (stat[,primary] > critical.primary)
	if(any(is.na(DE))) DE[is.na(DE)] <- FALSE
	ng <- sum(DE)
	if(sum(DE) == 0) {
		warning("Nothing significant to plot")
		return(invisible())
	}
	stat <- stat[DE,,drop=FALSE]
	coef <- coef[DE,,drop=FALSE]
	if(!is.null(names)) names <- names[DE]
	if(critical.other > critical.primary) warning("critical.other greater than critical.primary")
	otherDE <- (stat > critical.other)
	otherDE[,primary] <- TRUE
	coef[!otherDE] <- NA
	ord <- order(coef[,primary],decreasing=TRUE)

#  Check colours
	if(is.character(low)) low <- col2rgb(low)/255
	if(is.character(high)) high <- col2rgb(high)/255
	col <- rgb( seq(low[1],high[1],len=ncolors), seq(low[2],high[2],len=ncolors), seq(low[3],high[3],len=ncolors) )

#  Truncate coefficients if limit is preset
	if(!is.null(limit))
		if(limit > 0) {
			coef[coef < -limit] <- -limit
			coef[coef > limit] <- limit
		} else
			warning("limit ignored because not positive")

#  Heat plot
	coef <- coef[ord,,drop=FALSE]
	names <- names[ord]
	out <- data.frame(Name=names,coef)
	if(orientation=="portrait") {
		coef <- t(coef)
		coef <- coef[,ng:1,drop=FALSE]
	}
	on.exit(par(old.par))
	if(orientation=="portrait") {
		if(is.null(mar)) mar <- cex*c(1,1,4,3)
		old.par <- par(mar=mar)
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,8/nt)
		mtext(paste(" ",treatments,sep=""),side=3,las=2,at=(cext-1)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,40/ng)
		mtext(paste(" ",names,sep=""),side=4,las=2,at=(1-cex)*0.005+((ng-1):0)/(ng-1),cex=cex)
	} else {
		if(is.null(mar)) mar <- cex*c(5,6,1,1)
		old.par <- par(mar=mar)
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,12/nt)
		mtext(paste(treatments," ",sep=""),side=2,las=1,at=(1-cext)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,60/ng)
		mtext(paste(names," ",sep=""),side=1,las=2,at=(cex-1)*0.005+(0:(ng-1))/(ng-1),cex=cex)
	}
	invisible(out)
}

plotSA <- function(fit, xlab="Average log-expression", ylab="sqrt(sigma)", zero.weights=FALSE, pch=16, cex=0.3, col=c("black","red"),...)
#	Plot log-residual variance vs intensity
#	Gordon Smyth
#	Created 14 Jan 2009. Last modified 12 April 2017.
{
	if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM object")
	x <- fit$Amean
	y <- sqrt(fit$sigma)
	if(!(is.null(fit$weights) || zero.weights)) {
		allzero <- rowSums(fit$weights>0,na.rm=TRUE) == 0
		y[allzero] <- NA
	}
	colv <- rep_len(col[1],nrow(fit))

#	Check for outlier variances
	if(length(fit$df.prior)>1L) {
		df2 <- max(fit$df.prior)
		s2 <- fit$sigma^2 / fit$s2.prior
		pdn <- pf(s2, df1=fit$df.residual, df2=df2)
		pup <- pf(s2, df1=fit$df.residual, df2=df2, lower.tail=FALSE)
		FDR <- p.adjust(2*pmin(pdn,pup),method="BH")
		colv[FDR <= 0.5] <- col[2]
	}

	plot(x,y,xlab=xlab,ylab=ylab,pch=pch,cex=cex,col=colv,...)
#	if(anyNA(x) || anyNA(y)) {
#		ok <- !(is.na(x) | is.na(y))
#		lines(lowess(x[ok],y[ok],f=0.4),col="red")
#	} else {
#		lines(lowess(x,y,f=0.4),col="red")
#	}
	if(!is.null(fit$s2.prior)) {
		if(length(fit$s2.prior)==1L) {
			abline(h=sqrt(sqrt(fit$s2.prior)),col="blue")
		} else {
			o <- order(x)
			lines(x[o],sqrt(sqrt(fit$s2.prior[o])),col="blue")
#			legend("topright",legend=c("lowess","prior"),col=c("red","blue"),lty=1)
		}
	}

	if(length(fit$df.prior)>1L) legend("topright",legend=c("Normal","Outlier"),col=col,pch=pch)

	invisible()
}
