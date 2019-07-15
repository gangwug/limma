propexpr <- function(x,neg.x=NULL,status=x$genes$Status,labels=c("negative","regular"))
#	Estimate proportion of expressed probes on each array
#	Wei Shi and Gordon Smyth.
#	17 April 2009. Last modified 15 Oct 2009.
{
	if(is.null(neg.x)) {
		ineg <- grep(tolower(labels[1]),tolower(status))
		if(length(labels)>1) {
			ireg <- grep(tolower(labels[2]),tolower(status))
		} else {
			ireg <- -ineg
		}
		x <- as.matrix(x)
		neg.x <- x[ineg,,drop=FALSE]
		x <- x[ireg,,drop=FALSE]
	} else {
		x <- as.matrix(x)
		neg.x <- as.matrix(neg.x)
	}

	narrays <- ncol(x)
	p <- pb <- p1 <- rep(NA, narrays)

	for(i in 1:narrays) {
		b <- neg.x[, i]
		b <- b[!is.na(b)]
		nb <- length(b)
	
		r <- x[, i]
		r <- r[!is.na(r)]
		nr <- length(r)
	
		mu <- mean(b)
		alpha <- max(mean(r) - mu, 10)
	
		b1 <- median(b)
			
		p1[i] <- mean(pexp(b1-b, rate=1/alpha))	
		pb[i] <- (sum(b<b1) + sum(b==b1)/2)/nb
		p[i] <- (sum(r<b1) + sum(r==b1)/2)/nr
	}

	pi1 <- (pb - p)/(pb - p1)
	pi1[pi1 > 1] <- 1
	pi1[pi1 < 0] <- 0
	names(pi1) <- colnames(x)
	pi1
}
