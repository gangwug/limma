#  VENN DIAGRAM COUNTS AND PLOTS

vennCounts <- function(x,include="both") {
#	Venn diagram counts
#	Gordon Smyth  
#	4 July 2003.  Last modified 10 March 2013.

	x <- as.matrix(x)
	include <- match.arg(include,c("both","up","down"))
	x <- sign(switch(include,
		both = abs(x),
		up = x > 0,
		down = x < 0
	))
	nprobes <- nrow(x)
	ncontrasts <- ncol(x)
	names <- colnames(x)
	if(is.null(names)) names <- paste("Group",1:ncontrasts)
	noutcomes <- 2^ncontrasts
	outcomes <- matrix(0,noutcomes,ncontrasts)
	colnames(outcomes) <- names
	for (j in 1:ncontrasts)
		outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
	xlist <- list()
	for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
	counts <- as.vector(table(xlist))
	rownames(outcomes) <- 1:nrow(outcomes)
	structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

vennDiagram <- function(object,include="both",names=NULL,mar=rep(1,4),cex=c(1.5,1,0.7),lwd=1,circle.col=NULL,counts.col=NULL,show.include=NULL,...)
#	Plot Venn diagram
#	Gordon Smyth, James Wettenhall, Yifang Hu.
#	Capabilities for multiple counts and colors uses code by Francois Pepin.
#	4 July 2003.  Last modified 31 October 2014.
{
#	Check include
	include <- as.character(include)
	LenInc <- min(length(include),2)

#	Get counts
	if(is(object, "VennCounts")) {
		include <- include[1]
		LenInc <- 1
	} else {
		if(LenInc>1) z2 <- vennCounts(object, include=include[2])[,"Counts"]
		object <- vennCounts(object, include=include[1])
	}
	z <- object[,"Counts"]
	nsets <- ncol(object)-1
	if(nsets > 5) stop("Can't plot Venn diagram for more than 5 sets")

#	Attach zone names to counts
	VennZone <- object[,1:nsets,drop=FALSE]
	VennZone <- apply(VennZone,1,function(x) paste(x,sep="",collapse=""))
	names(z) <- VennZone
	if(length(include)==2) names(z2) <- VennZone

#	Set names
	if(is.null(names)) names <- colnames(object)[1:nsets]

#	Set colors
	FILL.COL <- TRUE
	if(is.null(circle.col)) { 
		circle.col <- par('col')
		FILL.COL <- FALSE
	}
	if(length(circle.col)<nsets) circle.col <- rep(circle.col,length.out=nsets)
	if(is.null(counts.col)) counts.col <- par('col')
	if(length(counts.col)<LenInc) counts.col <- rep(counts.col,length.out=LenInc)

#	Set show.include
	if(is.null(show.include)) show.include <- as.logical(LenInc-1)

#	Set plot margins
	old.par <- par()$mar
	on.exit(par(mar=old.par))
	par(mar=mar)

	##############################################
	# 1--3 sets
	##############################################

	if(nsets <= 3) {

#	Open plot
	plot(x=0,y=0,type="n",xlim=c(-4,4),ylim=c(-4,4),xlab="",ylab="",axes=FALSE,...)

#	Plot circles
	theta <- 2*pi*(0:360)/360
	xcentres <- switch(nsets, 0, c(-1,1), c(-1,1,0) )
	ycentres <- switch(nsets, 0, c(0,0), c(1,1,-2)/sqrt(3) )
	r <- 1.5
	xtext <- switch(nsets, -1.2, c(-1.2,1.2), c(-1.2,1.2,0) )
	ytext <- switch(nsets, 1.8, c(1.8,1.8), c(2.4,2.4,-3) )
#	Note that lines() better than symbols() to make circles follow aspect ratio of plot
	for(circle in 1:nsets) {

		if(!FILL.COL) lines(xcentres[circle]+r*cos(theta),ycentres[circle]+r*sin(theta),lwd=lwd,col=circle.col[circle])
		if(FILL.COL) {
			RGB <- col2rgb(circle.col[circle])/255
			ALPHA <- 0.06
			RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], alpha = ALPHA)
			polygon(xcentres[circle]+r*cos(theta),ycentres[circle]+r*sin(theta),border=circle.col[circle],lwd=lwd,col=RGB.ALP)
		}

		text(xtext[circle],ytext[circle],names[circle],cex=cex)

	}

#	Plot rectangles
	switch(nsets,
		rect(-3,-2.5,3,2.5),
		rect(-3,-2.5,3,2.5),
		rect(-3,-3.5,3,3.3)
	)

#	Plot counts
	showCounts <- switch(nsets,
		function(counts, cex, adj,col,leg) {
			text(2.3,-2.1,counts[1],cex=cex,col=col,adj=adj)
			text(0,0,counts[2],cex=cex,col=col,adj=adj)
			if(show.include) text(-2.3,-2.1,leg,cex=cex,col=col,adj=adj)
		},
		function(counts, cex, adj,col,leg) {
			text(2.3,-2.1,counts[1],cex=cex,col=col,adj=adj)
			text(1.5,0.1,counts[2],cex=cex,col=col,adj=adj)
			text(-1.5,0.1,counts[3],cex=cex,col=col,adj=adj)
			text(0,0.1,counts[4],cex=cex,col=col,adj=adj)
			if(show.include) text(-2.3,-2.1,leg,cex=cex,col=col,adj=adj)
		},
		function(counts, cex, adj,col,leg) {
			text(2.5,-3,counts[1],cex=cex,col=col,adj=adj)
			text(0,-1.7,counts[2],cex=cex,col=col,adj=adj)
			text(1.5,1,counts[3],cex=cex,col=col,adj=adj)
			text(.75,-.35,counts[4],cex=cex,col=col,adj=adj)
			text(-1.5,1,counts[5],cex=cex,col=col,adj=adj)
			text(-.75,-.35,counts[6],cex=cex,col=col,adj=adj)
			text(0,.9,counts[7],cex=cex,col=col,adj=adj)
			text(0,0,counts[8],cex=cex,col=col,adj=adj)
			if(show.include) text(-2.5,-3,leg,cex=cex,col=col,adj=adj)
		}
	)
	if(LenInc==1) adj <- c(0.5,0.5) else adj <- c(0.5,0)
	showCounts(counts=z,cex=cex[1],adj=adj,col=counts.col[1],leg=include[1])
	if(LenInc==2) showCounts(counts=z2,cex=cex[1],adj=c(0.5,1),col=counts.col[2],leg=include[2])

	return(invisible())
	}

	##############################################
	# > 3 sets
	##############################################
	
#	Open plot
	plot(c(-20, 420), c(-20, 420), type="n", axes=FALSE, ylab="", xlab="", ...)

#	Function to turn and move ellipses
	relocate_elp <- function(e, alpha, x, y)
	{
		phi <- (alpha/180)*pi
		xr <- e[,1]*cos(phi)+e[,2]*sin(phi)
		yr <- -e[,1]*sin(phi)+e[,2]*cos(phi)
		xr <- x+xr
		yr <- y+yr
		cbind(xr, yr)
	}

	##############################################
	# 4 sets
	##############################################

   if (4 == nsets) {
		rect(-20, -20, 420, 400)
		elps <- cbind(162*cos(seq(0,2*pi,len=1000)), 108*sin(seq(0,2*pi,len=1000)))

		if(!FILL.COL){
			polygon(relocate_elp(elps,  45, 130, 170),border=circle.col[1],lwd=lwd)
			polygon(relocate_elp(elps,  45, 200, 200),border=circle.col[2],lwd=lwd)
			polygon(relocate_elp(elps, 135, 200, 200),border=circle.col[3],lwd=lwd)
			polygon(relocate_elp(elps, 135, 270, 170),border=circle.col[4],lwd=lwd)
		}
		if(FILL.COL){
			RGB <- col2rgb(circle.col)/255
			ALPHA <- 0.06
			RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], alpha = ALPHA)
			RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2], alpha = ALPHA)
			RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3], alpha = ALPHA)
			RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4], alpha = ALPHA)
			polygon(relocate_elp(elps,  45, 130, 170),border=circle.col[1],lwd=lwd,col=RGB.ALP1)
			polygon(relocate_elp(elps,  45, 200, 200),border=circle.col[2],lwd=lwd,col=RGB.ALP2)
			polygon(relocate_elp(elps, 135, 200, 200),border=circle.col[3],lwd=lwd,col=RGB.ALP3)
			polygon(relocate_elp(elps, 135, 270, 170),border=circle.col[4],lwd=lwd,col=RGB.ALP4)
		}

		text( 35, 315, names[1], cex=cex[1])
		text(138, 350, names[2], cex=cex[1])
		text(262, 347, names[3], cex=cex[1])
		text(365, 315, names[4], cex=cex[1])
		
		text( 35, 250, z["1000"], cex=cex[2], col=counts.col[1],)
		text(140, 315, z["0100"], cex=cex[2], col=counts.col[1])
		text(260, 315, z["0010"], cex=cex[2], col=counts.col[1])
		text(365, 250, z["0001"], cex=cex[2], col=counts.col[1])
		text( 90, 282, z["1100"], cex=cex[3], col=counts.col[1]) 
		text( 95, 110, z["1010"], cex=cex[2], col=counts.col[1]) 
		text(200,  52, z["1001"], cex=cex[3], col=counts.col[1]) 
		text(200, 292, z["0110"], cex=cex[2], col=counts.col[1]) 
		text(300, 110, z["0101"], cex=cex[2], col=counts.col[1]) 
		text(310, 282, z["0011"], cex=cex[3], col=counts.col[1]) 
		text(130, 230, z["1110"], cex=cex[2], col=counts.col[1])
		text(245,  81, z["1101"], cex=cex[3], col=counts.col[1]) 
		text(155,  81, z["1011"], cex=cex[3], col=counts.col[1])
		text(270, 230, z["0111"], cex=cex[2], col=counts.col[1]) 
		text(200, 152, z["1111"], cex=cex[2], col=counts.col[1]) 
		text(400,  15, z["0000"], cex=cex[1], col=counts.col[1]) 

		if(length(include)==2) {
			text( 35, 238, z2["1000"], cex=cex[2], col=counts.col[2])
			text(140, 304, z2["0100"], cex=cex[2], col=counts.col[2])
			text(260, 304, z2["0010"], cex=cex[2], col=counts.col[2])
			text(365, 238, z2["0001"], cex=cex[2], col=counts.col[2])
			text( 90, 274, z2["1100"], cex=cex[3], col=counts.col[2])
			text( 95, 100, z2["1010"], cex=cex[2], col=counts.col[2])
			text(200,  43, z2["1001"], cex=cex[3], col=counts.col[2])
			text(200, 280, z2["0110"], cex=cex[2], col=counts.col[2])
			text(300, 100, z2["0101"], cex=cex[2], col=counts.col[2])
			text(310, 274, z2["0011"], cex=cex[3], col=counts.col[2])
			text(130, 219, z2["1110"], cex=cex[2], col=counts.col[2])
			text(245,  71, z2["1101"], cex=cex[3], col=counts.col[2])
			text(155,  72, z2["1011"], cex=cex[3], col=counts.col[2])
			text(270, 219, z2["0111"], cex=cex[2], col=counts.col[2])
			text(200, 140, z2["1111"], cex=cex[2], col=counts.col[2])
			text(400,  -2, z2["0000"], cex=cex[1], col=counts.col[2])
				
			if(show.include) {
				text(10,15,include[1], cex=cex[1], col=counts.col[1])
				text(10,-2,include[2], cex=cex[1], col=counts.col[2])
			}
		}
		return(invisible())
 	}

	##############################################
	# 5 sets
	##############################################

	rect(-20, -30, 430, 430)

	elps <- cbind(150*cos(seq(0,2*pi,len=1000)), 60*sin(seq(0,2*pi,len=1000)))

	if(!FILL.COL){
		polygon(relocate_elp(elps,  90,200, 250),border=circle.col[1],lwd=lwd)
		polygon(relocate_elp(elps, 162,250, 220),border=circle.col[2],lwd=lwd)
		polygon(relocate_elp(elps, 234,250, 150),border=circle.col[3],lwd=lwd)
		polygon(relocate_elp(elps, 306,180, 125),border=circle.col[4],lwd=lwd)
		polygon(relocate_elp(elps, 378,145, 200),border=circle.col[5],lwd=lwd)
	}
	if(FILL.COL){
		RGB <- col2rgb(circle.col)/255
		ALPHA <- 0.06
		RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], alpha = ALPHA)
		RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2], alpha = ALPHA)
		RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3], alpha = ALPHA)
		RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4], alpha = ALPHA)
		RGB.ALP5 <- rgb(RGB[1, 5], RGB[2, 5], RGB[3, 5], alpha = ALPHA)
		polygon(relocate_elp(elps,  90,200, 250),border=circle.col[1],lwd=lwd,col=RGB.ALP1)
		polygon(relocate_elp(elps, 162,250, 220),border=circle.col[2],lwd=lwd,col=RGB.ALP2)
		polygon(relocate_elp(elps, 234,250, 150),border=circle.col[3],lwd=lwd,col=RGB.ALP3)
		polygon(relocate_elp(elps, 306,180, 125),border=circle.col[4],lwd=lwd,col=RGB.ALP4)
		polygon(relocate_elp(elps, 378,145, 200),border=circle.col[5],lwd=lwd,col=RGB.ALP5)
	}

	text( 50, 285, names[1],cex=cex[1])
	text(200, 415, names[2],cex=cex[1])
	text(350, 305, names[3],cex=cex[1])
	text(350,  20, names[4],cex=cex[1])
	text(100, -10, names[5],cex=cex[1])

	text( 61, 231, z["10000"], cex=cex[2], col=counts.col[1])
	text(200, 332, z["01000"], cex=cex[2], col=counts.col[1])
	text(321, 248, z["00100"], cex=cex[2], col=counts.col[1])
	text(290,  84, z["00010"], cex=cex[2], col=counts.col[1])
	text(132,  72, z["00001"], cex=cex[2], col=counts.col[1])
	text(146, 253, z["11000"], cex=cex[3], col=counts.col[1]) 
	text(123, 191, z["10100"], cex=cex[3], col=counts.col[1]) 
	text(275, 155, z["10010"], cex=cex[3], col=counts.col[1]) 
	text(137, 149, z["10001"], cex=cex[3], col=counts.col[1]) 
	text(243, 271, z["01100"], cex=cex[3], col=counts.col[1]) 
	text(175, 270, z["01010"], cex=cex[3], col=counts.col[1]) 
	text(187, 120, z["01001"], cex=cex[3], col=counts.col[1]) 
	text(286, 193, z["00110"], cex=cex[3], col=counts.col[1]) 
	text(267, 238, z["00101"], cex=cex[3], col=counts.col[1]) 
	text(228, 108, z["00011"], cex=cex[3], col=counts.col[1]) 
	text(148, 213, z["11100"], cex=cex[3], col=counts.col[1])
	text(159, 255, z["11010"], cex=cex[3], col=counts.col[1]) 
	text(171, 144, z["11001"], cex=cex[3], col=counts.col[1]) 
	text(281, 178, z["10110"], cex=cex[3], col=counts.col[1]) 
	text(143, 166, z["10101"], cex=cex[3], col=counts.col[1]) 
	text(252, 148, z["10011"], cex=cex[3], col=counts.col[1]) 
	text(205, 258, z["01110"], cex=cex[3], col=counts.col[1]) 
	text(254, 248, z["01101"], cex=cex[3], col=counts.col[1]) 
	text(211, 121, z["01011"], cex=cex[3], col=counts.col[1]) 
	text(267, 214, z["00111"], cex=cex[3], col=counts.col[1]) 
	text(170, 234, z["11110"], cex=cex[3], col=counts.col[1]) 
	text(158, 172, z["11101"], cex=cex[3], col=counts.col[1]) 
	text(212, 142, z["11011"], cex=cex[3], col=counts.col[1])
	text(263, 183, z["10111"], cex=cex[3], col=counts.col[1]) 
	text(239, 235, z["01111"], cex=cex[3], col=counts.col[1])
	text(204, 193, z["11111"], cex=cex[2], col=counts.col[1])
	text(400,   7, z["00000"], cex=cex[1], col=counts.col[1])

	if(length(include)==2) {
		text( 61, 220, z2["10000"], cex=cex[2], col=counts.col[2])
		text(200, 321, z2["01000"], cex=cex[2], col=counts.col[2])
		text(321, 237, z2["00100"], cex=cex[2], col=counts.col[2])
		text(290,  73, z2["00010"], cex=cex[2], col=counts.col[2])
		text(132,  61, z2["00001"], cex=cex[2], col=counts.col[2])
		text(146, 244, z2["11000"], cex=cex[3], col=counts.col[2]) 
		text(123, 180, z2["10100"], cex=cex[3], col=counts.col[2]) 
		text(275, 144, z2["10010"], cex=cex[3], col=counts.col[2]) 
		text(137, 143, z2["10001"], cex=cex[3], col=counts.col[2]) 
		text(243, 260, z2["01100"], cex=cex[3], col=counts.col[2]) 
		text(175, 259, z2["01010"], cex=cex[3], col=counts.col[2]) 
		text(187, 110, z2["01001"], cex=cex[3], col=counts.col[2]) 
		text(286, 186, z2["00110"], cex=cex[3], col=counts.col[2]) 
		text(267, 230, z2["00101"], cex=cex[3], col=counts.col[2]) 
		text(228,  97, z2["00011"], cex=cex[3], col=counts.col[2]) 
		text(148, 203, z2["11100"], cex=cex[3], col=counts.col[2])
		text(159, 249, z2["11010"], cex=cex[3], col=counts.col[2]) 
		text(171, 137, z2["11001"], cex=cex[3], col=counts.col[2]) 
		text(281, 171, z2["10110"], cex=cex[3], col=counts.col[2]) 
		text(143, 155, z2["10101"], cex=cex[3], col=counts.col[2]) 
		text(252, 137, z2["10011"], cex=cex[3], col=counts.col[2]) 
		text(205, 247, z2["01110"], cex=cex[3], col=counts.col[2]) 
		text(254, 242, z2["01101"], cex=cex[3], col=counts.col[2]) 
		text(211, 112, z2["01011"], cex=cex[3], col=counts.col[2]) 
		text(267, 207, z2["00111"], cex=cex[3], col=counts.col[2]) 
		text(170, 223, z2["11110"], cex=cex[3], col=counts.col[2]) 
		text(158, 162, z2["11101"], cex=cex[3], col=counts.col[2]) 
		text(212, 133, z2["11011"], cex=cex[3], col=counts.col[2])
		text(263, 172, z2["10111"], cex=cex[3], col=counts.col[2]) 
		text(239, 228, z2["01111"], cex=cex[3], col=counts.col[2])
		text(204, 182, z2["11111"], cex=cex[2], col=counts.col[2])
		text(400, -10, z2["00000"], cex=cex[1], col=counts.col[2])				

		if(show.include) {
			text(10,  7,include[1], cex=cex[1], col=counts.col[1]) 
			text(10,-10,include[2], cex=cex[1], col=counts.col[2]) 
		}
	}
	invisible()

}
