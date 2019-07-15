library(limma)

set.seed(0); u <- runif(100)

### strsplit2

x <- c("ab;cd;efg","abc;def","z","")
strsplit2(x,split=";")

### removeext

removeExt(c("slide1.spot","slide.2.spot"))
removeExt(c("slide1.spot","slide"))

### printorder

printorder(list(ngrid.r=4,ngrid.c=4,nspot.r=8,nspot.c=6),ndups=2,start="topright",npins=4)
printorder(list(ngrid.r=4,ngrid.c=4,nspot.r=8,nspot.c=6))

### merge.rglist

R <- G <- matrix(11:14,4,2)
rownames(R) <- rownames(G) <- c("a","a","b","c")
RG1 <- new("RGList",list(R=R,G=G))
R <- G <- matrix(21:24,4,2)
rownames(R) <- rownames(G) <- c("b","a","a","c")
RG2 <- new("RGList",list(R=R,G=G))
merge(RG1,RG2)
merge(RG2,RG1)

### background correction

RG <- new("RGList", list(R=c(1,2,3,4),G=c(1,2,3,4),Rb=c(2,2,2,2),Gb=c(2,2,2,2)))
backgroundCorrect(RG)
backgroundCorrect(RG, method="half")
backgroundCorrect(RG, method="minimum")
backgroundCorrect(RG, offset=5)

### loessFit

x <- 1:100
y <- rnorm(100)
out <- loessFit(y,x)
f1 <- quantile(out$fitted)
r1 <- quantile(out$residual)
w <- rep(1,100)
w[1:50] <- 0.5
out <- loessFit(y,x,weights=w,method="weightedLowess")
f2 <- quantile(out$fitted)
r2 <- quantile(out$residual)
out <- loessFit(y,x,weights=w,method="locfit")
f3 <- quantile(out$fitted)
r3 <- quantile(out$residual)
out <- loessFit(y,x,weights=w,method="loess")
f4 <- quantile(out$fitted)
r4 <- quantile(out$residual)
w <- rep(1,100)
w[2*(1:50)] <- 0
out <- loessFit(y,x,weights=w,method="weightedLowess")
f5 <- quantile(out$fitted)
r5 <- quantile(out$residual)
data.frame(f1,f2,f3,f4,f5)
data.frame(r1,r2,r3,r4,r5)

### normalizeWithinArrays

RG <- new("RGList",list())
RG$R <- matrix(rexp(100*2),100,2)
RG$G <- matrix(rexp(100*2),100,2)
RG$Rb <- matrix(rnorm(100*2,sd=0.02),100,2)
RG$Gb <- matrix(rnorm(100*2,sd=0.02),100,2)
RGb <- backgroundCorrect(RG,method="normexp",normexp.method="saddle")
summary(cbind(RGb$R,RGb$G))
RGb <- backgroundCorrect(RG,method="normexp",normexp.method="mle")
summary(cbind(RGb$R,RGb$G))
MA <- normalizeWithinArrays(RGb,method="loess")
summary(MA$M)
#MA <- normalizeWithinArrays(RG[,1:2], mouse.setup, method="robustspline")
#MA$M[1:5,]
#MA <- normalizeWithinArrays(mouse.data, mouse.setup)
#MA$M[1:5,]

### normalizeBetweenArrays

MA2 <- normalizeBetweenArrays(MA,method="scale")
MA$M[1:5,]
MA$A[1:5,]
MA2 <- normalizeBetweenArrays(MA,method="quantile")
MA$M[1:5,]
MA$A[1:5,]

### unwrapdups

M <- matrix(1:12,6,2)
unwrapdups(M,ndups=1)
unwrapdups(M,ndups=2)
unwrapdups(M,ndups=3)
unwrapdups(M,ndups=2,spacing=3)

### trigammaInverse

trigammaInverse(c(1e-6,NA,5,1e6))

### lmFit, eBayes, topTable

M <- matrix(rnorm(10*6,sd=0.3),10,6)
rownames(M) <- LETTERS[1:10]
M[1,1:3] <- M[1,1:3] + 2
design <- cbind(First3Arrays=c(1,1,1,0,0,0),Last3Arrays=c(0,0,0,1,1,1))
contrast.matrix <- cbind(First3=c(1,0),Last3=c(0,1),"Last3-First3"=c(-1,1))
fit <- lmFit(M,design)
fit2 <- eBayes(contrasts.fit(fit,contrasts=contrast.matrix))
topTable(fit2)
topTable(fit2,coef=3,resort.by="logFC")
topTable(fit2,coef=3,resort.by="p")
topTable(fit2,coef=3,sort="logFC",resort.by="t")
topTable(fit2,coef=3,resort.by="B")
topTable(fit2,coef=3,lfc=1)
topTable(fit2,coef=3,p=0.2)
topTable(fit2,coef=3,p=0.2,lfc=0.5)
topTable(fit2,coef=3,p=0.2,lfc=0.5,sort="none")

designlist <- list(Null=matrix(1,6,1),Two=design,Three=cbind(1,c(0,0,1,1,0,0),c(0,0,0,0,1,1)))
out <- selectModel(M,designlist)
table(out$pref)

### marray object

#suppressMessages(suppressWarnings(gotmarray <- require(marray,quietly=TRUE)))
#if(gotmarray) {
#	data(swirl)
#	snorm = maNorm(swirl)
#	fit <- lmFit(snorm, design = c(1,-1,-1,1))
#	fit <- eBayes(fit)
#	topTable(fit,resort.by="AveExpr")
#}

### duplicateCorrelation

cor.out <- duplicateCorrelation(M)
cor.out$consensus.correlation
cor.out$atanh.correlations

### gls.series

fit <- gls.series(M,design,correlation=cor.out$cor)
fit$coefficients
fit$stdev.unscaled
fit$sigma
fit$df.residual

### mrlm

fit <- mrlm(M,design)
fit$coef
fit$stdev.unscaled
fit$sigma
fit$df.residual

# Similar to Mette Langaas 19 May 2004
set.seed(123)
narrays <- 9
ngenes <- 5
mu <- 0
alpha <- 2
beta <- -2
epsilon <- matrix(rnorm(narrays*ngenes,0,1),ncol=narrays)
X <- cbind(rep(1,9),c(0,0,0,1,1,1,0,0,0),c(0,0,0,0,0,0,1,1,1))
dimnames(X) <- list(1:9,c("mu","alpha","beta"))
yvec <- mu*X[,1]+alpha*X[,2]+beta*X[,3]
ymat <- matrix(rep(yvec,ngenes),ncol=narrays,byrow=T)+epsilon
ymat[5,1:2] <- NA
fit <- lmFit(ymat,design=X)
test.contr <- cbind(c(0,1,-1),c(1,1,0),c(1,0,1))
dimnames(test.contr) <- list(c("mu","alpha","beta"),c("alpha-beta","mu+alpha","mu+beta"))
fit2 <- contrasts.fit(fit,contrasts=test.contr)
eBayes(fit2)

### uniquegenelist

uniquegenelist(letters[1:8],ndups=2)
uniquegenelist(letters[1:8],ndups=2,spacing=2)

### classifyTests

tstat <- matrix(c(0,5,0, 0,2.5,0, -2,-2,2, 1,1,1), 4, 3, byrow=TRUE)
classifyTestsF(tstat)
classifyTestsF(tstat,fstat.only=TRUE)
limma:::.classifyTestsP(tstat)

### avereps

x <- matrix(rnorm(8*3),8,3)
colnames(x) <- c("S1","S2","S3")
rownames(x) <- c("b","a","a","c","c","b","b","b")
avereps(x)

### roast

y <- matrix(rnorm(100*4),100,4)
sigma <- sqrt(2/rchisq(100,df=7))
y <- y*sigma
design <- cbind(Intercept=1,Group=c(0,0,1,1))
iset1 <- 1:5
y[iset1,3:4] <- y[iset1,3:4]+3
iset2 <- 6:10
roast(y=y,iset1,design,contrast=2)
roast(y=y,iset1,design,contrast=2,array.weights=c(0.5,1,0.5,1))
w <- matrix(runif(100*4),100,4)
roast(y=y,iset1,design,contrast=2,weights=w)
mroast(y=y,list(set1=iset1,set2=iset2),design,contrast=2,gene.weights=runif(100))
mroast(y=y,list(set1=iset1,set2=iset2),design,contrast=2,array.weights=c(0.5,1,0.5,1))
mroast(y=y,list(set1=iset1,set2=iset2),design,contrast=2,weights=w)
mroast(y=y,list(set1=iset1,set2=iset2),design,contrast=2,weights=w,array.weights=c(0.5,1,0.5,1))
fry(y=y,list(set1=iset1,set2=iset2),design,contrast=2,weights=w,array.weights=c(0.5,1,0.5,1))
rownames(y) <- paste0("Gene",1:100)
iset1A <- rownames(y)[1:5]
fry(y=y,index=iset1A,design,contrast=2,weights=w,array.weights=c(0.5,1,0.5,1))

### camera

camera(y=y,iset1,design,contrast=2,weights=c(0.5,1,0.5,1),allow.neg.cor=TRUE,inter.gene.cor=NA)
camera(y=y,list(set1=iset1,set2=iset2),design,contrast=2,allow.neg.cor=TRUE,inter.gene.cor=NA)
camera(y=y,iset1,design,contrast=2,weights=c(0.5,1,0.5,1))
camera(y=y,list(set1=iset1,set2=iset2),design,contrast=2)
camera(y=y,iset1A,design,contrast=2)

### with EList arg

y <- new("EList",list(E=y))
roast(y=y,iset1,design,contrast=2)
camera(y=y,iset1,design,contrast=2,allow.neg.cor=TRUE,inter.gene.cor=NA)
camera(y=y,iset1,design,contrast=2)

### eBayes with trend

fit <- lmFit(y,design)
fit <- eBayes(fit,trend=TRUE)
topTable(fit,coef=2)
fit$df.prior
fit$s2.prior
summary(fit$s2.post)

y$E[1,1] <- NA
y$E[1,3] <- NA
fit <- lmFit(y,design)
fit <- eBayes(fit,trend=TRUE)
topTable(fit,coef=2)
fit$df.residual[1]
fit$df.prior
fit$s2.prior
summary(fit$s2.post)

### eBayes with robust

fitr <- lmFit(y,design)
fitr <- eBayes(fitr,robust=TRUE)
summary(fitr$df.prior)
topTable(fitr,coef=2)
fitr <- eBayes(fitr,trend=TRUE,robust=TRUE)
summary(fitr$df.prior)
topTable(fitr,coef=2)

### voom

y <- matrix(rpois(100*4,lambda=20),100,4)
design <- cbind(Int=1,x=c(0,0,1,1))
v <- voom(y,design)
names(v)
summary(v$E)
summary(v$weights)

### goana

EB <- c("133746","1339","134","1340","134083","134111","134147","134187","134218","134266",
"134353","134359","134391","134429","134430","1345","134510","134526","134549","1346",
"134637","1347","134701","134728","1348","134829","134860","134864","1349","134957",
"135","1350","1351","135112","135114","135138","135152","135154","1352","135228",
"135250","135293","135295","1353","135458","1355","1356","135644","135656","1357",
"1358","135892","1359","135924","135935","135941","135946","135948","136","1360",
"136051","1361","1362","136227","136242","136259","1363","136306","136319","136332",
"136371","1364","1365","136541","1366","136647","1368","136853","1369","136991",
"1370","137075","1371","137209","1373","137362","1374","137492","1375","1376",
"137682","137695","137735","1378","137814","137868","137872","137886","137902","137964")
go <- goana(fit,FDR=0.8,geneid=EB)
topGO(go,n=10,truncate.term=30)
topGO(go,n=10,truncate.term=30,sort="down")
