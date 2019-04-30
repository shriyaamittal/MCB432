setwd('./GSE63310_RAW/')

files <- c("GSM1545535_10_6_5_11.txt","GSM1545536_9_6_5_11.txt",   "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt",   "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",   "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt",   "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")

library(limma)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")

library(edgeR)

x <- readDGE(files, columns=c(1,3))
class(x)

dim(x)
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

source("https://bioconductor.org/biocLite.R")
biocLite("Mus.musculus")
library("Mus.musculus")

geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
##'select()' returned 1:many mapping between keys and columns
head(genes)
mat <- match(geneid, genes$ENTREZID)
genes <- genes[mat,]
genes[genes$ENTREZID %in% dup,][1:5,]
## Error in genes$ENTREZID %in% dup : object 'dup' not found

genes[1:5,]
mat <- match(geneid, genes$ENTREZID)
genes <- genes[mat,]
genes[genes$ENTREZID dup,][1:5,]

genes[5360,]

x$genes <- genes
x

genes[genes$ENTREZID %in% dup,][1:5,]
## Error in genes$ENTREZID %in% dup : object 'dup' not found

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
head(cpm)
head(lcpm)

table(rowSums(x$counts==0)==9)
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

library(RColorBrewer)

## Figure 1

lcpm.cutoff<- log2(10/M +/L)
nsamples<-ncol(x)
col<-brewer.pal(nsamples,"Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]),col=col[1],lwd=2,ylim=c(0,0.26),las=2,main="",xlab="")
title(main="A. Raw data",xlab="Log-cpm")
abline(v=lcpm.cutoff,lty=3)
for (i in 2:nsamples){  
        den<-density(lcpm[,i])  
        lines(den$x,den$y,col=col[i],lwd=2)
}
legend("topright", samplenames,text.col=col,bty="n")
lcpm<-cpm(x,log=TRUE)
plot(density(lcpm[,1]),col=col[1],lwd=2,ylim=c(0,0.26),las=2,main="",xlab="")
title(main="B. Filtered data",xlab="Log-cpm")
abline(v=lcpm.cutoff,lty=3)
for(i in 2:nsamples){  
        den<-density(lcpm[,i])  
        lines(den$x, den$y,col=col[i],lwd=2)
}
legend("topright", samplenames,text.col=col,bty="n")

## Figure 2 

x<-calcNormFactors(x,method="TMM")x$samples$norm.factors
x2<-x
x2$samples$norm.factors<-1x2$counts[,1]<-ceiling(x2$counts[,1]*0.05)
x2$counts[,2]<-x2$counts[,2]*5

par(mfrow=c(1,2))

lcpm<-cpm(x2,log=TRUE)

boxplot(lcpm,las=2,col=col,main="")

title(main="A. Example: Unnormalised data",ylab="Log-cpm")

x2<-calcNormFactors(x2)

x2$samples$norm.factors

lcpm<-cpm(x2,log=TRUE)
boxplot(lcpm,las=2,col=col,main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")


## Figure 3

cpm<-cpm(x,log=TRUE)
par(mfrow=c(1,2))
col.group<-group
levels(col.group)<-brewer.pal(nlevels(col.group),"Set1")
col.group<-as.character(col.group)
col.lane<-lane
levels(col.lane)<-brewer.pal(nlevels(col.lane),"Set2")
col.lane<-as.character(col.lane)
plotMDS(lcpm,labels=group,col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm,labels=lane,col=col.lane,dim=c(3,4))
title(main="B. Sequencing lanes")


design<-model.matrix(~0+group+lane)
colnames(design)<-gsub("group","",colnames(design))
design

contr.matrix<-makeContrasts(
	BasalvsLP= Basal-LP,
	BasalvsML= Basal - ML,
	LPvsML= LP - ML,
	levels=colnames(design))
contr.matrix

## Figure 4
v<-voom(x, design,plot=TRUE)
vfit<-lmFit(v, design)
vfit<-contrasts.fit(vfit,contrasts=contr.matrix)
efit<-eBayes(vfit)
plotSA(efit)


summary(decideTests(efit))
tfit<-treat(vfit,lfc=1)
dt<-decideTests(tfit)
summary(dt)

de.common<-which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common],n=20)

## Figure 5
vennDiagram(dt[,1:2],circle.col=c("turquoise","salmon"))
write.fit(tfit, dt,file="results.txt")

basal.vs.lp<-topTreat(tfit,coef=1,n=Inf)
basal.vs.ml<-topTreat(tfit,coef=2,n=Inf)

## Figure 6
plotMD(tfit,column=1,status=dt[,1],main=colnames(tfit)[1],xlim=c(-8,13))
## glMDPlot(tfit,coef=1,status=dt[,1],main=colnames(tfit)[1],counts=x$counts,samples=colnames(x),anno=x$genes,groups=group,id.column="ENTREZID",display.columns=c("SYMBOL","ENTREZID"),search.by="SYMBOL",launch=FALSE)


library(gplots)

## Figure 7
basal.vs.lp.topgenes<-basal.vs.lp$ENTREZID[1:100]
i<-which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol<-colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,],scale="row",labRow=v$genes$SYMBOL[i],labCol=group,   col=mycol,trace="none",density.info="none",   margin=c(8,6),lhei=c(2,10),dendrogram="column")
