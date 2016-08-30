#####################################################################
## For each synthetic dataset and eac method
## visualisation of the raw pvalues histogram and
## of the raw pvalue histogram of the genes coming from H1-rich dataset in grey 
## visualisation of the quantile distribution
## dashed lines is an estimation on the set [0.8 ; 1] of the quantile
## of H1-rich raw pvalues as a linear function of  the quantile
## of full H0 raw pvalues
## DESeq results on the synthetic dataset 0.4_5 were used to draw the Figure S1
#####################################################################

file<- system("ls ../synthetic_datasets/synt*/*Com*",T)
## remove filtered DESeq2
file<-file[-grep("filter-DESeq2_CompleteListe.txt",file)]

pdf("../results/pvalue-distribution.pdf")
for (i in 1:length(file))
{
    par(mfrow=c(1,2),oma = c(0, 0, 2, 0))
    
    print(i)
    data<-read.table(file[i],h=T)
    id<-rownames(data)
    column.pvalue<-grep("pvalue",names(data))
    hist(data[,column.pvalue],100,main="",proba=TRUE,xlab="pvalue")
    hist(data[-grep(".Hnull",id),column.pvalue],100,add=TRUE,col="grey",proba=TRUE)
    h0<-data[grep(".Hnull",id),column.pvalue]
    h1<-data[-grep(".Hnull",id),column.pvalue]
    q.h0<-quantile(h0,probs=seq(0.8,1,0.001))
    q.h1<-quantile(h1,probs=seq(0.8,1,0.001))

    qqplot(h0,h1,xlab=expression(paste("quantiles of full ",H[0]," raw pvalues",sep="")),ylab=expression(paste("quantiles of ",H[1],"-rich raw pvalues",sep="")),type="l")
    abline(c(-1, 1)*coefficients(lm(I(q.h1-1)~I(q.h0-1)-1)) + c(1, 0),lty=2)
    tmp<-strsplit(file[i],"../synthetic_datasets/synthetic_")
    strsplit(tmp[[1]][2],"/")
    mtext(paste("synthetic_datasets", strsplit(tmp[[1]][2],"/")[[1]][1]),outer=TRUE)
    tmp1=strsplit(tmp[[1]][2],"/")[[1]][2]
    mtext(strsplit(tmp1,"_")[[1]][1],outer=TRUE,line=-1)

    }
dev.off()

## fig S1
tiff("../results/S1.tiff",width=2250,height=2250,pointsize=36)
i=282
par(mfrow=c(1,2))

print(i)
data<-read.table(file[i],h=T)
id<-rownames(data)
column.pvalue<-grep("pvalue",names(data))
hist(data[,column.pvalue],100,main="",proba=TRUE,xlab="pvalue")
hist(data[-grep(".Hnull",id),column.pvalue],100,add=TRUE,col="grey",proba=TRUE)
h0<-data[grep(".Hnull",id),column.pvalue]
h1<-data[-grep(".Hnull",id),column.pvalue]
q.h0<-quantile(h0,probs=seq(0.8,1,0.001))
q.h1<-quantile(h1,probs=seq(0.8,1,0.001))

qqplot(h0,h1,xlab=expression(paste("quantiles of full ",H[0]," raw pvalues",sep="")),ylab=expression(paste("quantiles of ",H[1],"-rich raw pvalues",sep="")),type="l")
abline(c(-1, 1)*coefficients(lm(I(q.h1-1)~I(q.h0-1)-1)) + c(1, 0),lty=2)
dev.off()



