###########################################################################
# objective: Kolmogorov-Smirnov tests on genes of the full H0 dataset
# 100 resampling of 1000 genes for each synthetic dataset  
#
# output: a RData object with the pvalues and the test statistic values
# saved in the result folder
# output: graph of the median of the Kolmogorov-Smirnov test statistics for various proportions of genes of the full H0 dataset saved in the result folder  
#######################################################################
file<- system("ls ../synthetic_datasets/synt*/*Com*",T)
## remove filtered DESeq2
file<-file[-grep("filter-DESeq2_CompleteListe.txt",file)]
## remove  filtered DESeq2 from methods
methods.KS<-c("DESeq2","DESeq","edgeR","edgeR filtered", "glm edgeR filtered","limma-voom filtered","glm edgeR","limma-voom")
method.nb<-length(methods.KS)
type.line.KS<-c(1,1,1,2,2,2,1,1)

pch.index.KS<-c(1,2,3,3,4,5,4,5)
color.KS<-rep(1,9)
##color.KS<-c("cyan","blue","purple","purple","magenta","orange","magenta","orange")
## 
n=1000 # number of genes 
JJ=100 # number of sampling
pvalue.ks<-list()
stat.ks<-list()
for (i in 1:length(file))
  {
    print(i)
    data<-read.table(file[i],h=T)
    column.pvalue<-grep("pvalue",names(data))
    id<-rownames(data)
    set.seed(43)
    tmp1<-numeric()
    tmp2<-numeric()
    for (jj in 1:JJ)
    {
	ks.res<-ks.test(data[sample(grep("Hnull",id),n,replace=FALSE),column.pvalue],"punif")
        tmp1[jj]=ks.res$p.value
        tmp2[jj]=ks.res$statistic
     }
     pvalue.ks[[i]]=tmp1
     stat.ks[[i]]=tmp2
}

cat("Percentage of rejected tests after Bonferroni adjustment \n")
print(100*sum(p.adjust(unlist(pvalue.ks),"bonferroni")<=alpha)/(length(file) *JJ))
# 69.34 %

save(pvalue.ks,stat.ks,file="../results/KS_test.RData")

pvalue.ks.res<-matrix(unlist(lapply(pvalue.ks,mean)),nrow=length(file)/method.nb,byrow=T,ncol=method.nb)
colnames(pvalue.ks.res)<-methods.KS

## Work on the mean calculated over the JJ resamplings 
stat.ks.res<-matrix(unlist(lapply(stat.ks,mean)),nrow=length(file)/method.nb,byrow=T,ncol=method.nb)
colnames(stat.ks.res)<-methods.KS

## main graphic 
jpeg("../results/ks.jpg")
par(mfrow=c(1,1))
plot(x,boxplot(stat.ks.res[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index.KS[1],xlab=titre.x,ylab="statistics of Kolmogorov-Smirnov",ylim=range(stat.ks),type="b",xaxt="n",lwd=4,lty=type.line.KS[1],col=color.KS[1])
axis(1,at=x,labels=x)

for (j in 2:method.nb)
    lines(x,boxplot(stat.ks.res[,j]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index.KS[j],type="b",lwd=4,lty=type.line.KS[j],col=color.KS[j])

legend(x="topleft",legend=methods[which(type.line==1)],lty=1,pch=pch.index.KS[which(type.line==1)],bty="n",lwd=4),col=color.KS[which(type.line==1)])
dev.off()

## supplementary graph to show the variability of limma-voom
load("../results/KS_test.RData")

stat.ks.res<-matrix(unlist(lapply(stat.ks,mean)),nrow=length(file)/method.nb,byrow=T,ncol=method.nb)
colnames(stat.ks.res)<-methods.KS

tiff("../results/varibility-stat.ks.tiff",width=2250,height=2250,pointsize=36)

order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","glm edgeR filtered","limma-voom filtered")

stat.ks.res<-stat.ks.res[,match(order.methods,methods.KS)]
par(mfrow=c(3,3))
sapply(1:6,function(j) plot(rep(x,each=10),stat.ks.res[,j],xlab=titre.x,ylab="statistics of Kolmogorov-Smirnov",main=order.methods[j],ylim=range(stat.ks.res)))
plot(rep(x,each=10),axes=FALSE,xlab="",ylab="",type="n")
sapply(7:8,function(j) plot(rep(x,each=10),stat.ks.res[,j],xlab=titre.x,ylab="statistics of Kolmogorov-Smirnov",main=order.methods[j],ylim=range(stat.ks.res)))
dev.off()
