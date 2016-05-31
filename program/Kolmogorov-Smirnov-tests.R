###########################################################################
# objective: Kolmogorov-Smirnov tests on genes of the full H0 dataset
# 100 resampling of 1000 genes for each synthetic dataset  
#
# output: a RData object with the pvalues and the test statistic values
# saved in the result folder
# output: graph of the median of the Kolmogorov-Smirnov test statistics for various proportions of genes of the full H0 dataset saved in the result folder  
#######################################################################
file<- system("ls ../synthetic_datasets/synt*/*Com*",T)

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
# 67.07%

save(pvalue.ks,stat.ks,file="../results/KS_test.RData")

pvalue.ks.res<-matrix(unlist(lapply(pvalue.ks,mean)),nrow=length(file)/length(methods),byrow=T,ncol=length(methods))
colnames(pvalue.ks.res)<-methods

## Work on the mean calculated over the JJ resamplings 
stat.ks.res<-matrix(unlist(lapply(stat.ks,mean)),nrow=length(file)/length(methods),byrow=T,ncol=length(methods))
colnames(stat.ks.res)<-methods

## main graphic 
tiff(paste("../results/stat.ks.n=",n,"_nb_reech=",JJ,"_median-BW.tiff",sep=""),width=2250,height=2250,pointsize=36)
par(mfrow=c(1,1))
plot(x,boxplot(stat.ks.res[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[1],xlab=titre.x,ylab="statistics of Kolmogorov-Smirnov",ylim=c(0,.9),type="b",xaxt="n",lwd=4,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
    lines(x,boxplot(stat.ks.res[,j]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[j],type="b",lwd=4,lty=type.line[j])

legend(x="topleft",legend=methods[which(type.line==1)],lty=1,pch=pch.index[which(type.line==1)],bty="n",lwd=4)
dev.off()

## supplementary graph to show the variability of limma-voom
load("../results/KS_test.RData")

stat.ks.res<-matrix(unlist(lapply(stat.ks,mean)),nrow=length(file)/length(methods),byrow=T,ncol=length(methods))
colnames(stat.ks.res)<-methods

tiff("../results/varibility-stat.ks.tiff",width=2250,height=2250,pointsize=36)

order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","DESeq2 filtered","glm edgeR filtered","limma-voom filtered")

stat.ks.res<-stat.ks.res[,match(order.methods,methods)]
par(mfrow=c(3,3))
sapply(1:9,function(j) plot(rep(x,each=10),stat.ks.res[,j],xlab=titre.x,ylab="statistics of Kolmogorov-Smirnov",main=order.methods[j],ylim=range(stat.ks.res)))
dev.off()
