load("../results/TPR.with.fixed.fdr.inter.RData")
dat<-read.table("../real_datasets/Arabidopsis_dataset.txt",h=T)
id.DE=dat$ID[grep("DE",dat$ID,fixed=TRUE)]
nb.DE <- length(id.DE)

res.TPR.with.fixed.fdr<-res.TPR.with.fixed.fdr/nb.DE
tiff("../results/TPR-with-fixed-fdr-inter.tiff",width=2250,height=2250,pointsize=36)
plot(x,boxplot(res.TPR.with.fixed.fdr[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[1],xlab=titre.x,ylab="TPR",ylim=c(0,1),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
{
  tmp<-boxplot(res.TPR.with.fixed.fdr[,j]~rep(x,each=10),plot=F)$stats[3,]
  print(j)
  lines(x[1:length(tmp)],tmp,pch=pch.index[j],type="b",lwd=2,lty=type.line[j])
}
legend(x="bottomleft",legend=methods[type.line==1],lty=1,pch=pch.index[type.line==1],bty="n",lwd=2)
dev.off()

tiff(paste("../results/variability-TPR-with-fixed-fdr-inter.tiff",sep=""),width=2250,height=2250,pointsize=36)
order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","DESeq2 filtered","glm edgeR filtered","limma-voom filtered")
res.TPR.with.fixed.fdr<-res.TPR.with.fixed.fdr[,match(order.methods,methods)]
par(mfrow=c(3,3))
sapply(1:9,function(j) plot(rep(x,each=10),res.TPR.with.fixed.fdr[,j],xlab=titre.x,ylab="FDR",main=order.methods[j],ylim=range(res.TPR.with.fixed.fdr)))
dev.off()


######
tiff("../results/TPR-with-fixed-fdr-union.tiff",width=2250,height=2250,pointsize=36)
load("../results/TPR.with.fixed.fdr.union.RData")
res.TPR.with.fixed.fdr<-res.TPR.with.fixed.fdr/nb.DE

plot(x,boxplot(res.TPR.with.fixed.fdr[,1]~rep(x,each=10),plot=F)$stats[4,],pch=pch.index[1],xlab=titre.x,ylab="TPR",ylim=c(0,1),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
{
  tmp<-boxplot(res.TPR.with.fixed.fdr[,j]~rep(x,each=10),plot=F)$stats[4,]
  print(j)
  lines(x[1:length(tmp)],tmp,pch=pch.index[j],type="b",lwd=2,lty=type.line[j])
}
legend(x="bottomleft",legend=methods[type.line==1],lty=1,pch=pch.index[type.line==1],bty="n",lwd=2)
dev.off()


tiff(paste("../results/variability-TPR-with-fixed-fdr-union.tiff",sep=""),width=2250,height=2250,pointsize=36)
order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","DESeq2 filtered","glm edgeR filtered","limma-voom filtered")
res.TPR.with.fixed.fdr<-res.TPR.with.fixed.fdr[,match(order.methods,methods)]
par(mfrow=c(3,3))
sapply(1:9,function(j) plot(rep(x,each=10),res.TPR.with.fixed.fdr[,j],xlab=titre.x,ylab="FDR",main=order.methods[j],ylim=range(res.TPR.with.fixed.fdr)))
dev.off()
