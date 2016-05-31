load("../results/FDR.inter.RData")

tiff("../results/FDR.inter.tiff",width=2250,height=2250,pointsize=36)
plot(x,boxplot(FDR.res[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[1],xlab=titre.x,ylab="FDR",ylim=c(0,.25),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
{
  tmp<-boxplot(FDR.res[,j]~rep(x,each=10),plot=F)$stats[3,]
  print(j)
  lines(x[1:length(tmp)],tmp,pch=pch.index[j],type="b",lwd=2,lty=type.line[j])
}
abline(h=alpha,lty=4)
legend(x="topleft",legend=methods[type.line==1],lty=1,pch=pch.index[type.line==1],bty="n",lwd=2)
dev.off()

tiff(paste("../results/variability-FDR.inter.tiff",sep=""),width=2250,height=2250,pointsize=36)
order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","DESeq2 filtered","glm edgeR filtered","limma-voom filtered")
FDR.res<-FDR.res[,match(order.methods,methods)]
par(mfrow=c(3,3))
sapply(1:9,function(j) plot(rep(x,each=10),FDR.res[,j],xlab=titre.x,ylab="FDR",main=order.methods[j],ylim=range(FDR.res)))
dev.off()

######
load("../results/FDR.union.RData")

tiff("../results/FDR.union.tiff",width=2250,height=2250,pointsize=36)
plot(x,boxplot(FDR.res[,1]~rep(x,each=10),plot=F)$stats[4,],pch=pch.index[1],xlab=titre.x,ylab="FDR",ylim=c(0,.25),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
{
  tmp<-boxplot(FDR.res[,j]~rep(x,each=10),plot=F)$stats[4,]
  print(j)
  lines(x[1:length(tmp)],tmp,pch=pch.index[j],type="b",lwd=2,lty=type.line[j])
}
legend(x="topleft",legend=methods[type.line==1],lty=1,pch=pch.index[type.line==1],bty="n",lwd=2)
abline(h=alpha,lty=4)
dev.off()

tiff(paste("../results/variability-FDR.union.tiff",sep=""),width=2250,height=2250,pointsize=36)
order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","DESeq2 filtered","glm edgeR filtered","limma-voom filtered")
FDR.res<-FDR.res[,match(order.methods,methods)]
par(mfrow=c(3,3))
sapply(1:9,function(j) plot(rep(x,each=10),FDR.res[,j],xlab=titre.x,ylab="FDR",main=order.methods[j],ylim=range(FDR.res)))
dev.off()
