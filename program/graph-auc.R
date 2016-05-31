## AUC graph based on NDE.inter
load("../results/auc.inter.RData")

tiff("../results/auc-inter.tiff",width=2250,height=2250,pointsize=36)
plot(x,boxplot(res.auc[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[1],ylab="AUC",xlab=titre.x,ylim=c(0,1),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
{
  tmp<-boxplot(res.auc[,j]~rep(x,each=10),plot=F)$stats[3,]
  print(j)
  lines(x[1:length(tmp)],tmp,pch=pch.index[j],type="b",lwd=2,lty=type.line[j])
}
legend(x="bottomright",legend=methods[type.line==1],lty=1,pch=pch.index[type.line==1],bty="n",lwd=2)
dev.off()

## AUC graph based on NDE.union
load("../results/auc.union.RData")

tiff("../results/auc-union.tiff",width=2250,height=2250,pointsize=36)
plot(x,boxplot(res.auc[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[1],ylab="AUC",xlab=titre.x,ylim=c(0,1),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
{
  tmp<-boxplot(res.auc[,j]~rep(x,each=10),plot=F)$stats[3,]
  print(j)
  lines(x[1:length(tmp)],tmp,pch=pch.index[j],type="b",lwd=2,lty=type.line[j])
}
legend(x="bottomright",legend=methods[type.line==1],lty=1,pch=pch.index[type.line==1],bty="n",lwd=2)
dev.off()




