#####################################################################
## TPR evaluation based on the truly DE genes
## for each differential analysis and each synthetic dataset
## calculation of the number of DE genes among the genes declared DE (FDR 5%)
## output : a RData object and a graph
######################################################################

## DE set : identical for all the methods
dat<-read.table("../real_datasets/Arabidopsis_dataset.txt",h=T)
id.DE=dat$ID[grep("DE",dat$ID,fixed=TRUE)]
nb.DE <- length(id.DE)

## TPR evaluation for each synthetic dataset
file<- system("ls ../synthetic_datasets/synt*/*BH*",T)
TPR.test<-numeric()
for (i in 1:length(file))
  {
    print(i)
    data<-read.table(file[i],h=T)
    id<-rownames(data)
    DE.found=length(grep(".DE",id,fixed=TRUE))
    TPR.test[i]=DE.found/nb.DE
  }

## results store in a matrix and save in the folder results 
TPR.res=matrix(TPR.test,nrow=length(file)/length(methods),byrow=T,ncol=length(methods))
colnames(TPR.res)<-methods
save(TPR.res,file="../results/TPR.RData")

## graph of the mean value of the TPR for each proportion of full H0 dataset
tiff("../results/TPR.tiff",width=2250,height=2250,pointsize=36)

plot(x,boxplot(TPR.res[,1]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[1],xlab=titre.x,ylab="TPR",ylim=c(0,1),type="b",xaxt="n",lwd=2,lty=type.line[1])
axis(1,at=x,labels=x)

for (j in 2:length(methods))
  lines(x,boxplot(TPR.res[,j]~rep(x,each=10),plot=F)$stats[3,],pch=pch.index[j],type="b",lwd=2,lty=type.line[j])

legend(x="bottomleft",legend=methods[which(type.line==1)],lty=1,pch=pch.index[which(type.line==1)],bty="n",lwd=2)
dev.off()
## graph about the variability
tiff(paste("../results/variability-TPR.tiff",sep=""),width=2250,height=2250,pointsize=36)

order.methods<-c("DESeq","edgeR","edgeR filtered","DESeq2","glm edgeR","limma-voom","DESeq2 filtered","glm edgeR filtered","limma-voom filtered")

TPR.res<-TPR.res[,match(order.methods,methods)]


par(mfrow=c(3,3))
sapply(1:9,function(j) plot(rep(x,each=10),TPR.res[,j],xlab=titre.x,ylab="TPR",main=order.methods[j],ylim=c(0,1)))
dev.off()
