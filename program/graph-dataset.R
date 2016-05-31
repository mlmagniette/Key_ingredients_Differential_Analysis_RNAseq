library(FactoMineR)
graphDataset<-function(countTable)
{
    pdf("graphSyntheticDataset.pdf")
    barplot(apply(countTable,2,sum))
    
    plot(density(log2(countTable[,1]+1)),main="density plot of the raw counts after log2 transformation",ylim=c(0,0.3))
	    for (i in 2:ncol(countTable)){
	        lines(density(log2(countTable[,i]+1)),lwd=2,col=i)
	    }
     legend("topright",legend=names(countTable),lty=1,col=1:4)

    
    data.PCA<-cbind.data.frame(c("FB","FB","L","L"),sqrt(t(countTable))) 
    res.pca = PCA(data.PCA, quali.sup=1:1, ncp=2, graph=F)
    plot.PCA(res.pca)
    
    x <- apply(countTable[, 1:2], 1, mean)
    y <- apply(countTable[, 3:4], 1, mean)
    plot(log2(x), log2(y), xlab = expression(paste(log[2]("FB"))), 
         ylab = expression(paste(log[2]("L"))))
    points(log2(x)[grep(".Hnull", rownames(countTable))], log2(y)[grep(".Hnull", 
                                                                  rownames(countTable))], pch = 20, col = "cyan", cex = 0.75)
    points(log2(x)[grep("DE", rownames(countTable))], log2(y)[grep("DE", 
                                                              rownames(countTable))], pch = 20, col = 2)
    
    legend("topleft", bty = "n", col = c("cyan", 2), 
           pch = c(20, 20), c(expression(paste(H[0])), "DE"), cex = 1.5, pt.cex = c(0.75, 1, 1))
    dev.off()
}
