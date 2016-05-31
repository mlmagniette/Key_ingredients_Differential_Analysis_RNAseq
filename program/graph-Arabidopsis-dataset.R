library(FactoMineR)
setwd("../real_datasets/")
pdf("graph_Arabidopsis_dataset.pdf")

countTable <- read.table("Arabidopsis_dataset.txt",h=T)
rownames(countTable)<-countTable[,1]
countTable <-countTable[,-1]
## library size
barplot(apply(countTable,2,sum),main="library size")

## density of the raw counts after log2 transformation
plot(density(log2(countTable[,1]+1)),main="density plot of the raw counts after log2 transformation",ylim=c(0,0.3))
for (i in 2:ncol(countTable)){
    lines(density(log2(countTable[,i]+1)),lwd=2,col=i)
}
legend("topright",legend=names(countTable),lty=1,col=1:6,lwd=2)

## PCA
data.PCA<-cbind.data.frame(c(rep("FB",2),rep("L",4)),sqrt(t(countTable))) 
res.pca = PCA(data.PCA, quali.sup=1:1, ncp=2, graph=F)
plot.PCA(res.pca)

## expression
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(log2(countTable+1), lower.panel = panel.cor,diag.panel=panel.hist,pch=".")

dev.off()
setwd("../program")
