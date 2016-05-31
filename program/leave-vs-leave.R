origin<-getwd()
source("DiffAnalysis.limmavoom.R")
source("DiffAnalysis.edgeR.R")
source("DiffAnalysis.glm.edgeR.R")
source("DiffAnalysis.DESeq.R")
source("DiffAnalysis.DESeq2.R")

## data
dat <- read.table("../real_datasets/Arabidopsis_dataset.txt",h=T)
countTable<- dat[, c(4, 6, 5, 7)]
rownames(countTable)<-dat$ID

## Differential analysis
setwd("../real_datasets/full_H0_dataset/Leave_vs_Leave_dataset")
DiffAnalysis.limmavoom(countTable)
DiffAnalysis.limmavoom(countTable,filter=FALSE)
DiffAnalysis.edgeR(countTable)
DiffAnalysis.edgeR(countTable,filter=FALSE)
DiffAnalysis.glm.edgeR(countTable)
DiffAnalysis.glm.edgeR(countTable,filter=FALSE)
DiffAnalysis.DESeq(countTable)	
DiffAnalysis.DESeq2(countTable)
DiffAnalysis.DESeq2(countTable,filter=FALSE)

file<-system("ls *Comple*",T)
alpha=0.05
DE.FDR<-character()
for (i in 1:length(file))
{
    tmp<-read.table(file[i],sep="\t",h=T)
    DE.FDR<-append(DE.FDR,rownames(tmp)[which(tmp$padj<=alpha)])
}

## 17 genes are removed
write.table(unique(DE.FDR),file="../../removeAGI.txt",sep="\n",quote=F,col.names=F,row.names=F)

setwd(origin)
