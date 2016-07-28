## parameters and functions
origin<-getwd()
source("DiffAnalysis.limmavoom.R")
source("DiffAnalysis.edgeR.R")
source("DiffAnalysis.glm.edgeR.R")
source("DiffAnalysis.DESeq.R")
source("DiffAnalysis.DESeq2.R")

## dataset
dat <- read.table("../real_datasets/Arabidopsis_dataset.txt",h=T)
removeID<-read.table("../real_datasets/removeAGI.txt")
dat<-dat[!is.element(dat[,1],removeID[,1]),]
countTable<- dat[, c(4, 6, 5, 7)]
rownames(countTable)<-dat$ID

## Differential analyses 
setwd("../real_datasets/full_H0_dataset")
DiffAnalysis.limmavoom(countTable)
DiffAnalysis.limmavoom(countTable,filter=FALSE)
DiffAnalysis.edgeR(countTable)
DiffAnalysis.edgeR(countTable,filter=FALSE)
DiffAnalysis.glm.edgeR(countTable)
DiffAnalysis.glm.edgeR(countTable,filter=FALSE)
DiffAnalysis.DESeq(countTable)	
DiffAnalysis.DESeq2(countTable)
DiffAnalysis.DESeq2(countTable,filter=FALSE)

## type-I error rate
file<-system("ls *Comple*",T)
typeI.error<-numeric()
alpha=0.05
for (i in 1:length(file))
{
    tmp<-read.table(file[i],sep="\t",h=T)
    typeI.error[i]=100*sum(tmp$pvalue<=alpha)/nrow(tmp)
}
   ## Type-I error (%) for each method
print(data.frame(file,round(typeI.error,2)))

###########################################################################
# Evaluation of the distribution of the raw pvalues
# Kolmogorov-Smirnov test 
###########################################################################
file<- system("ls *Com*",T)
pvalue.ks<-numeric()
for (i in 1:length(file))
{
    data<-read.table(file[i],h=T)
    column.pvalue<-grep("pvalue",names(data))
    ks.res<-ks.test(unique(data[,column.pvalue]),"punif")
    pvalue.ks[i]=ks.res$p.value
}

## The 9 pvalues are equal to 0
print(pvalue.ks)
## The nine tests are rejected, the distribution is not uniform
setwd(origin)
