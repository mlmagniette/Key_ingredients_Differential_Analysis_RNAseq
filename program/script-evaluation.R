##############################################################################
## global parameters for the graphs
##############################################################################
alpha=0.05

methods<-c("DESeq2","DESeq","edgeR","DESeq2 filtered","edgeR filtered", "glm edgeR filtered","limma-voom filtered","glm edgeR","limma-voom")

nb.methods=length(methods)

titre.x<-expression(paste("proportion of full ",H[0]," dataset",sep=""))

x<-seq(0.1,0.9,by=0.1)

type.line<-c(1,1,1,2,2,2,2,1,1)

pch.index<-c(1,2,3,1,3,4,5,4,5)

##############################################################################
## evaluation of the results of the differential analyses 
##############################################################################
source("Kolmogorov-Smirnov-tests.R")
source("NDE-sets.R")

source("auc.R")
load("../results/NDE.inter.RData")
auc(NDE.inter,"../results/auc.inter.RData")
load("../results/NDE.union.RData")
auc(NDE.union,"../results/auc.union.RData")
source("graph-auc.R")

source("FDR.R")
load("../results/NDE.inter.RData")
FDR(NDE.inter,alpha,"../results/FDR.inter.RData")
load("../results/NDE.union.RData")
FDR(NDE.union,alpha,"../results/FDR.union.RData")
source("graph-fdr.R")


source("TPR.R")



source("TPR-with-fixed-fdr.R")
load("../results/NDE.inter.RData")
TPR.with.fixed.fdr(NDE.inter,0.05,"../results/TPR.with.fixed.fdr.inter.RData")
load("../results/NDE.union.RData")
TPR.with.fixed.fdr(NDE.inter,0.05,"../results/TPR.with.fixed.fdr.union.RData")
source("graph-TPR-with-fixed-fdr.R")

