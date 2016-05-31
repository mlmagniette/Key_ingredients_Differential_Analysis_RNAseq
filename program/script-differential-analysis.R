#######################################################
## Main script for the synthetic datasets      
##
## for each proportion of full H0
## creation of ten synthetic datasets
## analysis differential with the 9 pipelines
##
## result evaluation 
## Kolmogorov-Smirnov tests 
## creation of the truly NDE sets
## AUC analysis
## TPR analysis
## FDR analysis
##########################################################
set.seed(22)

source("syntheticData.R")
source("graph-dataset.R")
source("DiffAnalysis.limmavoom.R")
source("DiffAnalysis.edgeR.R")
source("DiffAnalysis.glm.edgeR.R")
source("DiffAnalysis.DESeq.R")
source("DiffAnalysis.DESeq2.R")
sink("../results/sessionInfo.log")
sessionInfo()
sink()

# creation of the synthetic datasets 
nb.rep<-10
propHnull<-rep(seq(0.1,0.9,by=0.1),each=nb.rep)
j=1
dir.create("../synthetic_datasets")
setwd("../synthetic_datasets")
for (i in 1:length(propHnull))
{ 
    dir.create(paste("synthetic_",propHnull[i],"_",j,sep=""))
    setwd(paste("synthetic_",propHnull[i],"_",j,sep=""))
    countTable<-syntheticData(propHnull[i])
    write.table(countTable,file="synthetic.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    graphDataset(countTable)
    DiffAnalysis.limmavoom(countTable)
    DiffAnalysis.limmavoom(countTable,filter=FALSE)
    DiffAnalysis.edgeR(countTable)
    DiffAnalysis.edgeR(countTable,filter=FALSE)
    DiffAnalysis.glm.edgeR(countTable)
    DiffAnalysis.glm.edgeR(countTable,filter=FALSE)
    DiffAnalysis.DESeq(countTable)	
    DiffAnalysis.DESeq2(countTable)
    DiffAnalysis.DESeq2(countTable,filter=FALSE)
    setwd("..")
    j<-ifelse(j==nb.rep,1,j+1)
}
setwd("../program")
