source("DiffAnalysis.limmavoom.R")
source("DiffAnalysis.edgeR.R")
source("DiffAnalysis.glm.edgeR.R")
source("DiffAnalysis.DESeq.R")
source("DiffAnalysis.DESeq2.R")

dat <- read.table("../real_datasets/Arabidopsis_dataset.txt",h=T)

countTable <- dat[, 2:5]
rownames(countTable)<-dat$ID

setwd("../real_datasets/H1_rich_dataset")
DiffAnalysis.limmavoom(countTable)
DiffAnalysis.limmavoom(countTable,filter=FALSE)
DiffAnalysis.edgeR(countTable)
DiffAnalysis.edgeR(countTable,filter=FALSE)
DiffAnalysis.glm.edgeR(countTable)
DiffAnalysis.glm.edgeR(countTable,filter=FALSE)
DiffAnalysis.DESeq(countTable)	
DiffAnalysis.DESeq2(countTable)
DiffAnalysis.DESeq2(countTable,filter=FALSE)
setwd("../../program")
