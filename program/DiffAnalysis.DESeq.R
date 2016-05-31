library("DESeq")
DiffAnalysis.DESeq<-function(countTable,alpha=0.05,repbio=2)
{
    conds<-c("Bud","Leaf")
    pdf("graph-DESeq.pdf")
 
    tmp<-which(apply(countTable,1,sum)==0)
    if(length(tmp)!=0)
        countTable<-countTable[-tmp,]
    cat("Number of genes under study", nrow(countTable),"\n")
### Binomial negative test 
    cds<-newCountDataSet(countTable,factor(rep(conds,each=repbio)))
    cds<-estimateSizeFactors(cds)
    cds<-estimateDispersions(cds)
 # args par defaut : method="pooled",sharingMode="maximum",fitType="parametric"
    res<-nbinomTest(cds,condA=conds[1],condB=conds[2])
    hist(res$pval,150,main="Histogram of raw pvalues",xlab="raw pvalues")
    DE.BH<-which(res$padj<alpha)
    rownames(res)=res$id
    names(res)[grep("pval",names(res))]<-"pvalue"
    res<-res[,-1]
    res<-data.frame(res,dispersion=fitInfo(cds)$perGeneDispEsts)
    write.table(res,"DESeq_CompleteListe.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    write.table(res[DE.BH,],paste("DESeq_BH-",alpha,".txt",sep=""),row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    dev.off()
}
