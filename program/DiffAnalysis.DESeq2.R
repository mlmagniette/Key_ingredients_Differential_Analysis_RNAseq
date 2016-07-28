library("DESeq2")
DiffAnalysis.DESeq2<-function(countTable,filter=TRUE,alpha=0.05,repbio=2)
{
    conds<-c("Bud","Leaf")
    if(filter)
        pdf("graph-DESeq2-filtered.pdf")
    else
        pdf("graph-DESeq2.pdf")

    tmp<-which(apply(countTable,1,sum)==0)
    if(length(tmp)!=0)
        countTable<-countTable[-tmp,]
    
    cat("Number of genes under study", nrow(countTable),"\n")   
 
    colData<-data.frame(row.names=names(countTable), group = rep(conds,each=repbio), batch=rep(paste("batch",1:2,sep=""),2)) 
    cds<-DESeqDataSetFromMatrix(countData=countTable,
                                     colData=colData,design=~batch+group)
    cds<-DESeq(cds)
    res<-results(cds,independentFiltering=filter,alpha=0.05)
    hist(res$pvalue[!is.na(res$padj)],150,main="Histogram of raw pvalues",xlab="raw pvalues")
    DE.BH<-which(res$padj<alpha)
    res<-data.frame(res,dispersion=mcols(cds)$dispGeneEst)
    if(filter)
        {
            write.table(res,"filter-DESeq2_CompleteListe.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
            write.table(res[DE.BH,],paste("filter-DESeq2_BH",alpha,".txt",sep=""),row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
        }
    else
        {
            write.table(res,"DESeq2_CompleteListe.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
            write.table(res[DE.BH,],paste("DESeq2_BH",alpha,".txt",sep=""),row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
        }
    dev.off()
}
