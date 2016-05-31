library("edgeR")
DiffAnalysis.edgeR<-function(countTable,filter=TRUE,alpha=0.05,repbio=2)
{
    conds<-c("Bud","Leaf")
    if(filter)
        pdf("graph-edgeR_filtered.pdf")
    else
        pdf("graph-edgeR.pdf")
    
    tmp<-which(apply(countTable,1,sum)==0)
    if(length(tmp)!=0)
        countTable<-countTable[-tmp,]
    
    if (!filter)
        cat("Number of genes under study", nrow(countTable),"\n")
    
## identify treatment groups
    group<-factor(rep(conds,each=repbio))
    batch<-factor(c("1","2","1","2"))
  
## create data structures
    cds<- DGEList( countTable,group=group)
  
## filter uniformative genes
    if(filter)
        {
            ridx <- rowSums(cpm(cds) > 1) >= 2 
            cds <- cds[ridx,]
            cat("Number of genes after filtering:", nrow(cds$counts),"\n")
            cds$samples$lib.size<-colSums(cds$counts)
        }
    
    cds <- calcNormFactors(cds)
    
  ## comparison between groups
  ## common dispersion 
    design <- model.matrix( ~ group )
    cds <- estimateCommonDisp( cds, design)
  ## tagwise dispersion 
    cds <- estimateTagwiseDisp(cds)
    de.tgw <- exactTest(cds)
    res<- topTags(de.tgw,n=nrow(de.tgw$table),sort.by="none")$table
    res<-data.frame(res,dispersion=cds$tagwise.dispersion)
    hist(res$PValue,150,main="Exact test tagwise dispersion",xlab="raw pvalue")
    cat("Exact test of edgeR with tagwise dispersion \n")
    DE.BH<-which(res$FDR<alpha)
    nom<-"edgeR_CompleteListe.txt"
    if (filter)
        nom<-paste("filter-",nom,sep="")
    names(res)[grep("PValue",names(res))]<-"pvalue"
    names(res)[grep("FDR",names(res))]<-"padj"
    write.table(res,nom,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    nom<-paste("edgeR_BH-",alpha,".txt",sep="")
    if (filter)
        nom<-paste("filter-",nom,sep="")
    write.table(res[DE.BH,],nom,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    dev.off()
}
