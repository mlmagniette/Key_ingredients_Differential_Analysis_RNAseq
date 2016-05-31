library("limma")
library("edgeR")
DiffAnalysis.limmavoom<-function(countTable,filter=TRUE,alpha=0.05,repbio=2)
{
    conds<-c("Bud","Leaf")
    if(filter)
        pdf("graph-limma-voom-filter.pdf")
    else
        pdf("graph-limma-voom.pdf")	
   
    tmp<-which(apply(countTable,1,sum)==0)
    if(length(tmp)!=0)
        countTable<-countTable[-tmp,]
    
    if (!filter)
        cat("Number of genes under study", nrow(countTable),"\n")
  
## identify treatment groups
    group<-factor(rep(conds,each=repbio))
    batch<-factor(c("1","2","1","2"))
    
## create data structures
    cds<-DGEList(counts=countTable,group=paste(batch,group))
## filter uniformative genes
#### meme filtre que edgeR
  if(filter)
  {
      ridx <- rowSums(cpm(cds) > 1) >= 2 
      cds <- cds[ridx,]
      cds$samples$lib.size<-colSums(cds$counts)
  }
    cds <- calcNormFactors(cds)
    
  ## comparison between groups
  
    design <- model.matrix( ~0+ group+batch )
    colnames(design) <- c(conds,"batch2")
    y<-voom(cds, design,plot=TRUE)	
    fit <-lmFit(y,design)
    contrast.matrix<-makeContrasts(contrast="Leaf-Bud",levels=design)
    fit2<-contrasts.fit(fit,contrast.matrix)
    fit<-eBayes(fit2)
    res<- topTable(fit,n=nrow(countTable),sort.by="none")
    hist(res$P.Value,150,main="limma-voom",xlab="raw pvalue")
    DE.BH<-which(res$adj.P.Val<alpha)
    names(res)[grep("P.Value",names(res))]<-"pvalue"
    names(res)[grep("adj.P.Val",names(res))]<-"padj"
    res<-res[,-ncol(res)]
    res<-data.frame(res,dispersion=fit$stdev.unscaled*sqrt(fit$s2.post))
    names(res)[ncol(res)]="dispersion"
    nom<-"limma-voom_CompleteListe.txt"		  
  if (filter)
    nom<-paste("filter-",nom,sep="")
    write.table(res,nom,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    nom<-paste("limma-voom_BH-",alpha,".txt",sep="")
    if (filter)
        nom<-paste("filter-",nom,sep="")
    write.table(res[DE.BH,],nom,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    dev.off()  
}

