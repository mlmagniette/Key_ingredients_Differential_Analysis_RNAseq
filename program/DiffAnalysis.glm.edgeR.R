library("edgeR")
DiffAnalysis.glm.edgeR<-function(countTable,filter=TRUE,alpha=0.05,repbio=2)
{
  conds<-c("Bud","Leaf")
  if(filter)  
      pdf("graph-glm_edgeR-filtered.pdf")
  else
      pdf("graph-glm_edgeR.pdf")
  tmp<-which(apply(countTable,1,sum)==0)
  if(length(tmp)!=0)
      countTable<-countTable[-tmp,]
  
  if(!filter)
      cat("Number of genes under study", nrow(countTable),"\n")

## identify treatment groups
   group<-factor(rep(conds,each=repbio))
   batch<-factor(c("1","2","1","2"))

## create data structures

  cds <- DGEList( countTable , group = paste(batch,group))
 
## filter uninformative genes
  if(filter)
  {
    ridx <- rowSums(cpm(cds) > 1) >= 2 
    cds <- cds[ridx,]
    cat("Number of genes after filtering:", nrow(cds$counts),"\n")
    cds$samples$lib.size<-colSums(cds$counts)
  }
## normalisation
  cds <- calcNormFactors( cds ) 
  design <- model.matrix( ~ batch+group)
  cds <- estimateGLMCommonDisp( cds, design)
  cds<-estimateGLMTrendedDisp(cds,design)
  cds <- estimateGLMTagwiseDisp(cds,design)
  ##
  plotBCV(cds)
  fit<-glmFit(cds,design)
  fit<-glmLRT(fit) 
   N <- cds$samples$lib.size
  f <- cds$samples$norm.factors
  TMM <- N*f / mean(N*f)
  hist(fit$table$PValue,150,main="glm edgeR tagwise dispersion",xlab="raw pvalue")
  baseMean <- rowMeans(scale(cds$counts, center = FALSE, scale = TMM))
  baseMeanA <- rowMeans(scale(cds$counts[, which(group == unique(conds)[1])], 
                              center = FALSE, scale = TMM[which(group == unique(conds)[1])]))
  baseMeanB <- rowMeans(scale(cds$counts[, which(group == unique(conds)[2])], 
                              center = FALSE, scale = TMM[which(group == unique(conds)[2])]))
  foldChange <- baseMeanB/baseMeanA
  log2FoldChange <- log2(foldChange)
  
  norm<-data.frame(baseMean,baseMeanA,baseMeanB,log2FoldChange)
  
  final<-data.frame(norm,fit$coefficients,dispersion=fit$dispersion,pvalue=fit$table$PValue)
  rownames(final)=rownames(norm)
  final<-cbind(final,padj=p.adjust(final$pvalue,"BH"))
  DE.BH<-which(final$padj<alpha)
  nom<-"glm_edgeR_CompleteListe.txt"
  if(filter)
      nom<-paste("filter-",nom,sep="")
  write.table(final,nom,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
  nom<-paste("glm_edgeR_BH-",alpha,".txt",sep="")
  if (filter)
      nom<-paste("filter-",nom,sep="")
  write.table(final[DE.BH,],nom,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
  dev.off()
 
}
