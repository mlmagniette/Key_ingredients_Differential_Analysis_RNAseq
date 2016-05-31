TPR.with.fixed.fdr<-function(id.NDE,threshold,outputRData)
{
    file<- system("ls ../synthetic_datasets/synthetic*/*_CompleteListe.txt",T)
    nn<-length(file)
    res.TPR.with.fixed.fdr<-numeric()
    
    j=1
    id.NDE.tmp<-id.NDE[[j]]	
    
    for (i in 1:nn)
        {
            print(i)
            data<-read.table(file[i],h=T)
            if (is.element(i,seq(length(methods)+1,length(file),by=length(methods))))
            {
                j<-j+1
                id.NDE.tmp<-id.NDE[[j]]
            }
        ## order the dataset according to padj (increasing order)
            data<-data[order(data$padj),]
            index.NDE<-sort(match(id.NDE.tmp,rownames(data)))
            top.index.NDE<-index.NDE[floor(threshold*length(index.NDE))]
            res.TPR.with.fixed.fdr[i]<-length(grep(".DE",rownames(data)[1:top.index.NDE]))
        }
    res.TPR.with.fixed.fdr=matrix(res.TPR.with.fixed.fdr,nrow=length(file)/length(methods),byrow=T,ncol=length(methods))
    colnames(res.TPR.with.fixed.fdr)<-methods
    save(res.TPR.with.fixed.fdr,file=outputRData)   
}
######
