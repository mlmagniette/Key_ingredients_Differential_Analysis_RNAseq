#############################################################################
## FDR evaluation written as a function to use it with a varying set of NDE
## function parameter :
## id.NDE := a NDE set
## alpha := FDR threshold to declare a gene differentially expressed
## RData file to store the result as a matrix 
FDR<-function(id.NDE,alpha,outputRData)
{
    file<- system("ls ../synthetic_datasets/synthetic*/*_CompleteListe.txt",T)
    nn<-length(file)
    FDR<-numeric()
    
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
            id.declared.DE<-rownames(data)[data$padj<=alpha]
            FDR[i]<-length(intersect(id.declared.DE,id.NDE.tmp))/length(id.declared.DE)
        }
    FDR.res=matrix(FDR,nrow=length(file)/length(methods),byrow=T,ncol=length(methods))
    colnames(FDR.res)<-methods
    save(FDR.res,file=outputRData)   
}
######
