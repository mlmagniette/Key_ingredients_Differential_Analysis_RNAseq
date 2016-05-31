## AUC calculation 
## function parameters: a truly NDE set, a name of RData for the output
## output:  
auc<-function(id.NDE,outputRData)
{
    file<- system("ls ../synthetic_datasets/synthetic*/*_CompleteListe.txt",T)
    nn<-length(file)
    ## NDE set initialization
    j<-1
    id.NDE.tmp<-id.NDE[[j]]

    auc<-numeric()
    for (i in 1:nn)
        {
            print(i)
            data<-read.table(file[i],h=T)
            if (is.element(i,seq(nb.methods+1,length(file),by=nb.methods)))
            {
                j<-j+1
                id.NDE.tmp<-id.NDE[[j]]
            }
            DE=grep("DE",rownames(data),fixed=TRUE)
            NDE=na.omit(match(id.NDE.tmp,rownames(data)))
            subset<-c(NDE,DE)
            nb0<-length(NDE) 
            nb1<-length(DE)
            truth<-c(rep(0,nb0),rep(1,nb1))
            index<-order(data$padj[subset])
            truth.ranked<-truth[index]
            stack_x = cumsum(truth.ranked==0)/nb0
            stack_y = cumsum(truth.ranked==1)/nb1
            auc[i] = sum((stack_x[2:length(truth.ranked)]-stack_x[1:length(truth.ranked)-1])*stack_y[2:length(truth.ranked)])
        }
        res.auc=matrix(auc,byrow=T,ncol=length(methods))
        colnames(res.auc)=methods
    save(res.auc,file=outputRData)
    }
