##################################################################
## determination of the truly NDE for each synthetic dataset
## Differential analysis methods used are all the methods without data filtering
##################################################################
## for each  method and for each  synthetic dataset
## 1) correction of the raw pvalues
## 2) definition of the NDE genes as genes having a correted pvalue greater than alpha. All Hnull genes are included
## 3) remove DE genes
## 4) NDE sets compose a list which is saved in the result folder
##################################################################
## for each synthetic dataset
## NDE.inter is created by the genes found no DE by the five methods
## NDE.union is created by the genes found no DE by at least one method
## output : two lists  NDE.inter and NDE.union saved in the result folder
##################################################################


method.ref<-c("DESeq2","DESeq","edgeR","glm_edgeR","limma-voom")

for (ii in 1:length(method.ref))
{
    ref<-method.ref[ii]
    fileRData<-paste("../results/NDE-",ref,".RData",sep="")
    file.ref<- system(paste("ls ../synthetic_datasets/synthetic*/",ref,"_CompleteListe.txt",sep=""),T)
    
    id.NDE<-list()
    for (i in 1:length(file.ref))
        {
            print(i)
            data<-read.table(file.ref[i],h=T)
            pval<-as.matrix(data$pvalue,ncol=1)  
            index.Hnull<-grep("Hnull",rownames(data))
            nstar<-length(index.Hnull)
            pval.correct<-apply(pval,1,function(x) sum(pval[index.Hnull,1]<=x)/nstar)
            tmp<-union(rownames(data)[which(pval.correct>alpha)],rownames(data)[index.Hnull])
            if (length(grep("DE",tmp,fixed=TRUE))!=0)
                id.NDE[[i]]<-tmp[-grep("DE",tmp,fixed=TRUE)]
            else
                id.NDE[[i]]<-tmp
        }
            save(id.NDE,file=fileRData)
}
   
obj<-system("ls ../results/NDE*.RData",T)
NDE<-list()
for (i in 1:length(obj))
{
    load(obj[i])
    NDE[[i]]<-id.NDE

}

r1<-list()
r2<-list()
k=1
for (i in 1:4)
{
  for (j in (i+1):5)
    {
     r1[[k]]<-sapply(1:90,function(x) intersect(NDE[[i]][[x]],NDE[[j]][[x]]))
     r2[[k]]<-sapply(1:90,function(x) union(NDE[[i]][[x]],NDE[[j]][[x]]))
     k<-k+1
    }
}

NDE.inter<-sapply(1:90,function(x)  intersect(intersect(intersect(intersect(r1[[1]][[x]],r1[[2]][[x]]),r1[[3]][[x]]),r1[[4]][[x]]),r1[[5]][[x]]))

NDE.union<-sapply(1:90,function(x)  union(union(union(union(r1[[1]][[x]],r1[[2]][[x]]),r1[[3]][[x]]),r1[[4]][[x]]),r1[[5]][[x]]))


save(NDE.union,file="../results/NDE.union.RData")
save(NDE.inter,file="../results/NDE.inter.RData")

