syntheticData<-function (H0number) 
{
    dat <- read.table("../../real_datasets/Arabidopsis_dataset.txt",h=T)
    removeID<-read.table("../../real_datasets/removeAGI.txt")
    dat<-dat[!is.element(dat[,1],removeID[,1]),]
    fix <- grep(".DE", dat[, 1],fixed=TRUE)
    modify <- dat[-fix, 1]
    synth <- dat[, 1:5]
    if (H0number < 1 & H0number >= 0) {
        subset <- sample(modify, H0number * length(modify))
    }
    if (H0number >= 1) {
        subset <- sample(modify, H0number)
    }
    tochange <- is.element(synth[, 1], subset)
    synth[tochange, -1] <- dat[tochange, c(4, 6, 5, 7)]
    synth[, 1] <- as.character(synth[, 1])
    synth[tochange, 1] <- paste(synth[tochange, 1], ".Hnull", 
        sep = "")
    rownames(synth) <- synth[, 1]
    return(synth[, -1])
}
