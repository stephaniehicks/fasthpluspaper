#cellbench 5cl data
#https://github.com/LuyiTian/sc_mixology/tree/master/data/csv/sc_10x.count.csv.gz
#https://github.com/LuyiTian/sc_mixology/blob/master/data/csv/sc_10x.metadata.csv.gz
#data comes from here, store locally somewhere and then change variable below
#datdir <- "C:\\Users\\Nathan\\Downloads\\" #change this to yours
datdir <- "/Users/daviddyjack/Downloads/"
dsn <- read.csv(gzfile(paste0(datdir,'sc_10x.metadata.csv.gz')),header=T)
l <- as.numeric(as.factor(dsn$cell_line))
dat <- read.csv(gzfile(paste0(datdir,'sc_10x.count.csv.gz')),header=T)
library(scran)
library(fasthplus)
library(here)
sce <- SingleCellExperiment(assays=list(logcounts=log2(dat+1)))
fit <- modelGeneVar(sce)
hvg <- getTopHVGs(fit, n=1000)
dat <- t(logcounts(sce)[hvg,])
rm(sce)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

dis_test <- c("euclidean","maximum","manhattan","canberra","binary")
hp_vals <- sapply(dis_test, function(d){
  dis <- dist(x=dat, method=d)
  h <- hpe(D=dis, L=l, p=10001)
  return(h)
})
dis_test <- sapply(dis_test, firstup)
res <- cbind(dis_test, hp_vals)
colnames(res) <- c("Dissimilarity","H+")
write.table(x=res,file=here("figures", "stab02-dissimilarities.csv"),quote=F,col.names=T,row.names=F,sep=',')

