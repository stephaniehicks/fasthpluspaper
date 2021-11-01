#layout
# % % %%
# % % %%
# %%%%%%
# %%%%%%
#top half: two 2x2s (MDS) next to a bigger 1x1 (scatter plot ARI vs H+)
#bottom half: elbow plot (k=2:10), double y axis with left is ARI and right is H+


#cellbench 5cl data
#https://github.com/LuyiTian/sc_mixology/tree/master/data/csv/sc_10x.count.csv.gz
#https://github.com/LuyiTian/sc_mixology/blob/master/data/csv/sc_10x.metadata.csv.gz
#data comes from here, store locally somewhere and then change variable below
datdir <- "C:\\Users\\Nathan\\Downloads\\" #change this to yours
dsn <- read.csv(gzfile(paste0(datdir,'sc_10x.metadata.csv.gz')),header=T)
l <- as.numeric(as.factor(dsn$cell_line))
dat <- read.csv(gzfile(paste0(datdir,'sc_10x.count.csv.gz')),header=T)
library(mclust)
library(scran)
library(fasthplus)
library(cluster)
library(purrr)
sce <- SingleCellExperiment(assays=list(logcounts=log2(dat+1)))
fit <- modelGeneVar(sce)
hvg <- getTopHVGs(fit, n=1000)
dat <- logcounts(sce)[hvg,]
rm(sce)

colref <- c('#117733','#88CCEE','#CC6677')
filref <- paste0(colref,'64')
shpref <- c(21,22,23)

shp <- sapply(l, function(i) shpref[i])

dl2 <- dist(t(dat),method='euclidean')
md1 <- cmdscale(d=dl2,k = 2)
xlim1 <- c(floor(min(md1[,1])),ceiling(max(md1[,1])))
ylim1 <- c(floor(min(md1[,2])),ceiling(max(md1[,2]))) 

hc_test <- c("ward.D", "single", "complete", "average")
names(hc_test) <- c("ward.D", "single", "complete", "UPGMA")
labs <- lapply(hc_test, function(x) cutree(hclust(dl2,method=x),k=3))
hc_cols <- lapply(labs, function(x) sapply(x, function(y) colref[y]))
hc_fils <- lapply(labs, function(x) sapply(x, function(y) filref[y]))

#hvals <- sapply(labs, function(l) hpe(D=dl2,L=l,alphas=T,p=101)$h)
hs <- sapply(labs, function(l) hpe(D=dl2,L=l,alphas=T,p=101)$h)
as <- sapply(labs, function(y) adjustedRandIndex(l,y)) 

map_fxn <- function(x,out_min,out_max,in_min,in_max){
  (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

calc_wss <- function(d,l){
  nr_tmp <- as.numeric(table(l))
  k <- length(nr_tmp)
  drs_tmp <- sapply(1:k, function(r) {
    idx_tmp <- which(l==r)
    sbt_tmp <- as.matrix(d)[idx_tmp,idx_tmp]
    sbt_tmp <- sbt_tmp[upper.tri(sbt_tmp)]^2
    dr_tmp <- sum(sbt_tmp)/( 2 * nr_tmp[r] )
  })
  w_tmp <- sum(drs_tmp)
  return(w_tmp)
}

k_vec <- 2:10
stats_k <- sapply(k_vec, function(k) {
  p <- pam(x=dl2,k,diss=T)
  h <- 1 - hpe(D=dl2,L=p$clustering,alphas=F,p=101)$h
  s <- p$silinfo$avg.width
  #mean(silhouette(x=p$clustering,dist=dl2)[,3])
  #w <- calc_wss(dl2,p$clustering)
  o <- min(p$objective)
  c(h,s,o)
})

l_bg <- mean(log(ws_bg)) #apply(ws_bg,1,function(x) mean(log(x)))

stats_plot <- t(apply(stats_k,1, function(x) map_fxn(x,0,1,min(x),max(x))))

cols_s <- c("#332288","#44AA99","#AA4499")
ltys_s <- c(1,2,3)

xlim_s <- c(1.9,10.1)
xticks_s <- k_vec
xlabs_s <- as.character(k_vec)

ylim_s <- c(-0.05,1.05)
yticks_s <- c(0.01,0.99) #c(0.00,0.20,0.40,0.60,0.80,1.0)
ylabs_s <- c('min','max') #formatC(yticks_s,format='f',digits=1)		      

xlim2 <- c(-0.1,1.1)
xtck2 <- c(0,0.25,0.50,0.75,1.0)
xlab2 <- c('0.00','0.25','0.50','0.75','1.00')
ylim2 <- c(-0.05,0.35)
ytck2 <- c(0,0.10,0.20,0.30)
ylab2 <- c('0.0','0.1','0.2','0.3')

plotlocs <- list(
  c(0.06,0.26,0.77,0.97), #1,1
  c(0.29,0.49,0.77,0.97), #1,2
  c(0.06,0.26,0.53,0.73), # 2,1
  c(0.29,0.49,0.53,0.73), # 2,2
  c(0.55,0.96,0.53,0.97), # 1:2,3:4
  c(0.06,0.96,0.06,0.46)  # 3:4,1:4
)

pdf('application_plot_v4.pdf',width=8,height=8)
  plot.new()

  par(new = "TRUE",plt = plotlocs[[1]],las = 1,cex.axis = 1)
  plot(x=md1[,1],y=md1[,2],pch=shp, col=hc_cols[[1]], bg=hc_fils[[1]],main='',xlab='',ylab='',
       xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',cex=0.7,xlim=xlim1,ylim=ylim1)
  mtext(side=3,line=0.0,text=names(hc_test)[1])
  mtext(side=2,text='MDS 2',line=0.3,las=3)

  par(new = "TRUE",plt = plotlocs[[2]],las = 1,cex.axis = 1)
  plot(x=md1[,1],y=md1[,2],pch=shp, col=hc_cols[[2]], bg=hc_fils[[2]],main='',xlab='',ylab='',
       xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',cex=0.7,xlim=xlim1,ylim=ylim1)
  u <- which(hc_cols[[2]] == colref[3])
  points(md1[u,1],y=md1[u,2],pch=shp[u], col=hc_cols[[2]][u], bg=hc_fils[[2]][u])
  mtext(side=3,line=0.0,text=names(hc_test)[2])

  par(new = "TRUE",plt = plotlocs[[3]],las = 1,cex.axis = 1)
  plot(x=md1[,1],y=md1[,2],pch=shp, col=hc_cols[[3]], bg=hc_fils[[3]],main='',xlab='',ylab='',
       xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',cex=0.7,xlim=xlim1,ylim=ylim1)
  mtext(side=3,line=0.0,text=names(hc_test)[3])
  mtext(side=2,text='MDS 2',line=0.3,las=3)
  mtext(side=1,text='MDS 1',line=0.3,las=1)

  par(new = "TRUE",plt = plotlocs[[4]],las = 1,cex.axis = 1)
  plot(x=md1[,1],y=md1[,2],pch=shp, col=hc_cols[[4]], bg=hc_fils[[4]],main='',xlab='',ylab='',
       xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',cex=0.7,xlim=xlim1,ylim=ylim1)
  mtext(side=3,line=0.0,text=names(hc_test)[4])
  mtext(side=1,text='MDS 1',line=0.3,las=1)


  par(new = "TRUE",plt = plotlocs[[5]],las = 1,cex.axis = 1)
  plot(x=as,y=hs,pch=1,main='',xlab='',ylab='',xaxs = "i",yaxs = "i",xaxt='n',yaxt='n', cex=1.5,xlim=xlim2,ylim=ylim2)
  text(x=as,y=hs,labels=names(hc_test),pos=1)
  mtext(side=2,text='H+',line=0.3,las=1)
  mtext(side=1,text='ARI',line=1.3,las=1)
  axis(side=1,labels=xlab2,at=xtck2,cex=1.0,las=1,mgp=c(3, .5, 0))
  axis(side=2,labels=ylab2,at=ytck2,cex=1.0,las=1,mgp=c(3, .5, 0))

  par(new = "TRUE",plt = plotlocs[[6]],las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',main='',xlab='',ylab='',xaxs = "i",yaxs = "i",xaxt='n',yaxt='n', cex=1.5,xlim=xlim_s,ylim=ylim_s)
  axis(side=2,labels=ylabs_s,at=yticks_s,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text="k (number clusters)",cex=1.1,line=1.1,las=1)
  mtext(side=2,text="Scaled Performance",line=0.3,las=3)
  axis(side=1,labels=xticks_s,at=xticks_s,cex=1.0,las=1,mgp=c(3, .5, 0))
  for(i in 1:2){
    points(x=k_vec,y=stats_plot[i,],type='b',pch=16,col=cols_s[i],lty=ltys_s[i],lwd=1.5,cex=1.5)
  }
  legend('right',bty='n',legend=c("1-H+","Mean Sil."), pt.cex=1.5, pch=16, pt.lwd=1.5,cex = 1, col=cols_s,lty=ltys_s,lwd=1.5)

dev.off()

ylim_s <- c(30,50)
yticks_s <- c(30,40,50) #c(0.01,0.99) #c(0.00,0.20,0.40,0.60,0.80,1.0)
ylabs_s <-  c('30','40','50') #formatC(yticks_s,format='f',digits=1)		      

pdf('wss_plot.pdf',width=6,height=4)
  plot.new()
  par(new = "TRUE",plt = c(0.1,0.95,0.15,0.95),las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',main='',xlab='',ylab='',xaxs = "i",yaxs = "i",xaxt='n',yaxt='n', cex=1.5,xlim=xlim_s,ylim=ylim_s)
  axis(side=2,labels=ylabs_s,at=yticks_s,cex=1.0,las=2,mgp=c(3, .5, 0))
  mtext(side=1,text="k (number clusters)",cex=1.1,line=1.2,las=1)
  mtext(side=2,text="Mean within-cluster distance to each medoid",line=1.5,las=3)
  axis(side=1,labels=xticks_s,at=xticks_s,cex=1.0,las=1,mgp=c(3, .5, 0))
  points(x=k_vec,y=stats_k[3,],type='b',pch=16,col=cols_s[3],lty=ltys_s[3],lwd=1.5,cex=1.5)
dev.off()
