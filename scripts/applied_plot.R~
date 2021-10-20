#cellbench 5cl data
#https://github.com/LuyiTian/sc_mixology/tree/master/data/csv/sc_10x.count.csv.gz
#https://github.com/LuyiTian/sc_mixology/blob/master/data/csv/sc_10x.metadata.csv.gz
#data comes from here, store locally somewhere and then change variable below
datdir <- "C:\\Users\\Nathan\\Downloads\\" #change this to yours
dsn <- read.csv(gzfile(paste0(datdir,'sc_10x.metadata.csv.gz')),header=T)
l <- as.numeric(as.factor(dsn$cell_line))
dat <- read.csv(gzfile(paste0(datdir,'sc_10x.count.csv.gz')),header=T)
library(scran)
library(fasthplus)
sce <- SingleCellExperiment(assays=list(logcounts=log2(dat+1)))
fit <- modelGeneVar(sce)
hvg <- getTopHVGs(fit, n=1000)
dat <- logcounts(sce)[hvg,]
rm(sce)

#library(fasthplus)
#hl1 <- hpe(D=dl1,L=l,p=101)
#hl2 <- hpe(D=dl2,L=l,p=101)

#first row
#histograms for Dw/Db for L1 and L2
#second row
#histograms for two hclust methods?

plotlocs <- list(
  c(0.08,0.51,0.56,0.99), #topleft
  c(0.55,0.98,0.56,0.99), #topright
  c(0.08,0.51,0.08,0.51), #bottomleft 
  c(0.55,0.98,0.08,0.51) #bottomrigt
)

cl1 <- '#8338EC' #Purple
cl2 <- '#06D6A0' #Teal
clw <- '#0000ff64' #blue (transparent)
clb <- '#ff000064' #red (transparent)

ylims_hist <- c(0,56000)
yticks_hist <- c(0,10000,50000)#seq(ylims_hist[1],ylims_hist[2],length.out=3) #seq(0,1,length.out=5)
ytclbs_hist <-  c('0',expression('1e'^4),expression('5e'^4))
#formatC(yticks_hist*100,digits=1,format='f')
#ylims_hist <- c(0,0.076)
#yticks_hist <- c(0,10^4,10^5,10^6)
#ytclbs_hist <- c('0',expression('10'^4),expression('10'^5),expression('10'^6))
#ylims_hist <- c(0,103500)


pdf('application_plot.pdf',width=6,height=6)
  plot.new()

  par(new = "TRUE",plt = plotlocs[[1]],las = 1,cex.axis = 1)
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  ib <- which(!ind)
  dl1 <- dist(t(dat),method='manhattan')
  dis <- as.matrix(dl1)
  dis <- dis[upper.tri(dis)]
  bins <- seq(min(dis),max(dis),length.out=20)
  xticks_hist <- c(ceiling(min(dis)*1.01),floor(max(dis*0.99)))
  hist(x=dis[iw],breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',
    col=clw, freq=T,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',ylim=ylims_hist)
  hist(x=dis[ib],breaks=bins,add=T,border='red',col=clb,freq=T,ylim=ylims_hist)
  mtext(side=2,text='Frequency',line=1.5,las=3)
  mtext(side=1,text='L1',line=0.8,las=1)
  axis(side=1,at=xticks_hist,labels=c('min','max'),las=1,mgp=c(3, .2, 0),line=0.1,cex.axis=0.8)
  axis(side=2,at=yticks_hist,labels=ytclbs_hist,las=2,mgp=c(3, .5, 0),cex.axis=0.8)
  legend('topright',bty='n',legend= c(expression('D'[W]), expression('D'[B])),
   pch=c(22,22),col= c('blue','red'),pt.cex=1.5,pt.bg=c(clw,clb) )
  #den <- as.numeric(length(dis[iw]))*as.numeric(length(dis[ib]))
  #ht <- sapply(dis[iw], function(x) x > dis[ib])/ den
  hp <- hpe(D=dl1,L=l,p=251)
  hpl <- formatC(hp$h,format='f',digits=4)
  text(x=0.88*max(xticks_hist),y=yticks_hist[2]*3.0,labels= bquote("H+" ~ "=" ~ .(hpl)),cex=0.8)

  par(new = "TRUE",plt = plotlocs[[2]],las = 1,cex.axis = 1)
  dl2 <- dist(t(dat),method='euclidean')
  dis <- as.matrix(dl2) 
  dis <- dis[upper.tri(dis)]
  bins <- seq(min(dis),max(dis),length.out=20)
  xticks_hist <- c(ceiling(min(dis)*1.01),floor(max(dis*0.99)))
  hist(x=dis[iw],breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',
    col=clw, freq=T,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',ylim=ylims_hist)
  hist(x=dis[ib],breaks=bins,add=T,border='red',col=clb,freq=T,ylim=ylims_hist)
  mtext(side=1,text='L2',line=0.8,las=1)
  axis(side=1,at=xticks_hist,labels=c('min','max'),las=1,mgp=c(3, .2, 0),line=0.1,cex.axis=0.8)
  hp <- hpe(D=dl2,L=l,p=251)
  hpl <- formatC(hp$h,format='f',digits=4)
  text(x=0.88*max(xticks_hist),y=yticks_hist[2]*3.0,labels= bquote("H+" ~ "=" ~ .(hpl)),cex=0.8)


  par(new = "TRUE",plt = plotlocs[[3]],las = 1,cex.axis = 1)
  hcc <- cutree(hclust(dl2, method = "complete"),k=3)
  l <- hcc
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  ib <- which(!ind)
  hist(x=dis[iw],breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',
    col=clw, freq=T,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',ylim=ylims_hist)
  hist(x=dis[ib],breaks=bins,add=T,border='red',col=clb,freq=T,ylim=ylims_hist)
  mtext(side=2,text='Frequency',line=1.5,las=3)
  mtext(side=1,text='HC complete',line=0.8,las=1)
  axis(side=1,at=xticks_hist,labels=c('min','max'),las=1,mgp=c(3, .2, 0),line=0.1,cex.axis=0.8)
  axis(side=2,at=yticks_hist,labels=ytclbs_hist,las=2,mgp=c(3, .5, 0),cex.axis=0.8)    
  hp <- hpe(D=dl2,L=l,p=251)
  hpl <- formatC(hp$h,format='f',digits=4)
  text(x=0.88*max(xticks_hist),y=yticks_hist[2]*3.0,labels= bquote("H+" ~ "=" ~ .(hpl)),cex=0.8)
  alphas1 <- as.numeric(c(length(iw),length(ib)) / length(ind))


  par(new = "TRUE",plt = plotlocs[[4]],las = 1,cex.axis = 1)
  hca <- cutree(hclust(dl2, method = "average"),k=3)
  l <- hca
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  ib <- which(!ind)
  hist(x=dis[iw],breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',
    col=clw, freq=T,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',ylim=ylims_hist)
  hist(x=dis[ib],breaks=bins,add=T,border='red',col=clb,freq=T,ylim=ylims_hist)
  mtext(side=1,text='HC average',line=0.8,las=1)
  axis(side=1,at=xticks_hist,labels=c('min','max'),las=1,mgp=c(3, .2, 0),line=0.1,cex.axis=0.8)
  hp <- hpe(D=dl2,L=l,p=251)
  hpl <- formatC(hp$h,format='f',digits=4)
  text(x=0.88*max(xticks_hist),y=yticks_hist[2]*3.0,labels= bquote("H+" ~ "=" ~ .(hpl)),cex=0.8)
  alphas2 <- as.numeric(c(length(iw),length(ib)) / length(ind))
  ascale <- prod(alphas2)/prod(alphas1)
  hpp <- hp$h * ascale
  hpl <- formatC(hpp,format='f',digits=4)
  text(x=0.88*max(xticks_hist),y=yticks_hist[2]*2.7,labels= bquote("H'+" ~ "=" ~ .(hpl)),cex=0.8)

dev.off()
