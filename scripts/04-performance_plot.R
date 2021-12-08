#heatmaps for 2 difficulties. include PCA plots to demonstrate the 
set.seed(1234)
library(RColorBrewer)
library(fasthplus)
#cl1 <- '#8338EC' #Purple
#cl2 <- '#06D6A0' #Teal

set.seed(1234)
n_vec <- c(500,seq(1000,10000,by=1000))
m <- 500

calcf <- function(d,l,r){
  stmp <- table(l)/length(l)
  stmp <- as.vector(sapply(names(stmp), function(x) sample(which(l==x),round(r/length(stmp)))))
  labtmp <- l[stmp]
  distmp <- as.matrix(dist(d[stmp,]))
  distmp <- distmp[upper.tri(distmp)]
  indtmp <- sapply(labtmp, function(x) x==labtmp)
  indtmp <- indtmp[upper.tri(indtmp)]
  iwtmp <- which(indtmp)
  ibtmp <- which(!indtmp)
  dwtmp <- distmp[iwtmp]
  dbtmp <- distmp[ibtmp]
  sptmp <- sum(sapply(dwtmp, function(x) sum(x>dbtmp)))
  hptmp <- sptmp / (as.numeric(length(dwtmp))*as.numeric(length(dbtmp)))
 return(hptmp) 
}


res <- sapply(n_vec, function(n) {
  dat <- t(sapply(1:n, function(i) rnorm(n=m,mean=0,sd=1)))
  lab <- c(rep(0,n/2),rep(1,n/2))
  #time distance calculation
  ptm <- proc.time()
  dis <- dist(dat)
  t_dis <- unname(proc.time() - ptm)[3]
  #time hpe calculation
  ptm <- proc.time()
  hpe <- hpe(D=dis,L=lab,p=251)
  t_hpe <- unname(proc.time() - ptm)[3]
  #extra dis calcs (untimed)
  dis <- as.matrix(dis)
  dis <- dis[upper.tri(dis)]
  #time adjacency calculation
  ptm <- proc.time()
  ind <- sapply(lab, function(x) x==lab)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  ib <- which(!ind)
  t_adj <- unname(proc.time() - ptm)[3]
  #time s+ calculation
  if(n < 2000){
    ptm <- proc.time()
    sp <- sum(sapply(dis[iw], function(x) sum(x>dis[ib])))
    t_spl <- unname(proc.time() - ptm)[3]
  } else {
    t_spl <- NA
  }
  #time hpb calculation (fixed n+r)
  ptm <- proc.time()
  hpb <- mean(replicate(30,calcf(dat,lab,100),T)) 
  t_hpb <- unname(proc.time() - ptm)[3]
  c(t_dis,t_adj,t_spl,t_hpe,t_hpb)
})

colref <- c('#000000','#E69F00','#009E73','#0072B2','#CC79A7') #palette.colors(palette = "Okabe-Ito")[2:6]
pchref <- 21:25

ylim <- c(-50,1.05*max(res,na.rm=T))
yaxl <- pretty(c(0,max(res,na.rm=T)),n=3)
yaxt <- formatC(yaxl,digits=0,format='f')

xlim <- c(0,10500)
xaxl <- n_vec
xaxt <- as.character(n_vec)

pdf("04-performance_plot.pdf",width=6,height=4)
  par(mfrow=c(1,1),mar=c(3.1,3.1,0.6,0.6))
  plot(x=1,y=1,type='n',xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim)
  for(i in 1:nrow(res)){
    lines(x=n_vec,y=res[i,],col=colref[i],lwd=2.0,type='b',pch=pchref[i],cex=1.5,lty=i)
  }
  mtext(side=2,text='computation time (sec.)',line=1.5,las=3,cex=1.3)
  mtext(side=1,text='number of observations',line=1.5,las=1,cex=1.3)
  axis(side=1,labels=xaxt,at=xaxl,cex.axis=0.5,las=1,mgp=c(1.0, .1, 0))
  axis(side=2,labels=yaxt,at=yaxl,cex.axis=0.7,las=3,mgp=c(3.0, .4, 0))
  legend(x='topright',legend=c('Dissimilarity','Adjacency','S+ (Dw-Db comparisons)','HPE (pre-calculated D)','HPB (r=30,s=100)'),
    lty=1:5,pch=pchref,col=colref,lwd=2.0,bty='n')

dev.off()


