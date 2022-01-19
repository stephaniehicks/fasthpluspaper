#performance plots for different calculations in H+ and its estimation
set.seed(1234)
library(RColorBrewer)
library(fasthplus)

set.seed(1234)
n_vec <- c(100,500,1000,3000)#seq(1000,5000,by=1000))
m <- 500

res <- sapply(n_vec, function(n) {
  print(n)
  dat <- t(sapply(1:n, function(i) rnorm(n=m,mean=0,sd=1)))
  lab <- c(rep(0,n/2),rep(1,n/2))
  #time distance calculation
  ptm <- proc.time()
  dis <- dist(dat)
  #dis <- as.matrix(dis)
  #dis <- dis[upper.tri(dis)]
  t_dis <- unname(proc.time() - ptm)[3]
  #time hpe calculation
  ptm <- proc.time()
  hpe <- hpe(D=dis,L=lab,p=1001,alg='grid_search')
  t_hpe <- unname(proc.time() - ptm)[3]
  #time adjacency calculation
  ptm <- proc.time()
  ind <- sapply(lab, function(x) x==lab)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  ib <- which(!ind)
  t_adj <- unname(proc.time() - ptm)[3]
  #additional distance formatting
  ptm <- proc.time()
  dis <- as.matrix(dis)
  dis <- dis[upper.tri(dis)]
  t_di2 <- unname(proc.time() - ptm)[3]
  t_dis <- t_dis + t_di2 
  #time s+ calculation
  if(n < 1000){
  ptm <- proc.time()
  sp <- sum(sapply(dis[iw], function(x) sum(x>dis[ib])))
  t_sum <- unname(proc.time() - ptm)[3]
  } else {
    t_sum <- NA
  }
  #time hpb calculation (fixed n+r)
  ptm <- proc.time()
  #hpb <- mean(replicate(30,calcf(dat,lab,100),T))
  h <- hpb(D=dat,L=lab,t=.05*n,r=30)
  t_hpb <- unname(proc.time() - ptm)[3]
  c(t_dis,t_adj,t_sum,t_hpe,t_hpb)
})

rownames(res) <- c("dis","adj","sum","hpe","hpb")

colref <- c(
  "dis"="#648FFF",
  "adj"="#FFB000",
  "sum"="#DC267F",
  "hpe"="#FE6100",
  "hpb"="#785EF0"
)
filref <- paste0(colref, '51')
names(filref) <- names(colref)
titref <- c(
  "dis"="Dissimilarity",
  "adj"="Adjacency",
  "sum"="s (Dw-Db comparisons)",
  "hpe"="hpe (pre-calculated D)",
  "hpb"="hpb (all components)"
)

ylim <- c(-1,1.05*max(res,na.rm=T))
yaxl <- pretty(c(0,max(res,na.rm=T)),n=3)
yaxt <- formatC(yaxl,digits=0,format='f')

ym2 <- max(res[-which(rownames(res)=='sum'),])*1.05 #mean(c(max(res[-which(rownames(res)=='sum'),]),ylim[2]))
yl2m <- c(-1,ym2)
ya2l <- pretty(c(0,ym2),n=3)
ya2t <- formatC(ya2l, digits=0,format='f')

con <- ya2l[3]
c2n <- (con/ym2)*ylim[2]

xlim <- c(0.99*min(n_vec),1.05*max(n_vec))
xaxl <- n_vec
xaxt <- as.character(n_vec)

plotlocs <- list(
  c(0.06,0.36,0.16,0.91),
  c(0.38,0.68,0.16,0.91),
  c(0.70,1.00,0.16,0.91)
)

#just fill the rest with pink for s in brute force H+
#res[is.na(res)] <- ylim[2]
#make sure the max is right for 2nd/3rd plots
#add tickmark to the right side of first plot to connect re-scaling

pdf("04-performance_plot.pdf",width=8,height=3)
  plot.new()

  use <- c('sum','dis','adj')
  par(new = "TRUE",plt = plotlocs[[1]],las = 1,cex.axis = 1, bty = 'n')
  plot(x=1,y=1,type='n',xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlim,ylim=ylim)
  #awkward but we want the lines to come after the shading, sometimes beauty is pain
  for(u in use){
    if(!any(is.na(res[u,]))){
      polygon(x=c(min(n_vec),n_vec,max(n_vec)), y=c(0,res[u,],0), border=NA, col=filref[u])
    }else{
      tx <- !is.na(res[u,])
      polygon(x=c(min(n_vec),n_vec[tx],max(n_vec[tx])),  y=c(0,res[u,tx],0), border=NA, col=filref[u])
      polygon(x=c(rep(max(n_vec[tx]),2),n_vec[!tx],max(n_vec)),y=c(0,rep(ylim[2],length(which(!tx))+1),0),border=NA,col=filref[u])
    }
  }
  for(u in use){
    lines(x=n_vec, y=res[u,],type='l',lty=2,col='black',lwd=0.8)
  }
  mtext(side=2,text='computation time (sec.)',line=1.0,las=3,cex=1.2)
  mtext(side=1,text='number of observations',line=1.3,las=1,cex=1.2)
  axis(side=1,labels=xaxt,at=xaxl,cex.axis=0.5,las=1,mgp=c(1.0, .1, 0))
  axis(side=2,labels=yaxt,at=yaxl,cex.axis=0.5,las=3,mgp=c(3.0, .3, 0))
  mtext(text="A",side=3,line=0.0,at=0.2*n_vec[1],cex=1.5)
  mtext(text="H+ (full)", side=3, line=0.0,cex=1.2)
  legend(x='top',legend=titref[use],pch=15,col=filref[use],cex=0.7,bg='white') 
  #axis(side=4,at=con,cex.axis=0.5,tick = TRUE, labels = FALSE,cex.axis=0.5,las=3,mgp=c(3.0, .3, 0))

  #draw two lines to connect axes between 1st and 2nd/3rd plots
  par(xpd=T)
  segments(x0=max(n_vec),y0=con,x1=mean(c(rep(xlim[2],3),max(n_vec))),y1=c2n,lty=1,col='black',lwd=1)
  #axis(side=4,at=con,cex.axis=0.5,tick = TRUE, labels = FALSE,cex.axis=0.5,las=3,mgp=c(3.0, .3, 0))
  par(xpd=F)

  use <- c('dis','hpe')
  par(new = "TRUE",plt = plotlocs[[2]],las = 1,cex.axis = 1)
  plot(x=1,y=1,type='n',xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlim,ylim=yl2m)
  for(u in use){
   polygon(x=c(min(n_vec),n_vec,max(n_vec)), y=c(0,res[u,],0), border=NA, col=filref[u])
  }
  for(u in use){
   lines(x=n_vec, y=res[u,],type='l',lty=2,col='black',lwd=0.8)
  }
  axis(side=1,labels=xaxt,at=xaxl,cex.axis=0.5,las=1,mgp=c(1.0, .1, 0))
  axis(side=2,labels=ya2t,at=ya2l,cex.axis=0.5,las=3,mgp=c(3.0, .3, 0))
  mtext(text="B",side=3,line=0.0,at=0.2*n_vec[1],cex=1.5)
  mtext(text="hpe (gridsearch)", side=3, line=0.0,cex=1.2)
  legend(x='top',legend=titref[use],pch=15,col=filref[use],cex=0.7,bg='white')
  mtext(side=1,text='number of observations',line=1.3,las=1,cex=1.2)
  #axis(side=4,at=c2n,cex.axis=0.5,tick = TRUE, labels = FALSE,col.ticks='grey60')

  use <- 'hpb'
  par(new = "TRUE",plt = plotlocs[[3]],las = 1,cex.axis = 1)
  plot(x=1,y=1,type='n',xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',xlim=xlim,ylim=yl2m)
  for(u in use){
   polygon(x=c(min(n_vec),n_vec,max(n_vec)), y=c(0,res[u,],0), border=NA, col=filref[u])
  }
  for(u in use){
   lines(x=n_vec, y=res[u,],type='l',lty=2,col='black',lwd=0.8)
  }
  axis(side=1,labels=xaxt,at=xaxl,cex.axis=0.5,las=1,mgp=c(1.0, .1, 0))
  mtext(text="C",side=3,line=0.0,at=0.2*n_vec[1],cex=1.5)
  mtext(text="hpb (bootstrap)", side=3, line=0.0,cex=1.2)
  legend(x='top',legend=titref[use],pch=15,col=filref[use],cex=0.7,bg='white')
  mtext(side=1,text='number of observations',line=1.3,las=1,cex=1.2)

dev.off()


