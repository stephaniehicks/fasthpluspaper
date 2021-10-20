#library(pbapply)
library(plyr)
map_fxn <- function(x,out_min,out_max,in_min,in_max){
  (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}


set.seed(1234)
n <- 1000
m <- 500
Nd <- (n*(n-1))/2
Nz <- (Nd*(Nd-1))/2

pct_samps <- c(.000005,0.00005,.0005,.005,.05,.125) 

#one dataset, different amounts of dw and db to sample
#2x2, g+ h+ columns, value time rows

#note: time is going to represent time for s+ calculation because the devision is trivial by comparison.
#one time thing is the same for both with s+ though? 

dat <- cbind(
  sapply(1:(n/2), function(i) rnorm(n=m,mean=0.1,sd=1)),
  sapply(1:(n/2), function(i) rnorm(n=m,mean=-0.1,sd=1))
)
labs <- c(rep(0,n/2),rep(1,n/2))
ind <- sapply(labs, function(x) sapply(labs, function(y) x==y))
ind <- ind[upper.tri(ind)]
iw <- which(ind)
ib <- which(!ind) 
dis <- as.matrix(dist(t(dat))) 
dis <- dis[upper.tri(dis)]
ptm <- proc.time()
dw <- sort(dis[iw])
db <- sort(dis[ib]) #think this needs to encorporated into the timing
tsort <- proc.time() - ptm


#calculate true s+
ptm <- proc.time()
sp <- sum(sapply(dw, function(x) sum(x>db)))
tfull <- proc.time() - ptm

gp <- sp/Nz
hp <- sp/(as.numeric(length(dw))*as.numeric(length(db)))

spestm <- lapply(pct_samps, function(x) {
  n_samp <- round(x*length(ind))
  Nz_samp <- ((n_samp+n_samp)*(n_samp+n_samp-1))/2
  ow <- round(seq(1,length(iw),length.out=n_samp))
  ob <- round(seq(1,length(ib),length.out=n_samp))
  ptm <- proc.time()
  sp <- sum(sapply(dw[ow], function(x) sum(x > db[ob])))
  t <- proc.time() - ptm
  hpe <- sp / (as.numeric(length(ow))*as.numeric(length(ob)))
  gpe <- (2*sp) / ((n_samp+n_samp)*(n_samp+n_samp-1))
  return(list(gpe,hpe,t[3]+tsort[3]))
})


#make paired heatmaps (colorbar on top?) for G+ and H+
#plotting parameters

#first row, G+ and H+
#second row, time 
plotlocs <- rbind(
  c(0.07,0.51,0.54,0.96), #topleft g+
  c(0.54,0.99,0.54,0.96), #topright h+
  c(0.25,0.75,0.05,0.47), #bottom time
  c(0.75,1.00,0.05,0.47) #bottomright legend
)


xlab <- '% Dists. Sampled'
xlims <- c(-.05,1.05)
xticks <- c(0.0,0.15,0.30,0.50,0.75,1.0) #pct_samps*200
xtclbs <- c('.0001', '0.01', '0.1','1','10','25')

ylims_p12 <- c(-.02,0.52)
yticks_p12 <- seq(0.0,.5,by=0.1)
ytclbs_p12 <- formatC(yticks_p12,digits=1,format='f')
ylab_p1 <- 'G+'
ylab_p2 <- 'H+'


et <- sapply(spestm, function(x) x[[3]])
tmin <- round_any(et[5],5,f=floor)
tmax <- round_any(tfull[3],10,f=ceiling)
et[5:6] <- map_fxn(et[5:6],0.8,2,tmin,tmax)
tfull_scld <-  map_fxn(tfull[3],0.8,2,tmin,tmax)
prtbig <- round_any(seq(tmin,tmax,length.out=3),5)
ylims_p3 <- c(-0.1,2.1) #c(tmin-2.1,tmax+2.1)
yticks_p3 <- c(0.0,0.5,seq(from=0.8,to=2.0,length.out=length(prtbig)))
ytclbs_p3 <- c('0.0','0.5',as.character(prtbig))
ylab_p3 <- 'Time (s)'


pdf("C:/Users/Nathan/Documents/Hplus/performance_plot.pdf",width=8,height=8)
  plot.new()

  par(new = "TRUE",plt = plotlocs[1,],las = 1, cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims,ylim=ylims_p12,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks, y=sapply(spestm, function(x) x[[1]]),pch=1,col='blue')
  abline(h=gp,lty=2,col="black",lwd=1.1)
  axis(side=2,labels=ytclbs_p12,at=yticks_p12,cex=1.0,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xtclbs[1:3],at=xticks[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs[4:6],at=xticks[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab,line=1.0,font=2,las=1)
  mtext(side=2,text=ylab_p1,line=1.5,font=2,las=3)

  par(new = "TRUE",plt = plotlocs[2,],las = 1, cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims,ylim=ylims_p12,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks, y=sapply(spestm, function(x) x[[2]]),pch=1,col='blue')
  abline(h=hp,lty=2,col="black",lwd=1.1)
  axis(side=1,labels=xtclbs[1:3],at=xticks[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs[4:6],at=xticks[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab,line=1.0,font=2,las=1)
  mtext(side=2,text=ylab_p2,line=0.05,font=2,las=3)

  par(new = "TRUE",plt = plotlocs[3,],las = 1, cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=xlims,ylim=ylims_p3,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  points(x=xticks, y=et,pch=1,col='blue')
  abline(h=tfull_scld,lty=2,col="black",lwd=1.1)
  axis(side=2,labels=ytclbs_p3[1:2],at=yticks_p3[1:2],cex=1.0,las=2,mgp=c(3, .5, 0))
  axis(side=2,labels=ytclbs_p3[3:5],at=yticks_p3[3:5],cex.axis=1.0,las=2,mgp=c(3, .5, 0))
  axis(side=1,labels=xtclbs[1:3],at=xticks[1:3],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  axis(side=1,labels=xtclbs[4:6],at=xticks[4:6],cex.axis=0.7,las=1,mgp=c(3, .3, 0),line=0)
  mtext(side=1,text=xlab,line=1.0,font=2,las=1)
  mtext(side=2,text=ylab_p3,line=2.0,font=2,las=3)

  par(new = "TRUE",plt = plotlocs[4,],las = 1, cex.axis = 1)
  plot(x=0,y=0,type='n',xlim=c(0,1),ylim=c(0,1),xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',frame.plot=F)
  legend('center',legend=c('Full','Estimated'),col=c('black','blue'),pch=c(NA,1),lty=c(2,NA),pt.cex=c(NA,1.5),lwd=c(1.5,NA),bty=T)

dev.off()

