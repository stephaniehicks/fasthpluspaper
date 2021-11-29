#heatmaps for 2 difficulties. include PCA plots to demonstrate the 
set.seed(1234)
library(RColorBrewer)
library(viridis)
library(leaflet)
library(fasthplus)
#cl1 <- '#8338EC' #Purple
#cl2 <- '#06D6A0' #Teal
colref <- palette.colors(palette = "Okabe-Ito")[2:8]

set.seed(1234)
#q, number of boostraps
#r, per-bootstrap sample size
q_vec <- c(1,3,5,10,20,30)
r_vec <- seq(10,100,by=10) #c(10,2,50,100)

#full H+ calculation
n <- 1000
#n1 <- n2 <- 500
m <- 500

setups <-  list(
 list(b=c(0.5,0.5),m=c(-0.05,0.05),s=0.5),
 list(b=c(0.25,0.25,0.25,0.25),m=c(-0.15,-0.1,0.1,0.15),s=0.5)
)


calcf <- function(d,l,r){
  #is this a general solution to sampling from classes in a balanced way?
  stmp <- table(l)/length(l)
  stmp <- as.vector(sapply(names(stmp), function(x) sample(which(l==x),round(r/length(stmp)))))
  #
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

dat <- lapply(setups, function(x) {
  len <- 1:length(x$b)
  z <- lapply(len, function(j) sapply(1:round(x$b[j]*n), function(i) rnorm(n=m,mean=x$m[j],sd=x$s)))
  z <- t(do.call(cbind,z))
  lab <- as.vector(sapply(len, function(j) rep(j,round(x$b[j]*n))))
  col <- sapply(lab, function(j) colref[j])
  pc <- prcomp(z)$x[,1:2]
  dis <- dist(z)
  hp <- hpe(D=dis,L=lab,p=10001)$h
  hpb <- sapply(r_vec, function(r) {
    sapply(q_vec, function(q) {
      mean(replicate(q,calcf(z,lab,r),T))
    })
  })
  fin <- list(p=pc, h=hp, hb=hpb, res=t(abs(hp-hpb)), c=col,b=x$b)
  return(fin)
})


nr <- length(r_vec)
nc <- length(q_vec)
gres <- dat[[1]]$res
hres <- dat[[2]]$res

plotlocs <- rbind(
  c(0.01,0.51,0.01,0.51),
  c(0.49,0.99,0.01,0.51),
  c(0.30,0.70,0.56,0.62),
  c(0.05,0.48,0.67,0.97), 
  c(0.53,0.98,0.67,0.97) 

#  c(0.01,0.51,0.35,0.85), #left heatmap
#  c(0.49,0.99,0.35,0.85), #right heatmap
#  c(0.30,0.70,0.90,0.97), # colorbar
#  c(0.05,0.48,0.02,0.32), #left pcaplot
#  c(0.53,0.98,0.02,0.32)  #right pcaplot
)

zlim <- c(0,max(rbind(gres,hres)))
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)

xlim <- c(-0.5,nc+1)
ylim <-c(-0.5,nr+0.4)
zleglocs <- c(0.01,0.50,0.99)
zleglabs <- formatC(seq(zlim[1],zlim[2],length.out=3),format='f',digits=2)
txtcut <- 0.05
xlabs <- as.character(q_vec)
ylabs <- as.character(r_vec)
yvals <- sapply(1:nr, function(i) c(i-0.5,i+0.5))

pdf("supp03-bootstrap_accuracy.pdf",width=8,height=10)
  plot.new()
  par(new = "TRUE",plt = plotlocs[1,],las = 1, cex.axis = 1)
  #plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(k in 1:nc){
    for(j in 1:nr){
      rect(xleft=k-0.5,xright=k+0.5,ybot=yvals[1,j],ytop=yvals[2,j],border='black',col=col_pal_fxn(gres[j,k]))
      text(x=k,y=mean(yvals[,j]),labels=formatC(gres[j,k],digits=2,format='f'),cex=0.6,col=ifelse(gres[j,k]>txtcut,'black','white'))
    }
  }
  text(x=1:nc,y=0.3,labels=xlabs,cex=1.0)
  text(x=0.2,y=1:nr,labels=ylabs,cex=1.0,srt=0)
  mtext(side=3,text='C',cex=1.3,at=0.5)
  par(xpd = TRUE)
  text(x=-0.2,y=mean(yvals),labels= "per-bootstrap sample size",srt=90,cex=1.3)
  text(x=7,y=-0.2,labels='number of bootstraps',cex=1.3)
  par(xpd=FALSE)

  par(new = "TRUE",plt = plotlocs[2,],las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(k in 1:nc){
    for(j in 1:nr){
      rect(xleft=k-0.5,xright=k+0.5,ybot=yvals[1,j],ytop=yvals[2,j],border='black',col=col_pal_fxn(hres[j,k]))
      text(x=k,y=mean(yvals[,j]),labels=formatC(hres[j,k],digits=2,format='f'),cex=0.6,col=ifelse(hres[j,k]>txtcut,'black','white'))
    }
  }
  text(x=1:nc,y=0.3,labels=xlabs,cex=1.0)
  mtext(side=3,text='D',cex=1.3,at=0.5)

  par(new = "TRUE",plt = plotlocs[3,],las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))
  mtext(side=3,text="|Hb - H+|", cex=1.3)

  par(new = "TRUE",plt = plotlocs[4,],las = 1,cex.axis = 1)
  plot(x=dat[[1]]$p[,1],y=dat[[1]]$p[,2],pch=16,col=dat[[1]]$c,cex=0.7,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  mtext(side=2,text='PC2',line=0.5,las=3)
  mtext(side=1,text='PC1',line=0.0,las=1)
  mtext(text='A',side=3,at=min(dat[[1]]$p[,1]),cex=1.3)
  legend('bottomleft',legend=as.character(1:length(dat[[1]]$b)), pch=16,
    col = colref[1:length(dat[[1]]$b)],horiz=T,box.lwd = 1,box.col = "black",bg = "white",x.intersp=0.5)

  par(new = "TRUE",plt = plotlocs[5,],las = 1,cex.axis = 1)
  plot(x=dat[[2]]$p[,1],y=dat[[2]]$p[,2],pch=16,col=dat[[2]]$c,cex=0.7,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
  mtext(side=1,text='PC1',line=0.0,las=1)
  mtext(text='B',side=3,at=min(dat[[2]]$p[,1]),cex=1.3)
  legend('bottomleft',legend=as.character(1:length(dat[[2]]$b)), pch=16,
    col = colref[1:length(dat[[2]]$b)],horiz=T,box.lwd = 1, box.col = "black",bg = "white",x.intersp=0.5)


dev.off()

