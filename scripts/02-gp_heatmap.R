#library(pbapply)
library(doParallel)
library(RColorBrewer)
library(viridis)
library(leaflet)

set.seed(1234)
#set up parallel backend
cl <- makeCluster(12)
registerDoParallel(cl) #needs to be a socket [] on windows

calcs <- function(b){
  n <- 1000
  m <- 500
  Nd <- (n*(n-1))/2
  Nz <- (Nd*(Nd-1))/2
  sigs <- seq(0.00,0.5,by=0.05)
  res <- sapply(sigs, function(s) {
    n1 <- b*n
    n2 <- (1-b)*n
    x1 <- sapply(1:n1, function(i) rnorm(n=m,mean=s,sd=1))
    x2 <- sapply(1:n2, function(i) rnorm(n=m,mean=-1*s,sd=1))
    x <- cbind(x1,x2)
    l <- c(rep(0,n1),rep(1,n2))
    ind <- sapply(l, function(x) sapply(l, function(y) x==y))
    ind <- ind[upper.tri(ind)]
    iw <- which(ind)
    ib <- which(!ind) 
    dis <- as.matrix(dist(t(x))) 
    dis <- dis[upper.tri(dis)]
    dw <- sort(dis[iw])
    db <- sort(dis[ib])
    sp <- sum(sapply(dw, function(x) sum(x>db)))
    #gp <- (2*sp)/Nz
    gp <- sp / Nz
    hp <- sp / (as.numeric(length(dw))*as.numeric(length(db)))
    a <- as.numeric(length(iw)) / as.numeric(length(ind))
    c(gp,hp,a)
  })
  return(res)
}

bals <- seq(0.05,0.5,by=0.05)
sres <- parLapply(cl=cl, X=bals, fun=calcs)

gres <- sapply(sres, function(x) x[1,] )
hres <- sapply(sres, function(x) x[2,] )

#make paired heatmaps (colorbar on top?) for G+ and H+
#plotting parameters
nr <- nrow(gres)
nc <- ncol(gres)

plotlocs <- rbind(
  c(0.05,0.53,0.05,0.82), #left heatmap
  c(0.51,0.99,0.05,0.82), #right heatmap
  c(0.3,0.7,0.9,0.97) #top colorbar
)

zlim <- c(0,max(rbind(gres,hres)))
nlevs <- 100
levs <- seq(0,1,length.out=nlevs)
col_pal <- viridis(nlevs-1)
col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)

xlim <- c(-0.5,nc+1)
ylim <-c(-0.5,nr+0.4)
zleglocs <- c(0.01,0.50,0.99)
zleglabs <- c("0.0","0.25","0.5")
txtcut <- 0.23

xlabs <- sapply(bals, function(x) {
 b <- (1-x)*100
 a <- 100*x
 c <- paste0(a,":",b)
 return(c)
})

#xlabs2 <- 

yvals <- sapply(1:nr, function(i) c(i-0.5,i+0.5))
ylabs <- formatC(seq(0.00,0.5,by=0.05)*2,digits=1,format='f')

pdf("02-signal_heatmap.pdf",width=10,height=8)
  plot.new()
  par(new = "TRUE",plt = plotlocs[1,],las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim, ylim, main="", xaxs = 'i', yaxs = 'i', asp = NA)
  for(k in 1:nc){
    for(j in 1:nr){
      rect(xleft=k-0.5,xright=k+0.5,ybot=yvals[1,j],ytop=yvals[2,j],border='black',col=col_pal_fxn(gres[j,k]))
      text(x=k,y=mean(yvals[,j]),labels=formatC(gres[j,k],digits=2,format='f'),cex=0.6,col=ifelse(gres[j,k]>txtcut,'black','white'))
    }
  }
  text(x=1:nc,y=c(0,0.3),labels=xlabs,cex=0.8)
  text(x=0.1,y=1:nr,labels=ylabs,cex=1.0,srt=0)
  mtext(side=3,text='G+',cex=1.2)
  par(xpd = TRUE)
  #text(x=-1.0,y=mean(yvals),labels=expression("E[Cl"[1] * "] - E[Cl"[2] * "]"),srt=90,cex=1.5)
  text(x=-1.0,y=mean(yvals),labels= expression("Mean cluster difference: E[Cl"[1] * "] - E[Cl"[2] * "]"),srt=90,cex=1.5)
  #text(x=11,y=-0.6,labels=expression("Cl"[1] * ":Cl"[2]),cex=1.5)
  text(x=11,y=-0.6,labels=expression(" Cl"[1] * ":Cl"[2]),cex=1.5)
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
  text(x=1:nc,y=c(0,0.3),labels=xlabs,cex=0.8)
  mtext(side=3,text='H+',cex=1.2)

  par(new = "TRUE",plt = plotlocs[3,],las = 1,cex.axis = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(0,1), xaxs = "i",yaxs = "i")
  rect(xleft=levs[-length(levs)], ybottom=0.0, xright=levs[-1L], ytop=1.0, col=col_pal, border=NA)
  box()
  axis(side=1,labels=zleglabs,at=zleglocs,cex=0.7,las=1,mgp=c(1.0, .4, 0))
  mtext(side=3,text='G+ / H+', cex=1.3)

dev.off()

