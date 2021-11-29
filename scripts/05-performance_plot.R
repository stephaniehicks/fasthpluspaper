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
n_vec <- seq(1000,10000,by=1000)
m <- 500
#q, number of boostraps
#r, per-bootstrap sample size
#q_vec <- c(1,3,5,10,20,30)
#r_vec <- seq(10,100,by=10) #c(10,2,50,100)
#q <- 30
#r <- 100

#full H+ calculation
#n <- 1000
#n1 <- n2 <- 500
#m <- 500

setups <-  list(
 list(b=c(0.5,0.5),m=c(-0.1,0.1),s=0.5),
 list(b=c(0.25,0.25,0.25,0.25),m=c(-0.2,-0.1,0.1,0.2),s=0.5)
)

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


dat <- lapply(n_vec, function(n) {
  dat <- sapply(1:n, function(i) rnorm(n=m,mean=x$m[j],sd=x$s))) 
  #time distance calculation
  #time adjacency calculation
  #time s+ calculation
  #time hpe calculation
  #time hpb calculation (fixed n+r)
})

#dat <- lapply(setups, function(x) {
#  len <- 1:length(x$b)
#  z <- lapply(len, function(j) sapply(1:round(x$b[j]*n), function(i) rnorm(n=m,mean=x$m[j],sd=x$s)))
#  z <- t(do.call(cbind,z))
#  lab <- as.vector(sapply(len, function(j) rep(j,round(x$b[j]*n))))
#  col <- sapply(lab, function(j) colref[j])
#  pc <- prcomp(z)$x[,1:2]
#  dis <- dist(z)
#  hp <- hpe(D=dis,L=lab,p=10001)$h
#  hpb <- sapply(r_vec, function(r) {
#    sapply(q_vec, function(q) {
#      mean(replicate(q,calcf(z,lab,r),T))
#    })
#  })
#  fin <- list(p=pc, h=hp, hb=hpb, res=t(abs(hp-hpb)), c=col,b=x$b)
#  return(fin)
#})
#nr <- length(r_vec)
#nc <- length(q_vec)
#gres <- dat[[1]]$res
#hres <- dat[[2]]$res
#plotlocs <- rbind(
#  c(0.01,0.51,0.35,0.85), #left heatmap
#  c(0.49,0.99,0.35,0.85), #right heatmap
#  c(0.30,0.70,0.90,0.97), #top colorbar
#  c(0.05,0.48,0.02,0.32), #bottomleft pcaplot
#  c(0.53,0.98,0.02,0.32)  #bottomright pcaplot
#)
#zlim <- c(0,max(rbind(gres,hres)))
#nlevs <- 100
#levs <- seq(0,1,length.out=nlevs)
#col_pal <- viridis(nlevs-1)
#col_pal_fxn <- colorNumeric(col_pal, domain=zlim, na.color = "#808080",alpha = FALSE,reverse = FALSE)
#xlim <- c(-0.5,nc+1)
#ylim <-c(-0.5,nr+0.4)
#zleglocs <- c(0.01,0.50,0.99)
#zleglabs <- formatC(seq(zlim[1],zlim[2],length.out=3),format='f',digits=2)
#txtcut <- 0.05
#xlabs <- as.character(q_vec)
#ylabs <- as.character(r_vec)
#yvals <- sapply(1:nr, function(i) c(i-0.5,i+0.5))

pdf("05-bootstrap_accuracy.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(1,1,1,1))

dev.off()


