cl1 <- '#8338EC' #Purple
cl2 <- '#06D6A0' #Teal
clw <- '#0000ff64' #blue (transparent)
clb <- '#ff000064' #red (transparent)

library(pbapply)


set.seed(1234)
n <- 1000
m <- 500
Nd <- (n*(n-1))/2
Nz <- (Nd*(Nd-1))/2

#datasets
#no difference 50:50, no difference 10:90
#some overlap 50:50, some overlap 10:90
#no overlap 50:50, no overlap 10:90
#fix variance = 1

neven <- n*0.5
nshort <- n*0.1
nlong <- n*0.9

dat <- list(

  list(
    x=sapply(1:n, function(i) rnorm(n=m,mean=0,sd=1)),
    l=round(runif(n=n))
  ),

  list(
    x=sapply(1:n, function(i) rnorm(n=m,mean=0,sd=1)),
    l=c(rep(0,nshort),rep(1,nlong))
  ),

  list(
    x=cbind(sapply(1:neven, function(i) rnorm(n=m,mean=0.1,sd=1)),sapply(1:neven, function(i) rnorm(n=m,mean=-0.1,sd=1))),
    l=c(rep(0,neven),rep(1,neven))
  ),

  list(
    x=cbind(sapply(1:nshort, function(i) rnorm(n=m,mean=0.1,sd=1)),sapply(1:nlong, function(i) rnorm(n=m,mean=-0.1,sd=1))),
    l=c(rep(0,nshort),rep(1,nlong))
  ),

  list(
    x=cbind(sapply(1:neven, function(i) rnorm(n=m,mean=0.3,sd=1)),sapply(1:neven, function(i) rnorm(n=m,mean=-0.3,sd=1))),
    l=c(rep(0,neven),rep(1,neven))
  ),

  list(
    x=cbind(sapply(1:nshort, function(i) rnorm(n=m,mean=0.3,sd=1)),sapply(1:nlong, function(i) rnorm(n=m,mean=-0.3,sd=1))),
    l=c(rep(0,nshort),rep(1,nlong))
  )

)

l <- length(dat)

#compute distances and pca for each

dat <- pblapply(dat, function(d) {
  pc <- prcomp(t(d$x))$x[,1:2]

  cols <- ifelse(d$l==0,cl1,cl2)
  cols <- paste0(cols,'64')

  ind <- sapply(d$l, function(x) sapply(d$l, function(y) x==y))
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  ib <- which(!ind) 

  dis <- as.matrix(dist(t(d$x))) 
  dis <- dis[upper.tri(dis)]
  dw <- sort(dis[iw])
  db <- sort(dis[ib])
  alf <- as.numeric(length(iw)) / as.numeric(length(ind))

  fin <- list(p = pc, c = cols, w = dw, b = db, a=alf)
  return(fin)
})

binmin <- min(sapply(dat, function(x) min(c(x$w,x$b))))
binmax <- max(sapply(dat, function(x) max(c(x$w,x$b))))
bins <- seq(binmin,binmax,length.out=20)


#g+ and h+ values
perfs <- pbsapply(dat, function(d) {
  sp <- sum(sapply(d$w, function(x) sum(x>d$b)))
  #gp <- (2*sp) / Nz
  gp <- sp / Nz
  hp <- sp / (as.numeric(length(d$w))*as.numeric(length(d$b)))
  c(g=gp,h=hp)
})


#2 rows (r1=PCA, r2=histgram), 6 columns (datasets)
r1t <- 0.95
r1b <- 0.60
r2t <- 0.50
r2b <- 0.15
r3t <- 0.10
r3b <- 0.00

c1l <- 0.00
c1r <- 0.10
c2l <- 0.11
c2r <- 0.24
c3l <- 0.26
c3r <- 0.39
c4l <- 0.41
c4r <- 0.54
c5l <- 0.56
c5r <- 0.69
c6l <- 0.71
c6r <- 0.84
c7l <- 0.86
c7r <- 0.99


plotlocs <- rbind(
  c(c2l,c2r,r1b,r1t), #row1 col2
  c(c3l,c3r,r1b,r1t), #row1 col3
  c(c4l,c4r,r1b,r1t), #row1 col4
  c(c5l,c5r,r1b,r1t), #row1 col5
  c(c6l,c6r,r1b,r1t), #row1 col6
  c(c7l,c7r,r1b,r1t), #row1 col7

  c(c2l,c2r,r2b,r2t), #row2 col2
  c(c3l,c3r,r2b,r2t), #row2 col3
  c(c4l,c4r,r2b,r2t), #row2 col4
  c(c5l,c5r,r2b,r2t), #row2 col5
  c(c6l,c6r,r2b,r2t), #row2 col6
  c(c7l,c7r,r2b,r2t), #row2 col7

  c(c2l,c2r,r3b,r3t), #row3 col2
  c(c3l,c3r,r3b,r3t), #row3 col3
  c(c4l,c4r,r3b,r3t), #row3 col4
  c(c5l,c5r,r3b,r3t), #row3 col5
  c(c6l,c6r,r3b,r3t), #row3 col6
  c(c7l,c7r,r3b,r3t), #row3 col7
  
  c(c1l,c1r,r1b,r1t), #row 1 col 1
  c(c1l,c1r,r2b,r2t) #row 2 col 2
)

xticks_hist <- pretty(bins,3)
yticks_hist <- c(0,10^4,10^5,10^6)
ytclbs_hist <- c('0',expression('10'^4),expression('10'^5),expression('10'^6))
ylims_hist <- c(0,103500)

alphloc <- c(37,8e4)

pdf("01-motivating_figure.pdf",width=10,height=5)
  plot.new()

  #1st row (PCA plots)
  for(i in 1:l){
    par(new = "TRUE",plt = plotlocs[i,],las = 1,cex.axis = 1)
    plot(x=dat[[i]]$p[,1],y=dat[[i]]$p[,2],pch=16,col=dat[[i]]$c,cex=0.7,xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n')
    if(i==1){mtext(side=2,text='PC2',line=1.5,las=3) }
    mtext(side=1,text='PC1',line=0.5,las=1)
    mtext(text=LETTERS[i],side=3,at=min(dat[[i]]$p[,1]),cex=1.5)
    if(i == 2){
      par(xpd=T)
      abline(v=1.15*max(dat[[i]]$p[,1]),lty=2,col='grey60',lwd=2)
      par(xpd=F)
    }else if(i ==4){
      par(xpd=T)
      abline(v=1.10*max(dat[[i]]$p[,1]),lty=2,col='grey60',lwd=2)
      par(xpd=F)      
    }
  }

  #2nd row (distance histograms)
  for(i in 1:l){
    par(new = "TRUE",plt = plotlocs[l+i,],las = 1,cex.axis = 1)
    hist(x=dat[[i]]$w,breaks=bins,main='',xlab='',ylab='',plot=T,border='blue',col=clw, freq=T,xaxs = "i",yaxs = "i",xaxt='n',yaxt='n',ylim=ylims_hist)
    hist(x=dat[[i]]$b,breaks=bins,add=T,border='red',col=clb,freq=T,ylim=ylims_hist)
    if(i==1){mtext(side=2,text='Frequency',line=1.5,las=3)}
    mtext(side=1,text='Dist. (L2)',line=1.0,las=1)
    #mtext(text=letters[l+i],side=3,at=bins[1],cex=1.5,line=0.5)
    axis(side=1,at=xticks_hist,las=1,mgp=c(3, .2, 0),line=0.1,cex.axis=0.8)
    if(i==1){axis(side=2,at=yticks_hist,labels=ytclbs_hist,las=2,mgp=c(3, .5, 0),cex.axis=0.8)}
    alphtxt <- formatC(dat[[i]]$a,digits=2,format='f')
    text(x=alphloc[1],y=alphloc[2], labels=bquote(alpha == .(alphtxt)),cex=0.8)
  }

  #3rd row (h+ and g+ values)
  for(i in 1:l){
    par(new = "TRUE",plt = plotlocs[(2*l)+i,],las = 1,cex.axis = 1)
    plot(x=0,y=0,type='n',ylim=c(0,1),xlim=c(0,1),xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
    text(x=0.5,y=0.3,labels=paste0('G+=',formatC(perfs[1,i],digits=2,format='f'),' H+=',formatC(perfs[2,i],digits=2,format='f')),cex=0.8)
  }

  #far left columns (legends)
  
  par(new = "TRUE",plt = plotlocs[(3*l)+1,],las = 1,cex.axis = 1)
  plot(x=0,y=0,type='n',ylim=c(0,1),xlim=c(0,1),xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',bty='n') 
  #legend(x=0.0,y=0.5,legend=c('cl1','cl2'), pch=16, col = c(cl1,cl2),bty='n',horiz=F)
  points(x=c(0.15,0.15),y=c(0.4,0.6),pch=c(16,16),col=c(cl1,cl2),cex=1.5)
  text(x=c(0.35,0.35),y=c(0.6,0.4),labels=c(expression('Cl'[1]), expression('Cl'[2])),cex=1.2)

  par(new = "TRUE",plt = plotlocs[(3*l)+2,],las = 1,cex.axis = 1) 
  plot(x=0,y=0,type='n',ylim=c(0,1),xlim=c(0,1),xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
  #legend(x=0.0,y=0.5,legend=c('Dw','Db'), pch=22, col = c('blue','red'),pt.bg=c(clw,clb),bty='n',horiz=F)
  points(x=c(0.15,0.15),y=c(0.4,0.6),pch=c(22,22),col= c('blue','red'),cex=1.5,bg=c(clw,clb))
  text(x=c(0.35,0.35),y=c(0.4,0.6),labels= c(expression('D'[W]), expression('D'[B])),cex=1.2)
dev.off()

