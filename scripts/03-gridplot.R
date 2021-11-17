#revison D, 3x1 plot
set.seed(1234)
b <- rnorm(mean=-0.3,n=10000)
a <- rnorm(mean=0.3,n=10000)
ht <- sum(sapply(a, function(x) sum(x>b)))/(length(a)*length(b))

syvl <- sapply(sort(a), function(x) sum(x>b)/length(b))
sxvl <- seq(0,1,length.out=length(syvl))

pltlocs <- list(
  c(0.07,0.37,0.16,0.91),
  c(0.38,0.68,0.16,0.91),
  c(0.69,0.99,0.16,0.91)
)
pvec <- c(11,26,51)
lwds <- c(2.0,1.0,0.5)
lwds2 <- c(6,5,4)

axloc <- c(0.01,0.50,0.98)
#axtxt <- c('0','','p-1')

pdf('03-gridplot.pdf',width=8,height=3)
plot.new()

for(k in 1:3) {
  p <- pvec[k]
  ps <- seq(0,1,length.out=p)
  qa <- quantile(a, probs = ps)
  qb <- quantile(b, probs = ps)
  l <- length(ps)
  clz <- comp <- t(sapply(qa, function(x) x>qb))
  clz[which(comp)] <- '#0072B2'
  clz[which(!comp)] <- '#D41159'
  comp <- sapply(ps, function(x) x*ps)
  brds<- matrix('black',nrow=nrow(comp),ncol=ncol(comp))
  brds2 <- matrix(NA,nrow=nrow(comp),ncol=ncol(comp))
  cut <- 1/(p-1)
  ai <- bi <- 1
  inty <- rep(0,l)
  while(ai <= l & bi <= l){
    brds2[ai,bi] <- '#56B4E9'
    si <- qa[ai] > qb[bi]
    if(si){ #increase qB
      bi <- bi + 1
    }else{ #increase qA
      inty[ai] <- (bi-1)
      ai <- ai + 1
    }
  }
  brds2[abs(comp-ht)<=cut & !is.na(brds2) ] <- '#009E73'
  brds2[abs(comp-ht)<=cut &  is.na(brds2)] <- '#F0E442'  
  par(new = "TRUE",plt = pltlocs[[k]],las = 1, cex.axis = 1)
  plot.new()
  plot.window(xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),ylab='',xlab='',xaxt='n',yaxt='n',xaxs = "i",yaxs = "i",bty='n')
  #axis(side=1,at=seq(0,1,by=0.2),las=1,mgp=c(0, .3, 0),line=0.0,cex.axis=0.7,labels=c('0','0.2','0.4','0.6','0.8','1'))
  #axis(side=1,at=axloc,las=1,mgp=c(0, .3, 0),line=0.0,cex.axis=0.7,labels=c('0','i','p-1'))
  mtext(side=1,text=c('0','i','p'),at=axloc,cex=0.8)
  mtext(side=1,text='q(A)',line=1.2,las=1,cex=1.2)
  if(k==1){
    #axis(side=2,at=seq(0,1,by=0.20),las=1,mgp=c(0, .5, 0),line=0.0,cex.axis=0.7,c('0','0.2','0.4','0.6','0.8','1'))
    #axis(side=2,at=axloc,las=1,mgp=c(0, .5, 0),line=0.0,cex.axis=0.7,labels=c('0','j','p-1'))
    mtext(side=2,text=c('0','j','p'),at=axloc,cex=0.8,las=1,line=0.1)
    mtext(side=2,text='q(B)',line=1.5,las=3,cex=1.2)
  }

  for(i in 1:l){
    for(j in 1:l){
      rect(xleft=(i-1)/l,xright=i/l,ybottom=(j-1)/l,ytop=j/l, col=clz[i,j],border=brds[i,j],lwd=lwds[k])
    }
  }

  for(i in 1:l){
    for(j in 1:l){
      rect(xleft=(i-1)/l,xright=i/l,ybottom=(j-1)/l,ytop=j/l, col=NA,border=brds2[i,j],
        lwd=ifelse(brds2[i,j]=='#009E73',lwds2[k],lwds2[k]*0.5))
    }
  }  

  lines(x=sxvl,y=syvl,type='l',lwd=1.5,col='white')  
  mtext(text=LETTERS[k],side=3,line=0.0,at=0,cex=1.5)
}	
	
dev.off()
