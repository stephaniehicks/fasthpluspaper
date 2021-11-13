library(here)
set.seed(1234)
colref <- c("#1AFF1A","#4B0092")

pdf(here("figures", "supp01-balance_figure.pdf"), width = 10,height=8)
  par(mfrow=c(2,2),mar=c(0.6,1.1,2.1,1.1))

  n <- 1000 # observations
  bals <-  c(0.5,0.5)
  l <- unlist(lapply(1:length(bals), function(i) rep(i,n*bals[i])))
  ind <- sapply(rev(l), function(x) x==l)
  image(ind, main = "",xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',useRaster=T,col = colref)
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  alpha <- formatC(round(length(iw) / (length(ind)),digits=2),digits=2,format='f')
  mtext(text= bquote( '2 balanced classes (' ~ alpha == .(alpha) ~ ')'))
  abline(a=1,b=-1,lty=3,lwd=2,col='black')

  bals <-  c(0.9,0.1)
  l <- unlist(lapply(1:length(bals), function(i) rep(i,n*bals[i])))
  ind <- sapply(rev(l), function(x) x==l)
  image(ind, main = "",xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',useRaster=T,col = colref)
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  alpha <- formatC(round(length(iw) / (length(ind)),digits=2),digits=2,format='f')
  mtext(text= bquote( '2 imbalanced classes (' ~ alpha == .(alpha) ~ ')'))
  abline(a=1,b=-1,lty=3,lwd=2,col='black')

  bals <-  rep(0.1,10)
  l <- unlist(lapply(1:length(bals), function(i) rep(i,n*bals[i])))
  ind <- sapply(rev(l), function(x) x==l)
  image(ind, main = "",xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',useRaster=T,col = colref)
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  alpha <- formatC(round(length(iw) / (length(ind)),digits=2),digits=2,format='f')
  mtext(text= bquote( '10 balanced classes (' ~ alpha == .(alpha) ~ ')'))
  abline(a=1,b=-1,lty=3,lwd=2,col='black')

  #tmp <- sample(1:20000,10)
  #bals <- tmp/sum(tmp)
  bals <- c(0.40,0.18,0.14,0.09,0.06,0.05,0.04,0.02,0.01,0.01)
  l <- unlist(lapply(1:length(bals), function(i) rep(i,n*bals[i])))
  ind <- sapply(rev(l), function(x) x==l)
  image(ind, main = "",xaxs = "i",yaxs = "i",xlab='',ylab='',xaxt='n',yaxt='n',useRaster=T,col = colref)
  ind <- sapply(l, function(x) x==l)
  ind <- ind[upper.tri(ind)]
  iw <- which(ind)
  alpha <- formatC(round(length(iw) / (length(ind)),digits=2),digits=2,format='f')
  mtext(text= bquote( '10 imbalanced classes (' ~ alpha == .(alpha) ~ ')'))
  abline(a=1,b=-1,lty=3,lwd=2,col='black')

dev.off()
