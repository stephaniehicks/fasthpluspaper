library(here)
pdf(here("figures", "01-mot-fig-balance.pdf"), width = 10)

par(mfrow=c(1,2))

n <- 100 # observations
m <- 100 # features
cl1 <- sapply(1:(n), function(i) rnorm(n=m,mean=0.5,sd=1))
cl2 <- sapply(1:(n), function(i) rnorm(n=m,mean=-0.5,sd=1))
dat <- cbind(cl1,cl2)
d <- dist(t(dat))
dvec <- as.matrix(d)
dvec <- dvec[upper.tri(dvec)]
l <- c(rep(0,n),rep(1,n))
ind <- sapply(l, function(x) x==l)
image(ind, main = "balanced classes (alpha = 0.5)")
ind <- ind[upper.tri(ind)]
iw <- which(ind)
ib <- which(!ind)
dw <- dvec[iw]
db <- dvec[ib]
length(dw) / (length(dw)+length(db)) # should be 0.5

n <- 100 # observations
m <- 100 # features
cl1 <- sapply(1:(n-n*.8), function(i) rnorm(n=m,mean=0.5,sd=1))
cl2 <- sapply(1:(n+n*.8), function(i) rnorm(n=m,mean=-0.5,sd=1))
dat <- cbind(cl1,cl2)
d <- dist(t(dat))
dvec <- as.matrix(d)
dvec <- dvec[upper.tri(dvec)]

l <- c(rep(0,n-n*.8),rep(1,n+n*.8))
ind <- sapply(l, function(x) x==l)
image(ind, main = "imbalanced classes (alpha = 0.82)")
ind <- ind[upper.tri(ind)]
iw <- which(ind)
ib <- which(!ind)
dw <- dvec[iw]
db <- dvec[ib]
length(dw) / (length(dw)+length(db)) # should this be 0.82 or 0.81? 

dev.off()



