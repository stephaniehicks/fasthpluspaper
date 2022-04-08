
# install.packages('devtools')
# library(devtools)
# install_github(repo="ntdyjack/fasthplus", ref = "main")#load relevant packages
library(fasthplus)
library(here)
#code for simulating data 
#https://github.com/stephaniehicks/benchmark-hdf5-clustering/blob/2020-05-07-freeze/scripts/simulate_gauss_mix_k.R

#params to fix
# k 
# means for each group
# sd for each group
#proportion belonging to each group
#num observations and num genes

simulate_gauss_mix_k <- function(n_cells, n_genes,
                                 k, x_mus, x_sds, 
                                 y_mus, y_sds, prop1)
{ 
  
  if(k != length(x_mus)){stop("k is not same as length of x_mus")} 
  if(k != length(x_sds)){stop("k is not same as length of x_sds")} 
  if(k != length(y_mus)){stop("k is not same as length of y_mus")} 
  if(k != length(y_sds)){stop("k is not same as length of y_sds")} 
  if(k != length(prop1)){stop("k is not same as length of prop1")} 
  
  comp1 <- sample(seq_len(k), prob=prop1, size=n_cells, replace=TRUE)
  
  # Sampling locations for cells in each component
  samples1 <- cbind(rnorm(n=n_cells, mean=x_mus[comp1],sd=x_sds[comp1]),
                    rnorm(n=n_cells, mean=y_mus[comp1],sd=y_sds[comp1]))
  #plot(samples1[,1], samples1[,2], col=c(1:15)[comp1])
  
  # Random projection to D dimensional space, to mimic high-dimensional expression data.
  proj <- matrix(rnorm(n_genes*n_cells), nrow=n_genes, ncol=2)
  A1 <- samples1 %*% t(proj)
  
  # Add normally distributed noise.
  A1 <- A1 + rnorm(n_genes*n_cells)
  rownames(A1) <- paste0("Cell", seq_len(n_cells), "-1")
  colnames(A1) <- paste0("Gene", seq_len(n_genes))
  
  list("true_center" = cbind("x" = x_mus, "y" = y_mus),
       "true_cluster_id" = comp1,
       "true_data" = samples1, 
       "obs_data" = A1)
}


data <- simulate_gauss_mix_k(n_cells = 1000, n_genes = 500, 
                     k = 3, x_mus = c(1 ,5, 10) , y_mus = c(1,5,10), 
                     x_sds = c(0.1,0.1,0.1), y_sds = c(0.1,0.1,0.1), prop1 = c(0.3, 0.4, 0.3))

res <- kmeans(data$obs_data, centers = 3)

table(data$true_cluster_id, res$cluster)

fasthplus::hpb(D = data$obs_data, L= res$cluster, t=10,r=10)

##plot x -axis different k, y is a single r 
##plot x- axis different rs, y compare to true H+ 
##easy medium and hard data in terms of means and sds 
##how would balance affect? 



#returns "counts data" and true "cluster id" for each observation 



#take your obs data (500 genes) apply to k means. 
#what happens when r varies, does hpb always get it right 



#Explore what happens when you apply hpb()  varying the number of bootstrap samples (`r).
ptm <- proc.time()
fasthplus::hpb(D = data$obs_data, L= data$true_cluster_id, t=10,r=10)
#D = truedata or obsdata? 
proc.time() - ptm

# user  system elapsed 
# 0.003   0.002   0.006 
x_mus = c(1 ,5, 10) 
y_mus = c(1,5,10) 
n_cells = 1000
n_genes = 3
x_sds = c(0.1,0.1,0.1)
y_sds = c(0.1,0.1,0.1)
comp1 <- sample(seq_len(3), prob= c(0.1, 0.1, 0.8), size=n_cells, replace=TRUE)
samples_1 <- cbind(rnorm(n= n_cells, mean=x_mus[comp1],sd=x_sds[comp1]),
      rnorm(n=n_cells, mean=y_mus[comp1],sd=y_sds[comp1])) #simulate a bivariate distribution

plot(samples_1[,1], samples_1[,2])

#project 2D to higher dim space 
proj <- matrix(rnorm(n_genes*n_cells), nrow=n_genes, ncol=2) 
A1 <- samples_1 %*% t(proj) 
dim(A1)

# Add normally distributed noise - true info is only captured by the two features, 
#projecting into a higher dim space means embedding true info with some noise in higher dim. 

A1 <- A1 + rnorm(n_genes*n_cells)
rownames(A1) <- paste0("Cell", seq_len(n_cells), "-1")
colnames(A1) <- paste0("Gene", seq_len(n_genes))


list("true_center" = cbind("x" = x_mus, "y" = y_mus)
     "true_cluster_id" = comp1,
     "true_data" = samples_1, 
     "obs_data" = A1)

