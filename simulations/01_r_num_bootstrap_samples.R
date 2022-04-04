
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


data <- simulate_gauss_mix_k(n_cells = 1000, n_genes = 10, 
                     k = 3, x_mus = c(1 ,2, 3) , y_mus = c(1,2,3), 
                     x_sds = c(1,1,1), y_sds = c(1,1,1), prop1 = c(0.3, 0.4, 0.3))

#returns "counts data" and true "cluster id" for each observation 

#qsns
#what's an ideal num n_cells and n_genes 
#counts data has negative values
#how do you generate fake cluster data for each observation? 



#Explore what happens when you apply hpb()  varying the number of bootstrap samples (`r).
ptm <- proc.time()
fasthplus::hpb(D = data$obs_data, L= data$true_cluster_id, t=10,r=10)
#D = truedata or obsdata? 
proc.time() - ptm

# user  system elapsed 
# 0.003   0.002   0.006 



