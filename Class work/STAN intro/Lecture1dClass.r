library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# setwd("C:\\courses\\FISH 558_25\\Lectures\\Stan\\")



 Analyze <- function()
  {
   Npop <- 20
   Nsamp <- 15
   TheData <- read.table("Ex1d.DAT",header=F)
   Index <- TheData[,1]
   Lengths <- TheData[,2]
   Weights <- TheData[,3]
   
   the_data <- list(Index=Index,Lengths=Lengths,Weights=Weights,Npop=Npop,Nsamp=Nsamp,N=Nsamp*Npop)
   inits1 <- list(sigma = 0.05, sigma_a = 0.05, sigma_b = 0.05)
   inits <- list(inits1)

   fit <- stan(file = 'Ex1dclass.stan', data = the_data, 
               iter = 10000, chains = 1,init=inits,verbose=F)
   la <- extract(fit,permuted = TRUE)
   par(mfrow=c(2,2))
   hist(la$mean_alpha)
   hist(la$mean_beta)
   
 }
 
 Fit <- Analyze()
 
 print(Fit)
 
 
 library(tidyverse)
 
 y <- rnorm(100, 5, 2)
 data_list <- list(N = length(y), y = y)
 
 fit <- stan(file = "Ex1dclass.stan", data = data_list, 
             iter = 1000, chains = 4)

 
 pairs(fit)
 