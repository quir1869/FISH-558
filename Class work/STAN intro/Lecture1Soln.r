library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("C:\\courses\\FISH 558_25\\Lectures\\Stan\\")


GenerateData <- function()
 {
  set.seed(777)
  Npop <- 20
  Nsamp <- 15
  LogA <- log(1) + rnorm(Npop,0,0.05)
  LogB <- log(3) + rnorm(Npop,0,0.1)
  print(LogA)
  print(LogB)
  sigma <- 0.2

  YY <- NULL
  for (Ipop in 1:Npop)
   {
    Lengths <- runif(Nsamp,0,100)
    XX <- cbind(rep(Ipop,Nsamp),Lengths,exp(LogA[Ipop])*Lengths^exp(LogB[Ipop])*exp(rnorm(Nsamp,0,sigma)))
    YY <- rbind(YY,XX)
   }
  ZZ <- cbind(LogA,LogB,exp(LogA),exp(LogB))
  write(t(YY),"Ex1d.DAT",ncolumns=3)
  write(t(ZZ),"True.DAT",ncolumns=4)

 }

 #GenerateData()
 Npop <- 20
 Nsamp <- 15
 TheData <- read.table("Ex1d.DAT",header=F)
 Index <- TheData[,1]
 Lengths <- TheData[,2]
 Weights <- TheData[,3]
 
 the_data <- list(Index=Index,Lengths=Lengths,Weights=Weights,Npop=Npop,Nsamp=Nsamp,N=Nsamp*Npop)
 inits1 <- list(log_mean_alpha=log(0.01),log_sigma_alpha=log(0.05),
                log_mean_beta=log(3),log_sigma_beta=log(0.05),
                alphas=rep(log(0.01),Npop),betas=rep(log(3),Npop),log_sigma=0)
 inits <- list(inits1)
 print(inits1)
 
 fit <- stan(file = 'Ex1d.stan', data = the_data, 
             iter = 10000, chains = 3,init=inits,verbose=F)
 la <- extract(fit,permuted = TRUE)
 par(mfrow=c(2,2))
 hist(la$mean_alpha)
 hist(la$mean_beta)
 
  print(fit)
 