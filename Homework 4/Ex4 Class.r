library(mvtnorm)
library(stats4)
library(MASS)
library(coda)


# ================================================================================

NegLogLike2 <- function(pars,DataUsed,Print=F)
{
  b <- 15.5; S1 <- 0.7;SJ <- 0.8; SA <- 0.95; gamma <- 0.2
  Nint <- c(100,10,5,1)
  
  Nyear <- DataUsed$Nyear

  S0 <- pars[1]
  Lambda <- pars[2]
  q <- pars[3]
  EggVal <- pars[4:(4+Nyear-1-1)]
  
  Effort = DataUsed$Effort
  EggObs <- DataUsed$EggObs
  SDEgg <- DataUsed$SDEgg
  AdultObs <- DataUsed$AdultObs
  SDAdult <- DataUsed$SDAdult
  CatchObs <- DataUsed$CatchObs
  SDCatch <- DataUsed$SDCatch
  Ayears <- DataUsed$Ayears
  AgeDataUse <- DataUsed$AgeDataUse
    
  N <- array(0,dim=c(Nyear+1,4))

  # This gets retrurned
  return(TotNegLog)
  
}

# =================================================================================================================

DoMCMC<-function(Xinit,DataUsed,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1)
{
  Xcurr <- Xinit
  Fcurr <- -1*NegLogLike2(Xcurr,DataUsed)
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
  Ipnt <- 0; Icnt <- 0
  for (Isim in 1:Nsim)
  {
    repeat
     {
      Xnext <- rmvnorm(1, mean=Xcurr, sigma=covar)
      Fnext <- -1*NegLogLike2(Xnext,DataUsed)
      if (!is.na(Fnext)) break
     }
    Rand1 <- log(runif(1,0,1))
    if (Fnext > Fcurr+Rand1)
    {Fcurr <- Fnext; Xcurr <- Xnext }   
    if (Isim > Nburn & Isim %% Nthin == 0)
    {
      Icnt <- Icnt + 1; Outs[Icnt,] <- c(Xcurr,Fcurr); cat("saving",Icnt,"\n")     
    }
  } 
  xx <- seq(1:Icnt)
  par(mfrow=c(4,4),mar=c(3,4,2,1))
  for (II in 1:(Ndim+1))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- exp(yy)
    plot(xx,yy,xlab="Cycle number",ylab=lab1,type='b',pch=16,cex=0.02)
  }
  par(mfrow=c(4,4),mar=c(3,4,2,1))
  for (II in 1:(Ndim+1))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- exp(yy)
    hist(yy,ylab=lab1,main="")
  }
  return(Outs[1:Icnt,])
}

# ================================================================================

setwd("C:\\courses\\FISH 558_25\\Bayes Workshop FHL\\Assignments\\Assignment 4\\")

OFile <- "Ex4a.dat"
Nyear <- scan(OFile,skip=1,n=1,quiet=T)
Effort <- scan(OFile,skip=5,n=Nyear,quiet=T)
EggsObs <- scan(OFile,skip=7,n=Nyear,quiet=T)
SDEgg <- scan(OFile,skip=9,n=1,quiet=T)
AdultObs  <-scan(OFile,skip=11,n=Nyear,quiet=T)
SDAdult <- scan(OFile,skip=13,n=1,quiet=T)
CatchObs  <-scan(OFile,skip=15,n=Nyear,quiet=T)
SDCatch <- scan(OFile,skip=17,n=1,quiet=T)
print(SDCatch)
Nayears <- scan(OFile,skip=20,n=1,quiet=T)
Ayears <- scan(OFile,skip=22,n=Nayears,quiet=T)
AgeComp <- matrix(scan(OFile,skip=24,n=Nayears*4,quiet=T),ncol=4,byrow=T)

OFile <- "Ex4b.dat"
Npars <- 14
Vectors <- scan(OFile,skip=2,n=Npars,quiet=T)
SDVec <- scan(OFile,skip=4,n=Npars,quiet=T)
Varco <- matrix(scan(OFile,skip=6,Npars*Npars,quiet=T),ncol=Npars)

DataUsed <- list(EggObs=EggsObs,SDEgg=SDEgg,AdultObs=AdultObs,SDAdult=SDAdult,CatchObs=CatchObs,SDCatch=SDCatch,Ayears=Ayears,AgeDataUse=AgeComp,Effort=Effort,Nyear=Nyear)
#print(DataUsed)
print(Vectors)
NegLogLike2(Vectors,DataUsed=DataUsed,Print=T)

# use

# Convert to Coda
MCMC1 <- mcmc(Outs)

