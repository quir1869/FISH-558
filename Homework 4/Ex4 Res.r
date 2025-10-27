library(mvtnorm)
library(stats4)
library(MASS)
library(coda)
library(ggplot2)
library(ggmcmc)

#Nyear <- 12

#Effort <- c(1,1,2,4,5,6,7,7,4,4,2,2)


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
  
  # Intialize  
  N[1,] <- Nint
  
  for (Iy in 1:Nyear)
  {
    N[Iy+1,2] <- N[Iy,1] * S0
    N[Iy+1,3] <- N[Iy,2] * S1 + N[Iy,3] * (1-gamma) * SJ
    N[Iy+1,4] <- N[Iy,3] * gamma * SJ + N[Iy,4] * (SA - (1-Lambda)*q*Effort[Iy])
    if ((SA - (1-Lambda)*q*Effort[Iy]) < 0)
     {
#      print(SA - (1-Lambda)*q*Effort[Iy])
#      AAAA 
     }  
    if (Iy < Nyear)
      N[Iy+1,1] <- exp(EggVal[Iy])
    else
      N[Iy+1,1] <- N[Iy,4] * b
  }  
 
    # Likelihood function
  
  EggPred <- N[1:Nyear,1]
  Like1 <- sum(-1*dlnorm(EggObs,log(EggPred),SDEgg,log=T))

  AdultPred <- N[1:Nyear,4]
  Like2 <- sum(-1*dlnorm(AdultObs,log(AdultPred),SDAdult,log=T))
  
  CatchPred <- q*N[1:Nyear,4]*Effort
  Like3 <- sum(-1*dlnorm(CatchObs,log(CatchPred),SDCatch,log=T))
  
  Like4 <- 0
  for (Iyr in 1:length(Ayears))
  {
    Jyr <- Ayears[Iyr]
    LikeA <- dmultinom(AgeDataUse[Iyr,],sum(AgeDataUse[Iyr,]),prob=abs(N[Jyr,]),log=T)
    Like4 <- Like4 -1*LikeA
  } 
  
  TotNegLog <- Like1 + Like2 + Like3 + Like4
  #print(N)
  #cat(Like1,Like2,Like3,Like4,TotNegLog,"\n")

  if (Print==T) 
  {
    print(N)
    print(round(EggObs,2))
    print(round(EggPred,2))
    print(round(AdultObs,2))
    print(round(AdultPred,2))
    print(round(CatchObs,3))
    print(round(CatchPred,3))
    for (Iyr in 1:length(Ayears))
    {
      Jyr <- Ayears[Iyr]
      print(round(AgeDataUse[Iyr,]/sum(AgeDataUse[Iyr,]),4))
      print(N[Jyr,])
      print(round(N[Jyr,]/sum(N[Jyr,]),4))
    } 
    cat(Like1,Like2,Like3,Like4,TotNegLog,"\n")
  }
  #cat(Like1,Like2,Like3,Like4,TotNegLog,"\n")
  
  return(TotNegLog)
  
}

# =================================================================================================================

Project <- function(pars,Nproj,Print=F,EffFut)
{
  b <- 15.5; S1 <- 0.7;SJ <- 0.8; SA <- 0.95; gamma <- 0.2
  Nint <- c(100,10,5,1)
  
  S0 <- pars[1]
  Lambda <- pars[2]
  q <- pars[3]
  EggVal <- pars[4:(4+Nyear-1-1)]
  
  EggObs <- DataUsed$EggObs
  SDEgg <- DataUsed$SDEgg
  AdultObs <- DataUsed$AdultObs
  SDAdult <- DataUsed$SDAdult
  CatchObs <- DataUsed$CatchObs
  SDCatch <- DataUsed$SDCatch
  Ayears <- DataUsed$Ayears
  AgeDataUse <- DataUsed$AgeDataUse
  
  N <- array(0,dim=c(Nyear+Nproj+1,4))
  EffortP <- c(Effort,rep(EffFut,Nproj))
  
  # Intialize  
  N[1,] <- Nint
  
  for (Iy in 1:(Nyear+Nproj))
  {
    N[Iy+1,2] <- N[Iy,1] * S0
    N[Iy+1,3] <- N[Iy,2] * S1 + N[Iy,3] * (1-gamma) * SJ
    N[Iy+1,4] <- N[Iy,3] * gamma * SJ + N[Iy,4] * (SA - (1-Lambda)*q*EffortP[Iy])
    if (N[Iy+1,4] < 0) N[Iy+1,4] <- 0
    if (Iy < Nyear)
      N[Iy+1,1] <- exp(EggVal[Iy])
    else
      N[Iy+1,1] <- N[Iy,4] * b
  }  
  Adults <- N[,4]  
  if (Print==T) print(N)
  
  return(Adults)
  
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
    repeat {
      Xnext <- rmvnorm(1, mean=Xcurr, sigma=covar)
      Fnext <- -1*NegLogLike2(Xnext,DataUsed)
      if (!is.na(Fnext)) break; }
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
Nproj <- 20
Adults <- Project(Vectors,Nproj,Print=F,1)
print(Adults)

Cases <- c(0,0,1)
set.seed(666)
#Outs <- DoMCMC(Vectors,DataUsed,14,Varco,Nsim=30000,Nburn=10000,Nthin=10)
#AA

# Apply the mcmc algorithm
if (Cases[1] == 1)
{
 set.seed (123456)
 Outs <- DoMCMC(Vectors,DataUsed,14,Varco,Nsim=300000,Nburn=0,Nthin=10)
 write(t(Outs),"ParVec.out",ncolumn=15)
} 
 
# Plot the diagnostics
if (Cases[2] == 1)
{
 Outs <- matrix(scan("ParVec.Out"),ncol=15,byrow=T)  
 MCMC1 <- mcmc(Outs)
 ggmcmc(ggs(MCMC1),file="E:\\output.pdf")

 # Tabular Diagnostics
 print(geweke.diag(MCMC1))
 print(heidel.diag(MCMC1))
 print(autocorr.diag(MCMC1))

 # Graphical diagnostics
 par(mfrow=c(4,4),mar=c(3,4,2,1))
 traceplot(MCMC1)
 autocorr.plot(MCMC1)
 crosscorr.plot(MCMC1)
} 

# Create the decision table
if (Cases[3] == 1)
{  
  Outs <- matrix(scan("ParVec.Out"),ncol=15,byrow=T)  

  # projections (further thinning()
  thin = 10
  Ivec <- seq(from=1,to=length(Outs[,1]),by=thin)

  Lambda <- Outs[Ivec,2]
  Nvec <- length(Lambda)
  SortLamb <- sort(Lambda)
  cat(SortLamb[Nvec/3],SortLamb[2*Nvec/3],"\n")

  EffFuts <- c(1,5,10)
  NEffFut <- length(EffFuts)

  
  Resu <- matrix(0,ncol=3,nrow=NEffFut)
  colnames(Resu) <- c("low 33%","mid 33%","upper 33%")
  
  Ncnt <- matrix(0,ncol=3,nrow=NEffFut)

  Nproj <- 20

  for (IE in 1:NEffFut)
   {  
    EffFut <- EffFuts[IE]
    for (II in Ivec)
     {  
      Vector <- Outs[II,]
      Adults <- Project(Vector,Nproj,Print=F,EffFut)
      if (Outs[II,2] <= SortLamb[Nvec/3])
       Itype <- 1
      else
       if (Outs[II,2] >= SortLamb[2*Nvec/3])
        Itype <- 3
       else
        Itype <- 2
      Ncnt[IE,Itype] <- Ncnt[IE,Itype] + 1
      Resu[IE,Itype] <- Resu[IE,Itype] + Adults[Nyear+Nproj+1]
     }
   }  

  print(Ncnt)
  print("Decision table")
  print(Resu/Ncnt)
}
  



