Ex3Soln <- function()
{
 # for case B set remove=2
 ProbData1 <- DoSir(Nout=200,Model=F,Remove=0) 
 ProbData2 <-DoSir(Nout=200,Model=T,Remove=0) 
 cat("prob Model 1 =",ProbData1/(ProbData1+ProbData2),ProbData1/(ProbData2),"\n")
 
}

# =================================================================================

DoSir <- function(Nout=1000,Model,Remove=0)
{

 set.seed(666)
 TheData <- ReadData()
 Yr1 <- TheData$CatchYr[1]
 Catch <- TheData$CatchVal
 Nsurv <- length(TheData$SurveyEst)-Remove
 SurveyEst <- TheData$SurveyEst[1:Nsurv]
 SurveyCV <- TheData$SurveyCV[1:Nsurv]
 SurveyYr <- TheData$SurveyYr[1:Nsurv]
 Nyears <- length(Catch)

 Threshold <- exp(0)
 Cumu <- 0
 Ndone <- 0
 AboveK <- 0
 Vals <- matrix(0,ncol=5,nrow=Nout)
 PopOut <- matrix(0,ncol=Nyears,nrow=Nout)
 AveLike <- 0
 Ntest <- 0
 while (Ndone < Nout)
  {
   r <- runif(1,0,0.15)
   K <- runif(1,20000,50000)
   Pop1965 <- runif(1,10000,15000)
   AddCV <- runif(1,0.1,0.2)
   #r <- 0.1
   #K <- 40000
   #Pop1965 <- 12000
   #AddCV <- 0.15
   #Model <- F
   #r <- 0.09579
   #K <- 22871.4
   #Pop1965 <- 10297.3
   #AddCV <- 0.08916
   #Model <- F
   Pop <- PopModel(Catch,r,K,Pop1965,Model)
   #print(Pop)
   #cat(r,K, Pop1965,AddCV,"\n")
   #AAA
   NegLogLike <- Likelihood(Pop,SurveyYr-Yr1+1,SurveyEst,SurveyCV,AddCV)
   #print(NegLogLike)
   TheLike <- exp(-1*NegLogLike-32.19)
   #print(TheLike)
   Cumu <- Cumu + TheLike
   AveLike <- AveLike + TheLike
   Ntest <- Ntest + 1
   while (Cumu > Threshold & Ndone < Nout)
    {
     cat("Saving",Ndone,"\n")
     Cumu <- Cumu - Threshold
     Ndone <- Ndone + 1
     Vals[Ndone,1] <- r
     Vals[Ndone,2] <- K
     Vals[Ndone,3] <- Pop1965
     Vals[Ndone,4] <- AddCV
     Vals[Ndone,5] <- Pop[Nyears]/K
     PopOut[Ndone,] <- Pop
     #print(Vals[Ndone,])
     #print(TheLike)
     #print(NegLogLike)
     #AAAA
     if (max(Pop/K) > 0.9) AboveK <- AboveK + 1/Nout
    }
  }
 cat(AboveK,AveLike/Ntest,Ntest,"\n")

 par(mfrow=c(3,2))
 hist(Vals[,1],main="",xlab="r")
 hist(Vals[,2],main="",xlab="K")
 hist(Vals[,3],main="",xlab="Population size 1965")
 hist(Vals[,4],main="",xlab="Addition CV")
 if (Model==F) hist(Vals[,5],main="",xlab="Depletion")

 quants <- matrix(0,nrow=5,ncol=Nyears)
 for (Iyear in 1:(Nyears))
 quants[,Iyear] <- quantile(PopOut[,Iyear],prob=c(0.05,0.25,0.5,0.75,0.95))
 ymax <- max(quants)*1.2
 plot(TheData$CatchYr,quants[3,],ylim=c(0,ymax),xlab="year",ylab="population size",type="l")
 xx <- c(TheData$CatchY,rev(TheData$CatchY))
 yy <- c(quants[1,],rev(quants[5,]))
 polygon(xx,yy,col="gray10")
 xx <- c(TheData$CatchY,rev(TheData$CatchY))
 yy <- c(quants[2,],rev(quants[4,]))
 polygon(xx,yy,col="gray80")
 lines(TheData$CatchYr,quants[3,],col="red",lwd=3)
 for (Iyear in 1:length( SurveyYr))
   points(SurveyYr[Iyear], SurveyEst[Iyear],pch=16,col="green") 
 return(AveLike/Ntest)
}

# =================================================================================

Likelihood <- function(Pop,SurveyYr,SurveyEst,SurveyCV,AddCV)
{
 UseCV <- sqrt(SurveyCV^2+AddCV^2)
 Preds <- Pop[SurveyYr]
 Residuals <- log(UseCV)+0.5*(log(Preds)-log(SurveyEst))^2/UseCV^2
 LogLike <- sum(Residuals)
}

# =================================================================================
PopModel <- function(Catch,r,K,InitPop,ExponModel)
{
 Nyears <- length(Catch)
 Pop <- rep(0,length=Nyears)
 
 Pop[1] <- InitPop
 for (Iyear in 1:(Nyears-1))
  if (ExponModel == T)
   Pop[Iyear+1] <- max(0.01,Pop[Iyear]*(1+r)-Catch[Iyear])
  else  
   Pop[Iyear+1] <- max(0.01,Pop[Iyear] + r*Pop[Iyear]*(1-Pop[Iyear]/K) - Catch[Iyear])
 return(Pop)  

}

# =================================================================================
ReadData <- function()
{
 FileName <- "C:\\courses\\FISH 558_25\\Bayes Workshop FHL\\Assignments\\Assignment 3\\Ex3a.csv"
 TheData1 <<- matrix(scan(FileName,skip=1,n=22*3,sep=','),ncol=3,byrow=T)
 FileName <- "C:\\courses\\FISH 558_25\\Bayes Workshop FHL\\Assignments\\Assignment 3\\Ex3b.csv"
 TheData2 <<- matrix(scan(FileName,skip=1,n=38*2,sep=','),ncol=2,byrow=T)
 
 Outs <- NULL
 Outs$SurveyYr <- TheData1[,1]
 Outs$SurveyEst <- TheData1[,2]
 Outs$SurveyCV <- TheData1[,3]
 Outs$CatchYr <- TheData2[,1]
 Outs$CatchVal <- TheData2[,2]
 return(Outs)
}

# =================================================================================
Ex3Soln()
