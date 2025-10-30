rm(list=ls())

# ===============================================================================================================================
GetPBR <- function(SurveyCV,SurveyEst,Rmax=0.04,Perc=0.5,Fr=1)
 {
  # Find the most recent index of abundance   
  Index <- max(which(SurveyCV>=0))
  Est <- SurveyEst[Index]; CV <- SurveyCV[Index]
  
  PBR = 0.5*Rmax * Est * exp(qnorm(Perc)*CV) * Fr
  #cat(Index,qnorm(Perc),Est,CV,Fr,PBR,"\n")
  return(PBR)
}

# ===============================================================================================================================

ResetBaseOMPars<- function(Species=1)
{
 # This is scenario 0
  
 if (Species==1) Rmax <- 0.04  
 if (Species==2) Rmax <- 0.12  

 OMPars <- NULL
 OMPars$Rmax <- Rmax
 OMPars$bc <- 1
 OMPars$bn <- 1
 OMPars$SigmaC <- 0.3
 OMPars$SigmaN <- 0.2
 OMPars$SurFreq <- 4
 OMPars$SigmaNmult <- 1
 OMPars$RmaxBias <- 1
 OMPars$CVmult <- 1
 return(OMPars)
}

# ===============================================================================================================================

GetScenario <- function(OMPars,Scenario)
 {
  # The specifications for the scenarios
  if (Scenario == 1) OMPars$SigmaC <- 1.2
  if (Scenario == 2) OMPars$SigmaN <- 0.8
  if (Scenario == 3) OMPars$SurFreq <- 8
  if (Scenario == 4) OMPars$bc <- 1.2
  if (Scenario == 5) OMPars$bn <- 1.2
  if (Scenario == 6) OMPars$RmaxBias <- 0.5
  if (Scenario == 7) OMPars$SigmaNmult <- 0.5
  return(OMPars)
 }

# ===============================================================================================================================
ProjectOM <- function(Species,Scenario,PBR.Perc=0.5,PBR.Fr=1,Nyear=100,Nsim=100)
{
 # Select the Rmax for use in the PBR formulae  
 if (Species==1) PBR.Rmax <- 0.04  
 if (Species==2) PBR.Rmax <- 0.12  
 
 # Set the value of K (does not matter)
 K <- 1000
 
 # set the OM parameters
 OMPars <- ResetBaseOMPars(Species)
 OMPars <- GetScenario(OMPars,Scenario)
   
 # Extract parameters
 Rmax <- OMPars$Rmax
 SigmaC <- sqrt(log(OMPars$SigmaC^2+1))
 SigmaN <- sqrt(log(OMPars$SigmaN^2+1))
 bc <- OMPars$bc
 bn <- OMPars$bn
 SurFreq <- OMPars$SurFreq
 SigmaNmult <- OMPars$SigmaNmult
 RmaxBias <- OMPars$RmaxBias
 Rmax <- Rmax * RmaxBias
 
 # Output vector
 Depletion <- rep(0,Nsim)
 
 set.seed(19102)
 for (Isim in 1:Nsim)
  {
   # Set up the vectors and the initial N
   N <- rep(0,101)
   N[1] <- 0.3*K
   SurveyEst <- rep(-1,Nyear+1)
   SurveyCV <- rep(-1,Nyear+1)
   
   # Generate a first abundance estimate
   SurveyCV[1] <- SigmaN*SigmaNmult
   SurveyEst[1] <- bn*N[1]*exp(rnorm(1,0,SigmaN))
 
   # Project forward
   for (Year in 1:Nyear)
    {
     # Extract a PBR value
     PBR <- GetPBR(SurveyCV,SurveyEst,RmaxBias*PBR.Rmax,PBR.Perc,PBR.Fr)
     # Convert to catch (and account for bias)
     Catch <- bc*PBR*exp(rnorm(1,0,SigmaC)-SigmaC^2/2)
     
     # Update numbers
     N[Year+1] <- N[Year] + Rmax*N[Year]*(1-N[Year]/K) - Catch
     if (N[Year+1]<1) N[Year+1] <- 1
     
     # Check if a survey estimate is to be generated
     if ((Year %% SurFreq==0)) 
      {
       SurveyCV[Year+1] <- SigmaN*SigmaNmult
       SurveyEst[Year+1] <- bn*N[Year+1]*exp(rnorm(1,0,SigmaN))
      }
    }
   Depletion[Isim] <- N[Nyear+1]/K
   
  }
 # Calculate the probability
 Prob <- sum(Depletion>0.5)/Nsim
 MeanD <- mean(Depletion)
 #cat(PBR.Perc,PBR.Fr,Prob,MeanD,"\n")
 
 # Return stuff
 Return.V <- NULL
 Return.V$Prob <- Prob
 Return.V$MeanD <- MeanD
 return(Return.V)
}

# ===============================================================================================================================

uniroot1 <- function(x,PBR.Fr,Species,Scenario)
 {
  # Probability when changing "y" (the percentile of recent abundance)
  Prob <- ProjectOM(Species,Scenario,PBR.Perc=x,PBR.Fr=PBR.Fr)$Prob - 0.8
  return(Prob)
 }

# ===============================================================================================================================

uniroot2 <- function(x,PBR.Perc,Species,Scenario)
 {
  # Probability when changing Fr
  Prob <- ProjectOM(Species,Scenario,PBR.Perc=PBR.Perc,PBR.Fr=x)$Prob - 0.8
  return(Prob)
 }

# ===============================================================================================================================

# Task 1: Find Fr
# ===============
cat("\nDoing Task 1\n")
# First species
# Uniroot code here
PBR.perc1 <- uni1$root
# Uniroot code herePBR.perc2 <- uni2$root

# Find the lower value
PBR.perc <- min(PBR.perc1,PBR.perc2)
cat("Percentile for task 1 is ",PBR.perc,"\n")

# ----------------------------------------------------------------------------------------------------------------------------

# Task 2, Part 1: Run through all scenarios
# =========================================
cat("\nDoing Part 1 of Task 2\n")
for (Species in 1:2)
 for (Scenario in 0:7)
  {
   cat("Task 2, Part 1",Species,Scenario,SimOut$Prob,SimOut$MeanD,"\n")
  } 

# Task 2, Part 2: Solve for Fr
# =============================
cat("\nDoing Part 2 of Task 2\n")
# More uni-root use
PBR.Fr <- uni1$root
cat("Fr is",PBR.Fr,"\n")

# Task 2, Part 3: Run through all scenarios again
# ===============================================
cat("\nDoing Part 3 of Task 2\n")
for (Species in 1:2)
  for (Scenario in 0:7)
  {
    cat("Task 2, Part 3",Species,Scenario,SimOut$Prob,SimOut$MeanD,"\n")
  } 






