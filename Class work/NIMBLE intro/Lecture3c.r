library(nimble)

#setwd("E:\\")

TheCode <- nimbleCode({ 
  #time step [1] conditions (note: T for truncation)
  Pmed[1] <- 0
  #P[1] ~ T(dlnorm(Pmed[1], isigma2),0.05,1.6)
  P[1] ~ dlnorm(Pmed[1], isigma2)
  
  #time steps of model  
  for( t in 2 : N )
  {
    Pmed[t] <- log(max(P[t - 1] + (r * P[t - 1]) * (1 - P[t - 1]) - C[t - 1] / K, 0.001) )	
    #P[t] ~ T(dlnorm(Pmed[t],isigma2),0.05,1.5) 
    P[t] ~ dlnorm(Pmed[t],isigma2)
  }
  
  # Likelihood
  for( t in 1 : N )
  {
    Imed[t] <- log((q * K) * P[t])
    I[t] ~ dlnorm(Imed[t],itau2)
    
    #posterior predictions (hint, the parameterization of dlnorm is not the same as in R)
    index[t]<-(q*K*P[t])
    I.new[t]~dlnorm(log(index[t]), itau2)
  }
  
  #priors 
  r ~ T(dlnorm( -1.38, 3.845),0.01,1.2)
  isigma2 ~ dgamma(3.785,0.0102)
  itau2 ~ dgamma(1.709,0.00861)
  iq ~ T(dgamma(0.001,0.001),0.5,100)
  K ~ T(dlnorm(5.0429,3.7603),10,1000)
  
  sigma2 <- 1/isigma2
  tau2 <- 1/itau2
  q <- 1/iq
  
  #additional parameters and preditions
  MSP <-  r*K/4
  EMSP <-  r/(2*q)
  P1990 <-  P[N] + r*P[N]*(1-P[N]) - C[N]/K
  B1990 <-  P[N]
})

C<-c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3) 
I<-c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)
N<-23

data <- list(C=C,I=I,N=N)


TheConsts <- list("N" = N, "C" = C)

TheData <- list("I" = I)

#initial values
inits1 <- list(r=0.8, K=200, iq=0.5, isigma2=100, itau2=100, P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
inits2 <- list(r=0.5, K=300, iq=0.8, isigma2=200, itau2=200, P=c(0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99))
TheInits <- list(inits1,inits2)
print(str(inits1))

# Create the model
test <- nimbleModel(code = TheCode, name = "code", constants = TheConsts,
                    data = TheData, inits = inits1)

# Conduct MCMC sampling
mcmc.out <- nimbleMCMC(code = TheCode, constants = TheConsts,
                       data = TheData, inits = TheInits,
                       nchains = 2, niter = 50000,thin=10,nburnin=10000,
#                       nchains = 2, niter = 2000,
                       summary = TRUE, WAIC = TRUE,
                       samplesAsCodaMCMC=TRUE,
                       monitors = c('r','index','K','P','itau2','iq','q','I.new'))

# Diagnostics
print(names(mcmc.out))

# Tabular summary
print(mcmc.out$summary)

# Output statistics
mcmc <- rbind(mcmc.out$samples$chain1,mcmc.out$samples$chain2)
par(mfrow=c(2,2))
hist(mcmc[,which(colnames(mcmc)=="r")],xlab="r",main="")
hist(mcmc[,which(colnames(mcmc)=="K")],xlab="K",main="")
hist(mcmc[,which(colnames(mcmc)=="q")],xlab="q",main="")
Index <-which(colnames(mcmc)=="P[1]")
ResSum <- matrix(0,ncol=23,nrow=5)
for (Iyr in 1:23)
  ResSum[,Iyr] <- quantile(mcmc[,Iyr+Index-1],probs=c(0.05,0.25,0.5,0.75,0.95))  
xx <- seq(1:23)
plot(xx,ResSum[3,],xlab="Year",ylab="Depletion",type='l',lwd=3,ylim=c(0,1.3))
xx2 <- c(xx,rev(xx))
polygon(xx2,c(ResSum[1,],rev(ResSum[5,])),col="gray50")
polygon(xx2,c(ResSum[2,],rev(ResSum[4,])),col="gray95")
lines(xx,ResSum[3,],lwd=3,lty=1)

# Posterior predictive
Index <-which(colnames(mcmc)=="I.new[1]")
par(mfrow=c(2,2))
plot(xx,I,pch=16,ylim=c(0,100))
ResSum <- matrix(0,ncol=23,nrow=5)
for (Iyr in 1:23)
  ResSum[,Iyr] <- quantile(mcmc[,Iyr+Index-1],probs=c(0.05,0.25,0.5,0.75,0.95))  
lines(xx,ResSum[3,],lwd=3,lty=1) 
lines(xx,ResSum[1,],lwd=1,lty=2) 
lines(xx,ResSum[5,],lwd=1,lty=2) 


