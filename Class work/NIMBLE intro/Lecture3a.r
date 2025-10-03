library(nimble)
library(igraph)

pumpCode <- nimbleCode({
  
  for (i in 1:N) {
    theta[i] ~ dgamma(alpha,beta)
    lambda[i] <- theta[i]*t[i]
    x[i] ~ dpois(lambda[i])
  }
  
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

pumpConsts <- list(N = 10, t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4,  22))


pumpInits <- list(alpha = 1, beta = 1, theta = rep(0.1, pumpConsts$N))

pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts, data = pumpData, inits = pumpInits)

pump$getNodeNames()

pump$logProb_x
pump$lifted_d1_over_beta


pump$modelDef

pump$modelDef$BUGScode

# Plot the graphs
pump$plotGraph()

pump$getDependencies(c("alpha"))
pump$getDependencies(c("beta"))

# Check the way the samplers are applied.
pump$checkConjugacy()

# generate from the distribution for theta
set.seed(0)
simulate(pump,"theta")
print(pump$theta)

# calculate the log probabilities
pump$calculate(pump$getDependencies(c("theta")))



mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                       data = pumpData, inits = pumpInits, 
                       monitors = c("alpha","beta","theta"),
                       nchains = 2, niter = 10000,thin=1,nburnin=2000,
                       samplesAsCodaMCMC = TRUE,
                       summary = TRUE, WAIC = TRUE)

# Compile the model
Cpump <- compileNimble(pump,showCompilerOutput = TRUE)

mcmc.out <- nimbleMCMC(model=Cpump,
                       monitors = c("alpha","beta","theta"),
                       nchains = 2, niter = 10000, thin=1,nburnin=2000,
                       samplesAsCodaMCMC = TRUE,
                       summary = TRUE, WAIC = TRUE)

print(str(mcmc.out))
mcmc.out$summary


pumpConf <- configureMCMC(pump, print = TRUE)

