# Connor Quiroz
# Homework Assignment 4
# FISH 558 - Dr. Andre Punt

library(tidyverse) # Load in package for data cleaning
library(SimDesign) # Load in package for multivariate random sampling


##########################################
#####        DATA CLEANING           #####
##########################################

# Read in observed mosquito data
stream_full <- read.csv("Homework 4/HWK4Pars.csv", 
                                   skip = 1, 
                                   nrows = 20,
                                   row.names = 1,
                                   header = TRUE)

mosquito_obs <- as.matrix(stream_full[, 1:11])
  

# Read Variance-Covariance matrix (starts at row 24, 46 rows)
vcov_matrix <- as.matrix(read.csv("Homework 4/HWK4Pars.csv", 
                                     skip = 24, 
                                     nrows = 45,
                                     header = TRUE))

# Read MCMC starting parameters (starts at row 72)
params <- as.matrix(read.csv("Homework 4/HWK4Pars.csv", 
                          skip = 72, 
                          nrows = 1, 
                          header = TRUE))



##########################################
#####     CALCULATE LIKELIHOOD       #####
##########################################

calculate_NLL <- function(params, DataUsed) {
  
  # Assign site specific mus (Latent variables)
  mu_site <- params[6:25]
  
  # Assign site specific sigmas (Latent variables)
  alpha_site <- params[26:length(params)]
  
  # Generate intial predicted matrix 
  mosquito_pred <- matrix(0,
                          nrow = nrow(DataUsed),
                          ncol = ncol(DataUsed))
  
  # Generate intial predicted mosquito abundances
  for (site in 1:nrow(mosquito_pred)) {
    
    for (distance in 1:ncol(mosquito_pred)) {
      
      mosquito_pred[site,distance] <- mu_site[site] * exp(alpha_site[site] * ((distance - 1)  * 10))
    }
    
  }
  
  # Calculate negative log likelihood
  NLL_data <- sum(-log(dpois(DataUsed,mosquito_pred)))
  
  
  # Calculate NLL on latent variable priors (first being mu)
  mu_F <- params[1]
  mu_L <- params[2]
  sigma_mu <- params[3]
  
  likelihood_mu_site <- c()
  
  # Calculate likelihood
  for (site in 1:nrow(mosquito_pred)) {
    likelihood_mu_site <- append(likelihood_mu_site, dnorm(mu_site[site],mu_F+((mu_L-mu_F) * (site - 1)/(20-1)),sigma_mu))
  }
  
  # Convert to NLL
  NLL_mu_site <- sum(-log(likelihood_mu_site))
  
  
  # Next do alpha
  bar_alpha <- params[4]
  alpha_sigma <- params[5]
  
  
  NLL_alpha_site <- sum(-log(dnorm(alpha_site, bar_alpha, alpha_sigma)))
  
  
  sigma_mu_NLL <- -log(dnorm(sigma_mu, 1, 0.1)) # Calculate likelihood on sigma_mu hyperparameter priors
  sigma_alpha_NLL <- -log(dnorm(alpha_sigma, 0.01, 0.1)) # Calculate likelihood on sigma_alpha hyperparameter priors
  
  # Calculate total NLL
  NLL <- sum(c(NLL_data,NLL_mu_site,NLL_alpha_site, sigma_mu_NLL, sigma_alpha_NLL))
  
  if (is.na(NLL)) {
    NLL = 100000
  }
  
  return(NLL)
  
}

# Get initial likelihood
NLL <- calculate_NLL(params, DataUsed = mosquito_obs)

print(NLL) # 645.4353

##########################################
#####        A. MCMC SAMPLING        #####
##########################################

DoMCMC<-function(Xinit,DataUsed,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1, scalar)
{
  set.seed(52)
  Xcurr <- Xinit
  Fcurr <- -1 * calculate_NLL(params = Xcurr, DataUsed = mosquito_obs)
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
  Ipnt <- 0; Icnt <- 0
  accept = 0
  for (Isim in 1:Nsim)
  {
    
    # Print out simulation number to keep track of progress
    if (Isim %% 1000 == 0) {
      print(paste0("Simulation #: ", Isim, " out of ", Nsim))
    }
    
    repeat
    {
  
        varcovar <- scalar * covar
      
      Xnext <- rmvnorm(1, mean=Xcurr, sigma=varcovar)
      Fnext <- -1 *calculate_NLL(Xnext, DataUsed = mosquito_obs)
      
      if (!is.na(Fnext)) break
    }
    Rand1 <- log(runif(1,0,1))
    if (Fnext > Fcurr+Rand1)
    {Fcurr <- Fnext; Xcurr <- Xnext; print("accept"); accept <- accept + 1}   
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
  print(paste0("Acceptance rate: ", accept / Nsim))
  
  return(Outs[1:Icnt,])
}

# Run MCMC Sampling
MCMC_output <- DoMCMC(Xinit = params, DataUsed = mosquito_obs, Ndim = length(params), covar = vcov_matrix, Nsim = 505000, Nburn = 500, Nthin = 500,
            scalar = 0.05)

# Write to RDS
saveRDS(MCMC_output, file = "Homework 4/MCMC_output.rds")

# Read in RDS
MCMC_output <- readRDS("Homework 4/MCMC_output.rds")



##########################################
#####              B.                #####
##########################################

# posterior distribution for the expected count at a distance of 15km along stream 12. 
B <- MCMC_output[,22] * exp(MCMC_output[,43] * (15)) %>%
  data.frame()

B %>%
  ggplot(aes(x = `.`)) +
  geom_histogram(color = "black") +
  labs(x = "Number of mosquitos", y = "# Samples from posterior") +
  theme_light()

##########################################
#####              C.                #####
##########################################

# Get alpha posterior predictive
mu_7 <- MCMC_output[,12]

alpha_7 <- MCMC_output[,32]

posterior_predictive <- data.frame(matrix(ncol = length(colnames(mosquito_obs)),
                               nrow = length(mu_7)))
colnames(posterior_predictive) <- c("0m","10m","20m","30m","40m","50m","60m","70m","80m","90m", "100m")

for (distance in 1:ncol(posterior_predictive)) {
  
  predictions <- mu_7 * exp(alpha_7 * ((distance - 1)  * 10))
  
  # Get predictive distribution
  predictions <- rpois(n = length(predictions), lambda = predictions)
  
  # Assign predictive distribution to data frame
  posterior_predictive[, distance] <- predictions

}

observed_values <- data.frame(distance = c("0m","10m","20m","30m","40m","50m","60m","70m","80m","90m", "100m"),
                              observation = mosquito_obs[7,])

 (posterior_predictive_plot <- posterior_predictive %>%
  pivot_longer(cols = everything(),
               names_to = "distance",
               values_to = "prediction") %>%
     left_join(observed_values, by = "distance") %>%
     mutate(distance = factor(distance, levels = c("0m","10m","20m","30m","40m","50m","60m","70m","80m","90m", "100m"))) %>%
  ggplot(aes(x = prediction, fill = distance)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~distance, ncol = 1, scales = "free_y") +
     geom_vline(aes(xintercept = observation), linetype = "dashed", linewidth = 0.8) +
  theme_light() +
  theme(legend.position = "bottom") +
     labs(x = "# of mosquitos", y = "# of posterior predictive samples", fill = "Transect distance"))

ggsave("Homework 4/posterior_predictive2.jpg", posterior_predictive_plot, width = 4, height = 13, units = "in", device = "jpg")




