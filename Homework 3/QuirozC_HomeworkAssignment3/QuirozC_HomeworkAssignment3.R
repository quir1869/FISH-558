# Connor Quiroz
# FISH 558 - Dr. Andre Punt
# Homework Assignment 3 - (Bowhead assessment and decision analysis) 

library(tidyverse)

# Load in catch data
data <- read.csv("Homework 3/Hwk3.csv", header = FALSE) %>%
  rename(year = "V1",
         catches = "V2") %>%
    drop_na()

# Load in abundance data
abundance_obs <- read.csv("Homework 3/Hwk3B.csv", header = FALSE) %>%
  rename(year = "V1",
         estimate = "V2",
         CV = "V3")

# Extract data
Nobs <- abundance_obs$estimate
CV <- abundance_obs$CV

# Set initial values
Nyear <- length(1848:2002)
z = 2.39

# ----------------------------------------------------#
#                     Problem A2                      #
# ----------------------------------------------------#

popModel <- function(s0 = 0.9,s1 = 0.95, k1_plus = 15000, fmax = 0.29, z = 2.39, max_age = 13, catches = data$catches) {
  
  
  # Pre exploitation values
  n_years <- length(catches)
  pre_exp <- rep(NA, 14)
  f0 <- (1 - s1) / s0
  N_tilda <- k1_plus
  N0 <- N_tilda * f0
  
  for (age in 1:14)
  {
    
    if (age == 1) {
      
      pre_exp[age] <- N0 # age 0
      
    } 
    
    if (age %in% 2:13) {
      
      pre_exp[age] <- N0 * s0 * (s1^(age - 2)) # ages 1–12
      
    } 
    
    if (age == 14) {
      
      pre_exp[age] <- (N0 * s0 * s1^(max_age - 1)) / (1 - s1) # plus group
      
    }
  }
  
  # Now to project forward
  a_s <- matrix(NA, nrow = n_years, ncol = 14)
  a_s[1, ] <- pre_exp
  
  for (year in 2:n_years) {
    
    N_plus <- sum(a_s[year - 1, 2:14]) # total age 1+ animals 
    
    # Age 1
    a_s[year, 2] <- a_s[year - 1, 1] * s0 
    
    # Ages 2–12
    for (age in 3:13) {
      
      c_ya <- (catches[year - 1] * a_s[year - 1, age - 1]) / N_plus
      a_s[year, age] <- (a_s[year - 1, age - 1] - c_ya) * s1
      
    }
    
    # Plus group
    c_ya_13 <- (catches[year - 1] * a_s[year - 1, 13]) / N_plus
    c_ya_14 <- (catches[year - 1] * a_s[year - 1, 14]) / N_plus
    a_s[year, 14] <- (a_s[year - 1, 13] - c_ya_13) * s1 + (a_s[year - 1, 14] - c_ya_14) * s1 
    
    # Age 0
    N_1plus <- sum(a_s[year, 2:14])
    N_tilda <- a_s[year, 14]
    a_s[year, 1] <- max(N_tilda * (f0 + (fmax - f0) * (1 - (N_1plus / k1_plus)^z)), 0)
  } 
  
  return(a_s)
}




a_s <- popModel(s0 = 0.9, s1 = 0.95, k1_plus = 15000, fmax = 0.29, z = 2.39, catch = AllCatch)

# Clean up dataset
colnames(a_s) <- c("Age 0", "Age 1", "Age 2", "Age 3", "Age 4", "Age 5", "Age 6", "Age 7", "Age 8", "Age 9", "Age 10", "Age 11", "Age 12", "Age 13")
rownames(a_s) <- c(data$year)

plot(rowSums(a_s[,2:14]))


calcLikelihood <- function(s0 = 0.9, s1 = 0.95,
                           k1_plus = 15000, fmax = 0.29, z = 2.39, catches = data$catches) {
  
  # Run population model
  a_s <- popModel(s0 = s0, s1 = s1, k1_plus = k1_plus, fmax = fmax, z = z, catches = data$catches)
  
  # Total abundance predictions
  Npred <- rowSums(a_s[, 2:14])
  
  # Get observed data years
  years_to_check <- abundance_obs$year
  idx <- which(1848:2002 %in% years_to_check)
  
  # Check for invalid predictions
  if (any(is.na(Npred[idx])) || any(Npred[idx] <= 0)) {
    return(10000)  
  }
  
  # Compute log-scale difference safely
  log_diff <- log(pmax(Npred[idx], 1e-8)) - log(pmax(Nobs, 1e-8))
  
  # Compute NLL
  NLL <- sum((log_diff)^2 / (2 * CV^2))
  
  # Protect against invalid numeric output
  if (!is.finite(NLL) || is.na(NLL)) {
    NLL <- 1e7
  }
  
  return(NLL)
}

# Numbers at each age
popModel(s0 = 0.9, s1 = 0.95,
         k1_plus = 15000, fmax = 0.29, z = 2.39)[length(1848:2002),]
calcLikelihood() # Likelihood: 96.14

# ----------------------------------------------------#
#                     Problem A3                      #
# ----------------------------------------------------#

DoSir_predata <- function(Nout=1000){
  
  set.seed(52)
  Nyear <- length(data$catches)
  Cumu <- 0
  Vals <- matrix(0,ncol=4,nrow=Nout) # save parameter vectors 
  PopOut <- matrix(0,ncol=Nyear,nrow=Nout) # save population vectors for each run 
  Ntest <- 0 # counting variable 
  Ndone <- 0
  
  while (Cumu < Nout) {
    # sample from prior distributions
    s0 <- runif(1, 0.8, 1.0)
    s1 <- runif(1, 0.9, 1.0)
    k1_plus <- runif(1, 10000, 20000)
    fmax <- runif(1, 0.25, 0.33)
    
    Ntest <- Ntest + 1
    
    # Get population estimates
    N <- popModel(s0 = s0, s1 = s1, k1_plus = k1_plus, fmax = fmax, catches = data$catches)
    
    # Check if population goes extinct 
    if (any(is.na(N[,2:14])) || any(N[,2:14] < 0)) {
      
      next
      
    } else {
      Cumu <- Cumu + 1 # contribute to "likelihood" 
      Ndone <- Ndone + 1
      Vals[Ndone,1] <- s0
      Vals[Ndone,2] <- s1
      Vals[Ndone,3] <- k1_plus
      Vals[Ndone,4] <- fmax
      PopOut[Ndone,] <- rowSums(N[,2:14])
      
      if (Cumu %% 50 == 0) {
        print(paste0("Iteration: ", Cumu))
      }
      
      
      
    }
  } # while Cumu < Nout 
  
  par(mfrow=c(2,2))
  hist(Vals[,1],main="",xlab="S0")
  hist(Vals[,2],main="",xlab="S1")
  hist(Vals[,3],main="",xlab="K")
  hist(Vals[,4],main="",xlab="fmax")
  print(Ndone)
  print(Cumu)
  print(Ntest)
  
  
  return(list(Vals, PopOut))
}

DoSir_predata <- DoSir_predata(Nout=1000)

# Save data
write_rds(DoSir_predata, "Homework 3/DoSir_predata.rds")

# Results: Our values that are somewhat safe combinations of values that don't cause the population to go extinct was a uniform distirbution of S0,
# A slightly skewed right distribution of S1, a slightly skewed left distribution of K, and a uniform distribution of fmax.

# ----------------------------------------------------#
#                     Problem A4                      #
# ----------------------------------------------------#

DoSir <- function(Nout=1000)
{
  Threshold <- exp(-1*2.22)
  Nyear <- length(data$catches)
  set.seed(52)
  Cumu <- 0
  Vals <- matrix(0,ncol=4,nrow=Nout) # save parameter vectors 
  PopOut <- matrix(0,ncol=Nyear,nrow=Nout) # save population vectors for each run 
  Ntest <- 0 # counting variable 
  Ndone <- 0
  
  while (Ndone < Nout)
  {
    # Sample priors
    s0 <- runif(1, 0.8, 1.0)
    s1 <- runif(1, 0.9, 1.0)
    k1_plus <- runif(1, 10000, 20000)
    fmax <- runif(1, 0.25, 0.33)
    
    
    N <- popModel(s0 = s0, s1 = s1, k1_plus = k1_plus, fmax = fmax, z = z, catches = data$catches)
    # Calculate likelihood
    TheLike <- exp(-1 *(calcLikelihood(s0 = s0, s1 = s1, k1_plus = k1_plus, fmax = fmax, z = 2.39, catches = data$catches)))
    
    # Information needed for SIR algorithm
    Cumu <- Cumu + TheLike
    
    if(any(is.na(N[,2:dim(N)[2]])) || any(N[,2:dim(N)[2]] < 0)){ # if any age group in population output is NA or negative
      Ntest <- Ntest + 1 # move on 
    } else {
      while (Cumu > Threshold & Ndone < Nout)
      {
        cat("Saving",Ndone,"\n")
        Cumu <- Cumu - Threshold
        Ndone <- Ndone + 1
        Vals[Ndone,1] <- s0
        Vals[Ndone,2] <- s1
        Vals[Ndone,3] <- k1_plus
        Vals[Ndone,4] <- fmax
        PopOut[Ndone, ] <- rowSums(N[, 2:14])
      }
    }
    
  }
  
  return(list(Vals, PopOut))

}

# Run SIR algorithm
SirOutput <- DoSir(Nout = 200)

write_rds(SirOutput, "Homework 3/QuirozC_Problem4SIROutput.rds")

# ----------------------------------------------------#
#                     Problem A5                      #
# ----------------------------------------------------#


# Plot posteriors
hist(SirOutput[[1]][,1]) # s0
hist(SirOutput[[1]][,2]) # s1
hist(SirOutput[[1]][,3]) # K1+
hist(SirOutput[[1]][,4]) # fmax


# Look at predata fit
N_plus_pre <- DoSir_predata[[2]]
N_plus_pre_quants <- apply(N_plus_pre, 2 , quantile , probs = c(0.05,0.50,0.95))

# Convert to dataframe
N_plus_pre_quants <- t(N_plus_pre_quants)

# Clean rownames
rownames(N_plus_pre_quants) <- NULL

# Look at post data fit
data.frame(N_plus_pre_quants) %>%
  mutate(year = data$year) %>%
  left_join(abundance_obs, by =  "year") %>%
  ggplot() +
  geom_line(aes(x = year, y = X50.)) +
  geom_ribbon(aes(x = year, ymin = X5., ymax = X95.), fill = "black", alpha = 0.2) +
  geom_point(aes(x = year, y = estimate), color = "red") +
  theme_light() +
  labs(x = "Year", y = "Mature population size\n(Ages 2-13)")

#######

# Get quantiles for N+ values for each year
N_plus_post <- SirOutput[[2]]
N_plus_post_quants <- apply(N_plus_post, 2 , quantile , probs = c(0.05,0.50,0.95))

# Convert to dataframe
N_plus_post_quants <- t(N_plus_post_quants)

# Clean rownames
rownames(N_plus_post_quants) <- NULL

# Look at post data fit
data.frame(N_plus_post_quants) %>%
  mutate(year = data$year) %>%
  left_join(abundance_obs, by =  "year") %>%
  ggplot() +
  geom_line(aes(x = year, y = X50.)) +
  geom_ribbon(aes(x = year, ymin = X5., ymax = X95.), fill = "black", alpha = 0.2) +
  geom_point(aes(x = year, y = estimate), color = "red") +
  theme_light() +
  labs(x = "Year", y = "Mature population size\n(Ages 2-13)")


# ----------------------------------------------------#
#                     PART B                          #
# ----------------------------------------------------#
DoSir <- function(Nout, catch_data, harvest)
{
  set.seed(52)
  Catch <- catch_data$catches
  CatchYr <- catch_data$year
  
  # add new years / catch
  Nyear <- length(1848:2023)
  new_catch <- rep(harvest,length(2003:2023))
  AllCatch <- c(Catch, new_catch)
  
  
  Threshold <- exp(-1*2.22) # set Threshold at the supposed minimum of NLL 
  Cumu <- 0
  Ndone <- 0
  Vals <- matrix(0,ncol=4,nrow=Nout) # save parameter vectors
  Nplus_2003 <- rep(0, length = Nout)
  N_tilda_2003 <- rep(0, length = Nout)
  N_tilda_2023 <- rep(0, length = Nout)
  
  while (Ndone < Nout)
  {
    # sample  priors
    s0 <- runif(1, 0.8, 1.0)
    s1 <- runif(1, 0.9, 1.0)
    k1_plus <- runif(1, 10000, 20000)
    fmax <- runif(1, 0.25, 0.33)
    
    
    # Get population values per SIR sample
    N <- popModel(s0 = s0, s1 = s1, k1_plus = k1_plus, fmax = fmax, z = z, catches = AllCatch)
    # N <- PopModel(s0, s1, k1_plus, fmax, Catch = AllCatch)
    
    # deal with extinct pops 
    if(any(is.na(N[,2:dim(N)[2]])) || any(N[,2:dim(N)[2]] < 0)){ # if any age group in population output is NA or negative
      
      next
      
    }
    else{
      
      TheLike <- exp(-1 * calcLikelihood(s0 = s0, s1 = s1, k1_plus = k1_plus, fmax = fmax, z = z))
      Cumu <- Cumu + TheLike # cumulative likelihood
      
      while (Cumu > Threshold & Ndone < Nout) # while the Cumu likelihood is above the threshold
      {
        cat("Saving",Ndone,"\n")  # save this parameter vector 
        Cumu <- Cumu - Threshold # subtract the threshold from the cumulative likelihood (helps to not sample the same vector a bunch)
        Ndone <- Ndone + 1
        # save parameter vector
        Vals[Ndone,1] <- s0
        Vals[Ndone,2] <- s1
        Vals[Ndone,3] <- k1_plus
        Vals[Ndone,4] <- fmax
        Nplus_2003[Ndone] <- sum(N[(length(catch_data$catches)+1),2:14])
        N_tilda_2003[Ndone] <- N[(length(catch_data$catches)+1),14]
        N_tilda_2023[Ndone] <- N[length(1848:2023),14]
      } 
    }
  } 
  
  return(list(Vals, Nplus_2003, N_tilda_2003, N_tilda_2023))
}

# Set up harvest values + decision table matrix
harvest <- c(67,134,201)
d_t <- matrix(data = NA, nrow = 3, ncol = 4)
d_t[,1] <- harvest
# Row 1: annual harvest 2003-2022 = 67 
harvest_67 <- DoSir(Nout = 200, data, harvest = harvest[1])

# Write to / read in RDS
write_rds(harvest_67, "Homework 3/harvest_67.rds")
harvest_67 <- readRDS("Homework 3/harvest_67.rds")

# Convert list into data frame
harv_67_df <- data.frame("popn_size_2003_plusone" = harvest_67[[2]], "mature_2003" = harvest_67[[3]], "mature_2023" = harvest_67[[4]])

# Complete row 1 of decision table
d_t[1,2] <- harv_67_df %>%
  filter(popn_size_2003_plusone > 7000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

d_t[1,3] <- harv_67_df %>%
  filter(popn_size_2003_plusone >= 7000,
         popn_size_2003_plusone <= 8000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  View()
  distinct(prob) %>%
  pull(prob)

d_t[1,4] <- harv_67_df %>%
  filter(popn_size_2003_plusone > 8000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

harvest_134 <- DoSir(Nout = 200, data, harvest = harvest[2])

# Write to / read in RDS
write_rds(harvest_134, "Homework 3/harvest_134.rds")
harvest_134 <- readRDS("Homework 3/harvest_134.rds")


harv_134_df <- data.frame("popn_size_2003_plusone" = harvest_134[[2]], "mature_2003" = harvest_134[[3]], "mature_2023" = harvest_134[[4]])

# Complete row 2 of decision table
d_t[2,2] <- harv_134_df %>%
  filter(popn_size_2003_plusone > 7000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

d_t[2,3] <- harv_134_df %>%
  filter(popn_size_2003_plusone >= 7000,
         popn_size_2003_plusone <= 8000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

d_t[2,4] <- harv_134_df %>%
  filter(popn_size_2003_plusone > 8000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)



harvest_201 <- DoSir(Nout = 200, data, harvest = harvest[3])

# Write to / read in RDS
write_rds(harvest_201, "Homework 3/harvest_201.rds")
harvest_201 <- readRDS("Homework 3/harvest_201.rds")

harv_201_df <- data.frame("popn_size_2003_plusone" = harvest_201[[2]], "mature_2003" = harvest_201[[3]], "mature_2023" = harvest_201[[4]])

# Complete row 3 of decision table
d_t[3,2] <- harv_201_df %>%
  filter(popn_size_2003_plusone > 7000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

d_t[3,3] <- harv_201_df %>%
  filter(popn_size_2003_plusone >= 7000,
         popn_size_2003_plusone <= 8000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

d_t[3,4] <- harv_201_df %>%
  filter(popn_size_2003_plusone > 8000) %>%
  mutate(condition = mature_2023 > mature_2003,
         prob = mean(condition)) %>%
  distinct(prob) %>%
  pull(prob)

# Save decision table
write_rds(d_t, "Homework 3/decision_table.rds")

# This decision table looks at the relative change in mature individuals from 2003-2023 given the three different states of nature of being below 7000
# individuals, between 7k and 8k, and above 8k. I think this can be useful if we are aiming to have higher individuals than before. But What
# if this doesn't tell us how reliable their efficiency at producing babies is? What if having more mature individuals than 2003 isn't enough? What if
# it needs to be above a certain threshold since their rate at producing babies (e.g., Beverton Holt Stock recruitment relationship) is low
# (i.e., steepness is low)? I trying to identify a threshold for the number of mature individuals there needs to be in a population in order to sustain
# healthy numbers would be better to look at, since that tells you whether the fishery under X conditions will perform well, or that X scenario would put
# the stock at risk due to too few mature adults, resulting in the fishery potentially needing to close due to overfishing. Instead of N~2023 > N~2003, it
# could be N~2003 > some biologically determined threshold, and create a separate decision table for N~2023, to look at the projections based on how far
# a fishery would look out. There are more decisions here since with each year you would need a decision table, but if we are interested in seeing if
# a fishery needs to worry about closing, then I think comparing relative change in mature animals can't give us a clear picture since we arent checking
# the probability that the mature animals would fall below a certain threshold.


