# Connor Quiroz
# FISH 558 - Dr. Andre Punt

  library(tidyverse)
  
  ####################################################################
  
  # Problem 2
  
  # Read in Data
  data <- read.table("Homework 1/HWK1_data.txt", header = TRUE)
  data$Catch
  
  pop.model <- function(r, K) {
    
    # Create an array of biomass values for each year from 1922 to 2004
    # K is the same thing as the starting biomass in 1922
    N <- rep(NA, dim(data)[1])
    N[1] <- K # Set starting value to K
    z <- 2.39 # Predetermined Pella-Tomlinson parameter
    
    # Calcualte subsequent years
    for (t in 2:length(N)) {
      N[t] <- N[t-1] + r * N[t-1] * (1 - (N[t-1]/K)^z) - data$Catch[t-1]
    }
    
    # Return biomass across years
    return(N)
    
    
  }
  
  likelihood <- function(param.vec) {
  
    # Get parameter values from parameter vector
    r <- param.vec[1]
    K <- param.vec[2]
    tao <- abs(param.vec[3])
    
    # Get observation error standard deviation
    sigma_t <- sqrt(data$CV + tao)
    
    # Get the number of indivudals
    N <- pop.model(r = r, K = K)
    
    
    # TESTING (COMMENT BEFORE RUNNING LIKELIHOOD FUNCTION)
    # test_vector <- rep(NA, length(data$Abundance))
    # indexes <- which(!is.na(data$Abundance))
    # test_vector[indexes] <-  900000 # Some really large number to test in likelihood function
    
    # Calculate likelihood
    Likelihood <- prod((1/sigma_t * exp(-((log(N)-log(data$Abundance))^2)/(2*sigma_t^2))), na.rm = TRUE)
    # Likelihood <- prod((1/sigma_t * exp(-((log(test_vector)-log(data$Abundance))^2)/(2*sigma_t^2))), na.rm = TRUE)
    return(-log(Likelihood))
    
  }
  
  # Testing some stuff with likelihood
  
  
  # Set starting values
  param.vec <- c(0.05, 1000, 0.25)
  
  # Run optim to get MLE
  optim(par = param.vec, fn = likelihood, method = "Nelder-Mead")
  
  # Trying different starting values
  param.vec <- c(0.05, 1000, 0.25)
  result <- optim(par = param.vec, fn = likelihood, method = "Nelder-Mead")
  param.vec <- c(0.45, 2000, 0.05)
  result <- optim(par = param.vec, fn = likelihood, method = "Nelder-Mead")
  param.vec <- c(0.65, 5000, 0.5)
  result <- optim(par = param.vec, fn = likelihood, method = "Nelder-Mead")
  
  # These all yielded different results, and an r of 0 created null. I need to try across a grid of parameters 
  # to try numerous combinations of starting values across a parameter space I believe to be around the MLE
  # (prior beliefs gathered from scatterplot of observed catches) to obtain the MLE. By searching across ~76,000
  # different starting combinations, the algorithm will likely cross a suitable starting value that can fit the
  # three data points better than anyone else. I will test for this by only keeping the current lowest MLE with
  # the respective parameters as the R moves through the loop. I had to run the iterations multiple times because I
  # realized the MLE was at the maximum K value, meaning it's likely past that interval. I originally started with
  # these intervals:
  
  # for (r in seq(0.2, 1, by = 0.01)) {
  #   for (K in seq(1000, 1500, by = 50)) {
  #     for (tao in seq(0, 2.5, by = 0.05))
  
  
  # I'm going to switch narrower intervals that are closer to what tao and r parameter are at their MLE since I figured out
  # From the first iteration they were within their intervals, so I think their values would be relatively the same. Because
  # Shifting only our K value parameter would just be translating by K but with the same growth and error coefficient to
  # produce the same line, but shifted up due to K being higher in the hypothesized new MLE. This and likelihood profiles
  # allowed me to ensure a global minimum.
  
  # It is a generalized method across variables, but looking for any of the variable MLE estimates being the same as the
  # bounds I set in the for loop. If the MLE is at a bound, it is likely (I would argue almost certainly) that the MLE
  # is beyond that value. However, it seems that through a Bayesian analysis, we could traverse the minimums easier through
  # different sampling proceses like SIR or MCMC to traverse that parameter space differently to compute more holistic
  # probabilities.
  
  
  # I DONT THINK WE CAN FIND THE MLE. I think you tricked us. The "MLE" eventually converges around parameters when past a certain threshold, will
  # cause the model itself to diverge and not produce any estimate because the data points fall below the 0 population level.
  
  # I can estimate it as well as I can, but I will likely never find the *true exact* global minimum. Practically being close
  # on a parameter within 0.01 would be great even in a scientific context, but it mathematically isn't the same as
  # the global minimum. We are only numerically estimating it. Is this also why we are numerically estimating the posterior?
  # Because if we knew a closed form solution the shape of the parameter space, then we could analytically derive it,
  # but we're not, hence numerical estimation! As is, this is a four dimensional parameter space with the three parameters and the likelihood as the fourth variable.
  
  # I wasn't expecting it, but this homework assignment taught me about the difference between numerical and analytically solutions for the MLE,
  # with future thought settling on the applications of applying that in a Bayesian framework.
  
  # That being said, I can't actually find the true minimum, but I've gone ahead and minimized the model as much as I'm
  # patiently possible to test to each parameter to a specifc decimal place. Exceeding past either of those values (i.e.,
  # going lower which is what yields a higher MLE) will cause model divergence
  
  
  # This is basically a four dimensional negative log likelihood surface we are trying to traverse. There's no way we can
  # locate the global minimum with the current method.
  # The lowest likelihood value given our setup is 0 and that's
  # another reason we can't ever find the true minimum.
  # Because all the good low points have a likelihood of 0. I think we can do given our likelihood formula is 0. Not sure if that's how the likelihood
  # formula actually bounds it at 0, but I'm not sure that doesn't captures the analytical bounds of the MLE value.
  
  # THIS IS WHY IT IS SUCH A BIG ISSUE Because wherever you traverse the likelihood terrain with a model setup like this, different scientists
  # Will converge to different theories and they technicaly are both correct if their likelihoods are both "0" due to software + numerical
  # limitations for actually finding the MLE. I think the primary issue with this approach lies in the numerical limitation because that is how the likelihood is set up
  # and would be the similar likelihood we use to calculate likelihood in Bayesian, which is why finding the true value would be a different approach. But it would still
  # Give us a better idea of what we don't exactly know. 
  
  
  # Grid search / SIR / MCMC would also work much better because the sampling approach is also different then Nelder-Mead, helping us traverse uncertain parameter landscapes
  # much more efficiently. With that logic, wouldn't Bayesian be another form of loss function optimization (e.g., minimizing function instead of likelihood function)
  # as loss would be a different way to represent likelihood in a broader machine learning problem? Since both frequentist / Bayesian methods are both trying to find
  # the global minimum in one way or another.
  
  param.vec <- c(0.1, 600, 255)
  result <- optim(par = param.vec, fn = likelihood, method = "Nelder-Mead")
  
  MLE <- Inf   # Start with Inf, since we are minimizing
  params <- NULL
  iteration <- 0
  
  for (r in seq(0.07, 0.09, by = 0.01)) {
    for (K in seq(900, 1000, by = 50)) {
      for (tao in seq(0.05, 0.6, by = 0.05)) {
        
        param.vec <- c(r, K, tao)
        
        # Insert a catch so that if an error occurs from bad starting params (e.g., r = 0 --> error), 
        # it will continue to the next loop iteration and repeat until the optim() works properly,
        # then compare against the current lowest NLL.
        result <- tryCatch(
          {
            res <- optim(par = param.vec, fn = likelihood, method = "Nelder-Mead")
            if (is.finite(res$value)) res else NULL
          },
          error = function(e) {
            message("Skipping bad initial parameters: ", paste(param.vec, collapse = ", "))
            NULL
          }
        )
        
        # If result is valid, and better than current best, update
        if (!is.null(result)) {
          if (result$value < MLE) {
            MLE <- result$value
            params <- param.vec
            print(paste0("New best MLE found of ", signif(result$value, 4), "! Parameters used: ", param.vec))
          }
        }
        
        iteration <- iteration + 1
        if (iteration %% 500 == 0) {
          print(paste0("Iteration # ", iteration))
        }
      }
    }
  }

# For putting into plot for graphing observed vs expected data. Tao not needed because it was only used for nll,
# not producing an actual model abundance point.
r_MLE <- params[1]

K_MLE <- params[2]

tao_MLE <- params[3]

# Plot model to observed data
B_MLE <- pop.model(r = r_MLE, K = K_MLE) %>% # MLE predictions
  tibble(Abundance_MLE = .)

# Join MLE estimates to observed data
(MLE_plot <- data %>%
  bind_cols(B_MLE) %>%
  pivot_longer(cols = c("Abundance", "Abundance_MLE")) %>%
  mutate(name = case_when(
    name == "Abundance" ~ "Observed abundance",
    name == "Abundance_MLE" ~ "MLE abundance prediction"
  )) %>%
  ggplot(aes(x = Year, y = value, color = name)) +
  geom_line() +
  geom_point() +
  theme_light() +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Abundance (# whales)",
       color = "Data type") +
  scale_color_manual(values = c("blue", "red")))

ggsave("MLE_plot.jpg", plot = MLE_plot, device = "png", units = "in", width = 6.5, height = 4)



####################################################################

# Problem 3

# Redefine likelihood function so that the parameter vector only includes r, we are optimizing for R and will change
# K to get the MLE at that given K for which will give us the resulting r value too, producing a likelihood profile
likelihood <- function(param.vec, r = 0.078) {
  
  # Get parameter values from parameter vector
  # r <- param.vec[1] <----- No more need for r for likelihood profile <-----
  K <- param.vec[1] 
  tao <- abs(param.vec[2])
  
  # Get observation error standard deviation
  sigma_t <- sqrt(data$CV + tao)
  
  # Get the number of indivudals
  N <- pop.model(r = r, K = K)
  
  # Calculate likelihood
  Likelihood <- sum((1/sigma_t * exp(-((log(N)-log(data$Abundance))^2)/(2*sigma_t^2))), na.rm = TRUE)
  return(-log(Likelihood))
  
}


# make r.seq once
r.seq <- seq(0.07, 0.091, 0.001)
param.vec <- c(K_MLE, tao_MLE)

profile <- matrix(NA, nrow = length(r.seq), ncol = 2)
colnames(profile) <- c("r", "Likelihood")
profile[,1] <- r.seq   # first col = fixed r values

for (i in seq_along(r.seq)) {
  r <- r.seq[i]
  output <- optim(par = param.vec, fn = likelihood, r = r)
  
  # use i instead of which()
  profile[i, 2] <- output$value
}




# Visualize profile
(likelihood_profile <- data.frame(profile) %>%
    mutate(MSYR = 0.705 * r) %>%
  ggplot((aes(x = r, y = Likelihood))) +
  geom_point() +
    geom_line() +
    theme_light() +
    labs(y = "Negative log likelihood"))

# Save likelihood profile
ggsave("likelihood_profile.jpg", plot = likelihood_profile, device = "png", units = "in", width = 6.5, height = 4)
