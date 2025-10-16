# Connor Quiroz
# FISH 558 - Dr. Andre Punt
# Homework Assignment 3 - (Bowhead assessment and decision analysis) 

library(tidyverse)

data <- bind_rows(data.frame(year = NA,
                             catches = NA),
                  read.csv("Homework 3/Hwk3.csv", header = FALSE) %>%
  rename(year = "V1",
         catches = "V2") %>%
    drop_na()
  )
  

z <- 2.39

pop.model <- function(s0 = 0.9, s1 = 0.95,
                      k1_plus = 15000, fmax = 0.29) {
  
  initial_values <- rep(300, 14)
  ages <- rep(NA, 14)
  years <- nrow(data) # +1 to account for initialized values
  
  # Set up age structure matrix
  age_structure <- matrix(nrow = years, ncol = 14)
  
  # Set up initial values for age structure
  age_structure[1,] <- initial_values
  rownames(age_structure) <- data$year
  
  
  # First want to iterate through ages, then through years
  for (year in 2:years) {
    
    for (age in 0:13) {
      
      # Age 0
      if (age == 0) {
        # N_tilda <- age_structure[year-1,age+1] * (1 / s1)
        # N_tilda <- sum(age_structure[year-1,])
        N_tilda <- sum(age_structure[year-1,14])
        N_plus <- sum(age_structure[year-1,])
        # f0 <- (1 - s1) / s0
        f0 <- (1-s1) / (s0 * s1^12)
        age_structure[year,age+1] <- N_tilda * (f0 + (fmax - f0) * (1-(N_plus / k1_plus)^z))
      }
      
      # Age 1
      if (age == 1) {
        age_structure[year,age+1] <- age_structure[year,1] * s1
      }
      
      # ages 2 to 12
      if (age %in% 2:12) {
        c_ya <- (data$catches[year] * age_structure[year-1,age-1]) / N_plus
        age_structure[year,age+1] <- (age_structure[year,age-1] - c_ya) * s1
      }
      
      # age 13
      if (age == 13) {
        c_ya_x <- (data$catches[year] * age_structure[year-1,age]) / N_plus
        age_structure[year, 13+1] <- (age_structure[year,13] - c_ya) * s1 + (age_structure[year-1,13] - c_ya_x) * s1
      }
      
      
    }
  }
  
  return(age_structure)
  
}


age_structure <- pop.model()


nrow(age_structure)
