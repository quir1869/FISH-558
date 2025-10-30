calculate_risk <- function(data) {
  # Define all possible adaptive capacity variables and their weights
  ac_variables <- c("gdp", "gdp_trade", "sanitation", "supermarkets", 
                    "life_expectancy", "prop_labor", "hci", "gov_effectiveness",
                    "fsc", "rol")
  
  # Define weights for each variable
  ac_weights <- c(gdp = 1/3, gdp_trade = 1/3, sanitation = 1/3, 
                  supermarkets = 1/3, life_expectancy = 1/3, prop_labor = 1/3,
                  hci = 1, gov_effectiveness = 1/3, fsc = 1/3, rol = 1/3)
  
  # Identify which AC variables are actually present in the data
  available_vars <- intersect(ac_variables, names(data))
  
  if (length(available_vars) == 0) {
    stop("No adaptive capacity variables found in the data")
  }
  
  # Calculate theoretical max based on ALL possible variables
  theoretical_max <- sum(ac_weights)
  
  # Calculate adaptive capacity across available variables only
  adaptive_capacity_calcs <- data %>%
    # Scale only the variables that exist
    mutate(across(all_of(available_vars), 
                  ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE)), 
                  .names = "{.col}_scaled")) %>%
    distinct(consumer_iso3c, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(
      # Calculate actual max based on available (present AND non-NA) variables for this row
      actual_max = {
        max_val <- 0
        for (var in available_vars) {
          scaled_col <- paste0(var, "_scaled")
          if (!is.na(get(scaled_col))) {
            max_val <- max_val + ac_weights[var]
          }
        }
        max_val
      },
      
      # Sum the weighted scaled values for available variables
      ac_sum = {
        sum_val <- 0
        for (var in available_vars) {
          scaled_col <- paste0(var, "_scaled")
          val <- get(scaled_col)
          if (!is.na(val)) {
            sum_val <- sum_val + val * ac_weights[var]
          }
        }
        sum_val
      },
      
      # Adjust adaptive capacity to always scale to max of 2
      adaptive_capacity = if_else(
        actual_max > 0,
        (ac_sum / actual_max) * 2,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(adaptive_capacity)) %>%
    select(consumer_iso3c, adaptive_capacity)
  
  
  # Calculate vulnerability (Risk = Exposure + Sensitivity - Adaptive Capacity)
  risk <- left_join(data %>% filter(scenario == "ssp126"), adaptive_capacity_calcs, by = "consumer_iso3c") %>%
    filter(!is.na(adaptive_capacity)) %>%
    mutate(across(c(pct_change, aa_reliance_pct), 
                  ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE)), 
                  .names = "{.col}_scaled")) %>%
    group_by(consumer_iso3c) %>%
    mutate(aa_reliance_pct_scaled = aa_reliance_pct,
           sensitivity = aa_reliance_pct_scaled) %>%
    filter(!is.na(aa_reliance_pct_scaled)) %>%
    rename(exposure = "pct_change_scaled") %>%
    select(consumer_iso3c, exposure, sensitivity, adaptive_capacity) %>%
    mutate(vulnerability = (exposure + sensitivity) - adaptive_capacity)
  
  risk <- risk %>%
    add_region(col = "consumer_iso3c", region.col.name = "region") %>%
    arrange(region, .after = "consumer_iso3c")
  
  return(risk)
}
