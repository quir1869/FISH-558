data {
  int<lower=0> N;
  vector[N] Lengths;
  vector[N] Weights;
}

parameters {
  vector[N] log_a;
  vector[N] log_b;
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  vector[N] log_a_mean;
  vector[N] log_b_mean;
}

transformed parameters {
  vector[N] a = exp(log_a);
  vector[N] b = exp(log_b);
}

model {
  
  // Priors
  sigma   ~ cauchy(0, 2.5);
  sigma_a ~ cauchy(0, 2.5);
  sigma_b ~ cauchy(0, 2.5);
  
  // A and b samples given prior information
  log_a_mean ~ normal(0, 1000);
  log_b_mean ~ normal(0, 1000);
  
  log_a ~ normal(log_a_mean, sigma_a);
  log_b ~ normal(log_b_mean, sigma_b);

  // Likelihood
  for (t in 1:N) {
    Weights[t] ~ ;
  }
}
