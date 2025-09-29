data {
  int<lower=0> N; // number of data points
  int<lower=0> Npop; // number of populations
  int<lower=0> Nsamp; // number of data points per population
  
  int Index[N]; // pointers to populations
  real Lengths[N]; // lengths
  real Weights[N]; // weight
}
parameters {
  real log_mean_alpha;
  real log_sigma_alpha;
  real log_mean_beta;
  real log_sigma_beta;
  real alphas[Npop];
  real<lower=0> betas[Npop];
  real log_sigma;
}
transformed parameters {
 real sigma; real mean_alpha; real mean_beta; real sigma_alpha; real sigma_beta;
 sigma = exp(log_sigma);
 mean_alpha = exp(log_mean_alpha);
 sigma_alpha = exp(log_sigma_alpha);
 mean_beta = exp(log_mean_beta);
 sigma_beta = exp(log_sigma_beta);
 }
model {
  Missing...
 
 }

