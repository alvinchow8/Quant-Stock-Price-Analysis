data {
  int<lower=0> T;                       // Number of time intervals
  int<lower=2> N;                       // Number of assets           
  vector[N] ret[T];                     // Stock return data
  real<lower=0> sigma_init[N];          // Initial values of Sigma squared
}
transformed data{
  int<lower=1> N_choose_2;              // Number of strictly lower-than-diagonal elements of W
  N_choose_2 = (N*(N - 1)) / 2;
}
parameters {
  vector[N] mu;                         // Vector of mean returns
  real<lower=0> alpha[N];               // Alpha parameters for N assets
  real<lower=0,upper=1> beta[N];        // Beta parameters for N assets
  real<lower=0,upper=1> gamma[N];       // Gamma parameters for N assets
  vector[N_choose_2] W_lower;           // Vector of strictly lower-than-diagonal elements of W
}
transformed parameters {
  cholesky_factor_cov[N] W;             // Cholesky Factor Matrix of an NxN covariance matrix
  real<lower=0> sig[N,T];               // N sigma squared values for each time interval of T     
  for (n in 1:N) {
    W[n,n] = 1;
  }
  for (i in 1:N_choose_2) {
    for(j in 2:N) {
      for(k in 1:(j-1)) {
        W[j,k] = W_lower[i];
        W[k,j] = 0;
      }
    }
  }
  for (n in 1:N)
    sig[n,1] = sigma_init[n];
  for (t in 2:T) {
    for (n in 1:N) {
      sig[n,t] = alpha[n] + beta[n]*pow((to_array_1d(ret[t-1] - mu)[n]), 2) + gamma[n]*sig[n,t-1];
    }
  }
}
model {
  for (t in 1:T) {
    ret[t] ~ multi_normal_cholesky(mu, W*diag_matrix(sqrt(to_vector(sig[1:N,t]))));
  } 
}
generated quantities {
  vector[N] ret_out[T];
  for (t in 1:T) {
    ret_out[t] = multi_normal_cholesky_rng(mu, W*diag_matrix(sqrt(to_vector(sig[1:N,t]))));
  }
}
