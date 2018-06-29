// GARCH(1,1)

data {
  int<lower=2> T;
  real r[T];
  real<lower=0> sigma1;
  int<lower=0> T_out;
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0, upper=(1-alpha1)> beta1;
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0 +
                      + alpha1 * square(r[t - 1]-mu)
                    + beta1  * square(sigma[t - 1]));
}
model {
  r ~ normal(mu,sigma);
}
generated quantities {
  real r_out[T_out];
  r_out[1] = normal_rng(mu, sigma[1]);
  for (t in 2:T_out) {
    r_out[t] = normal_rng(mu, sigma[t]);
  }
}
