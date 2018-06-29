// GARCH(2,2)

data {
  int<lower=3> T;
  real r[T];
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  int<lower=0> T_out;
}
parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower = 0, upper = 1-alpha1> alpha2;
  real<lower=0, upper = 1-alpha1-alpha2> beta1;
  real<lower = 0, upper = 1-alpha1-alpha2-beta1> beta2;
}
transformed parameters {
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  sigma[2] = sigma2;
  for (t in 3:T)
    sigma[t] = sqrt(alpha0 +
                      + alpha1 * square(r[t - 1]-mu)
                      + alpha2 * square(r[t - 2]-mu)
                    + beta1  * square(sigma[t - 1])
                    + beta2 * square(sigma[t-2]));
}
model {
  r ~ normal(mu,sigma);
}
generated quantities {
  real r_out[T_out];
  r_out[1] = normal_rng(mu, sigma[1]);
  r_out[2] = normal_rng(mu, sigma[2]);
  for (t in 3:T_out) {
    r_out[t] = normal_rng(mu, sigma[t]);
  }
}
