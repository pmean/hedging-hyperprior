data {
  int<lower=0> J; // number of items
  int<lower=0> y[J]; // number of successes for j
  int<lower=0> n[J]; // number of trials for j
}
parameters {
  real<lower=0, upper=1> theta[J]; // chance of success for j
  real<lower=0, upper=1> lambda; // prior mean chance of success
  real<lower=0.1> kappa; // prior count
}
transformed parameters {
  real<lower=0> alpha = lambda * kappa; // prior success count
  real<lower=0> beta = (1 - lambda) * kappa; // prior failure count
}
model {
  lambda ~ uniform(0, 1); // hyperprior
  kappa ~ pareto(2, 1.5); // hyperprior
  theta ~ beta(alpha, beta); // prior
  y ~ binomial(n, theta); // likelihood
}
