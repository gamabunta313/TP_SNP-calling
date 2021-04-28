// zero inflated Weibull
data {
  int<lower=1> n_obs;
  real y_shape;
  real y_scale;
  real zeros;
  real<lower=0> y[n_obs];
}

parameters {
  real t;
  real sh;
  real sc;
}

model {
  t ~ normal(logit(zeros), 1);
  sh ~ normal(log(y_shape),1);
  sc ~ normal(log(y_scale),1);

 // likelihood: Weibull hurdle
  for (i in 1:n_obs) {
    if (y[i] == 0) {
      target += bernoulli_lpmf(1 | t);
    } else {
      target += bernoulli_lpmf(0 | t) +
                weibull_lpdf(y[i] | sh, sc);
    }
  }
}

