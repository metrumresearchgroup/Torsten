data {
  int<lower=0> J;
  int<lower=0> N;
  int<lower=1,upper=J> county[N];
  vector[N] y;
}
parameters {
  vector[J] eta;
  real mu_a;
  real<lower=0,upper=100> sigma_a;
  real<lower=0,upper=100> sigma_y;
}
transformed parameters {
  vector[J] a;
  vector[N] y_hat;

  a = mu_a + sigma_a * eta;

  for (i in 1:N)
    y_hat[i] = a[county[i]];
}
model {
  mu_a ~ normal(0, 1);

  eta ~ normal(0, 1);
  y ~ normal(y_hat, sigma_y);
}
generated quantities {
  vector[N] e_y;
  real<lower=0> s_a;
  real<lower=0> s_y;

  e_y = y - y_hat;
  s_a = sd(a);
  s_y = sd(e_y);
}