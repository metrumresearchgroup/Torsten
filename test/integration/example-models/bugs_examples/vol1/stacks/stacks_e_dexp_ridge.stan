# Stacks: robust regression and ridge regression 
# http://mathstat.helsinki.fi/openbugs/Examples/Stacks.html
# Model e) double exponential error ridge regression

data {
  int<lower=0> N;
  int<lower=0> p;
  real Y[N];
  matrix[N,p] x;
} 

// to standardize the x's 
transformed data {
  real z[N,p];
  real mean_x[p];
  real sd_x[p];
  for (j in 1:p) { 
    mean_x[j] <- mean(col(x,j)); 
    sd_x[j] <- sd(col(x,j)); 
    for (i in 1:N)  
      z[i,j] <- (x[i,j] - mean_x[j]) / sd_x[j]; 
  }
}

parameters {
  real beta0; 
  real beta[p]; 
  real<lower=0> sigmasq; 
  real<lower=0> phi;
} 

transformed parameters {
  real<lower=0> sigma;
  real mu[N];

  sigma <- sqrt(2) * sigmasq;
  for (n in 1:N)
    mu[n] <- beta0 + beta[1] * z[n, 1] + beta[2] * z[n, 2] + beta[3] * z[n, 3];
}

model {
  beta0 ~ normal(0, 316); 
  phi ~ gamma(0.01, 0.01);
  beta ~ normal(0, sqrt(phi));
  sigmasq ~ inv_gamma(.001, .001); 
  for (n in 1:N) 
    Y[n] ~ double_exponential(mu[n], sigmasq); 
} 

generated quantities {
  real b0;
  real b[p];
  real outlier_1;
  real outlier_3;
  real outlier_4;
  real outlier_21;

  for (j in 1:p)
    b[j] <- beta[j] / sd_x[j];
  b0 <- beta0 - b[1] * mean_x[1] - b[2] * mean_x[2] - b[3] * mean_x[3];

  outlier_1  <- step(fabs((Y[1] - mu[1]) / sigma) - 2.5);
  outlier_3  <- step(fabs((Y[3] - mu[3]) / sigma) - 2.5);
  outlier_4  <- step(fabs((Y[4] - mu[4]) / sigma) - 2.5);
  outlier_21 <- step(fabs((Y[21] - mu[21]) / sigma) - 2.5);
}
