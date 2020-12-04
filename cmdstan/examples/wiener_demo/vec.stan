data {
  int N;
  vector[N] y;  
}

parameters {
  vector<lower=0>[N] alpha;
  vector<lower=0>[N] tau;
  vector<lower=0>[N] beta;
  vector<lower=0>[N] delta;
}

model {
  alpha ~ lognormal(0, 0.5);
  tau ~ lognormal(0, 0.5);
  beta ~ lognormal(0, 0.5);
  delta ~ lognormal(0, 0.5);
  target += wiener_lpdf(y | alpha, tau, beta, delta);    
}
