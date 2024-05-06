data { 
  int d_int;
  real d_real;
}
transformed data {
  int transformed_data_int;
  real transformed_data_real;

  transformed_data_real = tan(d_int);
  transformed_data_real = tan(d_real);
}
parameters {
  real p_real;
  real y_p;
}
transformed parameters {
  real transformed_param_real;

  transformed_param_real = tan(d_int);
  transformed_param_real = tan(d_real);
  transformed_param_real = tan(p_real);
}
model {  
  y_p ~ normal(0,1);
}
