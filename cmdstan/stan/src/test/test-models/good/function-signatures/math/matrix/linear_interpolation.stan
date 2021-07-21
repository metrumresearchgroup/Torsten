data { 
  int d_int;
  real x_vec[d_int];
  real y_vec[d_int];
  real xout;
  int dout_int;
  real xout_vec[dout_int];
}
transformed data {
  real yout;
  real yout_vec[dout_int];

  yout = linear_interpolation(xout, x_vec, y_vec);
  yout_vec = linear_interpolation(xout_vec, x_vec, y_vec);
}
parameters {
  real px_vec[d_int];
  real py_vec[d_int];
  real pxout;
  real pxout_vec[dout_int];
  real y_p;
}
transformed parameters {
  real tp_yout;
  real tp_yout_vec[dout_int];

  tp_yout = linear_interpolation(xout, x_vec, y_vec);
  tp_yout = linear_interpolation(xout, px_vec, py_vec);
  tp_yout = linear_interpolation(pxout, x_vec, y_vec);
  tp_yout = linear_interpolation(pxout, px_vec, py_vec);

  tp_yout_vec = linear_interpolation(xout_vec, x_vec, y_vec);
  tp_yout_vec = linear_interpolation(xout_vec, px_vec, py_vec);
  tp_yout_vec = linear_interpolation(pxout_vec, x_vec, y_vec);
  tp_yout_vec = linear_interpolation(pxout_vec, px_vec, py_vec);
}

model {  
  y_p ~ normal(0,1);
}
