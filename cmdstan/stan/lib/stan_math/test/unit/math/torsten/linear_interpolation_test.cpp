#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <test/unit/util.hpp>

TEST(linear_interpolation, linear_example) {
  int nx = 5, nout = 3;
  std::vector<double> x(nx), y(nx), xout(nout), yout;
  double youtTrue;

  for(int i = 0; i < nx; i++){
    x[i] = i;
    y[i] = i;
  }

  xout[0] = 1.5;
  xout[1] = 2.5;
  xout[2] = 4.5;

  yout = torsten::linear_interpolation(xout, x, y);
  for(int i = 0; i < nout; i++){
      if(xout[i] <= x[0]){
	youtTrue = x[0];
      }else if(xout[i] >= x[nx - 1]){
	youtTrue = x[nx - 1];
      }else{
	youtTrue = xout[i];
      }
      EXPECT_FLOAT_EQ(youtTrue, yout[i]);
  }
}

TEST(linear_interpolation, xgradient){
  using stan::math::var;
  int nx = 5, nout = 3;
  std::vector<double> xdbl(nx), ydbl(nx), xoutdbl(nout), thisGrad(nx);
  Eigen::MatrixXd trueJac = Eigen::MatrixXd::Zero(3, 5); 

  for(int i = 0; i < nx; i++){
    xdbl[i] = i;
    ydbl[i] = i;
  }

  xoutdbl[0] = 1.5;
  xoutdbl[1] = 2.5;
  xoutdbl[2] = 4.5;

  trueJac(0, 1) = (ydbl[2] - ydbl[1]) * (xoutdbl[0] - xdbl[2]) / pow(xdbl[2] - xdbl[1], 2);
  trueJac(0, 2) = -(ydbl[2] - ydbl[1]) * (xoutdbl[0] - xdbl[1]) / pow(xdbl[2] - xdbl[1], 2);

  trueJac(1, 2) = (ydbl[3] - ydbl[2]) * (xoutdbl[1] - xdbl[3]) / pow(xdbl[3] - xdbl[2], 2);
  trueJac(1, 3) = -(ydbl[3] - ydbl[2]) * (xoutdbl[1] - xdbl[2]) / pow(xdbl[3] - xdbl[2], 2);

  for(int i = 0; i < nout; i++){
    std::vector<var> x(nx), y(nx), xout(nout), yout;
    for(int k = 0; k < nx; k++){
      x[k] = k;
      y[k] = k;
    }
    for(int k = 0; k < nout; k++){
      xout[k] = xoutdbl[k];
    }
    yout = torsten::linear_interpolation(xout, x, y);

    yout[i].grad(x, thisGrad);

    for(int j = 0; j < nx; j++){
      EXPECT_EQ(trueJac(i, j), thisGrad[j]);
    }
  }
}

TEST(linear_interpolation, ygradient){
  using stan::math::var;
  int nx = 5, nout = 3;
  std::vector<double> xdbl(nx), ydbl(nx), xoutdbl(nout), thisGrad(nx);
  Eigen::MatrixXd trueJac = Eigen::MatrixXd::Zero(3, 5); 

  for(int i = 0; i < nx; i++){
    xdbl[i] = i;
    ydbl[i] = i;
  }

  xoutdbl[0] = 1.5;
  xoutdbl[1] = 2.5;
  xoutdbl[2] = 4.5;

  trueJac(0, 2) = (xoutdbl[0] - xdbl[1]) / (xdbl[2] - xdbl[1]);
  trueJac(0, 1) = 1 - trueJac(0, 2);

  trueJac(1, 3) = (xoutdbl[1] - xdbl[2]) / (xdbl[3] - xdbl[2]);
  trueJac(1, 2) = 1 - trueJac(1, 3);

  trueJac(2, 4) = 1;

  for(int i = 0; i < nout; i++){
    std::vector<var> x(nx), y(nx), xout(nout), yout;
    for(int k = 0; k < nx; k++){
      x[k] = k;
      y[k] = k;
    }
    for(int k = 0; k < nout; k++){
      xout[k] = xoutdbl[k];
    }
    yout = torsten::linear_interpolation(xout, x, y);

    yout[i].grad(y, thisGrad);
    
    for(int j = 0; j < nx; j++){
      EXPECT_EQ(trueJac(i, j), thisGrad[j]);
    }
  }
}

TEST(linear_interpolation, xoutgradient){
  using stan::math::var;
  int nx = 5, nout = 3;
  std::vector<double> xdbl(nx), ydbl(nx), xoutdbl(nout), thisGrad(nx);
  Eigen::MatrixXd trueJac = Eigen::MatrixXd::Zero(3, 3); 

  for(int i = 0; i < nx; i++){
    xdbl[i] = i;
    ydbl[i] = i;
  }

  xoutdbl[0] = 1.5;
  xoutdbl[1] = 2.5;
  xoutdbl[2] = 4.5;

  trueJac(0, 0) = (ydbl[2] - ydbl[1]) / (xdbl[2] - xdbl[1]);
  trueJac(1, 1) = (ydbl[3] - ydbl[2]) / (xdbl[3] - xdbl[2]);

  for(int i = 0; i < nout; i++){
    std::vector<var> x(nx), y(nx), xout(nout), yout;
    for(int k = 0; k < nx; k++){
      x[k] = k;
      y[k] = k;
    }
    for(int k = 0; k < nout; k++){
      xout[k] = xoutdbl[k];
    }
    yout = torsten::linear_interpolation(xout, x, y);

    yout[i].grad(xout, thisGrad);

    for(int j = 0; j < nout; j++){
      EXPECT_EQ(trueJac(i, j), thisGrad[j]);
    }
  }
}

TEST(linear_interpolation, error_conditions) {
  using stan::math::var;
  int nx = 5, nout = 3;
  std::vector<double> x(nx), y(nx), xout(nout), yout;

  for(int i = 0; i < nx; i++){
    x[i] = i;
    y[i] = i;
  }

  xout[0] = 1.5;
  xout[1] = 2.5;
  xout[2] = 4.5;

  std::vector<double> xout_bad;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout_bad, x, y),
                   std::invalid_argument,
                   "xout has size 0");

  std::vector<double> x_bad;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x_bad, y),
                   std::invalid_argument,
                   "x has size 0");

  std::vector<double> y_bad;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x, y_bad),
                   std::invalid_argument,
                   "y has size 0");

  std::vector<double> x3_bad = x;
  x3_bad[2] = 0.0;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x3_bad, y),
                   std::domain_error,
                   "x is not a valid ordered vector");

  std::vector<double> x2_bad(nx - 1);  
  for(int i = 0; i < (nx - 1); i++) x2_bad[i] = x[i];
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x2_bad, y),
                   std::invalid_argument,
                   "size of x (4) and size of y (5) must match in size");
}

TEST(linear_interpolation, error_conditions_inf) {
  using stan::math::var;
  std::stringstream expected_is_inf;
  expected_is_inf << "is " << std::numeric_limits<double>::infinity();
  std::stringstream expected_is_neg_inf;
  expected_is_neg_inf << "is " << -std::numeric_limits<double>::infinity();
  double inf = std::numeric_limits<double>::infinity();
  int nx = 5, nout = 3;
  std::vector<double> x(nx), y(nx), xout(nout), yout;

  for(int i = 0; i < nx; i++){
    x[i] = i;
    y[i] = i;
  }

  xout[0] = 1.5;
  xout[1] = 2.5;
  xout[2] = 4.5;

  std::vector<double> xout_bad = xout;
  xout_bad[0] = inf;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout_bad, x, y),
                   std::domain_error,
                   "xout");
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout_bad, x, y),
                   std::domain_error,
                   expected_is_inf.str());

  xout_bad = xout;
  xout_bad[0] = -inf;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout_bad, x, y),
                   std::domain_error,
                   "xout");
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout_bad, x, y),
                   std::domain_error,
                   expected_is_neg_inf.str());

  std::vector<double> x_bad = x;
  x_bad[0] = inf;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x_bad, y),
                   std::domain_error,
                   "x");
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x_bad, y),
                   std::domain_error,
                   expected_is_inf.str());

  std::vector<double> y_bad = y;
  y_bad[0] = -inf;
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x, y_bad),
                   std::domain_error,
                   "y");
  EXPECT_THROW_MSG(torsten::linear_interpolation(xout, x, y_bad),
                   std::domain_error,
                   expected_is_neg_inf.str());
}
