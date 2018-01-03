#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <test/unit/util.hpp>

struct univar_integral_const_fun {
  template <typename T0, typename T1, typename T2>
  inline
  std::vector<typename stan::return_type<T1,T2>::type>
  operator()(const T0& t_in, // initial time
             const std::vector<T1>& y_in, //initial positions
             const std::vector<T2>& theta, // parameters
             const std::vector<double>& x, // double data
             const std::vector<int>& x_int,
             std::ostream* msgs) const { // integer data

    std::vector<typename stan::return_type<T1,T2>::type> res;
    res.push_back(theta.front());

    return res;
  }
};

struct univar_integral_linear_fun {
  template <typename T0, typename T1, typename T2>
  inline
  std::vector<typename stan::return_type<T1,T2>::type>
  operator()(const T0& t_in, // initial time
             const std::vector<T1>& y_in, //initial positions
             const std::vector<T2>& theta, // parameters
             const std::vector<double>& x, // double data
             const std::vector<int>& x_int,
             std::ostream* msgs) const { // integer data

    std::vector<typename stan::return_type<T1,T2>::type> res;
    res.push_back(theta.front()*t_in);

    return res;
  }
};

struct univar_integral_quad_fun {
  template <typename T0, typename T1, typename T2>
  inline
  std::vector<typename stan::return_type<T1,T2>::type>
  operator()(const T0& t_in, // initial time
             const std::vector<T1>& y_in, //initial positions
             const std::vector<T2>& theta, // parameters
             const std::vector<double>& x, // double data
             const std::vector<int>& x_int,
             std::ostream* msgs) const { // integer data

    std::vector<typename stan::return_type<T1,T2>::type> res;
    res.push_back(theta.at(0) + theta.at(1)*t_in + theta.at(2)*t_in*t_in);

    return res;
  }
};

TEST(univariate_integral, const_example) {

  univar_integral_const_fun f;
  double t0 {0.0};
  double t1 {2.5};
  int n {100};
  std::vector<double> theta;
  theta.push_back(2.0);

  EXPECT_NEAR(5.0, stan::math::univariate_integral(f, theta, t0, t1, n), 1e-5);
}

TEST(univariate_integral, linear_example) {

  univar_integral_linear_fun f;
  double t0 {0.0};
  double t1 {2.5};
  int n {100};
  std::vector<double> theta;
  theta.push_back(1.2);

  EXPECT_NEAR(3.75, stan::math::univariate_integral(f, theta, t0, t1, n), 1e-5);
}

TEST(univariate_integral, quad_example) {

  univar_integral_quad_fun f;
  double t0 {0.0};
  double t1 {0.4};
  int n {10};
  std::vector<double> theta {2.3, 2.0, 1.5};

  EXPECT_NEAR(1.112, stan::math::univariate_integral(f, theta, t0, t1, n), 1e-5);
}

