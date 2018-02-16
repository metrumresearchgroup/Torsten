#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <test/unit/util.hpp>

class univar_functor_ord0 {
  double val_;
public:
  explicit univar_functor_ord0(double& val) : val_(val) {}
  template <typename T0, typename T2, typename T3>
  inline
  typename stan::return_type<T0, T2, T3>::type
  operator()(const T0& t,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* msgs) const {
    typename stan::return_type<T0, T2, T3>::type res {val_};
    return res;
  }
};

class univar_functor_ord1 {
  double k_;
public:
  explicit univar_functor_ord1(double& k) : k_(k) {}
  template <typename T0, typename T2, typename T3>
  inline
  typename stan::return_type<T0, T2, T3>::type
  operator()(const T0& t,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* msgs) const {
    typename stan::return_type<T0, T2, T3>::type res {k_*t};
    return res;
  }
};

class univar_functor_ord2 {
  double a_, b_, c_;
public:
  univar_functor_ord2(double& a, double& b, double& c) :
    a_(a), b_(b), c_(c) {}

  template <typename T0, typename T2, typename T3>
  inline
  typename stan::return_type<T0, T2, T3>::type
  operator()(const T0& t,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* msgs) const {
    typename stan::return_type<T0, T2, T3>::type res {a_ + b_*t + c_*t*t};
    return res;
  }
};

TEST(univariate_integral, const_example) {
  double val {2.0};
  univar_functor_ord0 f0{val};
  double t0 = 0.0;
  double t1 = 2.5;
  std::vector<double> theta {999.9};
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::vector<stan::math::var> theta_var = stan::math::to_var(theta);
  using stan::math::univariate_integral_rk45;
  using stan::math::univariate_integral_bdf;
  
  {
    auto res { univariate_integral_rk45(f0, t0, t1, theta, x_r, x_i) };
    EXPECT_NEAR(5.0, res, 1e-8);
  }

  {
    auto res { univariate_integral_rk45(f0, t0, t1, theta_var, x_r, x_i) };
    EXPECT_NEAR(5.0, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0, t0, t1, theta, x_r, x_i) };
    EXPECT_NEAR(5.0, res, 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0, t0, t1, theta_var, x_r, x_i) };
    EXPECT_NEAR(5.0, stan::math::value_of(res), 1e-8);
  }
}

TEST(univariate_integral, linear_example) {
  double k {1.2};
  univar_functor_ord1 f0 {k};
  double t0 = 0.0;
  double t1 = 2.5;
  std::vector<double> theta {999.9};
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::vector<stan::math::var> theta_var = stan::math::to_var(theta);  
  using stan::math::univariate_integral_rk45;
  using stan::math::univariate_integral_bdf;

  {
    auto res { univariate_integral_rk45(f0, t0, t1, theta, x_r, x_i) };
    EXPECT_NEAR(3.75, res, 1e-8);
  }

  {
    auto res { univariate_integral_rk45(f0, t0, t1, theta_var, x_r, x_i) };
    EXPECT_NEAR(3.75, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0, t0, t1, theta, x_r, x_i) };
    EXPECT_NEAR(3.75, res, 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0, t0, t1, theta_var, x_r, x_i) };
    EXPECT_NEAR(3.75, stan::math::value_of(res), 1e-8);
  }
}

TEST(univariate_integral, quad_example) {
  double a{2.3}, b{2.0}, c{1.5};
  univar_functor_ord2 f0(a, b, c);
  double t0 = 0.0;
  double t1 = 0.4;
  std::vector<double> theta {999.9};
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::vector<stan::math::var> theta_var = stan::math::to_var(theta);  
  using stan::math::univariate_integral_rk45;
  using stan::math::univariate_integral_bdf;

  {
    auto res { univariate_integral_rk45(f0, t0, t1, theta, x_r, x_i) };
    EXPECT_NEAR(1.112, res, 1e-8);
  }

  {
    auto res { univariate_integral_rk45(f0, t0, t1, theta_var, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0, t0, t1, theta, x_r, x_i) };
    EXPECT_NEAR(1.112, res, 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0, t0, t1, theta_var, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }
}

TEST(univariate_integral, quad_example_var) {
  double a{2.3}, b{2.0}, c{1.5};
  univar_functor_ord2 f0(a, b, c);
  double t0 = 0.0;
  double t1 = 0.4;
  std::vector<double> theta {999.9};
  std::vector<double> x_r;
  std::vector<int> x_i;
  stan::math::var t0_var = stan::math::to_var(t0);  
  stan::math::var t1_var = stan::math::to_var(t1);  
  std::vector<stan::math::var> theta_var = stan::math::to_var(theta);  
  using stan::math::univariate_integral_rk45;
  using stan::math::univariate_integral_bdf;

  {
    auto res { univariate_integral_rk45(f0,
                                        t0_var, t1_var,
                                        theta_var, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_rk45(f0,
                                        t0_var, t1_var,
                                        theta, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0,
                                       t0_var, t1_var,
                                       theta_var, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0,
                                       t0_var, t1_var,
                                       theta, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }

  {
    auto res { univariate_integral_bdf(f0,
                                       t0, t1_var,
                                       theta, x_r, x_i) };
    EXPECT_NEAR(1.112, stan::math::value_of(res), 1e-8);
  }
}
