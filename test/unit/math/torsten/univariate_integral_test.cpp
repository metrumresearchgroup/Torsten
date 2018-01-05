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
  univar_functor_ord0(double& val) : val_(val) {}
  template <typename T0>
  inline
  std::vector<typename stan::return_type<T0>::type>
  operator()(const T0& t_in) const { // init time
    std::vector<typename stan::return_type<T0>::type> res {val_};
    return res;
  }
};

class univar_functor_ord1 {
  double k_;
public:
  univar_functor_ord1(double& k) : k_(k) {}
  template <typename T0>
  inline
  std::vector<typename stan::return_type<T0>::type>
  operator()(const T0& t_in) const {	// initial time
    std::vector<typename stan::return_type<T0>::type> res {k_*t_in};
    return res;
  }
};

class univar_functor_ord2 {
  double a_, b_, c_;
public:
  univar_functor_ord2(double& a, double& b, double& c) :
    a_(a), b_(b), c_(c) {}
  template <typename T0>
  inline
  std::vector<typename stan::return_type<T0>::type>
  operator()(const T0& t_in) const {
    std::vector<typename stan::return_type<T0>::type> res;
    res.push_back(a_ + b_*t_in + c_*t_in*t_in);
    return res;
  }
};

TEST(univariate_integral, const_example) {
  double val {2.0};
  univar_functor_ord0 f0{val};
  double t0 {0.0};
  double t1 {2.5};

  EXPECT_NEAR(5.0, stan::math::univariate_integral(f0, t0, t1), 1e-5);
}

TEST(univariate_integral, linear_example) {
  double k {1.2};
  univar_functor_ord1 f0 {k};
  double t0 {0.0};
  double t1 {2.5};

  EXPECT_NEAR(3.75, stan::math::univariate_integral(f0, t0, t1), 1e-5);
}

TEST(univariate_integral, quad_example) {
  double a{2.3}, b{2.0}, c{1.5};
  univar_functor_ord2 f0(a, b, c);
  double t0 {0.0};
  double t1 {0.4};

  EXPECT_NEAR(1.112, stan::math::univariate_integral(f0, t0, t1), 1e-5);
}

