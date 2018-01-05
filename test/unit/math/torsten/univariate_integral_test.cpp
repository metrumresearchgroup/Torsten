#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <test/unit/util.hpp>

struct univar_functor_ord0 {
  template <typename T0, typename T2>
  inline
  std::vector<typename stan::return_type<T0,T2>::type>
  operator()(const T0& t_in,	// initial time
             const std::vector<T2>& theta) const { // param
    std::vector<typename stan::return_type<T0,T2>::type> res {theta.front()};
    return res;
  }
};

struct univar_functor_ord1 {
  template <typename T0, typename T2>
  inline
  std::vector<typename stan::return_type<T0,T2>::type>
  operator()(const T0& t_in,	// initial time
             const std::vector<T2>& theta) const { // param
    std::vector<typename stan::return_type<T0,T2>::type> res {theta.front()*t_in};
    return res;
  }
};

struct univar_functor_ord2 {
  template <typename T0, typename T2>
  inline
  std::vector<typename stan::return_type<T0,T2>::type>
  operator()(const T0& t_in,	// initial time
             const std::vector<T2>& theta) const { // param
    std::vector<typename stan::return_type<T0,T2>::type> res;
    res.push_back(theta.at(0) + theta.at(1)*t_in + theta.at(2)*t_in*t_in);
    return res;
  }
};

TEST(univariate_integral, const_example) {
  univar_functor_ord0 f0;
  double t0 {0.0};
  double t1 {2.5};
  std::vector<double> theta {2.0};

  EXPECT_NEAR(5.0, stan::math::univariate_integral(f0, theta, t0, t1), 1e-5);
}

TEST(univariate_integral, linear_example) {

  univar_functor_ord1 f0;
  double t0 {0.0};
  double t1 {2.5};
  std::vector<double> theta {1.2};

  EXPECT_NEAR(3.75, stan::math::univariate_integral(f0, theta, t0, t1), 1e-5);
}

TEST(univariate_integral, quad_example) {
  univar_functor_ord2 f0;
  double t0 {0.0};
  double t1 {0.4};
  std::vector<double> theta {2.3, 2.0, 1.5};

  EXPECT_NEAR(1.112, stan::math::univariate_integral(f0, theta, t0, t1), 1e-5);
}

