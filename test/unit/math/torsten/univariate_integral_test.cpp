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
  operator()(const T0& t_in,     // init time
             std::ostream* msgs) const {
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
  operator()(const T0& t_in,	// initial time
             std::ostream* msgs) const {
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
  operator()(const T0& t_in,	// initial time
             std::ostream* msgs) const {
    std::vector<typename stan::return_type<T0>::type> res;
    res.push_back(a_ + b_*t_in + c_*t_in*t_in);
    return res;
  }
};

class univar_functor_vec {
    univar_functor_ord0 f0_;
    univar_functor_ord1 f1_;
    univar_functor_ord2 f2_;
public:
  univar_functor_vec(double& val, double& k, double& a, double& b, double& c) :
    f0_(val), f1_(k), f2_(a, b, c)
  {}
  template <typename T0>
  inline
  std::vector<typename stan::return_type<T0>::type>
  operator()(const T0& t_in,	// initial time
             std::ostream* msgs) const {
    std::vector<typename stan::return_type<T0>::type> res {
      f0_(t_in, msgs)[0], f1_(t_in, msgs)[0], f2_(t_in, msgs)[0] };
    return res;
  }
};

TEST(univariate_integral, const_example) {
  double val {2.0};
  univar_functor_ord0 f0{val};
  std::vector<double> t {0.0, 2.5};
  std::vector<double> y {1.5};

  {
    auto res { stan::math::univariate_integral_rk45(f0, y, t) };
    EXPECT_NEAR(6.5, res[0], 1e-8);
  }

  {
    auto res { stan::math::univariate_integral_bdf(f0, y, t) };
    EXPECT_NEAR(6.5, res[0], 1e-8);
  }
}

TEST(univariate_integral, linear_example) {
  double k {1.2};
  univar_functor_ord1 f0 {k};
  std::vector<double> t {0.0, 2.5};
  std::vector<double> y {0.0};

  {
    auto res { stan::math::univariate_integral_rk45(f0, y, t) };
    EXPECT_NEAR(3.75, res[0], 1e-8);
  }

  {
    auto res { stan::math::univariate_integral_bdf(f0, y, t) };
    EXPECT_NEAR(3.75, res[0], 1e-8);
  }
}

TEST(univariate_integral, quad_example) {
  double a{2.3}, b{2.0}, c{1.5};
  univar_functor_ord2 f0(a, b, c);
  std::vector<double> t {0.0, 0.4};
  std::vector<double> y {0.0};

  {
    auto res { stan::math::univariate_integral_rk45(f0, y, t) };
    EXPECT_NEAR(1.112, res[0], 1e-8);
  }

  {
    auto res { stan::math::univariate_integral_bdf(f0, y, t) };
    EXPECT_NEAR(1.112, res[0], 1e-8);
  }

}

TEST(univariate_integral, vector_example) {
  double val {2.0}, k {1.2}, a{2.3}, b{2.0}, c{1.5};
  univar_functor_vec f(val, k, a, b, c);
  std::vector<double> t {0.0, 2.5};
  std::vector<double> y {1.5, 0.0, 0.0};

  {
    auto res { stan::math::univariate_integral_rk45(f, y, t) };
    EXPECT_NEAR(6.50, res[0], 1e-8);
    EXPECT_NEAR(3.75, res[1], 1e-8);
  }

  {
    auto res { stan::math::univariate_integral_bdf(f, y, t) };
    EXPECT_NEAR(6.50, res[0], 1e-8);
    EXPECT_NEAR(3.75, res[1], 1e-8);
  }

}
