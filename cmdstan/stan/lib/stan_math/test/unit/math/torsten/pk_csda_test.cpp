#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <stan/math/torsten/pk_csda.hpp>
#include <test/unit/util.hpp>
#include <boost/range/numeric.hpp>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

struct ComplexStepDerivativeScalTest : public ::testing::Test {
  struct Fexp {
    template <typename T>
    inline T operator()(const T &theta, const std::vector<double> &x_r,
                        const std::vector<int> &x_i, std::ostream *msgs) const {
      return exp(theta) / sqrt(theta) - 0.5 * exp(theta) * pow(theta, -1.5);
    }
  };

  struct Fv {
    template <typename T>
    inline std::vector<T> operator()(const std::vector<T> &y,
                                     const std::vector<double> &x_r,
                                     const std::vector<int> &x_i,
                                     std::ostream *msgs) const {
      T x = boost::inner_product(y, y, y[0]);
      return std::vector<T>(y.size(), x);
    }
  };

  void SetUp() {}
  Fexp f;
  Fv fv;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  std::ostream *msgs = nullptr;
};

TEST_F(ComplexStepDerivativeScalTest, func_exp_sqrt) {
  using torsten::pk_csda;
  using stan::math::var;

  /* f near x = 0 has very large derivative */
  var x = 0.01;
  var y = pk_csda(f, x, x_r, x_i, 1.e-20, msgs);

  ASSERT_FLOAT_EQ(y.val(), f(x.val(), x_r, x_i, msgs));

  std::vector<stan::math::var> xv{x};
  std::vector<double> g1, g;
  var y1 = f(x, x_r, x_i, msgs);
  stan::math::set_zero_all_adjoints();
  y1.grad(xv, g1);
  stan::math::set_zero_all_adjoints();
  y.grad(xv, g);
  ASSERT_FLOAT_EQ(g[0], g1[0]);
}

TEST_F(ComplexStepDerivativeScalTest, func_exp_sqrt_default_h) {
  using torsten::pk_csda;
  using stan::math::var;

  /* f near x = 0 has very large derivative */
  var x = 0.01;
  var y = pk_csda(f, x, x_r, x_i, msgs);

  ASSERT_FLOAT_EQ(y.val(), f(x.val(), x_r, x_i, msgs));

  std::vector<stan::math::var> xv{x};
  std::vector<double> g1, g;
  var y1 = f(x, x_r, x_i, msgs);
  stan::math::set_zero_all_adjoints();
  y1.grad(xv, g1);
  stan::math::set_zero_all_adjoints();
  y.grad(xv, g);
  ASSERT_FLOAT_EQ(g[0], g1[0]);
}

TEST_F(ComplexStepDerivativeScalTest, directional_derivative) {
  using torsten::pk_csda;
  using namespace std::placeholders;
  const double h = 1.e-8;
  std::vector<double> y{1.0, 0.1, 2.0};
  std::vector<double> dy{0.0, -1.0, -1.0};
  std::vector<double> g = pk_csda(fv, y, dy, x_r, x_i, msgs);
  auto fvv = std::bind(fv, _1, x_r, x_i, msgs);
  std::vector<double> g1 = pk_csda(fvv, y, dy);
  std::vector<double> f1 = fv(y, x_r, x_i, msgs);  
  std::transform(y.begin(), y.end(), dy.begin(), y.begin(),
                 [&h](const double &x1, const double &x2) -> double {
                    return x1 + h * x2; });
  std::vector<double> f2 = fv(y, x_r, x_i, msgs);  
  for (size_t i = 0 ; i < y.size(); ++i) {
    EXPECT_FLOAT_EQ(g[i], (f2[i] - f1[i]) / h);
    EXPECT_FLOAT_EQ(g1[i], g[i]);
  }


}
