#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

TEST(Torsten, fTwoCpt) {
  double dt = 0.25;

  int nParameters = 5;
  Eigen::VectorXd parameters(nParameters);
  parameters << 5, 8, 20, 70, 1.2;  // CL, Q, VC, VP, ka

  int nCmt = 3;
  Eigen::VectorXd init(nCmt);
  init << 1000, 0, 0;  // initial dose in the gut

  std::vector<double> rate(3, 0);  // no rate

  Eigen::VectorXd pred;
  pred = fTwoCpt(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred(0));
  EXPECT_FLOAT_EQ(238.3713, pred(1));
  EXPECT_FLOAT_EQ(12.75775, pred(2));
}
