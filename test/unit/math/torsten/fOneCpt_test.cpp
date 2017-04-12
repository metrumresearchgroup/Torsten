#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_near_matrix_eq.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>

TEST(Torsten, fOneCpt) {
  double dt = 0.25;

  double nParameters = 3;
  Eigen::VectorXd parameters(nParameters);
  parameters << 10, 80, 1.2;  // CL, V2, ka

  int nCmt = 2;
  Eigen::VectorXd init(nCmt);
  init << 1000, 0;  // initial dose in the gut

  std::vector<double> rate(2, 0);  // no rate

  Eigen::VectorXd pred;
  pred = fOneCpt(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred(0));
  EXPECT_FLOAT_EQ(254.97490, pred(1));
}
