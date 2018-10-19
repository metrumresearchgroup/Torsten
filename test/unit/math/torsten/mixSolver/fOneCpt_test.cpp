#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/PKModel/Pred/fOneCpt.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/torsten/expect_near_matrix_eq.hpp>
#include <test/unit/math/torsten/expect_matrix_eq.hpp>
#include <test/unit/math/torsten/util_generalOdeModel.hpp>

TEST(Torsten, fOneCpt_dbl) {
  double dt = 0.25;

  int nParameters = 3;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // V2
  parameters[2] = 1.2; // ka

  int nCmt = 2;
  std::vector<double> init(nCmt, 0);
  init[0] = 1000;  // initial dose in the gut

  std::vector<double> rate(2, 0);  // no rate

  std::vector<double> pred;
  pred = torsten::fOneCpt(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred[0]);
  EXPECT_FLOAT_EQ(254.97490, pred[1]);
}

TEST(Torsten, fOneCpt_var) {
  using stan::math::var;

  double dt = 0.25;
  int nParameters = 3;
  std::vector<var> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // V2
  parameters[2] = 1.2; // ka

  int nCmt = 2;
  std::vector<var> init(nCmt, 0);
  init[0] = 1000;  // initial dose in the gut

  std::vector<double> rate(2, 0);  // no rate;

  std::vector<var> pred;
  pred = torsten::fOneCpt(dt, parameters, init, rate);

  EXPECT_FLOAT_EQ(740.8182, pred[0].val());
  EXPECT_FLOAT_EQ(254.97490, pred[1].val());
}
