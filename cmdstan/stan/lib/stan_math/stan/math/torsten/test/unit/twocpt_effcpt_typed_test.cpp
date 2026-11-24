#include <stan/math/torsten/test/unit/test_fixture_twocpt_effcpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_twocpt_effcpt);
TYPED_TEST_P(test_twocpt_effcpt, multiple_bolus) {
  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 1.e-8, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 1.e-8, "lag time");
}

TYPED_TEST_P(test_twocpt_effcpt, multiple_bolus_addl) {
  this -> ii[0] = 1.55;
  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 5.e-6, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 1.e-8, "lag time");
}

REGISTER_TYPED_TEST_SUITE_P(test_twocpt_effcpt,
                            multiple_bolus,
                            multiple_bolus_addl);

using test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_twocpt_effcpt_functor>, // solver 1
  ::testing::Types<pmx_solve_linode_functor,
                   pmx_solve_rk45_functor>, // solver 2
  ::testing::Types<double>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double> , // BIOVAR
  ::testing::Types<double> , // TLAG
  ::testing::Types<torsten::PMXTwocptEffCptODE>             // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_twocpt_effcpt, test_types);
