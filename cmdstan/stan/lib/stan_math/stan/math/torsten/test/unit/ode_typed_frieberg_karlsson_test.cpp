#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/test/unit/test_fixture_ode_non_autonomous.hpp>
#include <stan/math/torsten/test/unit/test_fixture_friberg_karlsson.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>


TYPED_TEST_SUITE_P(test_fk);

TYPED_TEST_P(test_fk, multiple_bolus) {
  this -> rtol = 1e-12;
  this -> atol = 1e-12;
  this -> addl[0] = 2;
  this -> amt[0] = 8000;
  
  this -> test_finite_diff_amt(1.e-5, 2.e-3);
  this -> test_finite_diff_theta(1.e-5, 2.5e-3);
}

TYPED_TEST_P(test_fk, multiple_bolus_ss) {
  this -> rtol = 1e-6;
  this -> atol = 1e-6;
  this -> max_num_steps = 1000000;
  this -> as_rtol = 1.e-10;
  this -> as_atol = 1.e-5;
  this -> as_max_num_steps = 100;
  this -> ss[0] = 1;
  this -> amt[0] = 80000;
  
  Eigen::MatrixXd res(10, 8);
  res << 8.000000e+04, 11996.63, 55694.35, -3.636308, -3.653620,  -3.653933, -3.653748, -3.653622,
    6.566800e+03, 53123.67, 70649.28, -3.650990, -3.653172, -3.653910, -3.653755, -3.653627,
    5.390358e+02, 34202.00, 80161.15, -3.662446, -3.653349, -3.653883, -3.653761, -3.653632,
    4.424675e+01, 23849.69, 80884.40, -3.665321, -3.653782, -3.653870, -3.653765, -3.653637,
    3.631995e+00, 19166.83, 78031.24, -3.664114, -3.654219, -3.653876, -3.653769, -3.653642,
    2.981323e-01, 16799.55, 74020.00, -3.660988, -3.654550, -3.653896, -3.653774, -3.653647,
    2.447219e-02, 15333.26, 69764.65, -3.656791, -3.654722, -3.653926, -3.653779, -3.653653,
    2.008801e-03, 14233.96, 65591.05, -3.651854, -3.654708, -3.653957, -3.653786, -3.653658,
    1.648918e-04, 13303.26, 61607.92, -3.646317, -3.654488, -3.653983, -3.653793, -3.653663,
    1.353552e-05, 12466.56, 57845.10, -3.640244, -3.654050, -3.653995, -3.653801, -3.653668;
  Eigen::MatrixXd x = res.transpose();
  this -> compare_rel_val(x, 1.e-3);

  this -> amt[0] = 8000;
  this -> test_finite_diff_amt(1.e-5, 2.e-3);
  this -> test_finite_diff_theta(1.e-5, 5.0e-1);
}

REGISTER_TYPED_TEST_SUITE_P(test_fk,
                            multiple_bolus,
                            multiple_bolus_ss);

using fk_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_rk45_functor, pmx_solve_bdf_functor>, // solver 1
  ::testing::Types<pmx_solve_rk45_functor>, // solver 2
  ::testing::Types<stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double> ,                                // BIOVAR
  ::testing::Types<double> ,                                // TLAG
  ::testing::Types<FribergKarlssonFunc>             // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_fk, fk_test_types);
