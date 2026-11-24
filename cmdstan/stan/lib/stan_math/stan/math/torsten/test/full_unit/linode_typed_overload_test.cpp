#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_onecpt);

TYPED_TEST_P(test_onecpt, multiple_bolus) {
  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_infusion) {
  this -> cmt[0] = 2;
  this -> rate[0] = 350;
  this -> addl[0] = 2;

  this -> biovar[0] = {0.8, 0.9};
  this -> tlag[0] = {0.4, 0.8};

  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_onecpt, multiple_infusion_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> rate[0] = 190;        // amt/rate = even time will cause test failure
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> compare_param_overload(1.e-10);
}

REGISTER_TYPED_TEST_SUITE_P(test_onecpt,
                            multiple_bolus,
                            multiple_bolus_ss,
                            multiple_infusion,
                            multiple_infusion_ss);

using linode_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_linode_functor>, // solver 1
  ::testing::Types<pmx_solve_linode_functor>, // solver 2
  ::testing::Types<double, stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double, stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double, stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<double, stan::math::var_value<double>> , // TLAG
  ::testing::Types<torsten::PMXOneCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_onecpt, linode_test_types);
