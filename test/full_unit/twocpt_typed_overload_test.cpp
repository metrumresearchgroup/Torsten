#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_twocpt);

TYPED_TEST_P(test_twocpt, multiple_bolus) {
  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_twocpt, multiple_bolus_with_tlag) {
  this -> tlag[0][0] = 1.7;  // FIXME: tlag-generated event coincidents another event the test would fail
  this -> time[0] = 0;
  for(int i = 1; i < this -> nt; ++i) this -> time[i] = this -> time[i - 1] + 1;
  this -> amt[0] = 1200;

  this -> compare_param_overload(1.e-10);
}

TYPED_TEST_P(test_twocpt, multiple_bolus_ss) {
  this -> reset_events(3);
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 5;

  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_twocpt, multiple_infusion) {
  this -> reset_events(3);
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> addl[0] = 1;
  this -> amt[0] = 1200;
  this -> rate[0] = 150;
  this -> addl[0] = 10;

  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_twocpt, multiple_infusion_ss) {
  this -> reset_events(3);
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> addl[0] = 1;
  this -> amt[0] = 1200;
  this -> rate[0] = 150;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> compare_overload(1.e-10);
}

TYPED_TEST_P(test_twocpt, single_infusion_cent) {
  this -> reset_events(3);
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> theta[0][2] = 35;  // Vc
  this -> theta[0][3] = 105;  // Vp
  this -> amt[0] = 1200;
  this -> cmt[0] = 2;
  this -> rate[0] = 780;

  this -> compare_overload(1.e-10);
}

REGISTER_TYPED_TEST_SUITE_P(test_twocpt,
                            multiple_bolus,
                            multiple_bolus_with_tlag,
                            multiple_bolus_ss,
                            multiple_infusion,
                            multiple_infusion_ss,
                            single_infusion_cent);

using twocpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_twocpt_functor>, // solver 1
  ::testing::Types<pmx_solve_twocpt_functor>, // solver 2
  ::testing::Types<double, stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double, stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double, stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<double, stan::math::var_value<double>> , // TLAG
  ::testing::Types<torsten::PMXTwoCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_twocpt, twocpt_test_types);
