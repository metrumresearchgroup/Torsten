#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_onecpt);

TYPED_TEST_P(test_onecpt, multiple_bolus) {
  this -> test_finite_diff_amt(1.e-3, 1.e-6);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_with_tlag) {
  this -> tlag[0][0] = 5.7;  // FIXME: tlag-generated event coincidents another event the test would fail
  this -> tlag[0][1] = 0;  // tlag2
  this -> time[0] = 0;
  for(int i = 1; i < 10; ++i) this -> time[i] = this -> time[i - 1] + 1;
  this -> amt[0] = 1200;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_tlag(1.e-5, 1.e-4);
}

TYPED_TEST_P(test_onecpt, multiple_bolus_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_biovar(1.e-3, 1.e-4);
  this -> test_finite_diff_ii(1.e-3, 1.e-3);
}

TYPED_TEST_P(test_onecpt, multiple_infusion) {
  this -> cmt[0] = 2;
  this -> rate[0] = 350;
  this -> addl[0] = 2;

  this -> biovar[0] = {0.8, 0.9};
  this -> tlag[0] = {0.4, 0.8};

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_rate(1.e-3, 1.e-6);
  this -> test_finite_diff_tlag(1.e-3, 1.e-6);
}

TYPED_TEST_P(test_onecpt, multiple_infusion_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> rate[0] = 190;        // amt/rate = even time will cause test failure
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_rate(1.e-3, 1.e-6);
}

TYPED_TEST_P(test_onecpt, multiple_infusion_2) {
  this -> reset_events(3);
  this -> amt[0] = 1200;
  this -> rate[0] = 100;
  this -> addl[0] = 0;
  this -> ii[0] = 17.0;
  this -> ss[0] = 1;
  this -> time[0] = 0.0;
  this -> time[1] = 17.0 * 0.5;
  this -> time[2] = 17.0;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_rate(1.e-3, 1.e-6);
  this -> test_finite_diff_ii(1.e-3, 1.e-6);
}

REGISTER_TYPED_TEST_SUITE_P(test_onecpt,
                            multiple_bolus,           
                            multiple_bolus_with_tlag, 
                            multiple_bolus_ss,        
                            multiple_infusion_ss,     
                            multiple_infusion,
                            multiple_infusion_2);

using onecpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_linode_functor>, // solver 1
  ::testing::Types<pmx_solve_linode_functor>, // solver 2
  ::testing::Types<double>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double> , // BIOVAR
  ::testing::Types<double> , // TLAG
  ::testing::Types<torsten::PMXOneCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_onecpt, onecpt_test_types);

// twocpt
TYPED_TEST_SUITE_P(test_twocpt);

TYPED_TEST_P(test_twocpt, multiple_bolus) {
  this -> biovar[0] = {0.8, 0.9, 0.9};
  this -> tlag[0] = {0.4, 0.8, 0.8};
  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_biovar(1.e-3, 1.e-6);
  this -> test_finite_diff_tlag(1.e-3, 5.e-4);
}

TYPED_TEST_P(test_twocpt, multiple_bolus_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 5;
  
  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_biovar(1.e-3, 1.e-6);
  this -> test_finite_diff_ii(1.e-3, 1.5e-3);
}

TYPED_TEST_P(test_twocpt, multiple_infusion) {
  this -> amt[0] = 1000;
  this -> amt[1] = 1000;
  this -> amt[2] = 700;
  this -> rate[0] = 300.0;
  this -> rate[2] = 300.0;
  this -> rate[5] = 400.0;
  this -> addl[0] = 0;
  this -> time[this -> nt - 1] = 13.0;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_rate(1.e-3, 1.e-6);
  this -> test_finite_diff_biovar(1.e-3, 1.e-6);
}

TYPED_TEST_P(test_twocpt, multiple_infusion_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; ++i) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> rate[0] = 170;        // amt/rate == event time would cause failure
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_rate(1.e-3, 1.e-6);
  this -> test_finite_diff_biovar(1.e-3, 5.e-4);
  this -> test_finite_diff_ii(1.e-3, 5.e-4);
}

TYPED_TEST_P(test_twocpt, time_dependent_theta) {
  this -> reset_events(11);
  this -> theta.resize(this -> nt);
  for (int i = 0; i < this -> nt; i++) {
    this -> theta[i].resize(5);
    if (i < 6) {
      this -> theta[i][0] = 5; // CL
    } else {
      this -> theta[i][0] = 50;  // CL is piece-wise constant      
    }
    this -> theta[i][1] = 8;  // Q
    this -> theta[i][2] = 20;  // Vc
    this -> theta[i][3] = 70;  // Vp
    this -> theta[i][4] = 1.2;  // ka
  }
	
  this -> time[0] = 0.0;
  for(int i = 1; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 2.5;

  this -> amt[0] = 1000;
  this -> cmt[0] = 1;
  this -> evid[0] = 1;
  this -> ii[0] = 12;
  this -> addl[0] = 1;

  this -> test_finite_diff_amt(1.e-3, 5.e-4);
}

TYPED_TEST_P(test_twocpt, multiple_infusion_2) {
  this -> theta[0][0] = 5;  // CL
  this -> theta[0][1] = 8;  // Q
  this -> theta[0][2] = 35;  // Vc
  this -> theta[0][3] = 105;  // Vp
  this -> theta[0][4] = 1.2;  // ka
  this -> amt[0] = 1200;
  this -> rate[0] = 1100;       // amt/rate == event time would cause failure

  this -> test_finite_diff_amt(1.e-3, 1.e-6);
  this -> test_finite_diff_rate(1.e-3, 1.e-6);
}

REGISTER_TYPED_TEST_SUITE_P(test_twocpt,
                            multiple_bolus,
                            multiple_bolus_ss,
                            multiple_infusion,
                            multiple_infusion_ss,
                            multiple_infusion_2,
                            time_dependent_theta);

using twocpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_linode_functor>, // solver 1
  ::testing::Types<pmx_solve_linode_functor>, // solver 2
  ::testing::Types<double>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<stan::math::var_value<double>> , // II
  ::testing::Types<stan::math::var_value<double>> , // PARAM
  ::testing::Types<stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<stan::math::var_value<double>> , // TLAG
  ::testing::Types<torsten::PMXTwoCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_twocpt, twocpt_test_types);
