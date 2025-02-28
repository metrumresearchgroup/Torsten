#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_coupled_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_coupled_twocpt);

TYPED_TEST_P(test_coupled_twocpt, single_bolus) {
  Eigen::MatrixXd amounts(this -> nt, this -> nOde);
  amounts << 10000, 0, 0, 0, 0, 0,
    6065.306597, 3579.304,  212.1623, -0.0006874417, -1.933282e-06, -3.990995e-09,
    3678.794412, 5177.749,  678.8210, -0.0022297559, -1.318812e-05, -5.608541e-08,
    2231.301599, 5678.265, 1233.3871, -0.0041121287, -3.824790e-05, -2.508047e-07,
    1353.352829, 5597.489, 1787.0134, -0.0060546255, -7.847821e-05, -7.039447e-07,
    820.849983, 5233.332, 2295.7780, -0.0079139199, -1.335979e-04, -1.533960e-06,
    497.870681, 4753.865, 2741.1870, -0.0096246889, -2.025267e-04, -2.852515e-06,
    301.973832, 4250.712, 3118.6808, -0.0111649478, -2.838628e-04, -4.760265e-06,
    183.156387, 3771.009, 3430.9355, -0.0125357580, -3.761421e-04, -7.345522e-06,
    3.354626, 1601.493, 4374.6747, -0.0192607813, -1.370742e-03, -5.951920e-05;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x, 5.e-4);

  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 1.5e-5, "theta");
  this -> compare_solvers_adj(this -> time, 1.e-2, "time");

  // FIXME: lag time sensitivity
}

TYPED_TEST_P(test_coupled_twocpt, multiple_bolus) {
  this -> addl[0] = 3;
  this -> ii[0] = 1.2;
  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 1.5e-5, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 5.e-6, "bioavailability");
}

TYPED_TEST_P(test_coupled_twocpt, multiple_infusion) {
  this -> rate[0] = 7000;
  this -> addl[0] = 3;
  this -> ii[0] = 1.2;
  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-5, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 1.5e-5, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 5.e-6, "bioavailability");
}

TYPED_TEST_P(test_coupled_twocpt, multiple_bolus_ss) {
  this -> amt[0] = 1000;
  this -> addl[0] = 1;
  for(int i = 0; i < this -> nt; i++) {
    this -> time[i] = i * 5.0; 
  }
  this -> ii[0] = 8;
  this -> ss[0] = 1;
  this -> compare_solvers_val(2.e-5);
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 5e-4, "theta");
  this -> compare_solvers_adj(this -> ii, 1.e-3, "ii");
}

TYPED_TEST_P(test_coupled_twocpt, multiple_infusion_ss) {
  this -> amt[0] = 1000;
  this -> rate[0] = 300;
  this -> addl[0] = 1;
  for(int i = 0; i < this -> nt; i++) {
    this -> time[i] = i * 5.0; 
  }
  this -> ii[0] = 8;
  this -> ss[0] = 1;
  this -> compare_solvers_val(2.e-5);
  this -> compare_solvers_adj(this -> amt, 1.e-5, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-5, "RATE");
  this -> compare_solvers_adj(this -> theta[0], 5e-4, "theta");
  this -> compare_solvers_adj(this -> ii, 1.e-3, "ii");
}

REGISTER_TYPED_TEST_SUITE_P(test_coupled_twocpt,
                            single_bolus,
                            multiple_bolus,
                            multiple_infusion,
                            multiple_bolus_ss,
                            multiple_infusion_ss);

using twocpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_twocpt_rk45_functor,
                   pmx_solve_twocpt_bdf_functor>, // solver 1
  ::testing::Types<pmx_solve_rk45_functor>, // solver 2
  ::testing::Types< stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<stan::math::var_value<double>> , // TLAG
  ::testing::Types<CoupledTwoCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_coupled_twocpt, twocpt_test_types);
