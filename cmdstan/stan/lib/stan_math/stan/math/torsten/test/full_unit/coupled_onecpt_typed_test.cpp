#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_coupled_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_coupled_onecpt);

TYPED_TEST_P(test_coupled_onecpt, single_bolus) {
  Eigen::MatrixXd amounts(this -> nt, this -> nOde);
  amounts << 10000, 0, 0, 0, 0,
    6065.306623, 3786.208, -0.0007126788, -1.985660e-06, -4.076891e-09,
    3678.794399, 5821.649, -0.0023972982, -1.390885e-05, -5.850529e-08,
    2231.301565, 6813.189, -0.0045843443, -4.140146e-05, -2.670477e-07,
    1353.352804, 7188.323, -0.0069950553, -8.713363e-05, -7.646749e-07,
    820.849975, 7205.188, -0.0094660436, -1.520317e-04, -1.698982e-06,
    497.870679, 7019.273, -0.0119035120, -2.360139e-04, -3.219375e-06,
    301.973833, 6723.888, -0.0142554906, -3.384324e-04, -5.470933e-06,
    183.156389, 6374.696, -0.0164950219, -4.583390e-04, -8.591089e-06,
    3.354626, 3716.663, -0.0300528403, -1.915294e-03, -7.813746e-05;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x, 5.e-4);

  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 1.e-5, "theta");
  this -> compare_solvers_adj(this -> time, 1.e-2, "time");

  // FIXME: lag time sensitivity
}

TYPED_TEST_P(test_coupled_onecpt, multiple_bolus) {
  this -> addl[0] = 3;
  this -> ii[0] = 1.2;
  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 1.5e-5, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 5.e-6, "bioavailability");
}

TYPED_TEST_P(test_coupled_onecpt, multiple_infusion) {
  this -> rate[0] = 7000;
  this -> addl[0] = 3;
  this -> ii[0] = 1.2;
  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-5, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 1.e-5, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 5.e-6, "bioavailability");
}

TYPED_TEST_P(test_coupled_onecpt, multiple_bolus_ss) {
  this -> amt[0] = 1000;
  this -> addl[0] = 1;
  for(int i = 0; i < this -> nt; i++) {
    this -> time[i] = i * 5.0; 
  }
  this -> ii[0] = 8;
  this -> ss[0] = 1;
  this -> compare_solvers_val(1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 5e-4, "theta");
  this -> compare_solvers_adj(this -> ii, 1.e-3, "ii");
}

TYPED_TEST_P(test_coupled_onecpt, multiple_infusion_ss) {
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

REGISTER_TYPED_TEST_SUITE_P(test_coupled_onecpt,
                            single_bolus,
                            multiple_bolus,
                            multiple_infusion,
                            multiple_bolus_ss,
                            multiple_infusion_ss);

using onecpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_onecpt_rk45_functor,
                   pmx_solve_onecpt_bdf_functor>, // solver 1
  ::testing::Types<pmx_solve_rk45_functor>, // solver 2
  ::testing::Types<stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double, stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double, stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<stan::math::var_value<double>> , // TLAG
  ::testing::Types<CoupledOneCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_coupled_onecpt, onecpt_test_types);
