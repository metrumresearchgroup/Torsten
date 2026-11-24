#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_twocpt);

TYPED_TEST_P(test_twocpt, single_bolus_tlag) {
  this -> reset_events(2);
  this -> evid[0] = 1;
  this -> cmt[0] = 1;
  this -> ii[0] = 0;
  this -> addl[0] = 0;
  this -> time[0] = 0.0;
  this -> tlag[0][0] = 1.5;
  this -> amt[0] = 1000;

  this -> time[0] = 0.0;
  this -> time[1] = 2.5;

  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 5.e-6, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_twocpt, multiple_bolus) {
  Eigen::MatrixXd amounts(10, 3);
  amounts << 1000.0, 0.0, 0.0,
    740.818221, 238.3713, 12.75775,
    548.811636, 379.8439, 43.55827,
    406.569660, 455.3096, 83.95657,
    301.194212, 486.6965, 128.32332,
    223.130160, 489.4507, 173.01118,
    165.298888, 474.3491, 215.75441,
    122.456428, 448.8192, 255.23842,
    90.717953, 417.9001, 290.79297,
    8.229747, 200.8720, 441.38985;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x);

  this -> biovar[0] = {0.8, 0.9, 0.9};
  this -> tlag[0] = {0.4, 0.8, 0.8};
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_twocpt, multiple_bolus_cent) {
  this -> cmt[0] = 2;

  this -> biovar[0] = {0.8, 0.9, 0.9};
  this -> tlag[0] = {0.4, 2.8, 2.8};
  this -> compare_solvers_adj(this -> theta[0], 1.e-6, "amt");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_twocpt, multiple_bolus_tlag) {
  this -> tlag[0][0] = 1.7;
  this -> ii[0] = 1.4;

  this -> compare_solvers_adj(this -> amt, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-6, "lag time");
}

TYPED_TEST_P(test_twocpt, multiple_bolus_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 5;
  
  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  Eigen::MatrixXd amounts(10, 3);
  amounts << 1.200001e+03, 224.5332, 1196.900,
    1.200001e+03, 224.5332, 1196.900, 
    2.974504e+00, 360.2587, 1533.000,
    7.373059e-03, 244.3668, 1294.140,
    3.278849e+01, 548.9136, 1534.479,
    8.127453e-02, 270.3431, 1396.353,
    3.614333e+02, 799.6771, 1304.769,
    8.959035e-01, 316.6314, 1494.540,
    2.220723e-03, 234.0179, 1244.718,
    9.875702e+00, 432.4698, 1552.220;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);
  this -> compare_solvers_adj(this -> amt, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_twocpt, multiple_infusion_addl) {
  this -> cmt[0] = 2;
  this -> rate[0] = 350;
  this -> addl[0] = 2;
  this -> ii[0] = 3.5;
  this -> time[this -> nt - 1] = 12.0;

  this -> biovar[0] = {0.8, 0.9, 0.9};
  this -> tlag[0] = {2.4, 2.8, 2.8};
  this -> compare_solvers_adj(this -> amt, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
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

  this -> compare_solvers_adj(this -> amt, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_twocpt, multiple_infusion_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < this -> nt; ++i) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> rate[0] = 150;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  Eigen::MatrixXd amounts(10, 3);
  amounts << 1.028649, 286.5656, 1391.610,
    1.028649, 286.5656, 1391.610,
    124.692706, 452.4021, 1377.667,
    11.338982, 367.1773, 1461.416,
    121.612641, 410.2024, 1340.203,
    124.991604, 477.3286, 1452.499,
    87.660547, 315.1768, 1352.746,
    124.907445, 463.2236, 1402.095,
    3.415236, 318.7214, 1432.451,
    123.979747, 436.1057, 1355.890;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);
  this -> compare_solvers_adj(this -> amt, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "amt");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_twocpt, time_dependent_theta) {
  this -> reset_events(11);
  this -> theta.resize(this -> nt);
  this -> pMatrix.resize(this -> nt);
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

    this -> pMatrix[i].resize(3, 3);
    auto CL = this -> theta[i][0];
    auto Q = this -> theta[i][1];
    auto V2 = this -> theta[i][2];
    auto V3 = this -> theta[i][3];
    auto ka = this -> theta[i][4];
    auto k10 = CL / V2;
    auto k12 = Q / V2;
    auto k21 = Q / V3;
    this -> pMatrix[i] << -ka, 0.0, 0.0, ka, -(k10 + k12), k21, 0.0, k12, -k21;
  }
	
  this -> time[0] = 0.0;
  for(int i = 1; i < this -> nt; i++) this -> time[i] = this -> time[i - 1] + 2.5;

  this -> amt[0] = 1000;
  this -> cmt[0] = 1;
  this -> evid[0] = 1;
  this -> ii[0] = 12;
  this -> addl[0] = 1;

  Eigen::MatrixXd amounts(this -> nt, 3);
  amounts << 1.000000e+03,   0.000000,   0.0000,
    4.978707e+01, 352.089056, 349.4148,
    2.478752e+00, 146.871246, 458.3010,
    1.234098e-01,  93.537648, 442.6420,
    6.144212e-03,  77.732083, 405.7800,
    5.488119e+02, 449.105589, 412.0337,
    2.732374e+01,  36.675537, 430.0023,
    1.360369e+00,  14.886990, 341.6754,
    6.772877e-02,  10.966107, 267.7033,
    3.372017e-03,   8.549649, 209.5604,
    1.678828e-04,   6.690631, 164.0364;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x);

  for (auto i = 0; i < this -> nt; ++i) {
    this -> compare_solvers_adj(this -> theta[i], 5.e-6, "theta");
  }
}

TYPED_TEST_P(test_twocpt, multiple_infusion_2) {
  this -> theta[0][0] = 5;  // CL
  this -> theta[0][1] = 8;  // Q
  this -> theta[0][2] = 35;  // Vc
  this -> theta[0][3] = 105;  // Vp
  this -> theta[0][4] = 1.2;  // ka
  
  auto CL = this -> theta[0][0];
  auto Q = this -> theta[0][1];
  auto V2 = this -> theta[0][2];
  auto V3 = this -> theta[0][3];
  auto ka = this -> theta[0][4];
  auto k10 = CL / V2;
  auto k12 = Q / V2;
  auto k21 = Q / V3;
  this -> pMatrix[0] << -ka, 0.0, 0.0, ka, -(k10 + k12), k21, 0.0, k12, -k21;

  this -> amt[0] = 1200;
  this -> rate[0] = 1200;

  Eigen::MatrixXd amounts(10, 3);
  amounts << 0.00000,   0.00000,   0.0000000,
             259.18178,  39.55748,   0.7743944,
             451.18836, 139.65573,   5.6130073,
             593.43034, 278.43884,  17.2109885,
             698.80579, 440.32663,  37.1629388,
             517.68806, 574.76950,  65.5141658,
             383.51275, 653.13596,  99.2568509,
             284.11323, 692.06145, 135.6122367,
             210.47626, 703.65965, 172.6607082,
             19.09398, 486.11014, 406.6342765;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);

  std::fill(this -> rate.begin(), this -> rate.end(), 340);
  std::fill(this -> amt.begin(), this -> amt.end(), 1000);
  std::fill(this -> evid.begin(), this -> evid.end(), 1);
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");
}

TYPED_TEST_P(test_twocpt, single_infusion_cent) {
  this -> theta[0][2] = 35;  // Vc
  this -> theta[0][3] = 105;  // Vp
  auto CL = this -> theta[0][0];
  auto Q = this -> theta[0][1];
  auto V2 = this -> theta[0][2];
  auto V3 = this -> theta[0][3];
  auto ka = this -> theta[0][4];
  auto k10 = CL / V2;
  auto k12 = Q / V2;
  auto k21 = Q / V3;
  this -> pMatrix[0] << -ka, 0.0, 0.0, ka, -(k10 + k12), k21, 0.0, k12, -k21;

  this -> amt[0] = 1200;
  this -> rate[0] = 700;
  this -> cmt[0] = 2;

  Eigen::MatrixXd amounts(10, 3);
  amounts << 0.00000,   0.00000,   0.0000000,
    0, 167.150925, 4.81832410,
    0, 319.651326,  18.583668,
    0, 458.954120, 40.3399483,
    0, 586.366896, 69.2271570,
    0, 703.066462,  104.47175,
    0, 810.111927, 145.377981,
    0, 883.621482, 191.218629,
    0, 809.159915, 235.475039,
    0, 425.085221, 450.714833;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);
  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");
}

REGISTER_TYPED_TEST_SUITE_P(test_twocpt,
                            single_bolus_tlag, multiple_bolus,
                            multiple_bolus_tlag,
                            multiple_bolus_cent,
                            multiple_bolus_ss,
                            multiple_infusion_addl,
                            multiple_infusion,
                            multiple_infusion_ss,
                            multiple_infusion_2,
                            single_infusion_cent,
                            time_dependent_theta);

using twocpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_twocpt_functor>, // solver 1
  ::testing::Types<pmx_solve_linode_functor,
                   pmx_solve_rk45_functor,
                   pmx_solve_bdf_functor>, // solver 2
  ::testing::Types<double>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double> , // BIOVAR
  ::testing::Types<double> , // TLAG
  ::testing::Types<torsten::PMXTwoCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_twocpt, twocpt_test_types);
