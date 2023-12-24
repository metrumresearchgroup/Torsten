#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_fixture_twocpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_onecpt);

TYPED_TEST_P(test_onecpt, multiple_bolus) {
  Eigen::MatrixXd amounts(10, 2);
  amounts << 1000.0, 0.0,
    740.8182, 254.97490,
    548.8116, 436.02020,
    406.5697, 562.53846,
    301.1942, 648.89603,
    223.1302, 705.72856,
    165.2989, 740.90816,
    122.4564, 760.25988,
    90.71795, 768.09246,
    8.229747, 667.87079;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x);

  this -> biovar[0] = {0.8, 0.9};
  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 1.e-8, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 1.e-8, "lag time");
}

TYPED_TEST_P(test_onecpt, multiple_bolus_cent) {
  this -> cmt[0] = 2;
  Eigen::MatrixXd amounts(10, 2);
  amounts << 0, 1000.0000,
    0,  969.2332,
    0,  939.4131,
    0 , 910.5104,
    0,  882.4969,
    0,  855.3453,
    0,  829.0291,
    0,  803.5226,
    0,  778.8008,
    0,  606.5307;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x);

  this -> biovar[0] = {0.8, 0.9};
  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 1.e-8, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 1.e-8, "lag time");
}

TYPED_TEST_P(test_onecpt, multiple_bolus_with_tlag) {
  this -> tlag[0][0] = 5;  // tlag1
  this -> tlag[0][1] = 0;  // tlag2
  this -> time[0] = 0;
  for(int i = 1; i < 10; ++i) this -> time[i] = this -> time[i - 1] + 1;
  this -> amt[0] = 1200;

  Eigen::MatrixXd amounts(10, 2);
  amounts << 0, 0,
             0, 0,
             0, 0,
             0, 0,
             0, 0,
             0, 0,
             361.433054, 778.6752,
             108.861544, 921.7110,
             32.788467, 884.0469,
             9.875696, 801.4449;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 6.e-7, "lag time");
}

TYPED_TEST_P(test_onecpt, multiple_bolus_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  Eigen::MatrixXd amounts(10, 2);
  amounts << 1200.0, 384.7363,
    1200.0, 384.7363,
    2.974504, 919.6159,
    7.373062e-3, 494.0040,
    3.278849e+1, 1148.4725,
    8.127454e-2, 634.2335,
    3.614333e+2, 1118.2043,
    8.959035e-1, 813.4883,
    2.220724e-3, 435.9617,
    9.875702, 1034.7998;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_rel_val(x, 1.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");

  this -> tlag[0][0] = 1.7;
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_onecpt, multiple_infusion_ss) {
  this -> time[0] = 0.0;
  this -> time[1] = 0.0;
  for(int i = 2; i < 10; i++) this -> time[i] = this -> time[i - 1] + 5;
  this -> amt[0] = 1200;
  this -> rate[0] = 150;
  this -> addl[0] = 10;
  this -> ss[0] = 1;

  Eigen::MatrixXd amounts(10, 2);
  amounts << 1.028649, 659.9385,
    1.028649, 659.9385,
    124.692706, 837.1959,
    11.338982, 836.1947,
    121.612641, 737.4911,
    124.991604, 950.4222,
    87.660547, 642.9529,
    124.907445, 879.6271,
    3.415236, 745.2971,
    123.979747, 789.6393;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);

  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");

  this -> biovar[0] = {0.8, 0.9};
  this -> tlag[0] = {2.4, 1.7};
  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

TYPED_TEST_P(test_onecpt, single_infusion) {
  this -> amt[0] = 1200;
  this -> rate[0] = 1200;
  this -> addl[0] = 0;

  Eigen::MatrixXd amounts(10, 2);
  amounts << 0.00000,   0.00000,
             259.18178,  40.38605,
             451.18836, 145.61440,
             593.43034, 296.56207,
             698.80579, 479.13371,
             517.68806, 642.57025,
             383.51275, 754.79790,
             284.11323, 829.36134,
             210.47626, 876.28631,
             19.09398, 844.11769;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x);

  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");
  this -> compare_solvers_adj(this -> theta[0], 1.5e-6, "theta");
}

TYPED_TEST_P(test_onecpt, single_infusion_cent) {
  this -> amt[0] = 1000;
  this -> rate[0] = 600;
  this -> evid[0] = 1;
  this -> cmt[0] = 2;

  Eigen::MatrixXd amounts(10, 2);
  amounts << 
    0,           0,
    0, 147.6804745,
    0, 290.8172984,
    0, 429.5502653,
    0, 564.0148675,
    0, 694.3424289,
    0, 820.6602327,
    0, 893.3511610,
    0,  865.865635,
    0, 674.3368348;
  Eigen::MatrixXd x = amounts.transpose();
  this -> compare_val(x, 5.e-7);

  this -> amt[0] = 1200;
  this -> rate[0] = 1200;
  this -> cmt[0] = 2;
  this -> addl[0] = 0;

  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");
  this -> compare_solvers_adj(this -> theta[0], 1.e-6, "theta");
}

TYPED_TEST_P(test_onecpt, time_dependent_theta) {
  this -> reset_events(11);
  int nt = this -> nt;

  this -> theta.resize(nt);
  this -> pMatrix.resize(nt);
  for (int i = 0; i < nt; i++) {
    this -> theta[i].resize(3);
    this -> pMatrix[i].resize(2, 2);
    if (i < 6) {
      this -> theta[i][0] = 10; // CL      
    } else {
      this -> theta[i][0] = 50; // CL is piece-wise contant      
    } 
    this -> theta[i][1] = 80; // Vc
    this -> theta[i][2] = 1.2; // ka
    // for linode
    this -> pMatrix[i] << - this -> theta[i][2], 0.0,
      this -> theta[i][2], - this -> theta[i][0]/ this -> theta[i][1];
  }

  this -> time[0] = 0.0;
  for(int i = 1; i < nt; i++) {
    this -> time[i] = this -> time[i - 1] + 2.5; 
  }

  this -> amt[0] = 1000;
  this -> cmt[0] = 1;
  this -> evid[0] = 1;
  this -> ii[0] = 12;
  this -> addl[0] = 1;

  Eigen::MatrixXd amounts(nt, 2);
  amounts << 1000.0, 0.0,
    4.978707e+01, 761.1109513,
    2.478752e+00, 594.7341503,
    1.234098e-01, 437.0034049,
    6.144212e-03, 319.8124495,
    5.488119e+02, 670.0046601,
    2.732374e+01, 323.4948561,
    1.360369e+00, 76.9219400,
    6.772877e-02, 16.5774607,
    3.372017e-03, 3.4974152,
    1.678828e-04, 0.7342228;
  Eigen::MatrixXd x = amounts.transpose();

  this -> compare_val(x);
  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
}

REGISTER_TYPED_TEST_SUITE_P(test_onecpt,
                            multiple_bolus, 
                            multiple_bolus_cent, 
                            multiple_bolus_with_tlag, 
                            multiple_bolus_ss, 
                            multiple_infusion_ss, 
                            single_infusion, 
                            single_infusion_cent, 
                            time_dependent_theta);

using onecpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_linode_functor>, // solver 1
  ::testing::Types<pmx_solve_onecpt_functor>, // solver 2
  ::testing::Types<double, stan::math::var_value<double>>,  // TIME
  ::testing::Types<double, stan::math::var_value<double>>,  // AMT
  ::testing::Types<double, stan::math::var_value<double>> , // RATE
  ::testing::Types<double, stan::math::var_value<double>> , // II
  ::testing::Types<double, stan::math::var_value<double>> , // PARAM
  ::testing::Types<double, stan::math::var_value<double>> , // BIOVAR
  ::testing::Types<double, stan::math::var_value<double>> , // TLAG
  ::testing::Types<torsten::PMXOneCptODE>                   // ODE
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(PMX, test_onecpt, onecpt_test_types);

// twocpt
TYPED_TEST_SUITE_P(test_twocpt);

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

  this -> compare_val(x, 2.e-6);
  this -> compare_solvers_adj(this -> amt, 1.e-6, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-6, "RATE");
}

REGISTER_TYPED_TEST_SUITE_P(test_twocpt,
                            multiple_bolus, 
                            multiple_bolus_ss, 
                            multiple_infusion_ss, 
                            time_dependent_theta, 
                            multiple_infusion_2, 
                            single_infusion_cent);

using twocpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_linode_functor>, // solver 1
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
