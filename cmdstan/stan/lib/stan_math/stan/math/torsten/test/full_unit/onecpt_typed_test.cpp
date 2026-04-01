#include <stan/math/torsten/test/unit/test_fixture_onecpt.hpp>
#include <stan/math/torsten/test/unit/test_functors.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <boost/mp11.hpp>

TYPED_TEST_SUITE_P(test_onecpt);

TYPED_TEST_P(test_onecpt, single_bolus_with_tlag) {
  this -> reset_events(2);      // only need two events
  this -> amt[0] = 1000;
  this -> evid[0] = 1;
  this -> cmt[0] = 1;
  this -> ii[0] = 0;
  this -> addl[0] = 0;
  this -> time[0] = 0.0;
  this -> tlag[0][0] = 1.5;
  this -> time[1] = 2.5;

  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 5.e-6, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 5.e-7, "lag time");
}

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

TYPED_TEST_P(test_onecpt, multiple_bolus_addl) {
  this -> ii[0] = 1.3;          // ensure test II + ADDL by end of time
  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 5.e-6, "II");
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

TYPED_TEST_P(test_onecpt, multiple_infusion) {
  this -> reset_events(3);
  this -> amt[0] = 1200;
  this -> rate[0] = 100;
  this -> addl[0] = 0;
  this -> ii[0] = 17.0;
  this -> ss[0] = 1;
  this -> time[0] = 0.0;
  this -> time[1] = 17.0 * 0.5;
  this -> time[2] = 17.0;

  this -> compare_solvers_val();
  this -> compare_solvers_adj(this -> amt, 1.e-8, "AMT");
  this -> compare_solvers_adj(this -> rate, 1.e-8, "RATE");
  this -> compare_solvers_adj(this -> ii, 5.e-6, "II");
  this -> compare_solvers_adj(this -> theta[0], 5.e-6, "theta");
  this -> compare_solvers_adj(this -> biovar[0], 1.e-6, "bioavailability");
  this -> compare_solvers_adj(this -> tlag[0], 1.e-8, "lag time");  
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

TYPED_TEST_P(test_onecpt, reset_cmt) {
  this -> reset_events(3);
  this -> evid[0] = 1;
  this -> evid[1] = 2;
  this -> evid[2] = 1;
  this -> cmt[0] = 1;
  this -> cmt[1] = -2;
  this -> amt[0] = 1000;
  this -> amt[2]= 800;
  this -> ii[0] = 0;
  this -> addl[0] = 0;
  this -> time[0] = 0.0;
  this -> time[1] = 0.25;
  this -> time[2] = this -> time[1];
  
  auto y = this -> sol1(this -> time, this -> amt, this -> rate, this -> ii, this -> evid, this -> cmt, this -> addl, this -> ss, this -> theta, this -> biovar, this -> tlag);
  EXPECT_FLOAT_EQ(stan::math::value_of(y(0, 1)), 740.8182206817178);
  EXPECT_FLOAT_EQ(stan::math::value_of(y(1, 1)), 0.0);
  EXPECT_FLOAT_EQ(stan::math::value_of(y(0, 2)), 740.8182206817178);
  EXPECT_FLOAT_EQ(stan::math::value_of(y(1, 2)), 800.0);
}

// TEST(Torsten, pmx_solve_onecptModel_SS_rate_2) {
//   // Test the special case where the infusion rate is longer than
//   // the interdose interval.
//   // THIS TEST FAILS.
//   using std::vector;

//   vector<vector<double> > pMatrix(1);
//   pMatrix[0].resize(3);
//   pMatrix[0][0] = 10;  // CL
//   pMatrix[0][1] = 80;  // Vc
//   pMatrix[0][2] = 1.2;  // ka

//   int nCmt = 2;
//   vector<vector<double> > biovar(1);
//   biovar[0].resize(nCmt);
//   biovar[0][0] = 1;  // F1
//   biovar[0][1] = 1;  // F2

//   vector<vector<double> > tlag(1);
//   tlag[0].resize(nCmt);
//   tlag[0][0] = 0;  // tlag1
//   tlag[0][1] = 0;  // tlag2

//   vector<double> time(10);
//   time[0] = 0;
//   for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
//   time[9] = 4.0;

//   vector<double> amt(10, 0);
//   amt[0] = 1200;

//   vector<double> rate(10, 0);
//   rate[0] = 75;

//   vector<int> cmt(10, 2);
//   cmt[0] = 1;

//   vector<int> evid(10, 0);
//   evid[0] = 1;

//   vector<double> ii(10, 0);
//   ii[0] = 12;

//   vector<int> addl(10, 0);
//   addl[0] = 14;

//   vector<int> ss(10, 0);

//   MatrixXd x;
//   x = pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
//                     pMatrix, biovar, tlag);

//   std::cout << x << std::endl;

//   MatrixXd amounts(10, 2);
//   amounts << 62.50420, 724.7889,
//              78.70197, 723.4747,
//              90.70158, 726.3310,
//              99.59110, 732.1591,
//              106.17663, 740.0744,
//              111.05530, 749.4253,
//              114.66951, 759.7325,
//              117.34699, 770.6441,
//              119.33051, 781.9027,
//              124.48568, 870.0308;

//   // expect_matrix_eq(amounts, x);

//   // Test AutoDiff against FiniteDiff
//   // test_pmx_solve_onecpt(time, amt, rate, ii, evid, cmt, addl, ss,
//   //                   pMatrix, biovar, tlag, 1e-8, 5e-4);
// }

REGISTER_TYPED_TEST_SUITE_P(test_onecpt,
                            multiple_bolus,           
                            multiple_bolus_cent,      
                            multiple_bolus_addl,      
                            single_bolus_with_tlag,   
                            multiple_bolus_with_tlag, 
                            multiple_bolus_ss,        
                            multiple_infusion_ss,     
                            single_infusion,          
                            single_infusion_cent,
                            multiple_infusion,
                            time_dependent_theta,
                            reset_cmt);

using onecpt_test_types = boost::mp11::mp_product<
  std::tuple,
  ::testing::Types<pmx_solve_onecpt_functor>, // solver 1
  ::testing::Types<pmx_solve_linode_functor,
                   pmx_solve_rk45_functor,
                   pmx_solve_bdf_functor,
                   pmx_solve_adams_functor>, // solver 2
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
