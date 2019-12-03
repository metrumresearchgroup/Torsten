#include <stan/math/rev/mat.hpp>  // FIX ME - more specific
#include <stan/math/torsten/PKModel/Pred/PredSS_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_general.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

// Currently test doesn't work if I do not include rev/mat.hpp.
// If I remove the header file, I get an odd bug with the
// integrator.

struct OneCpt_functor {
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;

    scalar
      CL = theta[0],
      VC = theta[1],
      ka = theta[2];

    std::vector<scalar> dydt(2);
    dydt[0] = -ka * y[0];
    dydt[1] = ka * y[0] - CL / VC * y[1];

    return dydt;
  }
};

TEST(Torsten, predSS_general_OneCpt_bolus) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double dt = 0;

  int nParameters = 3;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // VC
  parameters[2] = 1.2;  // ka

  int nCmt = 2;
  std::vector<double> biovar(nCmt, 0);
  std::vector<double> tlag(nCmt, 0);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> K(0, 0);

  // initialize Model Parameters object
  torsten::ModelParameters<double, double, double, double>
    parms(dt, parameters, biovar, tlag);

  // bolus dose
  double amt = 1200;
  double rate = 0;  // no rate
  double ii = 12;
  double cmt = 1;  // compartment number starts at 1

  // arguments for integrator constructor
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e+6;

  typedef torsten::general_functor<OneCpt_functor> F0;
  torsten::PredSS_general<F0> PredSS(F0(OneCpt_functor()), rel_tol, abs_tol, 
                            max_num_steps, 0, "rk45", nCmt);
  Matrix<double, 1, Dynamic> pred = PredSS(parms, amt, rate, ii, cmt);

  // Compare to results obtained with analytical solution
  torsten::PredSS_oneCpt PredSS_one;
  Matrix<double, 1, Dynamic>
    pred_an = PredSS_one(parms, amt, rate, ii, cmt);

  // relative error for 1st term determined empirically
  EXPECT_NEAR(pred_an(0), pred(0), pred_an(0) * 5e-2);
  EXPECT_FLOAT_EQ(pred_an(1), pred(1));
}

TEST(Torsten, predSS_general_OneCpt_truncated_infusion) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double dt = 0;

  int nParameters = 3;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // VC
  parameters[2] = 1.2;  // ka

  int nCmt = 2;
  std::vector<double> biovar(nCmt, 0);
  std::vector<double> tlag(nCmt, 0);
  torsten::ModelParameters<double, double, double, double>
    parms(dt, parameters, biovar, tlag);

  // multiple truncated infusion
  double amt = 1200;
  double rate = 1200;
  double ii = 12;
  double cmt = 1;  // compartment number starts at 1

  // arguments for integrator constructor
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e+6;

  typedef torsten::general_functor<OneCpt_functor> F0;
  torsten::PredSS_general<F0> PredSS(F0(OneCpt_functor()), rel_tol, abs_tol, 
                            max_num_steps, 0, "rk45", nCmt);
  Matrix<double, 1, Dynamic> pred = PredSS(parms, amt, rate, ii, cmt);

  // Compare to results obtained with analytical solution
  // (note: matrix exponential solution agrees with analytical solution).
  torsten::PredSS_oneCpt PredSS_one;
  Matrix<double, 1, Dynamic>
    pred_an = PredSS_one(parms, amt, rate, ii, cmt);

  // relative error for 1st term determined empirically
  double rel_err = 2e-2;
  EXPECT_NEAR(pred_an(0), pred(0), pred_an(0) * rel_err);
  EXPECT_FLOAT_EQ(pred_an(1), pred(1));
}

TEST(Torsten, predSS_general_OneCpt_constant_infusion) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  
  double dt = 0;
  
  int nParameters = 3;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // VC
  parameters[2] = 1.2;  // ka
  
  int nCmt = 2;
  std::vector<double> biovar(nCmt, 0);
  std::vector<double> tlag(nCmt, 0);
  
  torsten::ModelParameters<double, double, double, double>
    parms(dt, parameters, biovar, tlag);
  
  // constant infusion
  double amt = 1200;
  double rate = 1200;
  double ii = 0;
  double cmt = 1;  // compartment number starts at 1

  // arguments for integrator constructor
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e+6;
  
  typedef torsten::general_functor<OneCpt_functor> F0;
  torsten::PredSS_general<F0> PredSS(F0(OneCpt_functor()), rel_tol, abs_tol, 
                            max_num_steps, 0, "rk45", nCmt);
  Matrix<double, 1, Dynamic> pred = PredSS(parms, amt, rate, ii, cmt);
  
  // Compare to results obtained with analytical solution
  torsten::PredSS_oneCpt PredSS_one;
  Matrix<double, 1, Dynamic>
    pred_an = PredSS_one(parms, amt, rate, ii, cmt);
  
  // relative error for 1st term determined empirically
  EXPECT_FLOAT_EQ(pred_an(0), pred(0));
  EXPECT_FLOAT_EQ(pred_an(1), pred(1));
}

TEST(Torsten, predSS_general_exception) {
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double dt = 0;

  int nParameters = 3;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // VC
  parameters[2] = 1.2;  // ka

  int nCmt = 2;
  std::vector<double> biovar(nCmt, 0);
  std::vector<double> tlag(nCmt, 0);

  torsten::ModelParameters<double, double, double, double>
    parms(dt, parameters, biovar, tlag);

  // multiple truncated infusion
  double amt = 1200;
  double rate = 75;
  double ii = 12;
  double cmt = 1;  // compartment number starts at 1

  // arguments for integrator constructor
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e+6;

  typedef torsten::general_functor<OneCpt_functor> F0;
  torsten::PredSS_general<F0> PredSS(F0(OneCpt_functor()), rel_tol, abs_tol, 
                            max_num_steps, 0, "rk45", nCmt);
  torsten::PredSS_oneCpt PredSS_one;
  torsten::PredSS_linOde PredSS_lin;

  std::stringstream err_msg;
  err_msg << "Steady State Event: Infusion time (F * amt / rate) is 16"
          << " but must be less than the interdose interval (ii): 12!";
  std::string msg = err_msg.str();

  EXPECT_THROW_MSG(PredSS(parms, amt, rate, ii, cmt),
                   std::invalid_argument,
                   msg);

  EXPECT_THROW_MSG(PredSS_one(parms, amt, rate, ii, cmt),
                   std::invalid_argument,
                   msg);

  // FIXME: runtime error in debug mode (-g -O0)
  // EXPECT_THROW_MSG(PredSS_lin(parms, amt, rate, ii, cmt),
  //                  std::invalid_argument,
  //                  msg);
}

// Use this test for future versions
/*
TEST(Torsten, predSS_general_OneCpt_truncated_infusion_2) {
  // test the case where the duration of the infusion is longer
  // than the inter-dose interval.
  using Eigen::Matrix;
  using Eigen::Dynamic;

  double dt = 0;

  int nParameters = 3;
  std::vector<double> parameters(nParameters);
  parameters[0] = 10;  // CL
  parameters[1] = 80;  // VC
  parameters[2] = 1.2;  // ka

  int nCmt = 2;
  std::vector<double> biovar(nCmt, 0);
  std::vector<double> tlag(nCmt, 0);
  // Matrix<double, Dynamic, Dynamic> K(0, 0);
  Matrix<double, Dynamic, Dynamic> K(2, 2);
  K << -parameters[2], 0,
       parameters[2], - parameters[0] / parameters[1];

  ModelParameters<double, double, double, double, double>
    parms(dt, parameters, biovar, tlag, K);

  // multiple truncated infusion
  double amt = 1200;
  double rate = 75;
  double ii = 12;
  double cmt = 1;  // compartment number starts at 1

  // arguments for integrator constructor
  double rel_tol = 1e-6, abs_tol = 1e-6;
  long int max_num_steps = 1e+6;

  Matrix<double, 1, Dynamic>
    pred = PredSS_general_solver(parms, amt, rate, ii, cmt,
                                 general_functor<OneCpt_functor>(OneCpt_functor()),
                                 integrator_structure(rel_tol, abs_tol,
                                                      max_num_steps, 0,
                                                      "rk45"));
  std::cout << pred << std::endl;

  // Compare to results obtained with analytical solution
  // mrgsolve solution: 62.50420 724.7889
  // Neither PredSS_linOde nor PredSS_one are able to run -- the
  // call causes a failed assertion. It could be the rates are handled
  // in pred.
  
  
  // std::cout << PredSS_linOde(parms, amt, rate, ii, cmt) << std::endl;
  
  std::cout << "marker a" << std::endl;
  Matrix<double, 1, Dynamic>
    pred_an = PredSS_one(parms, amt, rate, ii, cmt);
  std::cout << "marker b" << std::endl;

  // relative error for 1st term determined empirically
  double rel_err = 2e-2;
  EXPECT_NEAR(pred_an(0), pred(0), pred_an(0) * rel_err);
  EXPECT_FLOAT_EQ(pred_an(1), pred(1));
} */
