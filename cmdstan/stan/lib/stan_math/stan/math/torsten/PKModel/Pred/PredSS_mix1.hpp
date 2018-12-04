#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_MIX1_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_MIX1_HPP

#include <stan/math/torsten/PKModel/integrator.hpp>
#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/torsten/PKModel/functors/SS_system.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_oneCpt.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_oneCpt.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <vector>
#include <iostream>
#include <string>

namespace torsten {

template <typename F>
struct PredSS_mix1 {
  F f_;
  integrator_structure integrator_;
  int nOde_;  // number of states in the reduced system

  PredSS_mix1(const F& f,
              const double& rel_tol,
              const double& abs_tol,
              const long int& max_num_steps,  // NOLINT
              std::ostream* msgs,
              const std::string& integratorType,
              const int& nOde)
    : f_(f),
      integrator_(rel_tol, abs_tol, max_num_steps, msgs, integratorType),
      nOde_(nOde) { }

  /**
   * Mix1 compartment model using built-in ODE solver, mixed
   * solving method, and root-finder.
   * Calculate amount in each compartment at the end of a
   * steady-state dosing interval or during a steady-state
   * constant input (if ii = 0). The function is overloaded
   * to address the cases where amt or rate may be fixed or
   * random variables (yielding a total of 4 cases).
   * 
   * Case 1 (dd): amt and rate are fixed.
   *
   *	 @tparam T_time type of scalar for time
   *	 @tparam T_ii type of scalar for interdose interval
   *	 @tparam T_parameters type of scalar for ODE parameters
   *   @tparam T_biovar type of scalar for bio-availability
   *	 @tparam F type of ODE system function
   *	 @param[in] parameter model parameters at current event
   *	 @param[in] rate
   *	 @param[in] ii interdose interval
   *	 @param[in] cmt compartment in which the event occurs
   *	 @param[in] f functor for base ordinary differential equation
   *              that defines compartment model
   *   @return an eigen vector that contains predicted amount in each
   *           compartment at the current event.
   */
  template<typename T_time,
           typename T_ii,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_ii,
    T_parameters>::type, Eigen::Dynamic, 1>
  operator() (const ModelParameters<T_time,
                                    T_parameters,
                                    T_biovar,
                                    T_tlag>& parameter,
              const double& amt,
              const double& rate,
              const T_ii& ii,
              const int& cmt) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using Eigen::VectorXd;
    using std::vector;
    using stan::math::algebra_solver;
    using stan::math::to_vector;
    using stan::math::to_array_1d;

    typedef typename boost::math::tools::promote_args<T_ii,
      T_parameters>::type scalar;

    double ii_dbl = unpromote(ii);

    // Compute solution for base 1cpt PK
    Matrix<T_parameters, Dynamic, 1> predPK;
    int nPK = 2;
    if (cmt <= 2) {  // check dosing occurs in a base state
      PredSS_oneCpt PredSS_one;
      int nParmsPK = 3;
      predPK = PredSS_one(parameter.truncate(nParmsPK),
                          amt, (cmt <= 2) ? rate : 0, ii_dbl, cmt);
    } else {
      predPK = Matrix<scalar, Dynamic, 1>::Zero(nPK);
    }

    // Arguments for ODE integrator (and initial guess)
    Matrix<double, 1, Dynamic> init_dbl
      = Matrix<double, 1, Dynamic>::Zero(nOde_);
    vector<double> x_r(nPK + nOde_, 0);  // rate for the full system
    x_r.push_back(0);  // include initial time (at SS, t0 = 0)
    vector<int> x_i;

    // Tuning parameters for algebraic solver
    double rel_tol = 1e-10;  // default
    double f_tol = 1e-4;  // empirical
    long int max_num_steps = 1e3;  // default // NOLINT

    // construct algebraic system functor: note we adjust cmt
    // such that 1 corresponds to the first state we compute
    // numerically.
    SS_system_dd<ode_rate_dbl_functor<F>, Pred1_oneCpt >
      system(ode_rate_dbl_functor<F>(f_), Pred1_oneCpt(), ii_dbl, cmt,
             integrator_, nPK);

    Matrix<double, Dynamic, 1> predPD_guess;
    Matrix<scalar, 1, Dynamic> predPD;

    if (rate == 0) {  // bolus dose
      if (cmt > 2) init_dbl(cmt - 1) = amt;
      else
        predPK(cmt - 1) += amt;

      // Construct augmented parameters
      ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
        theta = parameter.augment(predPK);
      theta.time(ii_dbl);

      predPD_guess = to_vector(integrator_(ode_rate_dbl_functor<F>(f_),
                                        to_array_1d(init_dbl),
                                        0.0, std::vector<double>(1, ii_dbl),
                                        unpromote(theta.get_RealParameters()),
                                        x_r, x_i)[0]);

      x_r.push_back(amt);
      predPD = algebra_solver(system, predPD_guess,
                              to_vector(theta.get_RealParameters()),
                              x_r, x_i,
                              0, rel_tol, f_tol, max_num_steps);

      // Remove dose input in dosing compartment. Pred will add it
      // later, so we want to avoid redundancy.
      if (cmt <= 2) predPK(cmt - 1) -= amt;

    } else if (ii > 0) {  // multiple truncated infusions
      x_r[cmt - 1] = rate;

      // Construct augmented parameters
      ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
        theta = parameter.augment(predPK);
      theta.time(ii_dbl);

      predPD_guess = to_vector(integrator_(ode_rate_dbl_functor<F>(f_),
                                         to_array_1d(init_dbl),
                                         0.0, std::vector<double>(1, ii_dbl),
                                         unpromote(theta.get_RealParameters()),
                                         x_r, x_i)[0]);

      x_r.push_back(amt);  // needed?
      predPD = algebra_solver(system, predPD_guess,
                            to_vector(theta.get_RealParameters()),
                            x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    } else {  // constant infusion
      x_r[cmt - 1] = rate;

      // Construct augmented parameters
      ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
        theta = parameter.augment(predPK);
      theta.time(ii_dbl);

      predPD_guess = to_vector(integrator_(ode_rate_dbl_functor<F>(f_),
                                         to_array_1d(init_dbl),
                                         0.0, std::vector<double>(1, 100),
                                         unpromote(theta.get_RealParameters()),
                                         x_r, x_i)[0]);

      x_r.push_back(amt);
      predPD = algebra_solver(system, predPD_guess,
                              to_vector(theta.get_RealParameters()),
                              x_r, x_i,
                              0, rel_tol, f_tol, max_num_steps);
    }

    Matrix<scalar, Dynamic, 1> pred(nPK + nOde_);
    for (int i = 0; i < nPK; i++) pred(i) = predPK(i);
    for (int i = 0; i < nOde_; i++) pred(nPK + i) = predPD(i);

    return pred;
  }

  /**
   * Case 2 (vd): amt is random, rate is fixed.
   */
  template<typename T_time,
           typename T_ii,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename T_amt>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_ii, T_amt,
    T_parameters>::type, Eigen::Dynamic, 1>
  operator() (const ModelParameters<T_time,
           T_parameters,
           T_biovar,
           T_tlag>& parameter,
           const T_amt& amt,
           const double& rate,
           const T_ii& ii,
           const int& cmt) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using Eigen::VectorXd;
    using std::vector;
    using stan::math::algebra_solver;
    using stan::math::to_vector;
    using stan::math::to_array_1d;
    using stan::math::invalid_argument;

    typedef typename boost::math::tools::promote_args<T_ii, T_amt,
      T_parameters>::type scalar;

    double ii_dbl = unpromote(ii);

    // Compute solution for base 1cpt PK
    Matrix<scalar, Dynamic, 1> predPK;
    int nPK = 2;
    if (cmt <= 2) {  // check dosing occurs in a base state
      PredSS_oneCpt PredSS_one;
      int nParmsPK = 3;
      predPK = PredSS_one(parameter.truncate(nParmsPK),
                          amt, rate, ii_dbl, cmt);
      predPK(cmt - 1) = predPK(cmt - 1) + amt;
    } else {
      predPK = Matrix<scalar, Dynamic, 1>::Zero(nPK);
    }

    // Construct augmented parameters
    ModelParameters<T_time, scalar, T_biovar, T_tlag>
      theta = parameter.augment(predPK);
    theta = theta.augment(vector<T_amt>(1, amt));  // add amt
    theta.time(ii_dbl);

    // Arguments for ODE integrator (and initial guess)
    Matrix<double, 1, Dynamic> init_dbl
      = Matrix<double, 1, Dynamic>::Zero(nOde_);
    vector<double> x_r(nPK + nOde_, 0);  // rate for the full system
    x_r.push_back(0);  // include initial time (at SS, t0 = 0)
    vector<int> x_i;

    // Tuning parameters for algebraic solver
    double rel_tol = 1e-10;  // default
    double f_tol = 1e-4;  // empirical
    long int max_num_steps = 1e3;  // default // NOLINT

    SS_system_vd<ode_rate_dbl_functor<F> >
      system(ode_rate_dbl_functor<F>(f_), ii_dbl, cmt, integrator_, nPK);

    Matrix<double, Dynamic, 1> predPD_guess;
    Matrix<scalar, 1, Dynamic> predPD;

    if (rate == 0) {  // bolus dose
      predPD_guess = to_vector(integrator_(ode_rate_dbl_functor<F>(f_),
                                        to_array_1d(init_dbl),
                                        0.0, std::vector<double>(1, ii_dbl),
                                        unpromote(theta.get_RealParameters()),
                                        x_r, x_i)[0]);

      predPD = algebra_solver(system, predPD_guess,
                              to_vector(theta.get_RealParameters()),
                              x_r, x_i,
                              0, rel_tol, f_tol, max_num_steps);

      if (cmt <= 2) predPK(cmt - 1) -= amt;
    } else if (ii > 0) {
      invalid_argument("Steady State Event",
                       "Current version does not handle the case of",
                       "", " multiple truncated infusions ",
                       "(i.e ii > 0 and rate > 0) when F * amt is a parameter.");  // NOLINT
    } else {
      x_r[cmt - 1] = rate;

      predPD_guess = to_vector(integrator_(ode_rate_dbl_functor<F>(f_),
                                         to_array_1d(init_dbl),
                                         0.0, std::vector<double>(1, 100),
                                         unpromote(theta.get_RealParameters()),
                                         x_r, x_i)[0]);

      predPD = algebra_solver(system, predPD_guess,
                              to_vector(theta.get_RealParameters()),
                              x_r, x_i,
                              0, rel_tol, f_tol, max_num_steps);
    }

    Matrix<scalar, Dynamic, 1> pred(nPK + nOde_);
    for (int i = 0; i < nPK; i++) pred(i) = predPK(i);
    for (int i = 0; i < nOde_; i++) pred(nPK + i) = predPD(i);

    return pred;
  }
};

}
#endif
