#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_GENERAL_SOLVER_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_GENERAL_SOLVER_HPP

#include <stan/math/torsten/PKModel/integrator.hpp>
#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/torsten/PKModel/functors/SS_system.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_general.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_void.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <vector>
#include <iostream>
#include <string>

namespace torsten {

template <typename F>
struct PredSS_general {
  F f_;
  integrator_structure integrator_;
  int nCmt_;

  PredSS_general(const F& f,
                 const double& rel_tol,
                 const double& abs_tol,
                 const long int& max_num_steps,  // NOLINT
                 std::ostream* msgs,
                 const std::string& integratorType,
                 const int& nCmt)
    : f_(f),
      integrator_(rel_tol, abs_tol, max_num_steps, msgs, integratorType),
      nCmt_(nCmt) { }

  /**
   * General compartment model using built-in ODE solver
   * and root-finder.
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

    typedef typename boost::math::tools::promote_args<T_ii,
      T_parameters>::type scalar;

    Matrix<scalar, Dynamic, 1> pred;

    // Arguments for ODE integrator (and initial guess)
    double ii_dbl = unpromote(ii);
    Matrix<double, 1, Dynamic> init_dbl(nCmt_);
    for (int i = 0; i < nCmt_; i++) init_dbl(i) = 0;
    vector<double> x_r(nCmt_, 0);
    vector<int> x_i(0);

    // Arguments for algebraic solver
    Matrix<double, Dynamic, 1> y;
    double rel_tol = 1e-10;  // default
    double f_tol = 1e-4;  // empirical
    long int max_num_steps = 1e3;  // default // NOLINT

    // construct algebraic function
    SS_system_dd<ode_rate_dbl_functor<F>, Pred1_void>
      system(ode_rate_dbl_functor<F>(f_), Pred1_void(),
             ii_dbl, cmt, integrator_);

    // Construct Pred1_general functor
    Pred1_general<F> Pred1(f_, integrator_);

    if (rate == 0) {  // bolus dose
      // compute initial guess
      init_dbl(cmt - 1) = amt;
      y = Pred1(ii_dbl, unpromote(parameter), init_dbl, x_r);
      x_r.push_back(amt);
      pred = algebra_solver(system, y,
                            to_vector(parameter.get_RealParameters()),
                            x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
      // DEV - what tuning parameters should we use for the algebra solver?
      // DEV - update initial guess or tuning parameters if result not good?
    }  else if (ii > 0) {  // multiple truncated infusions
      x_r[cmt - 1] = rate;
      // compute initial guess
      y = Pred1(ii_dbl, unpromote(parameter), init_dbl, x_r);
      x_r.push_back(amt);
      pred = algebra_solver(system, y,
                            to_vector(parameter.get_RealParameters()),
                            x_r, x_i,
                            0, rel_tol, 1e-3, max_num_steps);  // FIX ME
                                                               // use ftol
    } else {  // constant infusion
      x_r[cmt - 1] = rate;
      y = Pred1(100.0, unpromote(parameter), init_dbl, x_r);

      x_r.push_back(amt);
      pred = algebra_solver(system, y,
                            to_vector(parameter.get_RealParameters()),
                            x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    }

    return pred;
  }

  /**
   * Case 2 (vd): amt is random, rate is fixed.
   */
  template<typename T_time,
           typename T_amt,
           typename T_ii,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_ii, T_parameters,
    T_amt>::type, Eigen::Dynamic, 1>
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
    using stan::math::invalid_argument;

    typedef typename boost::math::tools::promote_args<T_time, T_ii,
      T_parameters, T_amt>::type scalar;

    Matrix<scalar, Dynamic, 1> pred;

    // Arguments for the ODE integrator
    double ii_dbl = unpromote(ii);
    Matrix<double, 1, Dynamic> init_dbl(nCmt_);
    for (int i = 0; i < nCmt_; i++) init_dbl(i) = 0;
    vector<double> x_r(nCmt_, 0);
    vector<int> x_i(0);

    // Arguments for algebraic solver
    Matrix<double, 1, Dynamic> y;
    double rel_tol = 1e-10;  // default
    double f_tol = 5e-4;  // empirical (note: differs from other function)
    long int max_num_steps = 1e4;  // default  // NOLINT

    // construct algebraic function
    SS_system_vd<ode_rate_dbl_functor<F> >
      system(ode_rate_dbl_functor<F>(f_), ii_dbl, cmt, integrator_);

    // Construct Pred1_general functor
    Pred1_general<F> Pred1(f_, integrator_);

    int nParameters = parameter.get_RealParameters().size();
    Matrix<scalar, Dynamic, 1> parms(nParameters + 1);
    for (int i = 0; i < nParameters; i++)
      parms(i) = parameter.get_RealParameters()[i];
    parms(nParameters) = amt;

    if (rate == 0) {  // bolus dose
      // compute initial guess
      init_dbl(cmt - 1) = unpromote(amt);
      y = Pred1(ii_dbl, unpromote(parameter), init_dbl, x_r);

      pred = algebra_solver(system, y, parms, x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    }  else if (ii > 0) {  // multiple truncated infusions
      // compute initial guess
      x_r[cmt - 1] = rate;
      y = Pred1(ii_dbl, unpromote(parameter), init_dbl, x_r);

      pred = algebra_solver(system, y, parms, x_r, x_i,
                            0, rel_tol, 1e-3, max_num_steps);  // use ftol
    } else {  // constant infusion
      x_r[cmt - 1] = rate;
      y = Pred1(100.0, unpromote(parameter), init_dbl, x_r);

      pred = algebra_solver(system, y, parms, x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    }

    return pred;
  }
};

}

#endif
