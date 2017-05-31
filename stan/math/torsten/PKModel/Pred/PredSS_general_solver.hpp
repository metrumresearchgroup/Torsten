#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_GENERAL_SOLVER_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_GENERAL_SOLVER_HPP

#include <stan/math/torsten/PKModel/integrator.hpp>
#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/torsten/PKModel/functors/SS_system.hpp>
#include <stan/math/torsten/PKModel/Pred/pred1_general_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <iostream>

/**
 * General compartment model using built-in ODE solver
 * and root-finder.
 * Calculate amount in each compartment at the end of a
 * steady-state dosing interval or during a steady-state
 * constant input (if ii = 0)
 *
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_ii type of scalar for interdose interval
 *	 @tparam T_parameters type of scalar for model parameters
 *   @tparam T_biovar type of scalar for bio-availability
 *   @tparam T_system type of scalar for linear ODE matrix
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
         typename T_tlag,
         typename T_system,
         typename F>
Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_ii,
  T_parameters, typename boost::math::tools::promote_args<T_biovar,
  T_tlag>::type>::type, 1, Eigen::Dynamic>
PredSS_general_solver(const ModelParameters<T_time,
                                            T_parameters,
                                            T_biovar,
                                            T_tlag,
                                            T_system>& parameter,
                      const double& amt,
                      const double& rate,
                      const T_ii& ii,
                      const int& cmt,
                      const F& f,
                      const integrator_structure& integrator) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using Eigen::VectorXd;
  using std::vector;
  using stan::math::algebra_solver;
  using stan::math::to_vector;

  typedef typename boost::math::tools::promote_args<T_time, T_ii,
    T_parameters>::type scalar;

  Matrix<scalar, Dynamic, 1> pred;

  // Arguments for ODE integrator
  double ii_dbl = unpromote(ii);
  int nCmt = 2;  // DEV - how do I get nCmt? Use f...
  Matrix<double, 1, Dynamic> init_dbl(nCmt);
  for (int i = 0; i < nCmt; i++) init_dbl(i) = 0;
  vector<double> x_r(nCmt, 0);
  vector<int> x_i(0);

  // Arguments for algebraic solver
  Matrix<double, 1, Dynamic> y;
  double rel_tol = 1e-10;  // default
  double f_tol = 1e-4;  // empirical
  long int max_num_steps = 1e3;  // default

  // construct algebraic function
  SS_system<ode_rate_dbl_functor<F> >
    system(ode_rate_dbl_functor<F>(f), ii_dbl, cmt, integrator);

  if (rate == 0) {  // bolus dose
    // compute initial guess
    init_dbl(cmt - 1) = amt;
    y = Pred1_general_solver(ii_dbl, unpromote(parameter),
                             init_dbl, x_r, f, integrator);

    x_r.push_back(amt);

    pred = algebra_solver(system, y, 
                          to_vector(parameter.get_RealParameters()),
                          x_r, x_i,
                          0, rel_tol, f_tol, max_num_steps);

    // DEV - what tuning parameters should we use for the algebra solver?
    // DEV - update initial guess or tuning parameters if result not good?

  } else if (ii > 0) {  // multiple truncated infusions
    // compute initial guess
    x_r[cmt - 1] = rate;
    y = Pred1_general_solver(ii_dbl,
                             unpromote(parameter),
                             init_dbl, x_r, f, integrator);
    y << 0.0012932115491361106, 409.81573;

    x_r.push_back(amt);

    pred = algebra_solver(system, y,
                          to_vector(parameter.get_RealParameters()),
                          x_r, x_i,
                          0, rel_tol, 1e-3, max_num_steps);  // should use ftol
  } else {  // constant infusion
    x_r[cmt - 1] = rate;
    y = Pred1_general_solver(100.0, unpromote(parameter),
                             init_dbl, x_r, f, integrator);

    pred = algebra_solver(system, y,
                          to_vector(parameter.get_RealParameters()),
                          x_r, x_i,
                          0, rel_tol, f_tol, max_num_steps);
  }

  return pred;
}

#endif
