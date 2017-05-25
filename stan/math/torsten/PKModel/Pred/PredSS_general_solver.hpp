#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_GENERAL_SOLVER_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PREDSS_GENERAL_SOLVER_HPP

#include <stan/math/torsten/PKModel/integrator.hpp>
#include <stan/math/torsten/PKModel/Pred/pred1_general_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/prim/mat/fun/to_vector.hpp>
#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <iostream>

/**
 * A structure to store the algebraic system
 * which gets solved when computing the steady
 * state solution.
 */
template <typename F >
struct algebra_system {
  F f_;
  double ii_;
  integrator_structure integrator_;

  algebra_system() { }

  algebra_system(const F& f,
                 const double& ii,
                 const integrator_structure& integrator) {
    f_ = f;
    ii_ = ii;
    integrator_ = integrator;
  }

  template <typename T0, typename T1>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             const std::vector<double>& dat,
             const std::vector<int>& dat_int,
             std::ostream* msgs) const {
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    using Eigen::Matrix;
    using Eigen::Dynamic;
    typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;

    std::vector<double> ts(1, ii_);
    double t0 = 0;

    // Bolus dose 
    std::vector<T0> x0 = to_array_1d(x);
    int cmt = 0;  // FIX ME - need to pass cmt to the solver
    double amt = 1200;  // FIX ME - need to pass amt to the solver
    x0[cmt] += amt;  // initial state for evolution operator

    Matrix<scalar, Dynamic, 1>
      pred = to_vector(integrator_(f_, x0, t0, ts, to_array_1d(y), 
                                   dat, dat_int)[0]);

    Matrix<scalar, Dynamic, 1> result(x.size());
    for (int i = 0; i < result.size(); i++) result(i) = x(i) - pred(i);
    pred(cmt) -= amt;

    return result;
  }
};


/**
 * General compartment model using built-in ODE solver
 * and root-finder.
 * Calculate amount in each compartment at the end of a
 * steady-state dosing interval or during a steady-state
 * constant input (if ii = 0)
 *
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_amt type of scalar for amount
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_ii type of scalar for interdose interval
 *	 @tparam T_parameters type of scalar for model parameters
 *	 @tparam T_tlag type of scalar for lag time
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
         typename T_amt,
         typename T_rate,
         typename T_ii,
         typename T_parameters,
         typename T_biovar,
         typename T_tlag,
         typename T_system,
         typename F>
Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_amt, T_rate,
  typename boost::math::tools::promote_args< T_ii, T_parameters, T_biovar,
  T_tlag>::type>::type, 1, Eigen::Dynamic>
PredSS_general_solver(const ModelParameters<T_time,
                                            T_parameters,
                                            T_biovar,
                                            T_tlag,
                                            T_system>& parameter,
                      const T_amt& amt,
                      const T_rate& rate,
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

  typedef typename boost::math::tools::promote_args< T_time, T_amt, T_rate,
    typename boost::math::tools::promote_args< T_ii, T_parameters
                                               >::type>::type scalar;

  Matrix<scalar, Dynamic, 1> pred;

  double rate_dbl = unpromote(rate);
  vector<double> dat_real(1, rate_dbl);
  vector<int> dat_int(0);
  double ii_dbl = unpromote(ii);
  int nCmt = 2;  // FIX ME - arbitrary for now - how to get nCmt?
  VectorXd init(nCmt);
  for (int i = 0; i < nCmt; i++) init(i) = 0;
  vector<double> rate_v(nCmt, 0);
  rate_v[cmt] = rate_dbl;
  Matrix<double, Dynamic, 1> x;

  // Algebraic Solver parameters
  double rel_tol = 1e-10;  // default
  double f_tol = 1e-4;  // empirical upper bound
  long int max_num_steps = 1e3;  // default

  if (rate == 0) {  // bolus dose
    init(cmt - 1) = amt;
    x = Pred1_general_solver(ii_dbl, unpromote(parameter), init,
                             rate_v, f, integrator);  // initial guess
    x << 0.000668869, 384.736;

    pred = algebra_solver(algebra_system<F>(f, ii, integrator),
                          x, to_vector(parameter.get_RealParameters()),
                          dat_real, dat_int,
                          0,
                          rel_tol, f_tol, max_num_steps);
    
    // DEV - what tuning parameters should we use for the algebra solver?
    // DEV - update initial guess if result not good?

  } /* else if (ii > 0) {  // multiple truncated infusions
    x = pred1(amt / rate, unpromote(parameter.get_RealParameters()), init, rate_dbl);
    pred = algebra_solver(algebra_system<F>(f, ii_dbl),
                          x, parameter.get_RealParameters(),
                          rate_dbl, dat_int);
  } else {  // constant infusion
    Eigen::Matrix<double, Eigen::Dynamic, 1> x;  // DEV - need an initial guess
    // DEV -- do I need to agument f with R, or is it already part of the equation,
    // via rate_d?
    pred = algebra_solver(f, x, parameter.get_RealParameters(), rate_dbl, dat_int);
  }
    */
  return pred;
}

#endif
