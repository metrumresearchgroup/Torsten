#ifndef STAN_MATH_TORSTEN_MIXODEONECPTMODEL_RK45_HPP
#define STAN_MATH_TORSTEN_MIXODEONECPTMODEL_RK45_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

/**
 * Compute the predicted amounts in each compartment at each event
 * of an ODEs model. The model contains a base 1 Compartment PK
 * component which gets solved analytically, while the other ODEs
 * are solved numerically using stan::math::integrate_ode_rk45. This
 * amounts to using the mixed solver method.
 *
 * <b>Warning:</b> This prototype does not handle steady state events.
 *
 * @tparam T0 type of scalar for time of events.
 * @tparam T1 type of scalar for amount at each event.
 * @tparam T2 type of scalar for rate at each event.
 * @tparam T3 type of scalar for inter-dose inteveral at each event.
 * @tparam T4 type of scalars for the model parameters.
 * @tparam T5 type of scalars for the bio-variability parameters.
 * @tparam T6 type of scalars for the model tlag parameters.
 * @tparam F type of ODE system function.
 * @param[in] f functor for base ordinary differential equation
 *            which gets solved numerically.
 * @param[in] nOde number of ODEs we solve numerically.
 * @param[in] time times of events
 * @param[in] amt amount at each event
 * @param[in] rate rate at each event
 * @param[in] ii inter-dose interval at each event
 * @param[in] evid event identity:
 *                    (0) observation
 *                    (1) dosing
 *                    (2) other
 *                    (3) reset
 *                    (4) reset AND dosing
 * @param[in] cmt compartment number at each event
 * @param[in] addl additional dosing at each event
 * @param[in] ss steady state approximation at each event (0: no, 1: yes)
 * @param[in] theta vector of ODE parameters
 * @param[in] biovar bio-availability in each compartment
 * @param[in] tlag lag time in each compartment
 * @param[in] rel_tol relative tolerance for the Boost ode solver
 * @param[in] abs_tol absolute tolerance for the Boost ode solver
 * @param[in] max_num_steps maximal number of steps to take within
 *            the Boost ode solver
 * @return a matrix with predicted amount in each compartment
 *         at each event.
 *
 * FIX ME: msg should be passed on to functor (allows use of
 * print statement inside ODE system).
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOdeOneCptModel_rk45(const F& f,
                       const int nOde,
                       const std::vector<T0>& time,
                       const std::vector<T1>& amt,
                       const std::vector<T2>& rate,
                       const std::vector<T3>& ii,
                       const std::vector<int>& evid,
                       const std::vector<int>& cmt,
                       const std::vector<int>& addl,
                       const std::vector<int>& ss,
                       const std::vector<std::vector<T4> >& theta,
                       const std::vector<std::vector<T5> >& biovar,
                       const std::vector<std::vector<T6> >& tlag,
                       std::ostream* msgs = 0,
                       double rel_tol = 1e-6,
                       double abs_tol = 1e-6,
                       long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  using std::vector;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using boost::math::tools::promote_args;

  int nPK = 2;
  PKModel model(theta[0].size(), nOde + nPK);

  // check arguments
  static const char* function("mixOdeOneCptModel_rk45");
  pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
                theta, biovar, tlag, function, model);

  // define functors used in Pred()
  Pred1_structure new_Pred1("mixOdeOneCptModel");
  Pred1 = new_Pred1;
  PredSS_structure new_PredSS("error");
  PredSS = new_PredSS;  // WARNING: PredSS returns an error
  pmetrics_solver_structure new_pmetrics_solver(rel_tol, abs_tol,
    max_num_steps, "rk45");
  pmetrics_solver = new_pmetrics_solver;

  // Construct dummy matrix for last argument of pred
  Matrix<double, Dynamic, Dynamic> dummy_system;
  vector<Matrix<double, Dynamic, Dynamic> >
    dummy_systems(1, dummy_system);

 return  pred = Pred(time, amt, rate, ii, evid, cmt, addl, ss,
                     theta, biovar, tlag, model,
                     f, dummy_systems);
}

#endif
