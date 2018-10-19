#ifndef STAN_MATH_TORSTEN_REFACTOR_MIXODEONECPTMODEL_BDF_HPP
#define STAN_MATH_TORSTEN_REFACTOR_MIXODEONECPTMODEL_BDF_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/PKModel/functors/mix1_functor.hpp>
#include <stan/math/torsten/PKModel/Pred/Pred1_mix1.hpp>
#include <stan/math/torsten/PKModel/Pred/PredSS_mix1.hpp>
#include <stan/math/torsten/Pred2.hpp>
#include <stan/math/torsten/pk_coupled_model.hpp>
#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace torsten {

/**
 * Compute the predicted amounts in each compartment at each event
 * of an ODEs model. The model contains a base 1 Compartment PK
 * component which gets solved analytically, while the other ODEs
 * are solved numerically using stan::math::integrate_ode_bdf. This
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
mixOde1CptModel_bdf(const F& f,
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

  // check arguments
  static const char* function("mixOde1CptModel_bdf");
  torsten::pmetricsCheck(time, amt, rate, ii, evid, cmt, addl, ss,
                theta, biovar, tlag, function);

  // Construct dummy array of matrix for last argument of pred
  Matrix<T4, Dynamic, Dynamic> dummy_system;
  vector<Matrix<T4, Dynamic, Dynamic> >
    dummy_systems(1, dummy_system);

  typedef mix1_functor<F> F0;

  const int &nPK = refactor::PKOneCptModel<double, double, double, double>::Ncmt;
  PredWrapper<refactor::PkOneCptOdeModel> pr;
  PkOdeIntegrator<StanBdf> integrator(rel_tol, abs_tol, max_num_steps, msgs);

  Pred1_mix1<F0> pred1(F0(f), rel_tol, abs_tol, max_num_steps, msgs,
                       "bdf");
  PredSS_mix1<F0> predss(F0(f), rel_tol, abs_tol, max_num_steps, msgs,
                         "bdf", nOde);

#ifdef OLD_TORSTEN
  return Pred(time, amt, rate, ii, evid, cmt, addl, ss,
              theta, biovar, tlag, nPK + nOde, dummy_systems,
              pred1, predss);

#else
  return pr.Pred2(time, amt, rate, ii, evid, cmt, addl, ss,
                  theta, biovar, tlag, nPK + nOde, dummy_systems,
                  pred1, predss,
                  integrator,
                  f, nOde);
#endif
}

/**
 * Overload function to allow user to pass an std::vector for 
 * theta.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
  typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
                     const int nOde,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& theta,
                     const std::vector<std::vector<T5> >& biovar,
                     const std::vector<std::vector<T6> >& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_theta(1, theta);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_theta, biovar, tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

/**
 * Overload function to allow user to pass an std::vector for 
 * theta and biovar.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
                     const int nOde,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& theta,
                     const std::vector<T5>& biovar,
                     const std::vector<std::vector<T6> >& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_theta(1, theta);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_theta, vec_biovar, tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

/**
 * Overload function to allow user to pass an std::vector for 
 * theta, biovar, and tlag.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
                     const int nOde,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& theta,
                     const std::vector<T5>& biovar,
                     const std::vector<T6>& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_theta(1, theta);
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_theta, vec_biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

/**
 * Overload function to allow user to pass an std::vector for 
 * theta and tlag.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
                     const int nOde,
                     const std::vector<T0>& time,
                     const std::vector<T1>& amt,
                     const std::vector<T2>& rate,
                     const std::vector<T3>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<T4>& theta,
                     const std::vector<std::vector<T5> >& biovar,
                     const std::vector<T6>& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T4> > vec_theta(1, theta);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              vec_theta, biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

/**
 * Overload function to allow user to pass an std::vector for 
 * biovar.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
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
                     const std::vector<T5>& biovar,
                     const std::vector<std::vector<T6> >& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_biovar(1, biovar);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              theta, vec_biovar, tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

/**
 * Overload function to allow user to pass an std::vector for 
 * biovar and tlag.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
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
                     const std::vector<T5>& biovar,
                     const std::vector<T6>& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T5> > vec_biovar(1, biovar);
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              theta, vec_biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}


/**
 * Overload function to allow user to pass an std::vector for 
 * tlag.
 */
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename F>
Eigen::Matrix <typename boost::math::tools::promote_args<T0, T1, T2, T3,
  typename boost::math::tools::promote_args<T4, T5, T6>::type>::type,
  Eigen::Dynamic, Eigen::Dynamic>
mixOde1CptModel_bdf(const F& f,
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
                     const std::vector<T6>& tlag,
                     std::ostream* msgs = 0,
                     double rel_tol = 1e-6,
                     double abs_tol = 1e-6,
                     long int max_num_steps = 1e6) {  // NOLINT(runtime/int)
  std::vector<std::vector<T6> > vec_tlag(1, tlag);

  return mixOde1CptModel_bdf(f, nOde,
                              time, amt, rate, ii, evid, cmt, addl, ss,
                              theta, biovar, vec_tlag,
                              msgs, rel_tol, abs_tol, max_num_steps);
}

}
#endif
