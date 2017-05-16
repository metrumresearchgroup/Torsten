#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_MIX1_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_MIX1_HPP

#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <stan/math/torsten/PKModel/Pred/fOneCpt.hpp>
#include <iostream>
#include <vector>

/**
 *  Functor for mix solver with base
 *  one compartment model.
 */
template <typename F0>
struct mix1_functor {
  F0 f0_;

  mix1_functor(const F0 f0) : f0_(f0) { }

  /**
   *  Returns the derivative of the base ODE system. The base 1 PK
   *  component is calculated analytically. We use the last two
   *  elements of theta to store the inital PK states (init_PK) and
   *  the last element of x_r to store the initial time.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;
    typedef typename boost::math::tools::promote_args<T0, T2, T3>::type
      T_pk;  // return object of fOneCpt  doesn't depend on T1

    // Get PK parameters
    int nParmsPK = 3;
    std::vector<T2> thetaPK(nParmsPK);
    thetaPK[0] = theta[0];  // CL
    thetaPK[1] = theta[1];  // VC
    thetaPK[2] = theta[2];  // ka

    // Get initial PK states
    int nPK = 2;
    std::vector<T2> init_pk(nPK);
    size_t nTheta = theta.size();
    // The last two components of theta
    // should contain the initial PK states
    init_pk[0] = theta[nTheta - 2];
    init_pk[1] = theta[nTheta - 1];
    // Last element of x_r contains the initial time
    T0 dt = t - x_r[x_r.size() - 1];

    std::vector<T_pk> y_pk = fOneCpt(dt, thetaPK, init_pk, x_r);

    return f0_(dt, y, y_pk, theta, x_r, x_i, pstream_);
  }
};

/**
 *	ODE model with base 1 compartment PK. Use mix solver.
 *	Calculates the amount in each compartment at dt time units after the time
 *	of the initial state.
 *
 *	If the initial time equals the time of the event, the code does
 *	not run the mix solver, and sets the predicted amount equal to the
 *	initial condition. This can happen when we are dealing with events that
 *	occur simultaneously. The change to the predicted amount caused by bolus
 *	dosing events is handled later in the main Pred function.
 *
 *	 @tparam T_time type of scalar for time
 *	 @tparam T_rate type of scalar for rate
 *	 @tparam T_parameters type of scalar for model parameters
 *	 @tparam F type of ODE system function
 *	 @param[in] dt time between previous and current event
 *	 @param[in] parameter model parameters at current event
 *	 @param[in] init amount in each compartment at previous event
 *	 @param[in] rate rate in each compartment
 *	 @param[in] f functor for base ordinary differential equation that define
 *              the compartment model, passed the base 1 compartment PK.
 *   @return an eigen vector that contains predicted amount in each compartment
 *           at the current event.
 */
template<typename T_time,
         typename T_rate,
         typename T_parameters,
         typename T_biovar,
         typename T_tlag,
         typename T_system,
         typename F>
Eigen::Matrix<typename boost::math::tools::promote_args< T_time, T_rate,
  T_parameters>::type, 1, Eigen::Dynamic>
Pred1_mix1(const T_time& dt,
           const ModelParameters<T_time, T_parameters, T_biovar,
                                 T_tlag, T_system>& parameter,
           const Eigen::Matrix<typename boost::math::tools::promote_args<T_time,
                       T_rate, T_parameters>::type, 1, Eigen::Dynamic>& init,
           const std::vector<T_rate>& rate,
           const F& f) {
  using std::vector;
  using stan::math::to_array_1d;

  // FIX ME - should I revise the scalar time for T_biovar and T_lag
  // FIX ME - create other typedef to distinguish T_init and T_parms?
  typedef typename boost::math::tools::promote_args<T_time, T_rate,
    T_parameters>::type scalar;

  assert((size_t) init.cols() == rate.size());

  T_time t = parameter.get_time();  // time of current event
  T_time t0 = t - dt;  // time of previous event

  // Convert time parameters to fixed data for ODE integrator
  vector<double> t_dbl(1);
  t_dbl[0] = unpromote(t);
  double t0_dbl = unpromote(t0);
  vector<double> x_r(rate.size());  // FIX ME: currently store rate in x_r
  for (size_t i = 0; i < rate.size(); i++) x_r[i] = unpromote(rate[i]);
  x_r.push_back(t0_dbl);  // need to pass the initial time!

  vector<T_parameters> theta = parameter.get_RealParameters();
  vector<scalar> y0 = to_array_1d(init);

  Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
  if (t_dbl[0] == t0_dbl) { pred = init;
  } else {
    size_t nPK = 2;  // two states for 1Cpt with absorption
    // create vector with only PK initial states
    // (want to minmize the number of parameters that
    // get passed to a function)
    vector<scalar> y0_PK(nPK);
    y0_PK[0] = y0[0];
    y0_PK[1] = y0[1];

    // create vector with only PK parameters
    int nParmPK = 3;
    vector<scalar> thetaPK(nParmPK);
    thetaPK[0] = theta[0];  // CL
    thetaPK[1] = theta[1];  // VC
    thetaPK[2] = theta[2];  // ka

    vector<scalar> xPK = fOneCpt(dt, thetaPK, y0_PK, rate);
    for (size_t i = 0; i < nPK; i++) theta.push_back(init(i));

    // create vector with PD initial states
    size_t nPD = init.size() - nPK;
    vector<scalar> y0_PD(nPD);
    for (size_t i = 0; i < nPD; i++) y0_PD[i] = y0[nPK + i];
    vector<int> idummy;

    vector<vector<scalar> > pred_V = pmetrics_solver(mix1_functor<F>(f),
                                                     y0_PD, t0_dbl, t_dbl,
                                                     theta, x_r, idummy);
    size_t nOde = pred_V[0].size();

    pred.resize(nPK + nOde);
    for (size_t i = 0; i < nPK; i++) pred(i) = xPK[i];
    for (size_t i = 0; i < nOde; i++) pred(nPK + i) = pred_V[0][i];
  }
  return pred;
}

#endif
