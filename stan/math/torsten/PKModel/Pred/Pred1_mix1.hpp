#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_MIX1_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_MIX1_HPP

#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <stan/math/torsten/PKModel/Pred/pred1_one.hpp>
#include <iostream>
#include <vector>

/**
 *  Functor for mix solver.
 */
template <typename F0>
struct mix_functor {
  F0 f0_;

  mix1_functor(const F0 f0) : f0_(f0) { }

  /**
   *  Returns the derivative of the base ODE system. The base 1 PK
   *  component is calculated analytically. We use the last two
   *  elements of theta to store the inital PK states (init_PK) and
   *  the last element of x_r to store the initial time.
   */
  template <typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& theta,
             const std::vector<T3>& rate,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    // PK variables
    int nPK = 2;
    std::vector<T2> init_pk(nPK);
    size_t nTheta = theta.size();
    init_pk[0] = theta(nTheta - 1);
    init_pk[1] = theta(nTheta);
    double dt = t - x_r[x_r.size() - 1];
    std::vector<T2> x_pk(nPK) = stan::math::to_array_1d(pred1_one(dt, theta, init_pk, rate));

    return f0_(dt, x, theta, x_r, x_i, pstream_);
  }
};

/**
 *	Ode model with base 1 compartment PK. Use mix solver.
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

  typedef typename boost::math::tools::promote_args<T_time, T_rate,
    T_parameters>::type scalar;
  assert((size_t) init.cols() == rate.size());

  T_time t = parameter.get_time();  // time of current event
  T_time t0 = EventTime - dt;  // time of previous event

  // Convert time parameters to fixed data for ODE integrator
  vector<double> t_dbl(1);
  t_dbl[0] = unpromote(EventTime);
  double t0_dbl = unpromote(t0);
  vector<double> x_r(rate.size);  // FIX ME: currently store rate in x_r
  for (size_t i = 0; i < rate.size(); i++) x_r[i] = unpromote(rate[i]);
  x_r.push_back(t0_dbl);  // need to pass the initial time!

  vector<T_parameters> theta = parameter.get_RealParameters();
  vector<scalar> x0 = to_array_1d(init);

  Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
  if (EventTime_d[0] == InitTime_d) { pred = init;
  } else {
    vector<scalar> xPK = to_array_1d(pred1_one(dt, parameter, init, rate));
    size_t nPK = xPK.size();
    for (size_t i = 0; i < nPK; i++) {
      pred(i) = xPK[i];
      theta.push_back(xPK[i]);
    }

    mix1_functor mix1_f(f);

    vector<int> idummy;
    vector<vector<scalar> > pred_V = pmetrics_solver(mix1_f(), x0, t0_dbl, t_dbl, theta,
                                                     rate_dbl, idummy);

    size_t nOde = pred_V[0].size();
    pred.resize(nPK + nOde);
    for (size_t i = 0; i < nOde; i++) pred(nPK + i) = pred_V[0][i];
  }
  return pred;
}


#endif
