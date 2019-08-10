#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_MIX2_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_MIX2_HPP

#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <stan/math/torsten/PKModel/Pred/fTwoCpt.hpp>
#include <stan/math/torsten/PKModel/integrator.hpp>
#include <iostream>
#include <vector>
#include<string>

namespace torsten {

template <typename F>
struct Pred1_mix2 {
  F f_;
  integrator_structure integrator_;

  Pred1_mix2(const F& f,
             const double& rel_tol,
             const double& abs_tol,
             const long int& max_num_steps,  // NOLINT
             std::ostream* msgs,
             const std::string& integratorType)
    : f_(f),
      integrator_(rel_tol, abs_tol, max_num_steps, msgs, integratorType) { }

  Pred1_mix2(const F& f,
             const integrator_structure& integrator)
    : f_(f),
      integrator_(integrator) { }

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
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename T_init>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_time,
    T_parameters, T_init>::type, 1, Eigen::Dynamic>
  operator() (const T_time& dt,
              const ModelParameters<T_time, T_parameters, T_biovar,
                                    T_tlag>& parameter,
              const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& init,
              const std::vector<double>& rate) const {
    using std::vector;
    using stan::math::to_array_1d;
    using boost::math::tools::promote_args;

    typedef typename promote_args<T_time, T_parameters, T_init>::type scalar;
    typedef typename promote_args<T_parameters, T_init>::type T_theta;

    assert((size_t) init.cols() == rate.size());

    // pass fixed times to the integrator. FIX ME - see issue #30
    T_time t = parameter.get_time();  // time of current event
    T_time t0 = t - dt;  // time of previous event
    vector<double> t_dbl(1, unpromote(t));
    double t0_dbl = unpromote(t0);
    vector<double> x_r = rate;
    x_r.push_back(t0_dbl);  // need to pass the initial time!

    vector<scalar> y0 = to_array_1d(init);

    size_t nParm = parameter.get_RealParameters(false).size();
    vector<T_theta> theta(nParm);
    for (size_t i = 0; i < nParm; i++)
      theta[i] = parameter.get_RealParameters(false)[i];

    Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
    if (t_dbl[0] == t0_dbl) { pred = init;
    } else {
      size_t nPK = 3;  // three states for 1Cpt with absorption
      vector<scalar> y0_PK(nPK);
      y0_PK[0] = y0[0];
      y0_PK[1] = y0[1];
      y0_PK[2] = y0[2];

      int nParmPK = 5;
      vector<scalar> thetaPK(nParmPK);
      thetaPK[0] = theta[0];  // CL
      thetaPK[1] = theta[1];  // Q
      thetaPK[2] = theta[2];  // VC
      thetaPK[3] = theta[3];  // VP
      thetaPK[4] = theta[4];  // ka

      vector<scalar> xPK = fTwoCpt(dt, thetaPK, y0_PK, rate);

      // Add PK inits to theta
      for (size_t i = 0; i < nPK; i++) theta.push_back(init(i));

      // create vector with PD initial states
      size_t nPD = init.size() - nPK;
      vector<scalar> y0_PD(nPD);
      for (size_t i = 0; i < nPD; i++) y0_PD[i] = y0[nPK + i];
      vector<int> idummy;

      vector<vector<scalar> >
        pred_V = integrator_(ode_rate_dbl_functor<F>(f_),
                             y0_PD, t0_dbl, t_dbl,
                             theta, x_r, idummy);
      size_t nOde = pred_V[0].size();

      pred.resize(nPK + nOde);
      for (size_t i = 0; i < nPK; i++) pred(i) = xPK[i];
      for (size_t i = 0; i < nOde; i++) pred(nPK + i) = pred_V[0][i];
    }
    return pred;
  }

  /**
  * Overload function for case rate is a vector of var.
  * The parameters are stored in theta in this order:
  *   (1) ODE parameters
  *   (2) rate
  *   (3) initial states for base PK
  *
  * Unlike the general solver (pred1_general_solver), we'll
  * call the mix_rate_var_functor, which will know where rate
  * is located inside theta.
  */
  template<typename T_time,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename T_init,
           typename T_rate>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_time,
    T_rate, T_parameters, T_init>::type, 1, Eigen::Dynamic>
  operator() (const T_time& dt,
              const ModelParameters<T_time, T_parameters, T_biovar,
                                    T_tlag>& parameter,
              const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& init,
              const std::vector<T_rate>& rate) const {
    using std::vector;
    using stan::math::to_array_1d;
    using boost::math::tools::promote_args;

    typedef typename promote_args<T_time, T_rate,
      T_parameters, T_init>::type scalar;

    typedef typename promote_args<T_parameters, T_init, T_rate>::type T_theta;

    assert((size_t) init.cols() == rate.size());

    // pass fixed times to the integrator. FIX ME - see issue #30
    T_time t = parameter.get_time();
    T_time t0 = t - dt;
    vector<double> t_dbl(1, unpromote(t));
    double t0_dbl = unpromote(t0);

    vector<scalar> y0 = to_array_1d(init);

    size_t nParm = parameter.get_RealParameters(false).size();
    vector<T_theta> theta(nParm);
    for (size_t i = 0; i < nParm; i++)
      theta[i] = parameter.get_RealParameters(false)[i];

    Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
    if (t_dbl[0] == t0_dbl) { pred = init;
    } else {
      size_t nPK = 3;  // three states for 2Cpt with absorption
      vector<scalar> y0_PK(nPK);
      y0_PK[0] = y0[0];
      y0_PK[1] = y0[1];
      y0_PK[2] = y0[2];

      int nParmPK = 5;
      vector<scalar> thetaPK(nParmPK);
      thetaPK[0] = theta[0];  // CL
      thetaPK[1] = theta[1];  // Q
      thetaPK[2] = theta[2];  // VC
      thetaPK[3] = theta[3];  // VP
      thetaPK[4] = theta[4];  // ka

      vector<scalar> xPK = fTwoCpt(dt, thetaPK, y0_PK, rate);

      // Add rate and PK inits IN THIS ORDER to theta
      for (size_t i = 0; i < rate.size(); i++) theta.push_back(rate[i]);
      for (size_t i = 0; i < nPK; i++) theta.push_back(init(i));

      // create vector with PD initial states
      size_t nPD = init.size() - nPK;
      vector<scalar> y0_PD(nPD);
      for (size_t i = 0; i < nPD; i++) y0_PD[i] = y0[nPK + i];
      vector<int> idummy;

      vector<double> x_r(2);
      x_r[0] = nPK;
      x_r[1] = t0_dbl;

      vector<vector<scalar> >
        pred_V = integrator_(ode_rate_var_functor<F>(f_),
                             y0_PD, t0_dbl, t_dbl,
                             theta, x_r, idummy);
      size_t nOde = pred_V[0].size();

      pred.resize(nPK + nOde);
      for (size_t i = 0; i < nPK; i++) pred(i) = xPK[i];
      for (size_t i = 0; i < nOde; i++) pred(nPK + i) = pred_V[0][i];
    }
    return pred;
  }
};

}

#endif
