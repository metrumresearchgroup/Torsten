#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_GENERAL_SOLVER_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_PRED1_GENERAL_SOLVER_HPP

#include <stan/math/torsten/PKModel/integrator.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>
#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <stan/math/torsten/univariate_integral.hpp>

namespace torsten{

template <typename F>
struct Pred1_general {
  F f_;
  integrator_structure integrator_;

  Pred1_general(const F& f,
                const double& rel_tol,
                const double& abs_tol,
                const long int& max_num_steps,  // NOLINT
                std::ostream* msgs,
                const std::string& integratorType)
    : f_(f),
      integrator_(rel_tol, abs_tol, max_num_steps, msgs, integratorType) { }

  Pred1_general(const F& f,
                const integrator_structure& integrator)
    : f_(f), integrator_(integrator) { }

  /**
   *	General compartment model using the built-in ODE solver.
   *	Calculates the amount in each compartment at dt time units after the time
   *	of the initial condition.
   *
   *	If the initial time equals the time of the event, than the code does
   *	not run the ode integrator, and sets the predicted amount equal to the
   *	initial condition. This can happen when we are dealing with events that
   *	occur simultaneously. The change to the predicted amount caused by bolus
   *	dosing events is handled later in the main Pred function.
   *
   *  The function is overloaded for the cases where rate is a vector of double
   *  or var.
   *
   *	 @tparam T_time type of scalar for time
   *	 @tparam T_parameters type of scalar for Ode parameters in ModelParameters.
   *   @tparam T_biovar type of scalar of biovar in ModelParameters.
   *   @tparam T_tlag type of scalar of lag times in ModelParameters.
   *   @tparam T_init type of scalar for the initial state
   *	 @param[in] dt time between current and previous event
   *	 @param[in] parameter model parameters at current event
   *	 @param[in] init amount in each compartment at previous event
   *	 @param[in] rate rate in each compartment (here as fixed data)
   *	 @param[in] f functor for base ordinary differential equation that defines
   *              compartment model.
   *   @return an eigen vector that contains predicted amount in each compartment
   *           at the current event.
   */
  template<typename T_time,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename T_init>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_init,
    T_parameters>::type, Eigen::Dynamic, 1>
  operator() (const T_time& dt,
              const ModelParameters<T_time, T_parameters, T_biovar,
                                    T_tlag>& parameter,
              const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& init,
              const std::vector<double>& rate) const {
    using stan::math::to_array_1d;
    using std::vector;
    using boost::math::tools::promote_args;

    typedef typename promote_args<T_time, T_init, T_parameters>::type scalar;

    assert((size_t) init.cols() == rate.size());

    T_time EventTime = parameter.get_time();  // time of current event
    T_time InitTime = EventTime - dt;  // time of previous event

    // Convert time parameters to fixed data for ODE integrator
    // FIX ME - see issue #30
    vector<double> EventTime_d(1, unpromote(EventTime));
    double InitTime_d = unpromote(InitTime);

    vector<T_parameters> theta = parameter.get_RealParameters();
    vector<scalar> init_vector = to_array_1d(init);

    Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
    if (EventTime_d[0] == InitTime_d) { pred = init;
    } else {
      vector<int> idummy;
      vector<vector<scalar> >
        pred_V = integrator_(ode_rate_dbl_functor<F>(f_),
                             init_vector, InitTime_d,
                             EventTime_d, theta, rate,
                             idummy);

      // Convert vector in row-major vector (eigen Matrix)
      pred.resize(pred_V[0].size());
      for (size_t i = 0; i < pred_V[0].size(); i++) pred(i) = pred_V[0][i];
    }
    return pred;
  }

  /**
  * Overload function for case where rate is a vector of var.
  * That occurs either rate is passed as a parameter, or F,
  * the bio-availibility factor is a parameter -- thus making
  * the rate vector that gets passed a latent parameter.
  */
  template<typename T_time,
           typename T_parameters,
           typename T_biovar,
           typename T_tlag,
           typename T_init,
           typename T_rate>
  Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_init,
    T_parameters, T_rate>::type, Eigen::Dynamic, 1>
  operator() (const T_time& dt,
              const ModelParameters<T_time, T_parameters, T_biovar,
                                    T_tlag>& parameter,
              const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& init,
              const std::vector<T_rate>& rate) const {
    using stan::math::to_array_1d;
    using std::vector;
    using boost::math::tools::promote_args;

    typedef typename promote_args<T_time, T_init, T_parameters,
                                  T_rate>::type scalar;

    assert((size_t) init.cols() == rate.size());

    T_time EventTime = parameter.get_time();  // time of current event
    T_time InitTime = EventTime - dt;  // time of previous event

    // Convert time parameters to fixed data for ODE integrator
    // FIX ME - see issue #30
    vector<double> EventTime_d(1, unpromote(EventTime));
    double InitTime_d = unpromote(InitTime);

    // Construct theta with ode parameters and rates.
    vector<T_parameters> odeParameters = parameter.get_RealParameters();
    size_t nOdeParm = odeParameters.size();
    vector<typename promote_args<T_parameters, T_rate, T_time>::type>
      theta(nOdeParm + rate.size());
    for (size_t i = 0; i < nOdeParm; i++) theta[i] = odeParameters[i];
    for (size_t i = 0; i < rate.size(); i++)
      theta[nOdeParm + i] = rate[i];
    theta.push_back(InitTime);
    theta.push_back(EventTime);

    vector<scalar> init_vector = to_array_1d(init);
    vector<double> x_r;

    Eigen::Matrix<scalar, 1, Eigen::Dynamic> pred;
    ode_rate_var_functor<F> f0(f_);
    normalized_integrand_adaptor<ode_rate_var_functor<F>> f1(f0);
    vector<int> idummy;

    // if (EventTime_d[0] == InitTime_d) { pred = init;
    if (dt == 0) {
      std::vector<scalar> rhs {f0(EventTime, init_vector, theta, x_r, idummy, nullptr)};
      Eigen::Matrix<scalar, 1, Eigen::Dynamic> mat_rhs = stan::math::to_matrix(rhs, 1, rhs.size());
      pred = init + dt * mat_rhs;
    } else {
      // using normalized functor
      InitTime_d = 0.0;
      EventTime_d[0] = 1.0;
      vector<vector<scalar> >
        pred_V = integrator_(f1,
                             init_vector, InitTime_d,
                             EventTime_d, theta, x_r,
                             idummy);

      // Convert vector in row-major vector (eigen Matrix)
      // FIX ME - want to return column-major vector to use Stan's
      // to_vector function.
      pred.resize(pred_V[0].size());
      for (size_t i = 0; i < pred_V[0].size(); i++) pred(0, i) = pred_V[0][i];
    }
    return pred;
  }
};

}

#endif
