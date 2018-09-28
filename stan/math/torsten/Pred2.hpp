#ifndef STAN_MATH_TORSTEN_PKMODEL_REFACTOR_PRED_HPP
#define STAN_MATH_TORSTEN_PKMODEL_REFACTOR_PRED_HPP

#include <Eigen/Dense>
#include <vector>
#include <stan/math/torsten/pk_twocpt_model.hpp>
#include <stan/math/torsten/pk_onecpt_model.hpp>
#include <stan/math/torsten/pk_linode_model.hpp>


namespace torsten{
  namespace internal {
    /** 
     * Different model's @c solve functon has different
     * signature, so we need a helper struct. 
     */
    template<PkOdeIntegratorId It>
    struct SolverDispatcher
    {
      const double rtol;
      const double atol;
      const int max_num_steps;
      std::ostream* msgs;

      SolverDispatcher(const double rtol0,
                       const double atol0,
                       const int max_num_steps0,
                       std::ostream* msgs0) :
        rtol(rtol0), atol(atol0),
        max_num_steps(max_num_steps0),
        msgs(msgs0)
      {}

      /*
       * dispatch @c PKODEModel solver
       */
      template<typename T0, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKODEModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKODEModel<Ts...>& model, const T0& dt) { //NOLINT
        return model.template solve<It>(dt, rtol, atol, max_num_steps, msgs);
      }
      
      /*
       * dispatch @c PKOneCptmodel solver
       */
      template<typename T0, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKOneCptModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKOneCptModel<Ts...>& model, const T0& dt) { //NOLINT
        return model.solve(dt);
      }

      /*
       * dispatch @c PKTwoptmodel solver
       */
      template<typename T0, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKTwoCptModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKTwoCptModel<Ts...>& model, const T0& dt) { //NOLINT
        return model.solve(dt);
      }

      /*
       * dispatch @c PKLinodemodel solver
       */
      template<typename T0, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKLinODEModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKLinODEModel<Ts...>& model, const T0& dt) { //NOLINT
        return model.solve(dt);
      }

      /*
       * dispatch @c PKOneCptmodel steady-state solver
       */
      template<typename T_amt, typename T_r, typename T_ii, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKOneCptModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKOneCptModel<Ts...>& model,
            const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) { //NOLINT
        return model.solve(amt, rate, ii, cmt);
      }
      
      /*
       * dispatch @c PKTwoCptmodel steady-state solver
       */
      template<typename T_amt, typename T_r, typename T_ii, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKTwoCptModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKTwoCptModel<Ts...>& model,
            const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) { //NOLINT
        return model.solve(amt, rate, ii, cmt);
      }

      /*
       * dispatch @c PKLinODEModel steady-state solver
       */
      template<typename T_amt, typename T_r, typename T_ii, typename... Ts> //NOLINT
      Eigen::Matrix<scalar_t<refactor::PKLinODEModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKLinODEModel<Ts...>& model,
            const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) { //NOLINT
        return model.solve(amt, rate, ii, cmt);
      }

      /*
       * dispatch @c PKODEmodel steady-state solver.
       *  @c amt and @c rate are data
       */
      template<typename T_ii, typename... Ts>
      Eigen::Matrix<scalar_t<refactor::PKODEModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKODEModel<Ts...>& model,
            const double& amt, const double& rate, const T_ii& ii, const int& cmt) { //NOLINT
        return model.template solve<It, T_ii>(amt, rate, ii, cmt, rtol, atol, max_num_steps, msgs); //NOLINT
      }

      /*
       * dispatch @c PKODEmodel steady-state solver.
       * @c amt is param and @c rate is data
       */
      template<typename T_amt, typename T_ii, typename... Ts>
      Eigen::Matrix<scalar_t<refactor::PKODEModel<Ts...>>, Eigen::Dynamic, 1> //NOLINT
      solve(const refactor::PKODEModel<Ts...>& model,
            const T_amt& amt, const double& rate, const T_ii& ii, const int& cmt) { //NOLINT
        return model.template solve<It, T_amt, T_ii>(amt, rate, ii, cmt, rtol, atol, max_num_steps, msgs); //NOLINT
      }

      /*
       * dispatch @c PKODEmodel steady-state solver.
       * @c amt is param and @c rate is data
       */


    };
  }

  /*
   * the wrapper is aware of @c T_model so it build model
   * accordingly.
   */
  template<template<typename...> class T_model, PkOdeIntegratorId It=StanRk45>
  struct PredWrapper{
    internal::SolverDispatcher<It> sol;

    PredWrapper() : sol(0.0, 0.0, 0, nullptr)
    {}

    PredWrapper(const double rtol0,
                const double atol0,
                const int max_num_steps0,
                std::ostream* msgs0) :
      sol(rtol0, atol0, max_num_steps0, msgs0)
    {}
    
    /** 
     * Different model's @c solve functon has different
     * signature, so we need a helper function. 
     */

    /** 
     * Different model's @c solve functon has different
     * signature, so we need a helper function. 
     */

    /**
     * Every Torsten function calls Pred.
     *
     * Predicts the amount in each compartment for each event,
     * given the event schedule and the parameters of the model.
     *
     * Proceeds in two steps. First, computes all the events that
     * are not included in the original data set but during which
     * amounts in the system get updated. Secondly, predicts
     * the amounts in each compartment sequentially by going
     * through the augmented schedule of events. The returned pred
     * Matrix only contains the amounts in the event originally
     * specified by the users.
     *
     * This function is valid for all models. What changes from one
     * model to the other are the Pred1 and PredSS functions, which
     * calculate the amount at an individual event.
     *
     * @tparam T_time type of scalar for time
     * @tparam T_amt type of scalar for amount
     * @tparam T_rate type of scalar for rate
     * @tparam T_ii type of scalar for interdose interval
     * @tparam T_parameters type of scalar for the ODE parameters
     * @tparam T_biovar type of scalar for bio-variability parameters
     * @tparam T_tlag type of scalar for lag times parameters
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
     * @param[in] cmt compartment number at each event (starts at 1)
     * @param[in] addl additional dosing at each event
     * @param[in] ss steady state approximation at each event
     * (0: no, 1: yes)
     * @param[in] pMatrix parameters at each event
     * @param[in] addParm additional parameters at each event
     * @parem[in] model basic info for ODE model and evolution operators
     * @param[in] SystemODE matrix describing linear ODE system that
     * defines compartment model. Used for matrix exponential solutions.
     * Included because it may get updated in modelParameters.
     * @return a matrix with predicted amount in each compartment
     * at each event.
     */
    template<typename T_time,
             typename T_amt,
             typename T_rate,
             typename T_ii,
             typename T_parameters,
             typename T_biovar,
             typename T_tlag,
             typename F_one,
             typename F_SS,
             typename... Ts>
    Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_amt, T_rate,
                                                            T_ii, typename boost::math::tools::promote_args<T_parameters, T_biovar,
                                                                                                            T_tlag>::type >::type, Eigen::Dynamic, Eigen::Dynamic>
    Pred2(const std::vector<T_time>& time,
          const std::vector<T_amt>& amt,
          const std::vector<T_rate>& rate,
          const std::vector<T_ii>& ii,
          const std::vector<int>& evid,
          const std::vector<int>& cmt,
          const std::vector<int>& addl,
          const std::vector<int>& ss,
          const std::vector<std::vector<T_parameters> >& pMatrix,
          const std::vector<std::vector<T_biovar> >& biovar,
          const std::vector<std::vector<T_tlag> >& tlag,
          const int& nCmt,
          const std::vector<Eigen::Matrix<T_parameters,
          Eigen::Dynamic, Eigen::Dynamic> >& system,
          const F_one& Pred1,
          const F_SS& PredSS,
          // T_sol sol,                // FIXME:reference
          // T_ssol ssol,              // FIXME:reference
          const Ts... pars) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using boost::math::tools::promote_args;
      using std::vector;
      using::stan::math::multiply;
      using refactor::PKRec;

      typedef typename promote_args<T_time, T_amt, T_rate, T_ii,
                                    typename promote_args<T_parameters, T_biovar, T_tlag>::type >::type scalar;
      typedef typename promote_args<T_time, T_amt, T_tlag, T_rate>::type T_tau;
      typedef typename promote_args<T_rate, T_biovar>::type T_rate2;

      // BOOK-KEEPING: UPDATE DATA SETS
      EventHistory<T_tau, T_amt, T_rate, T_ii>
        events(time, amt, rate, ii, evid, cmt, addl, ss);

      ModelParameterHistory<T_tau, T_parameters, T_biovar, T_tlag>
        parameters(time, pMatrix, biovar, tlag, system);
      RateHistory<T_tau, T_rate> rates;

      events.Sort();
      parameters.Sort();
      int nKeep = events.get_size();

      events.AddlDoseEvents();
      parameters.CompleteParameterHistory(events);

      events.AddLagTimes(parameters, nCmt);
      rates.MakeRates(events, nCmt);
      parameters.CompleteParameterHistory(events);

      PKRec<scalar> zeros = PKRec<scalar>::Zero(nCmt);
      PKRec<scalar> init = zeros;

      // COMPUTE PREDICTIONS
      Matrix<scalar, Dynamic, Dynamic>
        pred = Matrix<scalar, Dynamic, Dynamic>::Zero(nKeep, nCmt);

      T_tau dt, tprev = events.get_time(0);
      Matrix<scalar, Dynamic, 1> pred1;
      Event<T_tau, T_amt, T_rate, T_ii> event;
      ModelParameters<T_tau, T_parameters, T_biovar, T_tlag> parameter;
      int iRate = 0, ikeep = 0;

      for (int i = 0; i < events.get_size(); i++) {
        event = events.GetEvent(i);

        // Use index iRate instead of i to find rate at matching time, given there
        // is one rate per time, not per event.
        if (rates.get_time(iRate) != events.get_time(i)) iRate++;
        Rate<T_tau, T_rate2> rate2;
        rate2.copy(rates.GetRate(iRate));

        for (int j = 0; j < nCmt; j++)
          rate2.rate[j] *= parameters.GetValueBio(i, j);

        parameter = parameters.GetModelParameters(i);
        if ((event.get_evid() == 3) || (event.get_evid() == 4)) {  // reset events
          dt = 0;
          init = zeros;
        } else {
          dt = event.get_time() - tprev;
          typedef typename promote_args<T_rate, T_biovar>::type T_rate2;
          using model_type = T_model<T_tau, scalar, T_rate2, T_parameters, Ts...>;
          T_tau                     model_time = tprev;
          std::vector<T_rate2>      model_rate = rate2.get_rate();

          // std::vector<T_parameters> model_par = parameter.get_RealParameters();

          // FIX ME: we need a better way to relate model type to parameter type
          std::vector<T_parameters> model_par = model_type::get_param(parameter);
          model_type pkmodel {model_time, init, model_rate, model_par, pars...};

          pred1 = sol.solve(pkmodel, dt);
          // pred1 = Pred1(dt, parameter, init, rate2.get_rate());
          init = pred1;
        }

        if (((event.get_evid() == 1 || event.get_evid() == 4)
             && (event.get_ss() == 1 || event.get_ss() == 2)) ||
            event.get_ss() == 3) {  // steady state event
          using model_type = T_model<T_tau, scalar, T_rate2, T_parameters, Ts...>;
          T_tau model_time = event.get_time(); // FIXME: time is not t0 but for adjust within SS solver
          auto model_rate = rate2.get_rate();
          // auto model_par = parameter.get_RealParameters();
          // FIX ME: we need a better way to relate model type to parameter type
          std::vector<T_parameters> model_par = model_type::get_param(parameter);
          model_type pkmodel {model_time, init, model_rate, model_par, pars...};
          pred1 = multiply(sol.solve(pkmodel,
                                     parameters.GetValueBio(i, event.get_cmt() - 1) * event.get_amt(), //NOLINT
                                     event.get_rate(),
                                     event.get_ii(),
                                     event.get_cmt()),
                           scalar(1.0));
          // pred1 = multiply(PredSS(parameter,
          //                         parameters.GetValueBio(i, event.get_cmt() - 1)
          //                           * event.get_amt(),
          //                         event.get_rate(), event.get_ii(),
          //                         event.get_cmt()),
          //                  scalar(1.0)); 


          // the object PredSS returns doesn't always have a scalar type. For
          // instance, PredSS does not depend on tlag, but pred does. So if
          // tlag were a var, the code must promote PredSS to match the type
          // of pred1. This is done by multiplying predSS by a Scalar.

          if (event.get_ss() == 2) init += pred1;  // steady state without reset
          else
            init = pred1;  // steady state with reset (ss = 1)
        }

        if (((event.get_evid() == 1) || (event.get_evid() == 4)) &&
            (event.get_rate() == 0)) {  // bolus dose
          init(0, event.get_cmt() - 1)
            += parameters.GetValueBio(i, event.get_cmt() - 1) * event.get_amt();
        }

        if (event.get_keep()) {
          pred.row(ikeep) = init;
          ikeep++;
        }
        tprev = event.get_time();
      }

      return pred;
    }


  };

}

#endif
