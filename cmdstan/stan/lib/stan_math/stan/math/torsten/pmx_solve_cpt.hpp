#ifndef STAN_MATH_TORSTEN_CPT_HPP
#define STAN_MATH_TORSTEN_CPT_HPP

#include <Eigen/Dense>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/is_std_vector.hpp>
#include <stan/math/torsten/pmx_solve_ode.hpp>
#include <stan/math/prim/err/check_positive_finite.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <stan/math/torsten/ev_solver.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <string>
#include <vector>

namespace torsten {

  template<template<class> class model_t>
  struct PMXCheckCPTModelParam;

  /** 
   * check validity of one-cpt model param 2d array
   * 
   */
  template<>
  struct PMXCheckCPTModelParam<PMXOneCptModel> {
    template<typename T>
    void operator()(const std::vector<std::vector<T> >& pMatrix) {
      const char* fun = "PMXOneCptModel";
      for (size_t i = 0; i < pMatrix.size(); ++i) {
        stan::math::check_positive_finite(fun, "CL", pMatrix[i][0]);
        stan::math::check_positive_finite(fun, "V2", pMatrix[i][1]);
        stan::math::check_nonnegative(fun, "ka", pMatrix[i][2]);
        stan::math::check_finite(fun, "ka", pMatrix[i][2]);
      }
    }
  };

  /** 
   * check validity of two-cpt model param 2d array
   * 
   */
  template<>
  struct PMXCheckCPTModelParam<PMXTwoCptModel> {
    template<typename T>
    void operator()(const std::vector<std::vector<T> >& pMatrix) {
      const char* fun = "PMXTwoCptModel";
      for (size_t i = 0; i < pMatrix.size(); ++i) {
        stan::math::check_positive_finite(fun, "CL", pMatrix[i][0]);
        stan::math::check_positive_finite(fun, "Q", pMatrix[i][1]);
        stan::math::check_positive_finite(fun, "V2", pMatrix[i][2]);
        stan::math::check_positive_finite(fun, "V3", pMatrix[i][3]);
        stan::math::check_nonnegative(fun, "ka", pMatrix[i][4]);
        stan::math::check_finite(fun, "ka", pMatrix[i][4]);
      }
    }
  };

  template<template<class> class model_t>
  struct PMXSolveCPT {
    /**
     * Computes the predicted amounts in each compartment at each event
     * for a compartment model with analytical solution, currently
     * supporting 1 & 2 cpt models.
     *
     * @tparam T0 type of scalar for time of events.
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam T5 type of scalar for bio-variability F.
     * @tparam T6 type of scalar for lag times.
     * @param[in] pMatrix parameters at each event
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
     * @return a matrix with predicted amount in each compartment
     *         at each event.
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename... Ts>
    Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, Ts...>,
                   Eigen::Dynamic, Eigen::Dynamic>
    static solve(TORSTEN_PMX_FUNC_EVENTS_ARGS,
                 const std::vector<std::vector<T4> >& pMatrix,
                 const std::vector<std::vector<Ts> >&... event_array_2d_params) {
      using std::vector;
      using Eigen::Dynamic;
      using Eigen::Matrix;
      using boost::math::tools::promote_args;
      using stan::math::check_positive_finite;

      static const char* function("pmx_solve_cpt");
      torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, function);
      torsten::pmx_check(time.size(), model_t<double>::Ncmt, event_array_2d_params...);
      PMXCheckCPTModelParam<model_t> param_check;
      param_check(pMatrix);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<Ts...> >>;
      const ER events_rec(model_t<T4>::Ncmt, time, amt, rate, ii, evid, cmt, addl, ss);

      Matrix<typename EM::T_scalar, -1, -1> pred(EM::nCmt(events_rec), events_rec.num_event_times());

      using model_type = model_t<typename EM::T_par>;
      EventSolver<model_type, EM> pr;
      pr.pred(0, events_rec, pred, PMXOdeIntegrator<Analytical>(), pMatrix, event_array_2d_params...);
      return pred;
    }

    /**
     * Overload function to allow user to pass an std::vector for 
     * <code>pMatrix/bioavailability/tlag</code>. This indicates the
     * entire population share these parameters.
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par>::value && torsten::is_std_vector<T_biovar>::value && torsten::is_std_vector<T_tlag>::value)>* = nullptr> //NOLINT
    Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                typename torsten::value_type<T_par>::type,
                                                typename torsten::value_type<T_biovar>::type,
                                                typename torsten::value_type<T_tlag>::type>,
                   Eigen::Dynamic, Eigen::Dynamic>
    static solve(TORSTEN_PMX_FUNC_EVENTS_ARGS,
                 const std::vector<T_par>& pMatrix,
                 const std::vector<T_biovar>& biovar,
                 const std::vector<T_tlag>& tlag) {
      return solve(time, amt, rate, ii, evid, cmt, addl, ss,
                   torsten::to_array_2d(pMatrix),
                   torsten::to_array_2d(biovar),
                   torsten::to_array_2d(tlag));
    }

    /**
     * Overload function to allow user to pass an std::vector for 
     * <code>pMatrix/bioavailability</code>. This indicates the
     * entire population share these parameters.
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par>::value && torsten::is_std_vector<T_biovar>::value)>* = nullptr> //NOLINT
    Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                typename torsten::value_type<T_par>::type,
                                                typename torsten::value_type<T_biovar>::type>,
                   Eigen::Dynamic, Eigen::Dynamic>
    static solve(TORSTEN_PMX_FUNC_EVENTS_ARGS,
                 const std::vector<T_par>& pMatrix,
                 const std::vector<T_biovar>& biovar) {
      return solve(time, amt, rate, ii, evid, cmt, addl, ss, torsten::to_array_2d(pMatrix), torsten::to_array_2d(biovar));
    }

    /**
     * Overload function to allow user to pass an std::vector for 
     * <code>pMatrix</code>. This indicates the
     * entire population share these parameters.
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par>::value)>* = nullptr> //NOLINT
    Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                typename torsten::value_type<T_par>::type>,
                   Eigen::Dynamic, Eigen::Dynamic>
    static solve(TORSTEN_PMX_FUNC_EVENTS_ARGS,
                 const std::vector<T_par>& pMatrix) {
      return solve(time, amt, rate, ii, evid, cmt, addl, ss, torsten::to_array_2d(pMatrix));
    }

  };
}

#endif
