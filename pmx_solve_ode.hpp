#ifndef STAN_MATH_TORSTEN_SOLVE_ODE_HPP
#define STAN_MATH_TORSTEN_SOLVE_ODE_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/pmx_population_check.hpp>
#include <stan/math/torsten/ev_solver.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <vector>

namespace torsten {

  /** 
   * helper type trait used for actual ODE solver functions to
   * detect if the last arg is <code>ostream*</code>
   * 
   * @tparam Ts all of the args for pmx solvers
   */
  template<typename... Ts>
  struct last_is_ostream_ptr;

  template<typename T>
  struct last_is_ostream_ptr<T> {
    static const bool value = std::is_same<T, std::ostream*>::value;
  };

  template<typename T, typename... Ts>
  struct last_is_ostream_ptr<T, Ts...> {
    static const bool value = last_is_ostream_ptr<Ts...>::value;
  };

/**
 * simpily pmx function declare
 * 
 */
#define TORSTEN_PMX_FUNC_EVENTS_ARGS const std::vector<T0>& time,\
    const std::vector<T1>& amt,                                  \
    const std::vector<T2>& rate,                                 \
    const std::vector<T3>& ii,                                   \
    const std::vector<int>& evid,                                \
    const std::vector<int>& cmt,                                 \
    const std::vector<int>& addl,                                \
    const std::vector<int>& ss

  template<typename integrator_type>
  struct PMXSolveODE {
    
    /// default tolerances & max steps for
    /// differential eq(DE) and algebra solver(AS)
    static constexpr double RTOL_DE = 1.e-6;
    static constexpr double ATOL_DE = 1.e-6;
    static constexpr int MAXSTEP_DE = 1e6;
    static constexpr double RTOL_AS = 1.e-6;
    static constexpr double ATOL_AS = 1.e-6;
    static constexpr int MAXSTEP_AS = 1e2;

    /**
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. 
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
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] as_rel_tol relative tolerance for the algebra solver
     * @param[in] as_abs_tol absolute tolerance for the algebra solver
     * @param[in] as_max_num_steps maximal number of steps to take within 
     *            the algebra solver
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     * FIX ME: currently have a dummy msgs argument. Makes it easier
     * to expose to stan grammar files, because I can follow more closely
     * what was done for the ODE integrator. Not ideal.
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
               const int nCmt,
               TORSTEN_PMX_FUNC_EVENTS_ARGS,
               const std::vector<std::vector<T4> >& pMatrix,
               const std::vector<std::vector<T5> >& biovar,
               const std::vector<std::vector<T6> >& tlag,
               double rel_tol,
               double abs_tol,
               long int max_num_steps,
               double as_rel_tol,
               double as_abs_tol,
               long int as_max_num_steps,
               std::ostream* msgs) {
      using std::vector;
      using Eigen::Dynamic;
      using Eigen::Matrix;

      // check arguments
      static const char* function("PMX SOLVE ODE");
      torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, function);
      torsten::pmx_check(time.size(), nCmt, biovar, tlag);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> >>;
      const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

      Matrix<typename EM::T_scalar, -1, -1> pred(EM::nCmt(events_rec), events_rec.num_event_times());

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;

      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      pr.pred(0, events_rec, pred, integrator, pMatrix, biovar, tlag, nCmt, f);
      return pred;
    }

    /*
     * Overload with default ODE & algebra solver controls 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
               const int nCmt,
               TORSTEN_PMX_FUNC_EVENTS_ARGS,
               const std::vector<std::vector<T4> >& pMatrix,
               const std::vector<std::vector<T5> >& biovar,
               const std::vector<std::vector<T6> >& tlag,
               std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * Overload with default algebra solver controls 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
               const int nCmt,
               TORSTEN_PMX_FUNC_EVENTS_ARGS,
               const std::vector<std::vector<T4> >& pMatrix,
               const std::vector<std::vector<T5> >& biovar,
               const std::vector<std::vector<T6> >& tlag,
               double rel_tol,
               double abs_tol,
               long int max_num_steps,
               std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /**
     * Overload function to allow user to pass an std::vector for 
     * pMatrix/bioavailability/tlag
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename F,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par, T_biovar, T_tlag>::value)>* = nullptr,
              typename... Ts,
              typename std::enable_if_t<torsten::none_std_vector<Ts...>::value>* = nullptr>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                       typename torsten::value_type<T_par>::type,
                                                       typename torsten::value_type<T_biovar>::type,
                                                       typename torsten::value_type<T_tlag>::type>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<T_par>& pMatrix,
          const std::vector<T_biovar>& biovar,
          const std::vector<T_tlag>& tlag,
          Ts... solver_ctrl) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   torsten::to_array_2d(pMatrix),
                   torsten::to_array_2d(biovar),
                   torsten::to_array_2d(tlag),
                   solver_ctrl...);
    }

    /** 
     * Overload: omitting lag time, with full tolerance spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
               const int nCmt,
               TORSTEN_PMX_FUNC_EVENTS_ARGS,
               const std::vector<std::vector<T4> >& pMatrix,
               const std::vector<std::vector<T5> >& biovar,
               double rel_tol,
               double abs_tol,
               long int max_num_steps,
               double as_rel_tol,
               double as_abs_tol,
               long int as_max_num_steps,
               std::ostream* msgs) {
      static const char* function("solve");
      torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, function);
      torsten::pmx_check(time.size(), nCmt, biovar);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5>>>;
      const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

      Matrix<typename EM::T_scalar, -1, -1> pred(EM::nCmt(events_rec), events_rec.num_event_times());

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;

      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      pr.pred(0, events_rec, pred, integrator, pMatrix, biovar, nCmt, f);
      return pred;
    }

    /** 
     * Overload: omitting lag time, with default ode & algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,                        
                   msgs);
    }

    /** 
     * Overload: omitting lag time, with default algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,                        
                   msgs);
    }

    /** 
     * Overload: omitting lag time, allow population-wise 1d array,
     * with full ode & algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename F,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par, T_biovar>::value)>* = nullptr,
              typename... Ts,
              typename std::enable_if_t<torsten::none_std_vector<Ts...>::value>* = nullptr>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                       typename torsten::value_type<T_par>::type,
                                                       typename torsten::value_type<T_biovar>::type>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<T_par>& pMatrix,
          const std::vector<T_biovar>& biovar,
          Ts... solver_ctrl) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   torsten::to_array_2d(pMatrix), torsten::to_array_2d(biovar),
                   solver_ctrl...);
    }

    /** 
     * Overload: omitting bioavailability & lag time,
     * with full ode & algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
               const int nCmt,
               TORSTEN_PMX_FUNC_EVENTS_ARGS,
               const std::vector<std::vector<T4> >& pMatrix,
               double rel_tol,
               double abs_tol,
               long int max_num_steps,
               double as_rel_tol,
               double as_abs_tol,
               long int as_max_num_steps,
               std::ostream* msgs) {
      // check arguments
      static const char* function("solve");
      torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, function);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<>>>;
      const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

      Matrix<typename EM::T_scalar, -1, -1> pred(EM::nCmt(events_rec), events_rec.num_event_times());

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;

      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      pr.pred(0, events_rec, pred, integrator, pMatrix, nCmt, f);
      return pred;
    }

    /** 
     * Overload: omitting bioavailability & lag time,
     * with default ode & algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,                        
                   msgs);
    }

    /** 
     * Overload: omitting bioavailability & lag time,
     * with default algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,                        
                   msgs);
    }

    /** 
     * Overload: omitting bioavailability & lag time, allow
     * population-wise 1d array,
     * with full algebra solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename F,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par>::value)>* = nullptr,
              typename... Ts,
              typename std::enable_if_t<torsten::none_std_vector<Ts...>::value>* = nullptr>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                       typename torsten::value_type<T_par>::type>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<T_par>& pMatrix,
          Ts... solver_ctrl) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   torsten::to_array_2d(pMatrix),
                   solver_ctrl...);
    }

    /** 
     * Overload: additional real data for ODE, with full ode & algebra
     * solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          double as_rel_tol,
          double as_abs_tol,
          long int as_max_num_steps,
          std::ostream* msgs) {
      // check arguments
      static const char* function("solve");
      torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, function);
      torsten::pmx_check(time.size(), nCmt, biovar, tlag);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6>, double>>;
      const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

      Matrix<typename EM::T_scalar, -1, -1> pred(EM::nCmt(events_rec), events_rec.num_event_times());

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;

      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      pr.pred(0, events_rec, pred, integrator, pMatrix, biovar, tlag, x_r, nCmt, f);
      return pred;
    }

    /** 
     * Overload: additional real data for ODE, with default ode & algebra 
     * solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                   Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /** 
     * Overload: additional real data for ODE, with default algebra 
     * solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    // Overload: with real data, PMX params can be either 1d or 2d arrays.
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename F,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par, T_biovar, T_tlag>::value)>* = nullptr,
              typename... Ts,
              typename std::enable_if_t<torsten::none_std_vector<Ts...>::value>* = nullptr>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                       typename torsten::value_type<T_par>::type,
                                                       typename torsten::value_type<T_biovar>::type,
                                                       typename torsten::value_type<T_tlag>::type>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<T_par>& pMatrix,
          const std::vector<T_biovar>& biovar,
          const std::vector<T_tlag>& tlag,
          const std::vector<std::vector<double> >& x_r,
          Ts... solver_ctrl) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   torsten::to_array_2d(pMatrix),
                   torsten::to_array_2d(biovar),
                   torsten::to_array_2d(tlag), x_r,
                   solver_ctrl...);
    }

    /** 
     * Overload: additional real & int data for ODE, with full ode & algebra
     * solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          const std::vector<std::vector<int> >& x_i,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          double as_rel_tol,
          double as_abs_tol,
          long int as_max_num_steps,
          std::ostream* msgs) {
      // check arguments
      static const char* function("solve");
      torsten::pmx_check(time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, function);
      torsten::pmx_check(time.size(), nCmt, biovar, tlag);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6>, double, int>>;
      const ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

      Matrix<typename EM::T_scalar, -1, -1> pred(EM::nCmt(events_rec), events_rec.num_event_times());

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;

      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      pr.pred(0, events_rec, pred, integrator, pMatrix, biovar, tlag, x_r, x_i, nCmt, f);
      return pred;
    }

    /** 
     * Overload: additional real & int data for ODE, with default ode & algebra
     * solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          const std::vector<std::vector<int> >& x_i,
          std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r, x_i,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /** 
     * Overload: additional real & int data for ODE, with default algebra 
     * solver spec
     * 
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3, T4, T5, T6>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
               const int nCmt,
               TORSTEN_PMX_FUNC_EVENTS_ARGS,
               const std::vector<std::vector<T4> >& pMatrix,
               const std::vector<std::vector<T5> >& biovar,
               const std::vector<std::vector<T6> >& tlag,
               const std::vector<std::vector<double> >& x_r,
               const std::vector<std::vector<int> >& x_i,
               double rel_tol,
               double abs_tol,
               long int max_num_steps,
               std::ostream* msgs) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r, x_i,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    // Overload: with real & int data, PMX params can be either 1d or 2d arrays.
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename F,
              typename std::enable_if_t<!(torsten::is_std_vector<T_par, T_biovar, T_tlag>::value)>* = nullptr,
              typename... Ts,
              typename std::enable_if_t<torsten::none_std_vector<Ts...>::value>* = nullptr>
    static Eigen::Matrix <typename stan::return_type_t<T0, T1, T2, T3,
                                                       typename torsten::value_type<T_par>::type,
                                                       typename torsten::value_type<T_biovar>::type,
                                                       typename torsten::value_type<T_tlag>::type>,
                          Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<T_par>& pMatrix,
          const std::vector<T_biovar>& biovar,
          const std::vector<T_tlag>& tlag,
          const std::vector<std::vector<double> >& x_r,
          const std::vector<std::vector<int> >& x_i,
          Ts... solver_ctrl) {
      return solve(f, nCmt,
                   time, amt, rate, ii, evid, cmt, addl, ss,
                   torsten::to_array_2d(pMatrix),
                   torsten::to_array_2d(biovar),
                   torsten::to_array_2d(tlag), x_r, x_i,
                   solver_ctrl...);
    }
  };
}  
#endif
