#ifndef STAN_MATH_TORSTEN_SOLVE_GROUP_ODE_HPP
#define STAN_MATH_TORSTEN_SOLVE_GROUP_ODE_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/to_array_2d.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/pmx_population_check.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <stan/math/torsten/pmx_solve_ode.hpp>
#include <vector>

namespace torsten {

  template<typename integrator_type>
  struct PMXSolveGroupODE {
    static constexpr double RTOL_DE = 1.e-6;
    static constexpr double ATOL_DE = 1.e-6;
    static constexpr int MAXSTEP_DE = 1e6;
    static constexpr double RTOL_AS = 1.e-6;
    static constexpr double ATOL_AS = 1.e-6;
    static constexpr int MAXSTEP_AS = 1e2;
    /*
     * For population models, more often we use ragged arrays
     * to describe the entire population, so in addition we need the arrays of
     * the length of each individual's data. The size of that
     * vector is the size of the population.
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
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
      static const char* caller("PMX SOLVE GROUP ODE");
      torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, caller);
      torsten::pmx_population_check(len, time, biovar, tlag, caller);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> >>;
      ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;
      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(nCmt, events_rec.total_num_event_times);

      pr.pred(events_rec, pred, integrator, pMatrix, biovar, tlag, nCmt, f);

      return pred;
    }

    /*
     * overload with default ode & algebra solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload with default algebra solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload for omitting <code>tlag</code>, with all the solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
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
      static const char* caller("PMX SOLVE GROUP ODE");
      torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, caller);
      torsten::pmx_population_check(len, time, biovar, caller);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5> >>;
      ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;
      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(nCmt, events_rec.total_num_event_times);

      pr.pred(events_rec, pred, integrator, pMatrix, biovar, nCmt, f);

      return pred;
    }

    /*
     * overload for omitting <code>tlag</code>, with default solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload for omitting <code>tlag</code>, with default algebra
     * solver controls and user-defined ode solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload for omitting <code>biuovar/tlag</code>, with all the solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          double as_rel_tol,
          double as_abs_tol,
          long int as_max_num_steps,
          std::ostream* msgs) {
      static const char* caller("PMX SOLVE GROUP ODE");
      torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, caller);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<> >>;
      ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;
      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(nCmt, events_rec.total_num_event_times);

      pr.pred(events_rec, pred, integrator, pMatrix, nCmt, f);

      return pred;
    }

    /*
     * overload for omitting <code>biovar/tlag</code>, with default solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload for omitting <code>tlag</code>, with default algebra
     * solver controls and user-defined ode solver controls
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * For population models, more often we use ragged arrays
     * to describe the entire population, so in addition we need the arrays of
     * the length of each individual's data. The size of that
     * vector is the size of the population. 
     * Overload with support of real data.
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
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
      static const char* caller("PMX SOLVE GROUP ODE");
      torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, caller);
      torsten::pmx_population_check(len, time, biovar, tlag, caller);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6>, double>>;
      ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;
      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(nCmt, events_rec.total_num_event_times);

      pr.pred(events_rec, pred, integrator, pMatrix, biovar, tlag, x_r, nCmt, f);

      return pred;
    }

    /*
     * overload with default ode & algebra solver controls, with support of real data
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload with default algebra solver controls, with support of real data
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          double rel_tol,
          double abs_tol,
          long int max_num_steps,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * For population models, more often we use ragged arrays
     * to describe the entire population, so in addition we need the arrays of
     * the length of each individual's data. The size of that
     * vector is the size of the population. 
     * Overload with support of real & integer data.
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
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
      static const char* caller("PMX SOLVE GROUP ODE");
      torsten::pmx_population_check(len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, caller);
      torsten::pmx_population_check(len, time, biovar, tlag, caller);

      using ER = NONMENEventsRecord<T0, T1, T2, T3>;
      using EM = EventsManager<ER, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6>, double, int>>;
      ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

      using model_type = torsten::PKODEModel<typename EM::T_par, F>;
      integrator_type integrator(rel_tol, abs_tol, max_num_steps, as_rel_tol, as_abs_tol, as_max_num_steps, msgs);
      EventSolver<model_type, EM> pr;

      Eigen::Matrix<typename EM::T_scalar, -1, -1> pred(nCmt, events_rec.total_num_event_times);

      pr.pred(events_rec, pred, integrator, pMatrix, biovar, tlag, x_r, x_i, nCmt, f);

      return pred;
    }

    /*
     * overload with default ode & algebra solver controls, with support of real & integer data
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
          TORSTEN_PMX_FUNC_EVENTS_ARGS,
          const std::vector<std::vector<T4> >& pMatrix,
          const std::vector<std::vector<T5> >& biovar,
          const std::vector<std::vector<T6> >& tlag,
          const std::vector<std::vector<double> >& x_r,
          const std::vector<std::vector<int> >& x_i,
          std::ostream* msgs) {
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r, x_i,
                   RTOL_DE, ATOL_DE, MAXSTEP_DE,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }

    /*
     * overload with default algebra solver controls, with support of real & integer data
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static Eigen::Matrix<typename  EventsManager<NONMENEventsRecord<T0, T1, T2, T3>, NonEventParameters<T0, T4, std::vector, std::tuple<T5, T6> > >::T_scalar, // NOLINT
                         Eigen::Dynamic, Eigen::Dynamic>
    solve(const F& f,
          const int nCmt,
          const std::vector<int>& len,
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
      return solve(f, nCmt, len, time, amt, rate,
                   ii, evid, cmt, addl, ss,
                   pMatrix, biovar, tlag, x_r, x_i,
                   rel_tol, abs_tol, max_num_steps,
                   RTOL_AS, ATOL_AS, MAXSTEP_AS,
                   msgs);
    }
  };

}
#endif
