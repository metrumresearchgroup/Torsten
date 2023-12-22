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
     *  
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] as_rel_tol relative tolerance for the algebra solver
     * @param[in] as_abs_tol absolute tolerance for the algebra solver
     * @param[in] as_max_num_steps maximal number of steps to take within 
     *            the algebra solver
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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

    /**
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE & algebra solver control parameter
     * take default values.
     *
     *  
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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

    /**
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. Algebra solver takes default control values.
     *
     *  
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. Overloaded @c pMatrix @c biovar @c tlag
     * arguments can be 1D arrays, indicating they are time-independent.
     * In the overloaded signature the three arguments can take
     * combination of 1D & 2D arrays.
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T_par type for model parameters (scalar or vector)
     * @tparam T_biovar type for the bio-variability parameters (scalar or vector)
     * @tparam T_tlag type for the model tlag parameters (scalar or vector)
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] solver_ctrl parameter pack for ODE and algebra solver controls,
     *            as well as I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename F,
              typename... Ts,
              typename = require_any_not_std_vector_t<T_par, T_biovar, T_tlag>,
              typename = require_all_not_std_vector_t<Ts...> >
    static stan::matrix_return_t<T0, T1, T2, T3, T_par, T_biovar, T_tlag>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam T5 type of scalars for the bio-variability parameters.
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] as_rel_tol relative tolerance for the algebra solver
     * @param[in] as_abs_tol absolute tolerance for the algebra solver
     * @param[in] as_max_num_steps maximal number of steps to take within 
     *            the algebra solver
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0, and ODE & algebra solver
     * controls also assume default values.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam T5 type of scalars for the bio-variability parameters.
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0, and algebra solver controls
     * also assume default values.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam T5 type of scalars for the bio-variability parameters.
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0. The @c pMatrix @c biovar @c tlag
     * arguments take combination of 1D & 2D arrays, with 1D array
     * indicating time-independence.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T_par type for model parameters (scalar or vector)
     * @tparam T_biovar type for the bio-variability parameters (scalar or vector)
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] solver_ctrl parameter pack for ODE & algebra solver controls,
     *            as well as I/O stream.
     *
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename F,
              typename... Ts,
              typename = require_any_not_std_vector_t<T_par, T_biovar>,
              typename = require_all_not_std_vector_t<Ts...> >
    static stan::matrix_return_t<T0, T1, T2, T3, T_par, T_biovar>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0, and @c biovar argument is
     * also omitted, with assumed value 1.0.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] as_rel_tol relative tolerance for the algebra solver
     * @param[in] as_abs_tol absolute tolerance for the algebra solver
     * @param[in] as_max_num_steps maximal number of steps to take within 
     *            the algebra solver
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0, and @c biovar argument is
     * also omitted, with assumed value 1.0, and ODE & algebra solver
     * controls assume default values.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0, and @c biovar argument is
     * also omitted, with assumed value 1.0. Algebra solver controls
     * assume default values.
     *
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T4 type of scalars for the model parameters.
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. In this signature the @c tlag argument
     * is omitted, with assumed value 0.0, and @c biovar argument is
     * also omitted, with assumed value 1.0. The @c pMatrix @c biovar @c tlag
     * arguments take combination of 1D & 2D arrays, with 1D array
     * indicating time-independence.
     *  
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T_par type for model parameters (scalar or vector)
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] solver_ctrl parameter pack for ODE & algebra solver controls,
     *            as well as I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename F,
              typename... Ts,
              typename = require_any_not_std_vector_t<T_par>,
              typename = require_all_not_std_vector_t<Ts...> >
    static stan::matrix_return_t<T0, T1, T2, T3, T_par>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional argument
     * for real data.
     *
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] as_rel_tol relative tolerance for the algebra solver
     * @param[in] as_abs_tol absolute tolerance for the algebra solver
     * @param[in] as_max_num_steps maximal number of steps to take within 
     *            the algebra solver
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional argument
     * for real data. ODE & algebra solver controls assume default values.
     *
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional argument
     * for real data. Algebra solver controls assume default values.
     *
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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

    /**
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional argument
     * for real data. Algebra solver controls assume default values.
     * The @c pMatrix @c biovar @c tlag
     * arguments take combination of 1D & 2D arrays, with 1D array
     * indicating time-independence.
     *
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T_par type for model parameters (scalar or vector)
     * @tparam T_biovar type for the bio-variability parameters (scalar or vector)
     * @tparam T_tlag type for the model tlag parameters (scalar or vector)
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] solver_ctrl parameter pack for ODE & algebra solver controls,
     *            as well as I/O stream.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename F,
              typename... Ts,
              typename = require_any_not_std_vector_t<T_par, T_biovar, T_tlag>,
              typename = require_all_not_std_vector_t<Ts...> >
    static stan::matrix_return_t<T0, T1, T2, T3, T_par, T_biovar, T_tlag>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional arguments
     * for real data & integer data.
     *
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] x_i integer data 2D array for ODE functor @c f
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] as_rel_tol relative tolerance for the algebra solver
     * @param[in] as_abs_tol absolute tolerance for the algebra solver
     * @param[in] as_max_num_steps maximal number of steps to take within 
     *            the algebra solver
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional arguments
     * for real data & integer data. ODE & algebra solver controls
     * assume default values.
     *
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] x_i integer data 2D array for ODE functor @c f
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional arguments
     * for real data & integer data. Algebra solver controls
     * assume default values.
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] x_i integer data 2D array for ODE functor @c f
     * @param[in] rel_tol relative tolerance for the Boost ode solver 
     * @param[in] abs_tol absolute tolerance for the Boost ode solver
     * @param[in] max_num_steps maximal number of steps to take within 
     *            the Boost ode solver 
     * @param[in] msg I/O stream for ODE & algebra solvers.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4,
              typename T5, typename T6, typename F>
    static stan::matrix_return_t<T0, T1, T2, T3, T4, T5, T6>
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

    /**
     * Computes the predicted amounts in each compartment at each event
     * for a general compartment model, defined by a system of ordinary
     * differential equations. ODE solver has an additional arguments
     * for real data & integer data. The @c pMatrix @c biovar @c tlag
     * arguments take combination of 1D & 2D arrays, with 1D array
     * indicating time-independence.
     *
     * @tparam T0 type of scalar for time of events. 
     * @tparam T1 type of scalar for amount at each event.
     * @tparam T2 type of scalar for rate at each event.
     * @tparam T3 type of scalar for inter-dose inteveral at each event.
     * @tparam T_par type for model parameters (scalar or vector)
     * @tparam T_biovar type for the bio-variability parameters (scalar or vector)
     * @tparam T_tlag type for the model tlag parameters (scalar or vector)
     * @tparam F type of ODE system function.
     * @param[in] f functor for base ordinary differential equation that defines 
     *            compartment model.
     * @param[in] nCmt number of compartments in model
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
     * @param[in] pMatrix 2D parameter array
     * @param[in] biovar bioavailability 2D array
     * @param[in] tlag lag time 2D array
     * @param[in] x_r real data 2D array for ODE functor @c f
     * @param[in] x_i integer data 2D array for ODE functor @c f
     * @param[in] solver_ctrl parameter pack for ODE & algebra solver controls,
     *            as well as I/O stream.
     * @return a matrix with predicted amount in each compartment 
     *         at each event. 
     *
     */
    template <typename T0, typename T1, typename T2, typename T3,
              typename T_par, typename T_biovar, typename T_tlag,
              typename F,
              typename... Ts,
              typename = require_any_not_std_vector_t<T_par, T_biovar, T_tlag>,
              typename = require_all_not_std_vector_t<Ts...> >
    static stan::matrix_return_t<T0, T1, T2, T3, T_par, T_biovar, T_tlag>
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
