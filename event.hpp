#ifndef STAN_MATH_TORSTEN_EVENT_HPP
#define STAN_MATH_TORSTEN_EVENT_HPP

#include<stan/math/torsten/torsten_def.hpp>
#include<stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include<vector>

namespace torsten {
  /** 
   * An event that affects dynamical system defined by the ODEs,
   * transplated from PKPD(NONMEM) inputs
   * 
   * @param id event id, translated from EVID and mrgsolve's def
   *    0. observation
   *    1. reset
   *    2. reset + evolve
   *    3. ovsteady state
   *    4. reset + steady state
   *    5. reset a single cmt
   *    6. overwrite a single cmt
   * @param t0 previous time
   * @param t1 current time(event time/solution time)
   * @param ii steady state solution interval
   * @param jump impulsive/discontinuous forcing corresponding to bolus dosing
   * @param force continous forcing corresponding to infusion dosing
   * @param force0 continous forcing corresponding to infusion dosing
   * @param cmt equation on which forcing is imposed, corresponding to dosing cmt
   * 
   */
  template<typename time_t, typename ii_t, typename jump_t, typename force_t, typename force0_t>
  struct Event {
    /// event id:
    const int id;
    time_t t0;
    time_t t1;
    ii_t ii;
    PKRec<jump_t> jump;
    std::vector<force_t> force;
    force0_t force0;
    const int cmt;

    Event(int id_, time_t t0_, time_t t1_, ii_t ii_,
          PKRec<jump_t> jump_,
          std::vector<force_t> force_, force0_t force0_,
          int cmt_) :
      id(id_), t0(t0_), t1(t1_), ii(ii_),
      jump(jump_),
      force(force_), force0(force0_), cmt(cmt_)
    {}

    /** 
     * Solve the event using given model and numerical integrator.
     *
     * @tparam T solution type
     * @tparam model_t PMX model type
     * @tparam integrator_type integrator type
     * @param[in, out] y initial condition and output solution
     * @param[in] model PMX model with model parameters
     * @param[in] integ numerical integrator with control parameters
     */
    template<typename T, typename model_t, typename integrator_type>
    inline void operator()(PKRec<T>& y,
                           const model_t& model,
                           const integrator_type& integ) {
      using stan::math::value_of;
      const double eps = 1.0E-12;
      const jump_t jp = force0 < eps ? jump(std::abs(cmt) - 1) : 0.0;
      switch(id) {
      case 1:                   // EVID=3: reset
        y.setZero();
        break;
      case 2:                   // EVID=4, transient
        y.setZero();
        // model.solve(y, t0, t1, force, integ);
        y(cmt - 1) += jp;
        break;
      case 3:                   // EVID=1, SS=2(overlay)
        y += T(1.0) * model.solve(value_of(t1), jump(cmt - 1), force0, ii, cmt, integ);
        y(cmt - 1) += jp;
        break;
      case 4:                   // EVID=1, SS=1(overwrite) or EVID=4, SS=1 or 2
        y = T(1.0) * model.solve(value_of(t1), jump(cmt - 1), force0, ii, cmt, integ);
        y(cmt - 1) += jp;
        break;
      case 5:                   // EVID=2, turn-off/reset a cmt
        model.solve(y, t0, t1, force, integ);
        y(-cmt - 1) = 0;        // NONMEN: negative "cmt" indicates turn-off/reset
        break;
      case 6:                   // mrgsolve's EVID=9, overwrite a cmt
        model.solve(y, t0, t1, force, integ);
        y(cmt - 1) = jp;
        break;
      default:                  // EVID=1, transient dosing
        model.solve(y, t0, t1, force, integ);
        y(cmt - 1) += jp;
      }
    }

    /** 
     * Solve the event using given model and numerical integrator and
     * save the results in both @<code>var</code> format and data-only
     * format for MPI communication.
     *
     * @tparam T solution type
     * @tparam model_t PMX model type
     * @tparam integrator_type integrator type
     * @tparam Ts scalar model parameters for <code>nCmt</code> and <code>f</code>.
     * @param[in, out] yd data-only solution
     * @param[in, out] y initial condition and output solution
     * @param[in] model PMX model
     * @param[in] integ numerical integrator with control parameters
     * @param[in] model_pars int & <code>F</code> functor type params for model construction
     */
    template<typename T, typename model_t, typename integrator_type, typename... Ts>
    inline void operator()(Eigen::VectorXd& yd,
                           PKRec<T>& y,
                           const model_t& model,
                           const integrator_type& integ,
                           const Ts... model_pars) {
      using stan::math::value_of;
      const double eps = 1.0E-12;
      const jump_t jp = force0 < eps ? jump(std::abs(cmt) - 1) : 0.0;
      std::vector<stan::math::var> vt;
      switch(id) {
      case 1:
        y.setZero();
        break;
      case 2:
        y.setZero();
        vt = pmx_model_vars<model_t>::vars(t1, y, force, model.par());
        if (t1 > t0) {
          yd = model_solve_d(model, y, t0, t1, force, integ, model_pars...);
          y = torsten::mpi::precomputed_gradients(yd, vt);
        }
        y(cmt - 1) += jp;
        break;
      case 3:
        yd = model_solve_d(model, value_of(t1), jump(cmt - 1), force0, ii, cmt, integ, model_pars...);
        vt = dsolve::pk_vars(jump(cmt - 1), force0, ii, model.par());
        y += torsten::mpi::precomputed_gradients(yd, vt);
        y(cmt - 1) += jp;
        break;
      case 4:
        yd = model_solve_d(model, value_of(t1), jump(cmt - 1), force0, ii, cmt, integ, model_pars...);
        vt = dsolve::pk_vars(jump(cmt - 1), force0, ii, model.par());
        y = torsten::mpi::precomputed_gradients(yd, vt);
        y(cmt - 1) += jp;
        break;
      case 5:
        vt = pmx_model_vars<model_t>::vars(t1, y, force, model.par());
        if (t1 > t0) {
          yd = model_solve_d(model, y, t0, t1, force, integ, model_pars...);
          y = torsten::mpi::precomputed_gradients(yd, vt);
        }
        y(-cmt - 1) = 0;        // NONMEM: negative "cmt" indicates turnoff/reset
        break;
      case 6:
        vt = pmx_model_vars<model_t>::vars(t1, y, force, model.par());
        if (t1 > t0) {
          yd = model_solve_d(model, y, t0, t1, force, integ, model_pars...);
          y = torsten::mpi::precomputed_gradients(yd, vt);
        }
        y(cmt - 1) = jp;
        break;
      default:
        vt = pmx_model_vars<model_t>::vars(t1, y, force, model.par());
        if (t1 > t0) {
          yd = model_solve_d(model, y, t0, t1, force, integ, model_pars...);
          y = torsten::mpi::precomputed_gradients(yd, vt);
        }
        y(cmt - 1) += jp;
      }
    }

  };
}


#endif
