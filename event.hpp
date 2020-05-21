#ifndef STAN_MATH_TORSTEN_EVENT_HPP
#define STAN_MATH_TORSTEN_EVENT_HPP

#include<stan/math/torsten/torsten_def.hpp>
#include<stan/math/torsten/pmx_ode_integrator.hpp>
#include <stan/math/torsten/mpi/precomputed_gradients.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include<vector>

namespace torsten {
  template<typename time_t, typename ii_t, typename jump_t, typename force_t, typename force0_t>
  struct Event {
    /// event id:
    /// 0. observation
    /// 1. reset
    /// 2. reset + evolve
    /// 3. ovsteady state
    /// 4. reset + steady state
    /// 5. reset a single cmt
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

    template<typename T, typename model_t, PMXOdeIntegratorId It>
    inline void operator()(PKRec<T>& y,
                           const model_t& model,
                           const PMXOdeIntegrator<It>& integ) {
      using stan::math::value_of;
      const double eps = 1.0E-12;
      const jump_t jp = force0 < eps ? jump(std::abs(cmt) - 1) : 0.0;
      switch(id) {
      case 1:
        y.setZero();
        break;
      case 2:
        y.setZero();
        model.solve(y, t0, t1, force, integ);
        y(cmt - 1) += jp;
        break;
      case 3:
        y += T(1.0) * model.solve(value_of(t1), jump(cmt - 1), force0, ii, cmt, integ);
        y(cmt - 1) += jp;
        break;
      case 4:
        y = T(1.0) * model.solve(value_of(t1), jump(cmt - 1), force0, ii, cmt, integ);
        y(cmt - 1) += jp;
        break;
      case 5:
        model.solve(y, t0, t1, force, integ);
        y(-cmt - 1) = 0;        // NONMEN: negative "cmt" indicates turn-off/reset
        break;
      case 6:
        model.solve(y, t0, t1, force, integ);
        y(cmt - 1) = jp;
        break;
      default:
        model.solve(y, t0, t1, force, integ);
        y(cmt - 1) += jp;
      }
    }

    template<typename T, typename model_t, PMXOdeIntegratorId It,
             typename... Ts>
    inline void operator()(Eigen::VectorXd& yd,
                           PKRec<T>& y,
                           const model_t& model,
                           const PMXOdeIntegrator<It>& integ,
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
