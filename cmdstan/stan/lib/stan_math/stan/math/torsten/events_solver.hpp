#ifndef STAN_MATH_TORSTEN_EVENTS_SOLVER_HPP
#define STAN_MATH_TORSTEN_EVENTS_SOLVER_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/events_manager.hpp>
#include <vector>

namespace torsten {

  /*
   * Expose PMX model information to events manager
   */
  template<typename T_events_record, typename T_model>
  struct EventsSolver;

  template<template<typename...> class NONMENEventsRecord, typename T_model, typename... Ts>
  struct EventsSolver<NONMENEventsRecord<Ts...>, T_model>
  {
    using ER = NONMENEventsRecord<Ts...>;
    using EM = EventsManager<ER>;

    static int system_size(int id, const ER& rec) {
      int ncmt = EM::nCmt(rec);
      int nvar = T_model::nvars(ncmt, EM::parameter_size(rec));
      int nvar_ss = T_model::template nvars<typename EM::T_amt, typename EM::T_par_rate, typename EM::T_par_ii>(EM::parameter_size(rec)); //NOLINT
      return EM::has_ss_dosing(id, rec) ? torsten::pk_nsys(ncmt, nvar, nvar_ss) : torsten::pk_nsys(ncmt, nvar);
    }
  };  
}


#endif
