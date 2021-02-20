#ifndef STAN_MATH_TORSTEN_MODEL_FACTORY_HPP
#define STAN_MATH_TORSTEN_MODEL_FACTORY_HPP

#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <vector>

namespace torsten{
  template<typename... Ts>
  using par_index_seq = std::make_integer_sequence<std::size_t, sizeof...(Ts)>;

  template<typename T_model, typename T_em, typename... scalar_pars_type>
  struct model_factory {
    template<typename std::size_t... Is>
    static T_model model(const T_em& em, int subj_id,
                         std::integer_sequence<std::size_t, Is...>,
                         const scalar_pars_type... scalar_pars) {
      return {em.theta(subj_id),
              em.template get_model_array_1d_param<Is>(subj_id)...,
              scalar_pars...};
    }
  };
}

#endif
