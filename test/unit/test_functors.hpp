#ifndef TORSTEN_PMX_TEST_FUNCTORS_HPP
#define TORSTEN_PMX_TEST_FUNCTORS_HPP

#include <stan/math/torsten/torsten.hpp>

#define PMX_ADD_FUNCTOR(name, func)             \
  struct name##_functor {                       \
    template <typename... Ts>                   \
    auto operator()(Ts... args) const {         \
      return func(args...);                     \
    }                                           \
  };

PMX_ADD_FUNCTOR(pmx_solve_onecpt, torsten::pmx_solve_onecpt);
PMX_ADD_FUNCTOR(pmx_solve_twocpt, torsten::pmx_solve_twocpt);
PMX_ADD_FUNCTOR(pmx_solve_linode, torsten::pmx_solve_linode);
PMX_ADD_FUNCTOR(pmx_solve_adams, torsten::pmx_solve_adams);
PMX_ADD_FUNCTOR(pmx_solve_bdf, torsten::pmx_solve_bdf);
PMX_ADD_FUNCTOR(pmx_solve_rk45, torsten::pmx_solve_rk45);
PMX_ADD_FUNCTOR(pmx_solve_onecpt_effcpt, torsten::pmx_solve_onecpt_effcpt);
PMX_ADD_FUNCTOR(pmx_solve_twocpt_effcpt, torsten::pmx_solve_twocpt_effcpt);
PMX_ADD_FUNCTOR(pmx_solve_onecpt_rk45, torsten::pmx_solve_onecpt_rk45);
PMX_ADD_FUNCTOR(pmx_solve_onecpt_bdf, torsten::pmx_solve_onecpt_bdf);
PMX_ADD_FUNCTOR(pmx_solve_twocpt_rk45, torsten::pmx_solve_twocpt_rk45);
PMX_ADD_FUNCTOR(pmx_solve_twocpt_bdf, torsten::pmx_solve_twocpt_bdf);

#endif
