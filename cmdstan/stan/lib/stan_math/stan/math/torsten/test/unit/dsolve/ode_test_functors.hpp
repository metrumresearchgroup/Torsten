#ifndef STAN_MATH_TEST_ODE_TORSTEN_TEST_FUNCTORS_HPP
#define STAN_MATH_TEST_ODE_TORSTEN_TEST_FUNCTORS_HPP

#include <test/unit/math/rev/functor/ode_test_functors.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_ckrk.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_rk45.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_bdf.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_rk45.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_bdf.hpp>

#define TORSTEN_DEF_ODE_SOLVER_FUNCTOR(solver_name, solver_func)                  \
  struct solver_name##_functor {                                               \
    const std::string functor_name = #solver_name;                             \
                                                                               \
    template <typename F, typename T_y0, typename T_t0, typename T_ts,         \
              typename... Args, stan::require_eigen_vector_t<T_y0>* = nullptr> \
    std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,  \
                              Eigen::Dynamic, 1>>                              \
    operator()(const F& f, const T_y0& y0, const T_t0& t0,                     \
               const std::vector<T_ts>& ts, std::ostream* msgs,                \
               const Args&... args) {                                          \
      return solver_func(f, y0, t0, ts, msgs, args...);                        \
    }                                                                          \
                                                                               \
    template <typename F, typename T_y0, typename T_t0, typename T_ts,         \
              typename... Args, stan::require_eigen_vector_t<T_y0>* = nullptr> \
    std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,  \
                              Eigen::Dynamic, 1>>                              \
    operator()(const F& f, const T_y0& y0_arg, const T_t0& t0,                 \
               const std::vector<T_ts>& ts, double rtol, double atol,          \
               size_t max_num_steps, std::ostream* msgs,                       \
               const Args&... args) {                                          \
      return solver_func##_ctrl(f, y0_arg, t0, ts, msgs, rtol, atol, max_num_steps, \
                               args...);                                 \
    }                                                                          \
  };

TORSTEN_DEF_ODE_SOLVER_FUNCTOR(pmx_ode_adams, torsten::pmx_ode_adams);
TORSTEN_DEF_ODE_SOLVER_FUNCTOR(pmx_ode_ckrk, torsten::pmx_ode_ckrk);
TORSTEN_DEF_ODE_SOLVER_FUNCTOR(pmx_ode_bdf, torsten::pmx_ode_bdf);
TORSTEN_DEF_ODE_SOLVER_FUNCTOR(pmx_ode_rk45, torsten::pmx_ode_rk45);

STAN_DEF_STD_ODE_SOLVER_FUNCTOR(pmx_integrate_ode_adams,
                                torsten::pmx_integrate_ode_adams);
STAN_DEF_STD_ODE_SOLVER_FUNCTOR(pmx_integrate_ode_bdf,
                                torsten::pmx_integrate_ode_bdf);
STAN_DEF_STD_ODE_SOLVER_FUNCTOR(pmx_integrate_ode_rk45,
                                torsten::pmx_integrate_ode_rk45);

#endif
