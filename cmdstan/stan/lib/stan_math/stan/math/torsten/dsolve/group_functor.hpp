#ifndef STAN_MATH_TORSTEN_DSOLVE_GROUP_FUNCTOR
#define STAN_MATH_TORSTEN_DSOLVE_GROUP_FUNCTOR

#include <stan/math/torsten/dsolve/pmx_ode_system.hpp>
#include <cvodes/cvodes.h>

namespace torsten {
  namespace dsolve {
    /*
     * ID for integrators, used for master passing integrator
     * information to slaves. Defaults(rk45) is of value 3.
     */
    template<typename T>
    struct integrator_id;

    template<typename scheme_t>
    struct integrator_id<dsolve::PMXOdeintIntegrator<scheme_t>> {
      static constexpr int value = 3;
    };

    /*
     * specialization for adams integrator.
     */
    template <int ism_type>
    struct integrator_id<dsolve::PMXCvodesIntegrator<CV_ADAMS, ism_type>> {
      static constexpr int value = 1;
    };

    /*
     * specialization for bdf integrator.
     */
    template <int ism_type>
    struct integrator_id<dsolve::PMXCvodesIntegrator<CV_BDF, ism_type>> {
      static constexpr int value = 2;
    };

  /*
   * Given a functor, return its id in
   * @c CALLED_FUNCTORS. We defer specilization in
   * @c stan::lang::generate_torsten_mpi.
   */
  template<typename F>
  struct pmx_ode_group_mpi_functor_id {
    static constexpr int value = -1;
  };

  /*
   * ode group functor can be used to replace @c F in class
   * template. Depends on its @c id, the functor dispatch to
   * different ODE functors. In master-slave MPI pattern,
   * master sends this integer id intead of actual functor to slaves.
   */
  struct pmx_ode_group_mpi_functor {
    int id;

    pmx_ode_group_mpi_functor(int i) : id(i) {}

    struct dummy_functor {
      template <typename T0, typename T1, typename T2>
      inline std::vector<typename stan::return_type_t<T0, T1, T2>>
      operator()(const T0& t_in, const std::vector<T1>& y_in,
                 const std::vector<T2>& theta, const std::vector<double>& x_r,
                 const std::vector<int>& x_i,
                 std::ostream* msgs = nullptr,
                 double rtol = 1e-6,
                 double atol = 1e-6,
                 long int max_num_step = 1e6) {
        return {};
      }
    };

    /*
     * Dispatch according to value of @c id, forward the args
     * to actual functors. Its definition is defered to
     * @c stan::lang::generate_torsten_mpi.
     */
    template<typename... Args>
    inline auto operator()(Args&&... args) const;

  };

  }
}

#endif
