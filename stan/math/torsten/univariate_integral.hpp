#ifndef STAN_MATH_PRIM_ARR_FUN_UNIVARIATE_INTEGRAL_HPP
#define STAN_MATH_PRIM_ARR_FUN_UNIVARIATE_INTEGRAL_HPP

#include <stan/math/prim/arr/functor/coupled_ode_observer.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>

#include <vector>
#include <algorithm>

namespace stan {
  namespace math {

    /**
     * Return the integral of a univariate function (provide
     * in the form of functor) at specifed integration interval.
     * Utilizing RK45 ODE solver, the integral of f(t) in
     * (t0, t1) can be obtained by solving u(t) from
     * dy/dt = f(t) with y(t0) = 0,
     * so that the seeked integral is
     * y(t1). To bypass the limit that t0 & t1 cannot be
     * parameters we put them in theta and perform a
     * change-of-variable over the integrand.
     *
     * @tparam F Type of functor that is to be integrated.
     * @param f function that is to be integrated
     * @param y0 init condition of ODE
     * @param theta integral limits as parameters
     * @return vector of integral value(as integrand is a * vector function).
     */

    template <typename F, typename T0, typename T1>
    inline
    std::vector<typename stan::return_type<T0, T1>::type>
    univariate_integral_bdf(const F &f0,    // integrand
                            const std::vector<T0>& y0,
                            const std::vector<T1>& theta,
                            const std::vector<double>& x_r,
                            const std::vector<int>& x_i,
                            std::ostream* msgs = 0) {
      double t0{0.0};
      std::vector<double> ts{1.0};  // integral in [0, 1]

      normalized_integrand_functor<F> f{f0};  // change-of-variable

      using scalar = typename stan::return_type<T0, T1>::type;
      std::vector<std::vector<scalar>> ode_res_vd =
        stan::math::integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i);

      return ode_res_vd.back();
    }


    /**
     * Return the integral of a univariate function (provide
     * in the form of functor) at specifed integration interval.
     * Utilizing RK45 ODE solver, the integral of f(t) in
     * (t0, t1) can be obtained by solving u(t) from
     * dy/dt = f(t) with y(t0) = 0,
     * so that the seeked integral is
     * y(t1). To bypass the limit that t0 & t1 cannot be
     * parameters we put them in theta and perform a
     * change-of-variable over the integrand.
     *
     * @tparam F Type of functor that is to be integrated.
     * @param f function that is to be integrated
     * @param y0 init condition of ODE
     * @param theta integral limits as parameters
     * @return vector of integral value(as integrand is a * vector function).
     */

    template <typename F, typename T0, typename T1>
    inline
    std::vector<typename stan::return_type<T0, T1>::type>
    univariate_integral_rk45(const F &f0,    // integrand
                             const std::vector<T0>& y0,
                             const std::vector<T1>& theta,
                             const std::vector<double>& x_r,
                             const std::vector<int>& x_i,
                             std::ostream* msgs = 0) {
      double t0{0.0};
      std::vector<double> ts{1.0};  // integral in [0, 1]

      normalized_integrand_functor<F> f{f0};  // change-of-variable

      using scalar = typename stan::return_type<T0, T1>::type;
      std::vector<std::vector<scalar>> ode_res_vd =
        stan::math::integrate_ode_rk45(f, y0, t0, ts, theta, x_r, x_i);

      return ode_res_vd.back();
    }
  }
}
#endif
