#ifndef STAN_MATH_PRIM_ARR_FUN_UNIVARIATE_INTEGRAL_HPP
#define STAN_MATH_PRIM_ARR_FUN_UNIVARIATE_INTEGRAL_HPP

#include <stan/math/prim/arr/functor/coupled_ode_observer.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>
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
     * y(t1)
     *
     * @tparam F Type of functor that is to be integrated.
     * @param f function that is to be integrated
     * @param theta parameters of f
     * @param t0 integration interval left end
     * @param t1 integration interval right end
     * @return Scalar of integral value.
     */

    template <typename F, typename T0, typename T1>
    inline
    typename stan::return_type<T0, T1>::type
    univariate_integral(const F &f0,  // function to be integrated
			       const T0 &t0,    // interval end: left
			       const T1 &t1) {  // interval end: right
      using std::vector;

      std::vector<double> y0{0.0};
      std::vector<double> ts{t1};
      std::vector<double> theta;
      std::vector<double> x;
      std::vector<int> x_int;

      ode_simple_functor<F> f{f0};

      using scalar = typename stan::return_type<T0, T1>::type;
      std::vector<std::vector<scalar>> ode_res_vd =
	stan::math::integrate_ode_rk45(f, y0, t0, ts, theta, x, x_int);

      double res{ode_res_vd.back().back()};

      return res;
    }
  }
}
#endif
