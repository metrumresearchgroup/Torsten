#ifndef STAN_MATH_TORSTEN_UNIVARIATE_INTEGRAL_HPP
#define STAN_MATH_TORSTEN_UNIVARIATE_INTEGRAL_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>
#include <stan/math/rev/arr/fun/to_var.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_observer.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>
#include <stan/math/rev/mat/functor/integrate_ode_bdf.hpp>
#include <stan/math/torsten/PKModel/functors/functor.hpp>
#include <vector>
#include <algorithm>

namespace torsten {
  template<typename T0, typename Tl, typename Tr>
  inline
  std::vector<typename stan::return_type<Tl, Tr, T0>::type>
  theta_x(const std::vector<T0>& theta,
          const Tl& t0, const Tr& t1){
    std::vector<stan::math::var> res{stan::math::to_var(theta)};
    res.push_back(stan::math::to_var(t0));
    res.push_back(stan::math::to_var(t1));
    return res;
  }

  inline
  std::vector<double> theta_x(const std::vector<double>&
                              theta,
                              const double& t0,
                              const double& t1){
    std::vector<double> res{theta};
    res.push_back(t0);
    res.push_back(t1);
    return res;
  }

  // for rate parameters
  template <typename F0>
  struct normalized_integrand_adaptor {
    const F0& f0_;
    normalized_integrand_adaptor() {}
    normalized_integrand_adaptor(const F0& f0) :
      f0_(f0)
    {}

    template <typename T0, typename T1, typename T2>
    inline
    std::vector<typename boost::math::tools::promote_args<
                  T0, T1, T2>::type >
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* pstream_) const {
      const T2& t1_{theta.rbegin()[0]};
      const T2& t0_{theta.rbegin()[1]};
      const T2& jacobian {t1_ - t0_};
      using scalar = typename boost::math::tools::promote_args<
        T0, T1, T2>::type;
      std::vector<scalar> res{f0_(t1_*t + t0_*(1.0-t), y, theta, x_r, x_i, pstream_)};
      std::transform(res.begin(), res.end(), res.begin(),
                     [&jacobian](scalar x){return jacobian * x;});
      return res;
    }
  };

}

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

    template <typename F, typename Tl, typename Tr, typename T0>
    inline
    typename stan::return_type<Tl, Tr, T0>::type
    univariate_integral_bdf(const F &f0,    // integrand
                            const Tl& t0,         // integral limit
                            const Tr& t1,         // integral limit
                            const std::vector<T0>& theta,
                            const std::vector<double>& x_r,
                            const std::vector<int>& x_i,
                            std::ostream* msgs = 0,
                             double relative_tolerance = 1e-10,
                            double absolute_tolerance = 1e-10,
                            int max_num_steps = 1E8) {
      using stan::math::to_var;
      static const double t{0.0};
      static const std::vector<double> ts{1.0};
      static const std::vector<double> y0{0.0};
      const normalized_integrand_functor<F> f{f0};
      using scalar = typename stan::return_type<T0,Tl,Tr>::type;
      std::vector<scalar> par = torsten::theta_x(theta, t0, t1);

      std::vector<std::vector<scalar>> ode_res_vd =
        stan::math::integrate_ode_bdf(f, y0, t, ts, par, x_r, x_i);

      return ode_res_vd.back().back();
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

    template <typename F, typename Tl, typename Tr, typename T0>
    inline
    typename stan::return_type<Tl, Tr, T0>::type
    univariate_integral_rk45(const F &f0,    // integrand
                            const Tl& t0,         // integral limit
                            const Tr& t1,         // integral limit
                            const std::vector<T0>& theta,
                            const std::vector<double>& x_r,
                            const std::vector<int>& x_i,
                            std::ostream* msgs = 0,
                            double relative_tolerance = 1e-6,
                            double absolute_tolerance = 1e-6,
                            int max_num_steps = 1E6) {
      using stan::math::to_var;
      static const double t{0.0};
      static const std::vector<double> ts{1.0};
      static const std::vector<double> y0{0.0};
      const normalized_integrand_functor<F> f{f0};
      using scalar = typename stan::return_type<T0,Tl,Tr>::type;
      std::vector<scalar> par = torsten::theta_x(theta, t0, t1);

      std::vector<std::vector<scalar>> ode_res_vd =
        stan::math::integrate_ode_rk45(f, y0, t, ts, par, x_r, x_i);

      return ode_res_vd.back().back();
    }
}
}
#endif
