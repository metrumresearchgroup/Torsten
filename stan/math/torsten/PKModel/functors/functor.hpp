#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_FUNCTOR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_FUNCTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>

namespace torsten {

/**
 * Convert a simple univariate function to the signature of
 * ODE system RHS.
 */
template <typename F0>
struct normalized_integrand_functor {
  F0 f0_;
  normalized_integrand_functor() {}
  explicit normalized_integrand_functor(const F0& f0) : f0_(f0) {}
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
    T2 t2 = theta.rbegin()[0];
    T2 t1 = theta.rbegin()[1];
    T2 jacobian {t2 - t1};
    std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
      res = f0_(t2*t + t1*(1.0-t), y, theta, x_r, x_i, pstream_);

    for(size_t i=0; i<res.size(); ++i ) res[i] *= jacobian;

    return res;
  }
};

/**
 * Functors for the general and the mix solver.
 * Returns the derivative of the ODE system describing
 * the pharmacometrics process, which is given by the
 * "natural" derivative + the rate.
 *
 * In this first structure, pass the rate as a standard
 * vector of real. In the second one, pass the rate as
 * as parameters.
 *
 * DEV NOTE: This is a non-mix approach, in the sense
 * I'm assuming all rates are parameters. An alternative
 * would be to only treat non-zero rates as parameters. While
 * this reduces the size of theta, it seems a bit risky. Until
 * I figure out if the latter method is reliable, I'll play it
 * on the safe side.
 */
template <typename F0>
struct ode_rate_dbl_functor {
  F0 f0_;

  ode_rate_dbl_functor() { }

  explicit ode_rate_dbl_functor(const F0& f0) : f0_(f0) { }

  /**
   * Assume x_r only contains rate.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
     return f0_.rate_dbl(t, y, theta, x_r, x_i, pstream_);
  }
};

template <typename F0>
struct ode_rate_var_functor {
  F0 f0_;

  explicit ode_rate_var_functor(const F0& f0) : f0_(f0) { }

  /**
   * Assume the last nOde elements of theta contain
   * the rates.
   */
   template <typename T0, typename T1, typename T2, typename T3>
   inline
   std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
   operator()(const T0& t,
              const std::vector<T1>& y,
              const std::vector<T2>& theta,
              const std::vector<T3>& x_r,
              const std::vector<int>& x_i,
              std::ostream* pstream_) const {
      return f0_.rate_var(t, y, theta, x_r, x_i, pstream_);
  }
};

}

#endif
