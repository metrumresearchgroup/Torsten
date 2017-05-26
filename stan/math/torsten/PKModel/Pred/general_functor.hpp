#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_GENERAL_FUNCTOR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_GENERAL_FUNCTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>

/**
 * Functor for the general solver.
 * Returns the derivative of the ODE system describing
 * the pharmacometrics process, which is given by the
 * "natural" derivative + the rate.
 *
 * In this first structure, pass the rate as a standard
 * vector of real. In the second one, pass the rate as
 * as parameters.
 */
template <typename F0>
struct general_rate_dbl_functor {
  F0 f0_;

  explicit general_rate_dbl_functor(const F0& f0) : f0_(f0) { }

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
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;
    assert(y.size() == x_r.size());

    std::vector<scalar> dydt = f0_(t, y, theta, x_r, x_i, pstream_);
    for (size_t i = 0; i < dydt.size(); i++)
      dydt[i] += x_r[i];

    return dydt;
  }
};

template <typename F0>
struct general_rate_var_functor {
  F0 f0_;

  explicit general_rate_var_functor(const F0& f0) : f0_(f0) { }

  /**
   * Assume the last nOde elements of theta contain
   * the rates in each compartment.
   *
   * DEV NOTE: This is a non-mix approach, in the sense 
   * I'm assuming all rates are parameters. An alternative
   * would be to only treat non-zero rates as parameters. While
   * this reduces the size of theta, it seems a bit risky. Until
   * I figure out if the latter method is reliable, I'll play it 
   * on the safe side.
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
      typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
        scalar;
      
      std::vector<scalar> dydt = f0_(t, y, theta, x_r, x_i, pstream_);

      size_t nTheta = theta.size();
      size_t nOde = dydt.size();
      int indexRate = nTheta - nOde;
      for (size_t i = 0; i < dydt.size(); i++)
        dydt[i] += theta[indexRate + i];

      return dydt;
  }
};

#endif
