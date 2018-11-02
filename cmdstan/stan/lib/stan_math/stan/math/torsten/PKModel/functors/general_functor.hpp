#ifndef STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_GENERAL_FUNCTOR_HPP
#define STAN_MATH_TORSTEN_PKMODEL_FUNCTORS_GENERAL_FUNCTOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>
#include <iostream>

namespace torsten {

/**
 *  Functor for general ODE solver.
 */
template <typename F0>
struct general_functor {
  F0 f0_;

  general_functor() { }

  explicit general_functor(const F0& f0) : f0_(f0) { }

  /**
   *  Returns the derivative of the base ODE system.
   *
   *  Case 1: rate is fixed data and is passed through x_r.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  rate_dbl(const T0& t,
           const std::vector<T1>& y,
           const std::vector<T2>& theta,
           const std::vector<T3>& x_r,
           const std::vector<int>& x_i,
           std::ostream* pstream_) const {
    using stan::math::to_array_1d;
    using stan::math::to_vector;
    typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
      scalar;

    std::vector<scalar> dydt = f0_(t, y, theta, x_r, x_i, pstream_);

    for (size_t i = 0; i < dydt.size(); i++)
     dydt[i] += x_r[i];

    return dydt;
  }

  /**
   *  Case 2: rate is a parameter, stored in theta.
   *  Theta contains in this order: ODE parameters, rates, and
   *  initial base PK states.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3>::type>
  rate_var(const T0& t,
           const std::vector<T1>& y,
           const std::vector<T2>& theta,
           const std::vector<T3>& x_r,
           const std::vector<int>& x_i,
           std::ostream* pstream_) const {
      typedef typename boost::math::tools::promote_args<T0, T1, T2, T3>::type
        scalar;

      std::vector<scalar> dydt = f0_(t, y, theta, x_r, x_i, pstream_);

      const size_t n_boundary = 2;    // integral boundary as parameter
      size_t nTheta = theta.size();
      size_t nOde = dydt.size();
      int indexRate = nTheta - nOde - n_boundary;
      for (size_t i = 0; i < dydt.size(); i++)
        dydt[i] += theta[indexRate + i];

    return dydt;
  }
};

}
#endif
