#ifndef STAN_MATH_TORSTEN_DSOLVE_ODE_CHECK_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODE_CHECK_HPP

#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/system_error.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>

namespace torsten {
  namespace dsolve {

    template <typename Tt, typename T_initial, typename T_param>
    void ode_check(const std::vector<T_initial>& y0,
                   double t0,
                   const std::vector<Tt>& ts,
                   const std::vector<T_param>& theta,
                   const std::vector<double>& x_r,
                   const std::vector<int>& x_i,
                   const char* caller) {
      stan::math::check_finite(caller, "initial time", t0);
      stan::math::check_finite(caller, "times", ts);
      stan::math::check_ordered(caller, "times", ts);
      stan::math::check_nonzero_size(caller, "times", ts);
      stan::math::check_less(caller, "initial time", t0, ts.front());
      stan::math::check_finite(caller, "initial state", y0);
      stan::math::check_finite(caller, "parameters", theta);
      stan::math::check_finite(caller, "continuous data", x_r);
      stan::math::check_nonzero_size(caller, "initial state", y0);
    }

    template <typename Tt, typename T_initial, typename T_param>
    void ode_group_check(const std::vector<std::vector<T_initial> >& y0,
                         double t0,
                         const std::vector<int>& len,
                         const std::vector<Tt>& ts,
                         const std::vector<std::vector<T_param> >& theta,
                         const std::vector<std::vector<double> >& x_r,
                         const std::vector<std::vector<int> >& x_i,
                         const char* caller) {
      stan::math::check_finite(caller, "initial time", t0);
      stan::math::check_finite(caller, "times", ts);
      stan::math::check_nonzero_size(caller, "times", ts);
      stan::math::check_nonzero_size(caller, "initial state", y0);
      stan::math::check_consistent_sizes(caller, "parameters"      , theta, "subject data length", len);
      stan::math::check_consistent_sizes(caller, "initial state"   , y0,    "subject data length", len);
      stan::math::check_consistent_sizes(caller, "continuous data" , x_r,   "subject data length", len);
      stan::math::check_consistent_sizes(caller, "discrete data"   , x_i,   "subject data length", len);

      int ibegin = 0;
      std::vector<double> ts_i(len[0]);
      for (size_t i = 0; i < len.size(); ++i) {
        ts_i.resize(len[i]);
        std::transform(ts.begin() + ibegin, ts.begin() + ibegin + len[i], ts_i.begin(),
                       [](const Tt& ti) { return stan::math::value_of(ti); });
        stan::math::check_ordered(caller, "times", ts_i);
        stan::math::check_less(caller, "initial time", t0, ts_i.front());
        stan::math::check_finite(caller, "initial state", y0[i]);
        stan::math::check_finite(caller, "parameters", theta[i]);
        stan::math::check_finite(caller, "continuous data", x_r[i]);
        ibegin += len[i];
      }
    }

  }
}

#endif
