#ifndef STAN_MATH_TORSTEN_PK_NVARS_HPP
#define STAN_MATH_TORSTEN_PK_NVARS_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace torsten {
  /*
   * Calculate the number of @c var for transient dosing event
   */
  template<typename T0, typename T1, typename T2, typename T3>
  inline int pk_nvars(const T0& t0,
                      const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                      const std::vector<T2> &rate,
                      const std::vector<T3> &par) {
      using stan::is_var;
      int res = 0;
      if (is_var<T0>::value) res += 1;
      if (is_var<T1>::value) res += y0.size();
      if (is_var<T2>::value) res += rate.size();
      if (is_var<T3>::value) res += par.size();
      return res;
  }

  /*
   * Calculate the number of @c var for steady-state dosing event
   */
  template<typename T_a, typename T_r, typename T_ii, typename T_par>
  inline int pk_nvars(const T_a& a, const T_r& r, const T_ii& ii, const std::vector<T_par>& par)
  {
      using stan::is_var;
      int res = 0;
      if (is_var<T_a>::value) res += 1;
      if (is_var<T_r>::value) res += 1;
      if (is_var<T_ii>::value) res += 1;
      if (is_var<T_par>::value) res += par.size();
      return res;
  }

}
#endif
