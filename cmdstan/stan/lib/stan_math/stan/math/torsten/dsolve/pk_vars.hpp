#ifndef STAN_MATH_TORSTEN_DSOLVE_PK_VARS_HPP
#define STAN_MATH_TORSTEN_DSOLVE_PK_VARS_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace torsten {
  namespace dsolve {
    /*
     * Given multiple vectors, concat @c var vectors among them.
     */
    template<typename T1, typename T2, typename T3>
    inline std::vector<stan::math::var>
    pk_vars(const std::vector<T1>& v1, const std::vector<T2>& v2, const std::vector<T3>& v3) // NOLINT
    {
      using stan::is_var;
      std::vector<stan::math::var> res;
      if (is_var<T1>::value) res.insert(res.end(), v1.begin(), v1.end());
      if (is_var<T2>::value) res.insert(res.end(), v2.begin(), v2.end());        
      if (is_var<T3>::value) res.insert(res.end(), v3.begin(), v3.end());        
      return res;
    }

    /*
     * For transient PK solvers such as the twpcpt solver,
     * at each time step the vars are from (t1, init, rate, param).
     */
    template<typename T0, typename T1, typename T2, typename T3>
    inline std::vector<stan::math::var>
    pk_vars(const T0& x0, const Eigen::Matrix<T1, 1, -1>& m1, const std::vector<T2>& v2, const std::vector<T3>& v3) // NOLINT
    {
      using stan::is_var;
      std::vector<stan::math::var> res;
      if (is_var<T0>::value) res.push_back(x0);
      if (is_var<T1>::value) {
        std::vector<stan::math::var> v1(m1.data(), m1.data() + m1.size());
        res.insert(res.end(), v1.begin(), v1.end());
      }
      if (is_var<T2>::value) res.insert(res.end(), v2.begin(), v2.end());
      if (is_var<T3>::value) res.insert(res.end(), v3.begin(), v3.end());
      return res;
    }

    template<typename T0, typename T1, typename T2, typename T3>
    inline std::vector<stan::math::var>
    pk_vars(const T0& x0, const Eigen::Matrix<T1, 1, -1>& m1, const std::vector<T2>& v2, const Eigen::Matrix<T3,-1,-1>& m3) // NOLINT
    {
      using stan::is_var;
      std::vector<stan::math::var> res;
      if (is_var<T0>::value) res.push_back(x0);
      if (is_var<T1>::value) {
        std::vector<stan::math::var> v1(m1.data(), m1.data() + m1.size());
        res.insert(res.end(), v1.begin(), v1.end());
      }
      if (is_var<T2>::value) res.insert(res.end(), v2.begin(), v2.end());
      if (is_var<T3>::value) res.insert(res.end(), m3.size(), m3(0));
      return res;
    }

    /*
     * For stead-state PK solvers such as the twpcpt solver,
     * at each time step the vars are from (amt, ev_rate, ii, param).
     */
    template<typename T_a, typename T_r, typename T_ii, typename T3> // NOLINT
    inline std::vector<stan::math::var>
    pk_vars(const T_a& a, const T_r& r, const T_ii& ii, const std::vector<T3>& v3) // NOLINT
    {
      using stan::is_var;
      std::vector<stan::math::var> res;
      if (is_var<T_a>::value) res.push_back(a);
      if (is_var<T_r>::value) res.push_back(r);
      if (is_var<T_ii>::value) res.push_back(ii);
      if (is_var<T3>::value) res.insert(res.end(), v3.begin(), v3.end());
      return res;
    }

    template<typename T_a, typename T_r, typename T_ii, typename T3> // NOLINT
    inline std::vector<stan::math::var>
    pk_vars(const T_a& a, const T_r& r, const T_ii& ii, const Eigen::Matrix<T3,-1,-1>& m3) // NOLINT
    {
      using stan::is_var;
      std::vector<stan::math::var> res;
      if (is_var<T_a>::value) res.push_back(a);
      if (is_var<T_r>::value) res.push_back(r);
      if (is_var<T_ii>::value) res.push_back(ii);
      if (is_var<T3>::value) res.insert(res.end(), m3.size(), m3(0));
      return res;
    }
  }
}

#endif
