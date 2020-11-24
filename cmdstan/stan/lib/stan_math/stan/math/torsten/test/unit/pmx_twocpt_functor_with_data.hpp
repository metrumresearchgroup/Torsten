#ifndef STAN_MATH_TORSTEN_TEST_TWOCPT_FUNCTOR_WITH_DATA_HPP
#define STAN_MATH_TORSTEN_TEST_TWOCPT_FUNCTOR_WITH_DATA_HPP

struct twocpt_ode_with_data {
  /**
   * standard two compartment PK ODE RHS function, with <code>CL</code>
   * added with real & integer data.
   * @tparam T0 t type
   * @tparam T1 initial condition type
   * @tparam T2 parameter type
   * @tparam T3 real data/rate type
   * @param t type
   * @param x initial condition type
   * @param parms parameters
   * @param rate dosing rate
   * @param dummy dummy
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename stan::return_type<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i, std::ostream* pstream__) const {
    typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

    scalar
      CL = parms.at(0) + x_r[0] + double(x_i[0]),
      Q = parms.at(1),
      V1 = parms.at(2),
      V2 = parms.at(3),
      ka = parms.at(4),
      k10 = CL / V1,
      k12 = Q / V1,
      k21 = Q / V2;

    std::vector<scalar> y(3, 0);
    y.at(0) = -ka * x.at(0);
    y.at(1) = ka * x.at(0) - (k10 + k12) * x.at(1) + k21 * x.at(2);
    y.at(2) = k12 * x.at(1) - k21 * x.at(2);

    return y;
  }
};

struct twocpt_ode_with_real_data {
  /**
   * standard two compartment PK ODE RHS function, with <code>CL</code>
   * added with real data.
   * @tparam T0 t type
   * @tparam T1 initial condition type
   * @tparam T2 parameter type
   * @tparam T3 real data/rate type
   * @param t type
   * @param x initial condition type
   * @param parms parameters
   * @param rate dosing rate
   * @param dummy dummy
   */
  template <typename T0, typename T1, typename T2, typename T3>
  inline
  std::vector<typename stan::return_type<T0, T1, T2, T3>::type>
  operator()(const T0& t,
             const std::vector<T1>& x,
             const std::vector<T2>& parms,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i, std::ostream* pstream__) const {
    typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

    scalar
      CL = parms.at(0) + x_r[0],
      Q = parms.at(1),
      V1 = parms.at(2),
      V2 = parms.at(3),
      ka = parms.at(4),
      k10 = CL / V1,
      k12 = Q / V1,
      k21 = Q / V2;

    std::vector<scalar> y(3, 0);
    y.at(0) = -ka * x.at(0);
    y.at(1) = ka * x.at(0) - (k10 + k12) * x.at(1) + k21 * x.at(2);
    y.at(2) = k12 * x.at(1) - k21 * x.at(2);

    return y;
  }
};

#endif
