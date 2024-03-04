#ifndef STAN_MATH_TORSTEN_ODE_RHS_OSTREAM_ADATPOR_HPP
#define STAN_MATH_TORSTEN_ODE_RHS_OSTREAM_ADATPOR_HPP

namespace torsten {
  template <typename F>
  struct pmx_ode_ostream_adapted {
    using yes = double;
    using no = bool;

    template <typename Functor>
    static double test(decltype(std::declval<Functor&>()(0.0, Eigen::VectorXd(), nullptr,
                                                         std::vector<double>(),
                                                         std::vector<double>(),
                                                         std::vector<int>()))*);

    template <typename Functor>
    static no test(...);

    static constexpr bool value = sizeof(test<F>(nullptr)) == sizeof(yes);
  };
  
  template <typename F>
  struct pmx_ode_ostream_not_adapted {
    using yes = double;
    using no = bool;

    template <typename Functor>
    static double test(decltype(std::declval<Functor&>()(0.0, Eigen::VectorXd(),
                                                         std::vector<double>(),
                                                         std::vector<double>(),
                                                         std::vector<int>(), nullptr))*);

    template <typename Functor>
    static no test(...);

    static constexpr bool value = sizeof(test<F>(nullptr)) == sizeof(yes);
  };

  /** 
   * This adaptor changes the location of ostream arg in functor call so
   * we don't need to make change in stanc compiler for variadic solver
   * backend.
   * 
   */
  template<typename F>
  struct ode_rhs_ostream_adaptor {
    F const& f;

    ode_rhs_ostream_adaptor(F const& f_) : f(f_) {}

    template <typename F0, typename T0, typename T1, typename T2,
              std::enable_if_t<pmx_ode_ostream_not_adapted<F0>::value>* = nullptr>
    Eigen::Matrix<stan::return_type_t<T0, T1, T2>, -1, 1>
    impl_(F0 const& f0,
          const T0& t, const Eigen::Matrix<T1, -1, 1> &x,
          std::ostream* pstream,
          const std::vector<T2>& parms,
          const std::vector<double>& x_r, const std::vector<int>& x_i) const {
      return f0(t, x, parms, x_r, x_i, pstream);
    }

    template <typename F0, typename T0, typename T1, typename T2,
              std::enable_if_t<pmx_ode_ostream_adapted<F0>::value>* = nullptr>
    Eigen::Matrix<stan::return_type_t<T0, T1, T2>, -1, 1>
    impl_(F0 const& f0,
          const T0& t, const Eigen::Matrix<T1, -1, 1> &x,
          std::ostream* pstream,
          const std::vector<T2>& parms,
          const std::vector<double>& x_r, const std::vector<int>& x_i) const {
      return f0(t, x, pstream, parms, x_r, x_i);
    }

    template <typename T0, typename T1, typename T2>
    Eigen::Matrix<stan::return_type_t<T0, T1, T2>, -1, 1>
    operator()(const T0& t, const Eigen::Matrix<T1, -1, 1> &x,
               std::ostream* pstream,
               const std::vector<T2>& parms,
               const std::vector<double>& x_r, const std::vector<int>& x_i) const {
      return impl_(f, t, x, pstream, parms, x_r, x_i);
    }
  };

}
#endif
