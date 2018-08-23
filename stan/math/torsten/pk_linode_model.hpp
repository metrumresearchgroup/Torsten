#ifndef PK_LINODE_MODEL_HPP
#define PK_LINODE_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using Eigen::Map;

  struct PKLinODE {
    /**
     * standard two compartment PK ODE RHS function
     * @tparam T0 t type
     * @tparam T1 initial condition type
     * @tparam T2 parameter type
     * @tparam T3 real data/rate type
     * @param t type
     * @param x initial condition type
     * @param parms parameters that form a n x n matrix, where
     * @c (n = x.size()). We assume the matrix is col-majored.
     * @param rate dosing rate
     * @param dummy dummy
     */
    template <typename T0, typename T1, typename T2>
    inline
    std::vector<typename stan::return_type<T0, T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<double>& rate,
               const std::vector<int>& dummy,
               std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2>::type scalar;
      
      size_t n = x.size();
      std::vector<scalar> res(n);
      Matrix<scalar, Dynamic, Dynamic> m(n, n);
      Matrix<scalar, Dynamic, 1> v(n);
      for (size_t i = 0; i < n * n; ++i) m(i) = parms[i];
      for (size_t i = 0; i < n; ++i) v(i) = x[i];
      PKLin<scalar> rv = m * v;
      for (size_t i = 0; i < n; ++i) res[i] = rv(i);

      return res;
    }
  };

  /**
   * Linear ODE model.
   *
   * @tparam T_time type of time
   * @tparam T_init type of init condition
   * @tparam T_rate type of dosing rate
   * @tparam T_par  type of parameters.
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  class PKLinODEModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const std::vector<T_par> & par_;

  public:
    static constexpr PKLinODE f_ = PKLinODE();

  /**
   * Constructor
   * FIXME need to remove parameter as this is for linode only.
   *
   * @tparam T_mp parameters class
   * @tparam Ts parameter types
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param parameter ModelParameter type
   */
    // template<template<typename...> class T_mp, typename... Ts>
    // PKLinODEModel(const T_time& t0,
    //               const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
    //               const std::vector<T_rate> &rate,
    //               const std::vector<T_par> & par,
    //               const T_mp<Ts...> &parameter) :
    //   t0_(t0),
    //   y0_(y0),
    //   rate_(rate),
    // {}

    PKLinODEModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par>& par) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      par_(par)
    {}

    const T_time              & t0()    const { return t0_; }
    const PKRec<T_init>       & y0()    const { return y0_; }
    const std::vector<T_rate> & rate()  const { return rate_; }
    const std::vector<T_par>  & par ()  const { return par_; }
    const PKLinODE            & f()     const { return f_; }

    const int ncmt () const {
      return y0_.size();
    }
    const PKLin<T_par> coef () const {
      return Eigen::Map<const PKLin<T_par> >(par_.data(), ncmt(), ncmt());
    }

  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PKLinODE PKLinODEModel<T_time, T_init, T_rate, T_par>::f_;

}

#endif
