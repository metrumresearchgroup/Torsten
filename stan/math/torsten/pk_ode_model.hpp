#ifndef STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  /**
   * ODE-based PKPD models.
   *
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   * @tparam F ODE functor
   * @tparam Ti ODE additional parameter type, usually the ODE size
   */
  template<typename T_time,
           typename T_init,
           typename T_rate,
           typename T_par,
           typename F,
           typename Ti>
  class PKODEModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const std::vector<T_par> &par_;
    const F &f_;
    const int ncmt_;
  public:
    using scalar_type = typename promote_args<T_time,
                                              T_rate, T_par, T_init>::type;
    using aug_par_type = typename promote_args<T_rate, T_par, T_init>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;
    using f_type      = F;

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
   * @param f ODE functor
   * @param ncmt the ODE size.
   */
    template<template<typename...> class T_mp, typename... Ts>
    PKODEModel(const T_time& t0,
               const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
               const std::vector<T_rate> &rate,
               const std::vector<T_par> &par,
               const T_mp<Ts...> &parameter,
               const F& f,
               const Ti &ncmt) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      par_(par),
      f_(f),
      ncmt_(ncmt)
    {}

  /**
   * Constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param f ODE functor
   * @param ncmt the ODE size.
   */
    PKODEModel(const T_time& t0,
               const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
               const std::vector<T_rate> &rate,
               const std::vector<T_par> &par,
               const F& f,
               const Ti &ncmt) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      par_(par),
      f_(f),
      ncmt_(ncmt)
    {}

  /**
   * Constructor from any other model type as long as it
   * provides enough information to build an ODE model.
   *
   * @tparam T_model model type that can be used as ODE model.
   * @tparam Ts type parameters for @c T_model
   * @param m model that provides ODE information.
   */
    template<template<typename...> class T_model, typename... Ts>    
    PKODEModel(const T_model<Ts...>& m) :
      t0_(m.t0()),
      y0_(m.y0()),
      rate_(m.rate()),
      par_(m.par()),
      f_(m.f()),
      ncmt_(m.ncmt())
    {}
    
    // // copy constructor
    // template<typename F1>
    // PKODEModel<T_time, T_init, T_rate, T_par, F1, Ti>    
    // with_f(const F1& f_new) const {
    //   return PKODEModel<T_time, T_init, T_rate, T_par, F1, Ti>
    //     (this -> t0_, this -> y0_, this -> rate_,
    //      this -> par_, f_new, this -> ncmt_);
    // }

    /**
     * @return initial time
     */
    const T_time              & t0()       const { return t0_; }
    /**
     * @return initial condition
     */
    const PKRec<T_init>    & y0()       const { return y0_; }
    /**
     * @return dosing rate
     */
    const std::vector<T_rate> & rate()     const { return rate_; }
    /**
     * @return model parameters
     */
    const std::vector<T_par>  & par()      const { return par_; }
    /**
     * @return RHS functor
     */
    const F                   & f ()      const { return f_; }
    /**
     * @return ODE size
     */
    const int                 & ncmt ()    const { return ncmt_; }
  };

}




#endif
