#ifndef STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/torsten/pk_ode_integrator.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::PkOdeIntegrator;
  using torsten::PkOdeIntegratorId;

  template<typename F, typename T_rate, typename T_par>
  struct PKOdeFunctorRateAdaptor;

  /*
   * Adaptor for ODE functor when rate is data. In this
   * case rate should be passed in by @c x_r.
   */
  template<typename F, typename T_par>
  struct PKOdeFunctorRateAdaptor<F, double, T_par> {
    const F& f;
    explicit PKOdeFunctorRateAdaptor(const F& f0) : f(f0) {}

    template <typename T0, typename T1, typename T2>
    inline std::vector<typename stan::return_type<T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      std::vector<typename stan::return_type<T1, T2>::type> res;
      res = f(t, y, theta, x_r, x_i, msgs);
      for (size_t i = 0; i < y.size(); i++) res.at(i) += x_r.at(i);
      return res;
    }
  };

  /*
   * Adaptor for ODE functor when rate is @c var. In this
   * case rate should be passed in by @c theta. When
   * constructed, this type stores an index @c index_rate of @c rate
   * within @c theta. Note that we only allow @c T_rate to
   * be @c var when @c T_par is @c var, so in this case @c
   * theta is always a @c var vector.
   */
  template<typename F>
  struct PKOdeFunctorRateAdaptor<F, stan::math::var, stan::math::var> {
    const F& f;
    const int index_rate;
    PKOdeFunctorRateAdaptor(const F& f0, int i) : f(f0), index_rate(i) {}

    template <typename T0, typename T1, typename T2>
    inline std::vector<typename stan::return_type<T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      std::vector<typename stan::return_type<T1, T2>::type> res;
      res = f(t, y, theta, x_r, x_i, msgs);
      for (size_t i = 0; i < y.size(); i++) res.at(i) += theta.at(i + index_rate);
      return res;
    }
  };


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

  private:
    /*
     * We overload @c integrate so that we can pass @c rate
     * with different types, due to limit of c++ of partial
     * spec of member functions.
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    integrate(const std::vector<double> &rate,
             const T_time& dt,
             const double rtol,
             const double atol,
             const int max_num_steps,
             std::ostream* msgs) {
      std::vector<T_time> ts{t0_ + dt};
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        PKOdeFunctorRateAdaptor<F, double, T_par> f(f_);
        std::vector<int> x_i;
        std::vector<std::vector<scalar_type> > res_v =
          PkOdeIntegrator<It>()(f, y, t0_, ts, par_,
                                rate, x_i, msgs, rtol, atol, max_num_steps);
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * We overload @c integrate so that we can pass @c rate
     * with different types, due to limit of c++ of partial
     * spec of member functions.
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    integrate(const std::vector<stan::math::var> &rate,
             const T_time& dt,
             const double rtol,
             const double atol,
             const int max_num_steps,
             std::ostream* msgs) {
      using stan::math::var;

      std::vector<T_time> ts{t0_ + dt};
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        PKOdeFunctorRateAdaptor<F, var, T_par> f(f_, par_.size());
        std::vector<double> x_r;
        std::vector<int> x_i;
        std::vector<stan::math::var> theta(par_.size() + rate.size());
        for (size_t i = 0; i < par_.size(); ++i) theta[i] = par_[i];
        for (size_t i = 0; i < rate.size(); ++i) theta[i + par_.size()] = rate[i];
        std::vector<std::vector<scalar_type> > res_v =
          PkOdeIntegrator<It>()(f, y, t0_, ts, theta,
                                x_r, x_i, msgs, rtol, atol, max_num_steps);
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

  public:
    /*
     * solve ODE system. The different cases when @c rate is
     * @c var or data are handled by private methods @c integrate.
     *
     * @parm[in] dt next time point when result is to be *
     *           solved. The actual time point will be @c t0+dt
     * @parm[in] rtol relative tolerance for ODE solver
     * @parm[in] atol absolute tolerance for ODE solver
     * @parm[in] max_num_steps max number of steps between @c t0 and @c t0+dt.
     * @parm[in] msgs output stream.
     * @return an Eigen matrix with each row being the
     *         solution at certain time step. Hence the returned
     *         matrix is of dim (numer of time steps) x (siez of ODE system).
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& dt,
             const double rtol,
             const double atol,
             const int max_num_steps,
             std::ostream* msgs) {
      return integrate<It>(rate_, dt, rtol, atol, max_num_steps, msgs);
    }
  };

}




#endif
