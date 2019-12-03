#ifndef STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_REFACTOR_ODE_MODEL_HPP

#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/torsten/PKModel/Pred/unpromote.hpp>

namespace refactor {

  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::PMXOdeIntegrator;
  using torsten::PMXOdeIntegratorId;

  /*
   * Adaptor for ODE functor when rate is @c var and
   * appended to @c theta. When
   * constructed, this type stores an index @c index_rate of @c rate
   * within @c theta. Note that we only allow @c T_rate to
   * be @c var when @c T_par is @c var, so in this case @c
   * theta is always a @c var vector.
   */
  template<typename F, typename T_rate>
  struct PMXOdeFunctorRateAdaptor {
    PMXOdeFunctorRateAdaptor() {}

    /*
     * Evaluate ODE functor, with @c theta contains original
     * parameter vector followed by @c rate params. So the
     * @c theta passed in must be modfified to reflect this
     * data arrangement.
     */
    template <typename T0, typename T1, typename T2>
    inline std::vector<typename stan::return_type<T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      std::vector<typename stan::return_type<T1, T2>::type> res;
      res = F()(t, y, theta, x_r, x_i, msgs);
      for (size_t i = 0; i < y.size(); i++) res.at(i) += theta.at(i + theta.size() - res.size()); // NOLINT

      return res;
    }

    /*
     * when solving ODE model with @c var rate, we
     * append rate to parameter vector.
     * FIXME: spurious @c var parameters will be generated
     * if the original parameters are data.
     */
    template<typename T>
    static std::vector<stan::math::var>
    adapted_param(const std::vector<T> &par, const std::vector<T_rate> &rate) {
      std::vector<stan::math::var> theta(par.size() + rate.size());
      for (size_t i = 0; i < par.size();  ++i) theta[i] = par[i];
      for (size_t i = 0; i < rate.size(); ++i) theta[i + par.size()] = rate[i];
      return theta;
    }

    /*
     * When @c rate is @c var, the @c x_r is empty
     */
    static const std::vector<double>
    adapted_x_r(const std::vector<T_rate> &rate) {
      return {};
    }
  };

  /*
   * Specification when rate is data. In this
   * case rate should be passed in by @c x_r. 
   */
  template<typename F>
  struct PMXOdeFunctorRateAdaptor<F, double> {

    PMXOdeFunctorRateAdaptor() {}

    /*
     * Evaluate ODE functor and add @c rate data afterwards.
     */
    template <typename T0, typename T1, typename T2>
    inline std::vector<typename stan::return_type<T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& y,
               const std::vector<T2>& theta,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      std::vector<typename stan::return_type<T1, T2>::type> res;
      res = F()(t, y, theta, x_r, x_i, msgs);
      for (size_t i = 0; i < y.size(); i++) res.at(i) += x_r.at(i);
      return res;
    }

    template<typename T>
    static const std::vector<T>&
    adapted_param(const std::vector<T> &par, const std::vector<double> &rate) {
      return par;
    }

    /*
     * When @c rate is data, it's passed in as @c x_r.
     */
    static const std::vector<double>&
    adapted_x_r(const std::vector<double> &rate) {
      return rate;
    }
  };

  namespace internal {
    /*
     * depends on the type to be retrieved, we retrieve from
     * one of the two vectors. This happens when unpacking
     * from @c theta and @c x_r in ODE adaptors.
     *
     * @tparam T type to be unpacked
     */
    template<typename T>
    struct VectorUnpacker {
      template<typename T1>
      static auto& get(const Eigen::Matrix<T1, -1, 1>& v1,
                       const std::vector<double>& v2, int i) {
        return v1[i];
      }
    };

    template<>
    struct VectorUnpacker<double> {
      template<typename T1>
      static auto& get(const Eigen::Matrix<T1, -1, 1>& v1,
                       const std::vector<double>& v2, int i) {
        return v2[i];
      }
    };
  }

  /**
   * A structure to pack & unpack the algebraic system
   * which gets solved when computing the steady
   * state solution for ODE models. Note that Stan algebra
   * solver demands default constructor for passed-in system
   * functor, so we need to put all relevant info into 
   * @c theta, @c x_r, and @c x_i.
   *
   * @tparam F ODE RHS functor type
   * @tparam T_amt @c amt type
   * @tparam T_rate @c rate type
   * @tparam T_ii @c dosing interval type
   */
  template <typename F, typename T_amt, typename T_rate, typename T_ii>
  struct PMXOdeFunctorSSAdaptorPacker {
    PMXOdeFunctorSSAdaptorPacker() {}

    static constexpr bool is_var_amt = stan::is_var<T_amt>::value;
    static constexpr bool is_var_ii = stan::is_var<T_ii>::value;

    template<typename T>
    inline Eigen::Matrix<typename stan::return_type<T, T_amt, T_rate, T_ii>::type, -1, 1>
    adapted_param(const std::vector<T> &par, const T_amt& amt, const T_rate& rate, const T_ii& ii,
                  const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int ncmt_ = x_i[1];
      PMXOdeFunctorRateAdaptor<F, T_rate> f_;
      using scalar_t = typename stan::return_type<T_amt, T_rate, T_ii>::type;
      std::vector<scalar_t> vec(ncmt_, 0.0);
      vec[cmt_ - 1] = rate;
      if (is_var_amt) {
        vec.push_back(amt);
      }
      if (is_var_ii) {
        vec.push_back(ii);
      }
      return stan::math::to_vector(f_.adapted_param(par, vec));
    }

    template<typename integrator_t>
    inline const std::vector<double>
    adapted_x_r(const T_amt& amt, const T_rate& rate, const T_ii& ii,
                const std::vector<int>& x_i,
                const integrator_t& integrator) const {
      using stan::math::value_of;
      std::vector<double> res;
      if (!is_var_amt) {
        res.push_back(value_of(amt));
      }
      if (!is_var_ii) {
        res.push_back(value_of(ii));
      }
      res.push_back(integrator.rtol);
      res.push_back(integrator.atol);
      return res;
    }
    
    template<typename T1>
    const T1& unpack_rate(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                          const std::vector<double>& x_r,
                          const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int npar_ = x_i[2];
      return y(npar_ + cmt_ - 1);
    }

    template<typename T1>
    const auto& unpack_amt(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i) const {
      int ncmt_ = x_i[1];
      int npar_ = x_i[2];
      int i = is_var_amt ? npar_ + ncmt_ : 0;
      return internal::VectorUnpacker<T_amt>::get(y, x_r, i);
    }

    template<typename T1>
    const auto& unpack_ii(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                          const std::vector<double>& x_r,
                          const std::vector<int>& x_i) const {
      int ncmt_ = x_i[1];
      int npar_ = x_i[2];
      int i;
      if (is_var_amt) {
        i = is_var_ii ? npar_ + ncmt_ + 1 : 0;
      } else {
        i = is_var_ii ? npar_ + ncmt_ : 1;
      }
      return internal::VectorUnpacker<T_ii>::get(y, x_r, i);
    }

    template<typename T1>
    std::vector<T1> unpack_ode_theta(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                                     const std::vector<int>& x_i) const {
      int ncmt_ = x_i[1];
      int npar_ = x_i[2];
      std::vector<T1> theta(npar_ + ncmt_);
      for (size_t i = 0; i < theta.size(); i++) theta[i] = y(i);
      return theta;
    }

    std::vector<double> unpack_ode_x_r(const std::vector<double>& x_r,
                                       const std::vector<int>& x_i) const {
      return {};
    }

    template<typename T1>
    inline void nullify_truncated_rate(std::vector<T1>& ode_theta, std::vector<double>& ode_x_r,
                                       const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int npar_ = x_i[2];
      ode_theta[npar_ + cmt_ - 1] = 0.0;
    }    
  };

  /**
   * Partial specialization of @c PMXOdeFunctorSSAdaptorPacker:
   * @c rate is data.
   *
   * @tparam F ODE RHS functor type
   * @tparam T_amt @c amt type
   * @tparam T_ii @c dosing interval type
   */
  template <typename F, typename T_amt, typename T_ii>
  struct PMXOdeFunctorSSAdaptorPacker<F, T_amt, double, T_ii> {

    PMXOdeFunctorSSAdaptorPacker() {}

    static constexpr bool is_var_amt = stan::is_var<T_amt>::value;
    static constexpr bool is_var_ii = stan::is_var<T_ii>::value;

    /*
     * Append @c amt parameter to original parameter vector
     */
    template<typename T>
    inline Eigen::Matrix<typename stan::return_type<T, T_amt, T_ii>::type, -1, 1>
    adapted_param(const std::vector<T> &par, const T_amt& amt, double rate, const T_ii& ii,
                  const std::vector<int>& x_i) const {
      int npar_ = x_i[2];
      using scalar_t = typename stan::return_type<T, T_amt, T_ii>::type;
      PMXOdeFunctorRateAdaptor<F, scalar_t> f_;
      std::vector<scalar_t> vec;
      if (is_var_amt) {
        vec.push_back(amt);
      }
      if (is_var_ii) {
        vec.push_back(ii);
      }
      return stan::math::to_vector(f_.adapted_param(par, vec));
    }

    template<typename integrator_t>
    inline const std::vector<double>
    adapted_x_r(const T_amt& amt, double rate, const T_ii& ii,
                const std::vector<int>& x_i,
                const integrator_t& integrator) const {
      using stan::math::value_of;
      int cmt_ = x_i[0];
      int ncmt_ = x_i[1];
      std::vector<double> res(ncmt_, 0.0);
      res[cmt_ - 1] = rate;
      if (!is_var_amt) {
        res.push_back(value_of(amt));
      }
      if (!is_var_ii) {
        res.push_back(value_of(ii));
      }
      res.push_back(integrator.rtol);
      res.push_back(integrator.atol);
      return res;
   }

    template<typename T1>
    const double& unpack_rate(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                              const std::vector<double>& x_r,
                              const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      return x_r[cmt_ - 1];
    }

    template<typename T1>
    const auto& unpack_amt(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i) const {
      int npar_ = x_i[2];
      int ncmt_ = x_i[1];
      int i = is_var_amt ? npar_ : ncmt_;
      return internal::VectorUnpacker<T_amt>::get(y, x_r, i);
    }

    template<typename T1>
    const auto& unpack_ii(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                            const std::vector<double>& x_r,
                        const std::vector<int>& x_i) const {
      int npar_ = x_i[2];
      int ncmt_ = x_i[1];
      int i;
      if (is_var_amt) {
        i = is_var_ii ? npar_ + 1 : ncmt_;
      } else {
        i = is_var_ii ? npar_ : ncmt_ + 1;
      }
      return internal::VectorUnpacker<T_ii>::get(y, x_r, i);
    }
    
    template<typename T1>
    std::vector<T1> unpack_ode_theta(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                                     const std::vector<int>& x_i) const {
      int npar_ = x_i[2];
      std::vector<T1> theta(npar_);
      for (size_t i = 0; i < theta.size(); i++) theta[i] = y(i);
      return theta;
    }

    std::vector<double> unpack_ode_x_r(const std::vector<double>& x_r,
                                       const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int ncmt_ = x_i[1];
      std::vector<double> rate(ncmt_, 0.0);
      rate[cmt_ - 1] = x_r[cmt_ - 1];
      return rate;
    }

    template<typename T1>
    inline void nullify_truncated_rate(std::vector<T1>& ode_theta,
                                       std::vector<double>& ode_x_r,
                                       const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      ode_x_r[cmt_ - 1] = 0.0;
    }
  };

  /**
   * A structure to store the algebraic system
   * which gets solved when computing the steady
   * state solution for ODE models.
   *
   * @tparam It integrator type
   * @tparam T_amt @c amt type
   * @tparam T_rate @c rate type
   * @tparam T_ii @c dosing interval type
   * @tparam F ODE RHS functor type
   */
  template <PMXOdeIntegratorId It, typename T_amt, typename T_rate, typename T_ii, typename F>
  struct PMXOdeFunctorSSAdaptor {

    PMXOdeFunctorSSAdaptor() {}

    /**
     * When rate is RV, it's attached to parameters vector.
     * IN this case parameter @c y consists of {theta, rate}
     */
    template <typename T0, typename T1>
    inline
    Eigen::Matrix<typename torsten::return_t<T0, T1>::type,
                  Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;

      typedef typename stan::return_type<T0, T1>::type scalar_t;
      const PMXOdeFunctorSSAdaptorPacker<F, T_amt, T_rate, T_ii> packer_;
      const PMXOdeIntegrator<It> integrator_(*(x_r.rbegin() + 1), x_r.back(), x_i.back(), msgs);
      const PMXOdeFunctorRateAdaptor<F, T_rate> f_;
      int cmt_ = x_i[0];
      int ncmt_ = x_i[1];
      int npar_ = x_i[2];

      const auto& ii_ = packer_.unpack_ii(y, x_r, x_i);
      const auto& rate = packer_.unpack_rate(y, x_r, x_i);
      const auto& amt = packer_.unpack_amt(y, x_r, x_i);

      double t0 = 0;

      std::vector<scalar_t> x0(x.data(), x.data() + x.size());
      std::vector<T1> ode_theta(packer_.unpack_ode_theta(y, x_i));
      std::vector<double> ode_x_r(packer_.unpack_ode_x_r(x_r, x_i));

      Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> result(x.size());

      static const char* function("Steady State Event");

      if (rate == 0) {  // bolus dose
        x0[cmt_ - 1] += amt;
        std::vector<scalar_t> pred = integrator_(f_, x0, t0, ii_, ode_theta, ode_x_r, x_i)[0];
        for (int i = 0; i < result.size(); i++) {
          result(i) = x(i) - pred[i];
        }
      } else if (ii_ > 0) {  // multiple truncated infusions
        T1 dt = amt / rate;
        torsten::check_mti(amt, dt, ii_, function);

        x0 = integrator_(f_, to_array_1d(x), t0, dt, ode_theta, ode_x_r, x_i)[0];

        dt = ii_ - dt;
        packer_.nullify_truncated_rate(ode_theta, ode_x_r, x_i);
        std::vector<scalar_t> pred = integrator_(f_, x0, t0, dt, ode_theta, ode_x_r, x_i)[0];
        for (int i = 0; i < result.size(); i++) result(i) = x(i) - pred[i];
      } else {  // constant infusion
        stan::math::check_less_or_equal(function, "AMT", amt, 0);
        std::vector<scalar_t> derivative = f_(0, to_array_1d(x), ode_theta, ode_x_r, x_i, 0);
        result = to_vector(derivative);
      }

      return result;
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
  template<typename T_time, typename T_init, typename T_rate, typename T_par, typename F> // NOLINT
  class PKODEModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const std::vector<T_par> &par_;
    const F &f_;
    const PMXOdeFunctorRateAdaptor<F, T_rate> f1;
    const int ncmt_;
  public:
    using scalar_type = typename torsten::return_t<T_time, T_rate, T_par, T_init>::type; // NOLINT
    using aug_par_type = typename torsten::return_t<T_rate, T_par, T_init>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;
    using f_type      = F;

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
               const F& f) :
      t0_(t0), y0_(y0), rate_(rate), par_(par), f_(f), f1(), ncmt_(y0.size()) // NOLINT
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
      f1(),
      ncmt_(m.ncmt())
    {}
    
    /*
     * The number of @c var that will be in the ODE
     * integrator. Note that not all @c var will be
     * explicitly in an integrator's call signature.
     * Since we only step one time-step, if @c t0_ is @c var
     * it only adds one(because @c ts will be of size one). Also note that regardless if @c
     * par_ is @c var or not, when @c rate_ is @c var we
     * generate a new @c theta of @c var vector to pass to ODE integrator that
     * contains both @c par_ and @c rate_.
     */

    /*
     * calculate number of @c vars for transient dosing.
     */
    static int nvars(int ncmt, int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_time>::value) n++; // t0
      if (is_var<T_init>::value) n += ncmt;
      if (is_var<T_rate>::value) {
        n += ncmt + npar;
      } else if (is_var<T_par>::value) {
        n += npar;
      }
      return n;
    }

    /*
     * calculate number of @c vars for steady-state dosing.
     */
    template<typename T_a, typename T_r, typename T_ii>
    static int nvars(int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_a>::value) n++; // amt
      if (is_var<T_r>::value) n++; // rate
      if (is_var<T_ii>::value) n++; // ii
      if (is_var<T_par>::value) n += npar;
      return n;
    }

    template<typename T0, typename T1, typename T2, typename T3>
    static int nvars(const T0& t0,
                     const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                     const std::vector<T2> &rate,
                     const std::vector<T3> &par) {
      using stan::is_var;
      int res = 0;
      if (is_var<T0>::value) res += 1;
      if (is_var<T1>::value) res += y0.size();
      if (is_var<T2>::value) {
        res += rate.size() + par.size();
      } else if (is_var<T3>::value) {
        res += par.size();
      }
      return res;
    }

    /*
     * generate @c var vector that consisting of all the
     * parameters used passed to an ODE integrator. The output
     * will be arranged in order (y0, par, rate, dt)
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static std::vector<stan::math::var> vars(const T0& t1,
                                             const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                                             const std::vector<T2> &rate,
                                             const std::vector<T3> &par) {
      using stan::math::to_var;      
      using stan::is_var;      
      using stan::math::var;
      std::vector<stan::math::var> res(nvars(t1, y0, rate, par));
      int i = 0;
      if (is_var<T1>::value) {
        for (int j = 0; j < y0.size(); ++j) {
          res[i] = y0(j);
          i++;
        }
      }
      if (is_var<T2>::value) {
        std::vector<stan::math::var> theta(to_var(PMXOdeFunctorRateAdaptor<F, T2>::adapted_param(par, rate)));
        for (size_t j = 0; j < theta.size(); ++j) {
          res[i] = theta[j];
          i++;
        }
      } else if (is_var<T3>::value) {
        for (size_t j = 0; j < par.size(); ++j) {
          res[i] = par[j];
          i++;
        }
      }

      if (is_var<T0>::value) {
        res[i] = t1;
      }
      return res;
    }

    /*
     * Calculate the size of the entire system, with ODe
     * solutions and their gradients w.r.t. the parameters.
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static int n_sys(const T0& t0,
                     const Eigen::Matrix<T1, 1, Eigen::Dynamic>& y0,
                     const std::vector<T2> &rate,
                     const std::vector<T3> &par) {
      using stan::is_var;
      int n = nvars(t0, y0, rate, par);
      return y0.size() * (n + 1);
    }

    /*
     * calculate number of vars with constructed data
     */
    template<typename T0>
    int nvars(const T0& t1) {
      return nvars(t1, y0_, rate_, par_);
    }

    /*
     * return the number @c var that will be the parameters
     * of the stead-state dosing event's solution
     */
    template<typename T_a, typename T_r, typename T_ii>
    int nvars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::pk_nvars(a, r, ii, par_);
    }

    /*
     * return @c vars that will be in integrator calls
     */
    template<typename T0>
    std::vector<stan::math::var> vars(const T0 t1) {
      return vars(t1, y0_, rate_, par_);
    }

    /*
     * return @c vars that will be steady-state
     * solution. For SS solution @c rate_ or @ y0_ will not
     * be in the solution.
     */
    template<typename T_a, typename T_r, typename T_ii>
    std::vector<stan::math::var> vars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::dsolve::pk_vars(a, r, ii, par_);
    }

    /*
     * calculate size of the entire system
     */
    int n_sys() const {
      return n_sys(t0_, y0_, rate_, par_);
    }

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
     * When time is param, the time step needs to
     * incorporate information of initial time.
     */
    std::vector<stan::math::var> time_step(const stan::math::var& t_next) const {
      return {stan::math::value_of(t0_) + t_next - t0_};
    }

    /*
     * When time is data, the time step is trivial
     */
    std::vector<double> time_step(const double t_next) const {
      return {t_next};
    }

    /*
     * For steady-state with data rate. 
     * Same as the non-steady version but with given @c init.
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<double, Eigen::Dynamic, 1>
    integrate(const std::vector<double> &rate,
              const Eigen::Matrix<double, 1, Eigen::Dynamic>& y0,
              const double& dt,
              const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;

      const double t0 = value_of(t0_) - dt;
      std::vector<double> ts{value_of(t0_)};
      Eigen::Matrix<double, Eigen::Dynamic, 1> res;
      if (ts[0] == t0) {
        res = y0;
      } else {
        auto y = stan::math::to_array_1d(y0);
        PMXOdeFunctorRateAdaptor<F, double> f;
        std::vector<int> x_i;
        const std::vector<double> pars{value_of(par_)};
        std::vector<std::vector<double> > res_v =
          integrator(f, y, t0, ts, pars, rate, x_i); // NOLINT
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

  public:
    /*
     * solve ODE system. The different cases when @c rate is
     * @c var or data are handled by private methods @c integrate.
     *
     * @parm[in] t_next next time point when result is to be
     *           solved. The actual time point will be @c t0+dt
     * @parm[in] rtol relative tolerance for ODE solver
     * @parm[in] atol absolute tolerance for ODE solver
     * @parm[in] max_num_steps max number of steps between @c t0 and @c t0+dt.
     * @parm[in] msgs output stream.
     * @return an Eigen matrix with each row being the
     *         solution at certain time step. Hence the returned
     *         matrix is of dim (numer of time steps) x (siez of ODE system).
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& t_next,
          const PMXOdeIntegrator<It>& integrator) const {
      const double t0 = stan::math::value_of(t0_);
      std::vector<T_time> ts(time_step(t_next));
      Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> res;
      if (ts[0] == t0_) {
        res = y0_;
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<int> x_i;
        std::vector<std::vector<scalar_type> > res_v =
          integrator(f1, y, t0, ts,
                     f1.adapted_param(par_, rate_),
                     f1.adapted_x_r(rate_),
                     x_i);
        res = stan::math::to_vector(res_v[0]);
      }
      return res;
    }

    /*
     * Solve the ODE but return the results in form of data.
     * The default behavior is defined in function template
     * @c model_solve_d(), using autodiff to recalculate gradients.
     * Some integrators, such as @c PkBdf, have their own implementation
     * that can return data directly, so we skip them.
     */
    template<PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::StanBdf || It == torsten::StanAdams || It == torsten::StanRk45>* = nullptr>
    Eigen::VectorXd solve_d(const T_time& t_next,
                            const PMXOdeIntegrator<It>& integrator) const
    {
      static const char* caller = "PKODEModel::solve";
      stan::math::check_greater(caller, "time step", t_next, t0_);

      return torsten::model_solve_d(*this, t_next, integrator);
    }

    /*
     * @c PkBdf can return results in form of data directly,
     * thanks to @c pk_cvodes_integrator implementation.
     */
    template<PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::PkBdf || It == torsten::PkAdams || It == torsten::PkRk45>* = nullptr>
    Eigen::VectorXd solve_d(const T_time& t_next,
                            const PMXOdeIntegrator<It>& integrator) const {
      static const char* caller = "PMXOdeModel::solve_d";
      stan::math::check_greater(caller, "next time", t_next, t0_);

      // return integrate_d(rate_, t_next, integrator);
      using stan::math::var;
      using stan::math::value_of;
      using stan::math::to_var;

      const double t0 = value_of(t0_);
      std::vector<T_time> ts(time_step(t_next));
      Eigen::VectorXd res(n_sys());
      if (ts[0] == t0_) {
        res = stan::math::value_of(y0_);
        res.resize(n_sys());
      } else {
        auto y = stan::math::to_array_1d(y0_);
        std::vector<int> x_i;
        res = integrator.solve_d(f1, y, t0, ts,
                                 f1.adapted_param(par_, rate_),
                                 f1.adapted_x_r(rate_),
                                 x_i).col(0);
      }
      return res;
    }

    /**
     * Solve steady state ODE-based PKPD model with @c amt
     * and @c rate both being data. Specifically for
     * constant rate we solve ODEs of the form
     *
     * x′(t)=F(x,t,θ)+ Rate,
     *
     * and solution x satisfies
     *
     * F(x,t,θ) + Rate = 0.
     *
     * For the periodic steady-state resulting from a
     * periodic input D into one or more cpts at a constant
     * interval τ, the root finder is combined with the ODE
     * solver to solve
     *
     * x(t0 + τ) − x(t0) = 0.
     *
     * Denote x(t0) by x0, it is equivalent to solve x0 in
     *
     * x(t0 + τ) = x0, where in the interval (t0, t0 + τ)
     * x(t) satisfies the ODE with x0 as init condition.
     *
     * @tparam It ODE numerical integrator ID.
     * @tparam T_ii dosing interval type
     * @param amt dosing amount
     * @param rate dosing rate
     * @param ii dosing interval
     * @param cmt compartment where the dosing occurs
     * @param rtol relative tolerance of ODE integrator.
     * @param atol absolute tolerance of ODE integrator.
     * @param max_num_steps max number of steps of ODE integrator.
     * @param msgs @c ostream output of ODE integrator.
     * @return col vector of the steady state ODE solution
     * Solve steady state ODE-based PKPD model with @c amt
     * being param and @c rate both being data.
     *
     */
    template<PMXOdeIntegratorId It, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename torsten::return_t<T_amt, T_r, T_par, T_ii>::type, Eigen::Dynamic, 1> // NOLINT
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;
      using stan::math::algebra_solver;

      typedef typename torsten::return_t<T_amt, T_r, T_par, T_ii>::type scalar;

      double ii_dbl = value_of(ii);
      Eigen::Matrix<double, 1, -1> init_dbl(Eigen::Matrix<double, 1, -1>::Zero(ncmt_));
      std::vector<double> rate_vec(ncmt_, 0);
      std::vector<int> x_i{cmt, ncmt_, int(par_.size()), int(integrator.max_num_step)};

      if (rate == 0) {                     // bolus dose
        init_dbl(cmt - 1) = value_of(amt); // bolus as initial condition
      } else {                             // infusion
        rate_vec[cmt - 1] = value_of(rate);
      }

      const double init_dt = (rate == 0.0 || ii > 0) ? ii_dbl : 24.0;
      PMXOdeFunctorSSAdaptor<It, T_amt, T_r, T_ii, F> fss;
      PMXOdeFunctorSSAdaptorPacker<F, T_amt, T_r, T_ii> packer;
      return algebra_solver(fss, integrate(rate_vec, init_dbl, init_dt, integrator),
                            packer.adapted_param(par_, amt, rate, ii, x_i),
                            packer.adapted_x_r(amt, rate, ii, x_i, integrator),
                            x_i, 0,
                            integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step);
    }

    /*
     * return steady state solution in form of data, use
     * default behavior, namely take gradients using autodiff.
     */
    template<PMXOdeIntegratorId It, typename T_amt, typename T_r, typename T_ii>
    Eigen::VectorXd solve_d(const T_amt& amt,
                            const T_r& rate,
                            const T_ii& ii,
                            const int& cmt,
                            const PMXOdeIntegrator<It>& integrator) const {
      return torsten::model_solve_d(*this, amt, rate, ii, cmt, integrator);
    }
  };
}

#endif
