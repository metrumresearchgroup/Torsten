#ifndef STAN_MATH_TORSTEN_ODE_MODEL_HPP
#define STAN_MATH_TORSTEN_ODE_MODEL_HPP

#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/rev/fun/to_var.hpp>
#include <stan/math/rev/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/functor/algebra_solver_newton.hpp>
#include <stan/math/rev/functor/algebra_solver_fp.hpp>
#include <stan/math/prim/err/check_less_or_equal.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>

namespace torsten {

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

  namespace {
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

    /** 
     * Formulate <code>var</code> vector for Steady state
     * system. The vector should include <code>var</code> variables in
     * ODE's params, along with SS system params. In particular, since
     * the current template takes <code>T_rate</code> as
     * <code>var</code>, the return vector contains SS dosing rate.
     * 
     * @param par ODE's params
     * @param amt SS dosing amount
     * @param rate SS dosing rate
     * @param ii SS dosing inverval
     * @param x_i integer vector that contains system & dosing information
     * 
     * @return <code>var</code> vector of params for SS system.
     */
    template<typename T>
    inline Eigen::Matrix<typename stan::return_type_t<T, T_amt, T_rate, T_ii>, -1, 1>
    adapted_param(const std::vector<T> &par, const T_amt& amt, const T_rate& rate, const T_ii& ii,
                  const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int ncmt_ = x_i[1];
      PMXOdeFunctorRateAdaptor<F, T_rate> f_;
      using scalar_t = typename stan::return_type_t<T_amt, T_rate, T_ii>;
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

    /** 
     * Formulate <code>double</code> vector for Steady state
     * system, to be passed into algebra solver. The vector should
     * include data on
     * ODE's params, along with SS system params. In particular, since
     * the current template takes <code>T_rate</code> as
     * <code>var</code>, the return vector does not contains SS dosing rate.
     * 
     * @param amt 
     * @param rate 
     * @param ii 
     * @param x_i 
     * @param integrator 
     * 
     * @return 
     */
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
    
    /** 
     * Get dosing rate out of param vector <code>y</code>
     * 
     * @param y param vector
     * @param x_r real data
     * @param x_i integer data
     * 
     * @return dosing rate
     */
    template<typename T1>
    const T1& unpack_rate(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                          const std::vector<double>& x_r,
                          const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int npar_ = x_i[2];
      return y(npar_ + cmt_ - 1);
    }

    /** 
     * Get dosing amt out of param vector <code>y</code> 
     * or <code>x_r</code>
     * 
     * @param y param vector
     * @param x_r real data
     * @param x_i integer data
     * 
     * @return dosing amt
     */
    template<typename T1>
    const auto& unpack_amt(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                         const std::vector<double>& x_r,
                         const std::vector<int>& x_i) const {
      int ncmt_ = x_i[1];
      int npar_ = x_i[2];
      int i = is_var_amt ? npar_ + ncmt_ : 0;
      return VectorUnpacker<T_amt>::get(y, x_r, i);
    }

    /** 
     * Get dosing interval out of param vector <code>y</code>
     * or <code>x_r</code>
     * 
     * @param y param vector
     * @param x_r real data
     * @param x_i integer data
     * 
     * @return dosing interval
     */
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
      return VectorUnpacker<T_ii>::get(y, x_r, i);
    }

    /** 
     * Get ODE param out of packed SS system param vector. Since here 
     * <code>rate</code> is <code>var</code>, the return vector is
     * original ODE param followed by dosing rates.
     * 
     * @param y packed SS system param vector generated 
     * from <code>adapted_param</code>
     * @param x_i integer data
     * 
     * @return ODE param vector including dosing rate
     */
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

    template<typename T1>
    void scale_rate(double n,
                    std::vector<T1>& ode_theta,
                    const std::vector<double>& ode_x_r,
                    const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      int npar_ = x_i[2];
      ode_theta[npar_ + cmt_ - 1] *= n;
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
      return VectorUnpacker<T_amt>::get(y, x_r, i);
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
      return VectorUnpacker<T_ii>::get(y, x_r, i);
    }
    
    /** 
     * Get ODE param out of packed SS system param vector. Since here 
     * <code>rate</code> is <code>double</code>, the return vector is
     * just original ODE params.
     * 
     * @param y packed SS system param vector generated 
     * from <code>adapted_param</code>
     * @param x_i integer data
     * 
     * @return ODE param vector
     */
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

    template<typename T1>
    void scale_rate(double n,
                    const std::vector<T1>& ode_theta,
                    std::vector<double>& ode_x_r,
                    const std::vector<int>& x_i) const {
      int cmt_ = x_i[0];
      ode_x_r[cmt_ - 1] *= n;
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
    Eigen::Matrix<typename stan::return_type_t<T0, T1>,
                  Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      using stan::math::value_of;

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

      double t0 = 0;            // FIXME: ODE explicitly depdends on time

      std::vector<scalar_t> x0(x.data(), x.data() + x.size());
      Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> result(x.size());

      static const char* func("Steady State Event");

      if (rate == 0) {  // bolus dose
        std::vector<T1> ode_theta(packer_.unpack_ode_theta(y, x_i));
        std::vector<double> ode_x_r(packer_.unpack_ode_x_r(x_r, x_i));

        x0[cmt_ - 1] += amt;
        std::vector<scalar_t> pred = integrator_(f_, x0, t0, ii_, ode_theta, ode_x_r, x_i)[0];
        for (int i = 0; i < result.size(); i++) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else if (ii_ > 0) {  // multiple truncated infusions
        T1 dt = amt / rate;
        int n = int(std::floor(value_of(dt) / value_of(ii_)) + 0.1);

        auto dt1 = dt - n * ii_;
        std::vector<T1> ode_theta1(packer_.unpack_ode_theta(y, x_i));
        std::vector<double> ode_x_r1(packer_.unpack_ode_x_r(x_r, x_i));
        packer_.scale_rate(n + 1, ode_theta1, ode_x_r1, x_i);
        x0 = integrator_(f_, to_array_1d(x), t0, dt1, ode_theta1, ode_x_r1, x_i)[0];

        auto dt2 = (n + 1) * ii_ - dt;
        std::vector<T1> ode_theta2(packer_.unpack_ode_theta(y, x_i));
        std::vector<double> ode_x_r2(packer_.unpack_ode_x_r(x_r, x_i));
        packer_.scale_rate(n, ode_theta2, ode_x_r2, x_i);
        std::vector<scalar_t> pred = integrator_(f_, x0, t0, dt2, ode_theta2, ode_x_r2, x_i)[0];
        for (int i = 0; i < result.size(); i++) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else {  // constant infusion
        std::vector<T1> ode_theta(packer_.unpack_ode_theta(y, x_i));
        std::vector<double> ode_x_r(packer_.unpack_ode_x_r(x_r, x_i));
        stan::math::check_less_or_equal(func, "AMT", amt, 0);
        std::vector<scalar_t> derivative = f_(0, to_array_1d(x), ode_theta, ode_x_r, x_i, 0);
#ifdef TORSTEN_AS_FP
        stan::math::throw_domain_error(func, "algebra_solver_fp used for ", 1, "constant infusion");
#else
        result = to_vector(derivative);
#endif
      }

      return result;
    }
  };

  /**
   * ODE-based PKPD models.
   *
   * @tparam T_time t type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   * @tparam F ODE functor
   * @tparam Ti ODE additional parameter type, usually the ODE size
   */
  template<typename T_par, typename F> // NOLINT
  class PKODEModel {
    const std::vector<T_par> &par_;
    const F &f_;
    const int ncmt_;
  public:
    using par_type    = T_par;
    using f_type      = F;

    /**
     * Constructor
     *
     * @param t0 initial time
     * @param rate dosing rate
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par,
               int ncmt,
               const F& f) :
      par_(par), f_(f), ncmt_(ncmt)
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
      par_(m.par()),
      f_(m.f()),
      ncmt_(m.ncmt())
    {}
    
    template<typename T0, typename T1, typename T2, typename T3>
    static int nvars(const T0& t0,
                     const PKRec<T1>& y0,
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
                                             const PKRec<T1>& y0,
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
    template<typename Tt0>
    std::vector<stan::math::var> time_step(const Tt0& t0, const stan::math::var& t1) const {
      return {stan::math::value_of(t0) + t1 - t0};
    }

    /*
     * When time is data, the time step is trivial
     */
    template<typename Tt0>
    std::vector<double> time_step(const Tt0& t0, const double t1) const {
      return {t1};
    }

    /*
     * For steady-state with data rate. 
     * Same as the non-steady version but with given @c init.
     */
    template<PMXOdeIntegratorId It>
    Eigen::Matrix<double, Eigen::Dynamic, 1>
    integrate(double t1,
              const std::vector<double> &rate,
              const Eigen::Matrix<double, 1, Eigen::Dynamic>& y0,
              const double& dt,
              const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;

      const double t0 = t1 - dt;
      std::vector<double> ts{t1};
      // const double t0 = value_of(t0_) - dt;
      // std::vector<double> ts{value_of(t0_)};
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
    template<typename Tt0, typename Tt1, typename T, typename T1, PMXOdeIntegratorId It>
    void solve(Eigen::Matrix<T, -1, 1>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const PMXOdeIntegrator<It>& integrator) const {
      const double t0_d = stan::math::value_of(t0);
      std::vector<Tt1> ts(time_step(t0, t1));
      PMXOdeFunctorRateAdaptor<F, T1> f_rate;
      if (t1 > t0) {
        auto y_vec = stan::math::to_array_1d(y);
        std::vector<int> x_i;
        std::vector<std::vector<T> > res_v =
          integrator(f_rate, y_vec, t0_d, ts,
                     f_rate.adapted_param(par_, rate),
                     f_rate.adapted_x_r(rate),
                     x_i);
        y = stan::math::to_vector(res_v[0]);
      }
    }

    /*
     * <code>PKBdf, PKAdams, PKRk45</code> integrators can return 
     * results in form of data directly,
     * thanks to @c pk_cvodes_integrator implementation.
     */
    template<typename T0, typename T, typename T1, PMXOdeIntegratorId It,
             typename std::enable_if_t<It == torsten::PkBdf || It == torsten::PkAdams || It == torsten::PkRk45>* = nullptr>
    void solve_d(Eigen::VectorXd& yd,
                 const PKRec<T>& y,
                 const T0& t0, const T0& t1,
                 const std::vector<T1>& rate,
                 const PMXOdeIntegrator<It>& integrator) const {
      // static const char* caller = "PMXOdeModel::solve_d";
      // stan::math::check_greater(caller, "next time", t1, t0);

      using stan::math::var;
      using stan::math::value_of;
      using stan::math::to_var;

      const double t0_d = value_of(t0);
      std::vector<T0> ts(time_step(t0, t1));
      PMXOdeFunctorRateAdaptor<F, T1> f_rate;
      if (t1 > t0) {
        auto y1d = stan::math::to_array_1d(y);
        std::vector<int> x_i;
        yd = integrator.solve_d(f_rate, y1d, t0_d, ts,
                                f_rate.adapted_param(par_, rate),
                                f_rate.adapted_x_r(rate),
                                x_i).col(0);
      }
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
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_par, T_ii>, Eigen::Dynamic, 1> // NOLINT
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const PMXOdeIntegrator<It>& integrator) const {
      using stan::math::value_of;
      using stan::math::algebra_solver_powell;
      using stan::math::algebra_solver_newton;
      using stan::math::algebra_solver_fp;

      typedef typename stan::return_type_t<T_amt, T_r, T_par, T_ii> scalar;

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
#ifdef TORSTEN_AS_POWELL
      return algebra_solver_powell(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                                   packer.adapted_param(par_, amt, rate, ii, x_i),
                                   packer.adapted_x_r(amt, rate, ii, x_i, integrator),
                                   x_i, 0,
                                   integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step);
#elif defined(TORSTEN_AS_FP)
      std::vector<double> u_scale(ncmt_, 1.0);
      std::vector<double> f_scale(ncmt_, 1.0);
      return algebra_solver_fp(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                               packer.adapted_param(par_, amt, rate, ii, x_i),
                               packer.adapted_x_r(amt, rate, ii, x_i, integrator),
                               x_i, u_scale, f_scale, 0,
                               integrator.as_atol, integrator.as_max_num_step);
#else
      return algebra_solver_newton(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                                   packer.adapted_param(par_, amt, rate, ii, x_i),
                                   packer.adapted_x_r(amt, rate, ii, x_i, integrator),
                                   x_i, 0,
                                   integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step);
#endif
    }
  };
}

#endif
