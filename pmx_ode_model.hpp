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
#include <stan/math/torsten/dsolve/pmx_algebra_solver_newton.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>

namespace torsten {

  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::dsolve::PMXOdeIntegrator;

  /** 
   * Rate adapter for ODE
   * 
   * @tparam T_theta type of original ODE param
   * @tparam T_rate  type of rate
   * 
   */
  template<typename T_theta, typename T_rate>
  struct PMXOdeFunctorRateAdaptorImpl {
    const std::vector<T_theta>& theta_;
    const std::vector<T_rate>& rate_;

    PMXOdeFunctorRateAdaptorImpl(const std::vector<T_theta>& theta,
                                     const std::vector<T_rate>& rate) :
      theta_(theta), rate_(rate) {}

    /** 
     * generate adapted ODE parameter vector for adatepd ODE
     * 
     * @return adapted parameter: theta_ & rate_
     */
    const std::vector<typename stan::return_type_t<T_theta, T_rate>>
    adapted_param() const {
      std::vector<typename stan::return_type_t<T_theta, T_rate>> res(theta_);
      res.insert(res.end(), rate_.begin(), rate_.end());
      return res;
    }

    /** 
     * original ode param
     * 
     * @return original ode param vector
     */
    template<typename t2>
    const std::vector<t2> ode_par(const std::vector<t2>& theta) const {
      return {theta.begin(), theta.begin() + theta_.size()};
    }
    
    /** 
     * add rate to original ode result
     * 
     * @param res original ode result
     * @param theta adapted ode param that contains theta & rate
     */
    template<typename t, typename t2>
    void apply_rate(std::vector<t>& res, const std::vector<t2>& theta) const {
      for (size_t i = 0; i < res.size(); i++) {
        res.at(i) += theta.at(i + theta.size() - res.size()); 
      }
    }
  };

  /** 
   * rate adapter specification when original ode param is data
   * 
   * @tparam t_rate  type of rate, must be <code>var</code>
   * 
   */
  template<>
  struct PMXOdeFunctorRateAdaptorImpl<double, stan::math::var> {
    
    const std::vector<double>& theta_;
    const std::vector<stan::math::var>& rate_;

    PMXOdeFunctorRateAdaptorImpl(const std::vector<double>& theta,
                                     const std::vector<stan::math::var>& rate) :
      theta_(theta), rate_(rate) {}

    /** 
     * generate adapted ode parameter vector for adatepd ode
     * 
     * @return adapted parameter: rate_
     */
    const std::vector<stan::math::var>& adapted_param() const {
      return rate_;
    }

    /** 
     * original ode param
     * 
     * 
     * @return rate vector
     */
    template<typename T2>
    const std::vector<double>& ode_par(const std::vector<T2>& theta) const {
      return theta_;
    }
    
    /** 
     * add rate to original ode result
     * 
     * @param res original ode result
     * @param theta adapted ode param that equals to rate
     */
    template<typename T, typename T2>
    void apply_rate(std::vector<T>& res, const std::vector<T2>& theta) const {
      for (size_t i = 0; i < res.size(); i++) {
        res.at(i) += theta.at(i + theta.size() - res.size());
      }
    }
  };

  /** 
   * rate adapter specification when rate is data
   * 
   * @tparam t_theta  type of original ode param
   * 
   */
  template<typename T_theta>
  struct PMXOdeFunctorRateAdaptorImpl<T_theta, double> {
    const std::vector<T_theta>& theta_;
    const std::vector<double>& rate_;

    PMXOdeFunctorRateAdaptorImpl(const std::vector<T_theta>& theta,
                                 const std::vector<double>& rate) :
      theta_(theta), rate_(rate) {}

    /** 
     * generate adapted ode parameter vector for adatepd ode
     * 
     * @return adapted parameter: rate_
     */
    const std::vector<T_theta>& adapted_param() const {
      return theta_;
    }

    /** 
     * original ode param
     * 
     * 
     * @return rate vector
     */
    template<typename T2>
    const std::vector<T2>& ode_par(const std::vector<T2>& theta)const {
      return theta;
    }
    
    /** 
     * add rate to original ODE result
     * 
     * @param res original ODE result
     * @param theta adapted ODE param that equals to rate
     */
    template<typename T, typename T2>
    void apply_rate(std::vector<T>& res, const std::vector<T2>& theta) const {
      for (size_t i = 0; i < res.size(); i++) {
        res.at(i) += rate_[i];
      }
    }
  };

  /** 
   * Functor class that solves ODE & apply infusion dosing.
   */
  template<typename F, typename T_theta, typename T_rate>
  struct PMXOdeFunctorRateAdaptor {

    PMXOdeFunctorRateAdaptorImpl<T_theta, T_rate> adaptor;

    PMXOdeFunctorRateAdaptor(const std::vector<T_theta>& theta,
                                const std::vector<T_rate>& rate) :
      adaptor(theta, rate) {}

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
      auto fy = F()(t, y, adaptor.ode_par(theta), x_r, x_i, msgs);
      std::vector<typename stan::return_type<T1, T2>::type>
        res(fy.begin(), fy.end());
      adaptor.apply_rate(res, theta);

      return res;
    }
  };

  template<typename T_theta>
  struct PMXOdeFunctorSSAdaptorComponentTheta {
    const std::vector<T_theta>& theta_;    
    
    PMXOdeFunctorSSAdaptorComponentTheta(const std::vector<T_theta>& theta)
      : theta_(theta)
    {}

    template<typename T>
    inline std::vector<T> operator()(const std::vector<T>& theta) const {
      return {theta.begin(), theta.begin() + theta_.size()};
    }
  };

  template<>
  struct PMXOdeFunctorSSAdaptorComponentTheta<double> {
    const std::vector<double>& theta_;    
    
    PMXOdeFunctorSSAdaptorComponentTheta(const std::vector<double>& theta)
      : theta_(theta)
    {}

    template<typename T>
    inline const std::vector<double>& operator()(const std::vector<T>& theta) const {
      return theta_;
    }
  };

  template<typename T_theta, typename T_amt>
  struct PMXOdeFunctorSSAdaptorComponentAmt {
    const T_amt& amt_;    
    const int npar_;
    
    PMXOdeFunctorSSAdaptorComponentAmt(const T_amt& amt, int npar)
      : amt_(amt), npar_(npar)
    {}

    template<typename T>
    inline const T& operator()(const std::vector<T>& theta) const {
      if (stan::is_var<T_theta>::value) {
        return theta[npar_];
      } else {
        return theta[0];
      }
    }
  };

  template<typename T_theta>
  struct PMXOdeFunctorSSAdaptorComponentAmt<T_theta, double> {
    const double amt_;    
    const int npar_;
    
    PMXOdeFunctorSSAdaptorComponentAmt(const double amt, int npar)
      : amt_(amt), npar_(npar)
    {}

    template<typename T>
    inline const double operator()(const std::vector<T>& theta) const {
      return amt_;
    }
  };

  template<typename T_theta, typename T_amt, typename T_rate>
  struct PMXOdeFunctorSSAdaptorComponentRate {
    const T_rate& rate_;    
    const int npar_;
    
    PMXOdeFunctorSSAdaptorComponentRate(const T_rate& rate, int npar)
      : rate_(rate), npar_(npar)
    {}

    template<typename T>
    inline const T& operator()(const std::vector<T>& theta) const {
      if (stan::is_var<T_theta>::value) {
        return (stan::is_var<T_amt>::value) ? theta[npar_ + 1] : theta[npar_];
      } else {
        return (stan::is_var<T_amt>::value) ? theta[1] : theta[0];
      }
    }
  };

  template<typename T_theta, typename T_amt>
  struct PMXOdeFunctorSSAdaptorComponentRate<T_theta, T_amt, double> {
    const double rate_;
    const int npar_;
    
    PMXOdeFunctorSSAdaptorComponentRate(const double rate, int npar)
      : rate_(rate), npar_(npar)
    {}

    template<typename T>
    inline const double operator()(const std::vector<T>& theta) const {
      return rate_;
    }
  };

  template<typename T_theta, typename T_amt, typename T_rate, typename T_ii>
  struct PMXOdeFunctorSSAdaptorComponentII {
    const T_ii& ii_;
    const int npar_;
    
    PMXOdeFunctorSSAdaptorComponentII(const T_ii& ii, int npar)
      : ii_(ii), npar_(npar)
    {}

    template<typename T>
    inline const T& operator()(const std::vector<T>& theta) const {
      if (stan::is_var<T_theta>::value) {
        return (stan::is_var<T_amt>::value) ?
          ((stan::is_var<T_rate>::value) ? theta[npar_ + 2] : theta[npar_ + 1]) :
          ((stan::is_var<T_rate>::value) ? theta[npar_ + 1] : theta[npar_]);
      } else {
        return (stan::is_var<T_amt>::value) ?
          ((stan::is_var<T_rate>::value) ? theta[2] : theta[1]) :
          ((stan::is_var<T_rate>::value) ? theta[1] : theta[0]);
      }
    }
  };

  template<typename T_theta, typename T_amt, typename T_rate>
  struct PMXOdeFunctorSSAdaptorComponentII<T_theta, T_amt, T_rate, double> {
    const double ii_;
    const int npar_;
    
    PMXOdeFunctorSSAdaptorComponentII(const double ii, int npar)
      : ii_(ii), npar_(npar)
    {}

    template<typename T>
    inline const double operator()(const std::vector<T>& theta) const {
      return ii_;
    }
  };

  template <typename integrator_type, typename T_theta,
            typename T_amt, typename T_rate, typename T_ii, typename F>
  struct PMXOdeFunctorSSAdaptor {
    PMXOdeFunctorSSAdaptorComponentTheta<T_theta> theta_;
    PMXOdeFunctorSSAdaptorComponentAmt<T_theta, T_amt> amt_;
    PMXOdeFunctorSSAdaptorComponentRate<T_theta, T_amt, T_rate> rate_;
    PMXOdeFunctorSSAdaptorComponentII<T_theta, T_amt, T_rate, T_ii> ii_;
    const int cmt_;
    const int ncmt_;    
    const integrator_type& integrator_;
    
    PMXOdeFunctorSSAdaptor(const std::vector<T_theta>& theta, 
                              const T_amt& amt, const T_rate& rate,
                              const T_ii& ii, int cmt, int ncmt,
                              const integrator_type& integrator)
      : theta_(theta),
        amt_(amt, theta.size()),
        rate_(rate, theta.size()),
        ii_(ii, theta.size()),
        cmt_(cmt), ncmt_(ncmt),
        integrator_(integrator)
    {}

    inline Eigen::Matrix<typename stan::return_type_t<T_theta, T_amt, T_rate, T_ii>, -1, 1>
    adapted_param() const {
      using scalar_t = typename stan::return_type_t<T_theta, T_amt, T_rate, T_ii>;
      std::vector<scalar_t> vec;
      vec.reserve(theta_.theta_.size() + 3);
      if (stan::is_var<T_theta>::value) {
        vec.insert(vec.end(), theta_.theta_.begin(), theta_.theta_.end());
      }
      if (stan::is_var<T_amt>::value) {
        vec.push_back(amt_.amt_);
      }
      if (stan::is_var<T_rate>::value) {
        vec.push_back(rate_.rate_);
      }
      if (stan::is_var<T_ii>::value) {
        vec.push_back(ii_.ii_);
      }
      return stan::math::to_vector(vec);
    }

    /**
     * When rate is RV, it's attached to parameters vector.
     * IN this case parameter @c y consists of {theta, rate}
     * This is used AFTER newton iteration to calculate <code>Jxy</code>.
     */
    template <typename T0>
    inline Eigen::Matrix<stan::math::var, -1, 1>
    operator()(const Eigen::Matrix<T0, -1, 1>& x,
               const Eigen::Matrix<stan::math::var, -1, 1>& y,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      using stan::math::value_of;

      static const char* func("Steady State Event");
      using scalar_t = stan::math::var;
      auto y_vec(stan::math::to_array_1d(y));

      double t0 = 0;            // FIXME: ODE explicitly depdends on time

      std::vector<scalar_t> x0(x.data(), x.data() + x.size());
      Eigen::Matrix<scalar_t, -1, 1> result(x.size());

      if (rate_(y_vec) == 0) {  // bolus dose
        x0[cmt_ - 1] += amt_(y_vec);
        std::vector<scalar_t> pred = integrator_(F(), x0, t0, ii_(y_vec),
                                                 theta_(y_vec), x_r, x_i)[0];
        for (int i = 0; i < result.size(); ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else if (ii_(y_vec) > 0) {  // multiple truncated infusions
        auto dt = amt_(y_vec) / rate_(y_vec);
        // consider dosing internal > ii
        int n = int(std::floor(value_of(dt) / value_of(ii_(y_vec))) + 0.1);

        std::vector<T_rate> rate_vec(ncmt_, 0.0);

        auto dt1 = dt - n * ii_(y_vec);
        rate_vec[cmt_ - 1] = (n + 1) * rate_(y_vec);
        auto ode_par = theta_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, T_theta, T_rate> f1(ode_par, rate_vec);
        x0 = integrator_(f1, to_array_1d(x), t0, dt1, f1.adaptor.adapted_param(), x_r, x_i)[0];

        auto dt2 = (n + 1) * ii_(y_vec) - dt;
        rate_vec[cmt_ - 1] = n * rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, T_theta, T_rate> f2(ode_par, rate_vec);
        std::vector<scalar_t> pred = integrator_(f2, x0, t0, dt2, f2.adaptor.adapted_param(), x_r, x_i)[0];
        for (int i = 0; i < result.size(); ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else {  // constant infusion
        stan::math::check_less_or_equal(func, "AMT", amt_(y_vec), 0); 
        std::vector<T_rate> rate_vec(ncmt_, 0.0);
        rate_vec[cmt_ - 1] = rate_(y_vec);

        auto ode_par = theta_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, T_theta, T_rate> f(ode_par, rate_vec);
#ifdef TORSTEN_AS_FP
        stan::math::throw_domain_error(func, "algebra_solver_fp used for ", 1, "constant infusion");
#else
        result = to_vector(f(0, to_array_1d(x), f.adaptor.adapted_param(), x_r, x_i, 0));
#endif
      }

      return result;
    }

    /** 
     * All-data version of the functor, used in KINSOL newton solver
     * 
     * 
     * @return functor evaluation
     */
    inline Eigen::VectorXd
    operator()(const Eigen::VectorXd& x,
               const Eigen::VectorXd& y,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      using stan::math::value_of;

      static const char* func("Steady State Event");
      auto y_vec(stan::math::to_array_1d(y));
      std::vector<double> theta = value_of(theta_.theta_);

      double amt = value_of(amt_.amt_);
      double rate = value_of(rate_.rate_);
      double ii = value_of(ii_.ii_);

      double t0 = 0;            // FIXME: ODE explicitly depdends on time

      std::vector<double> x0 = stan::math::to_array_1d(x);
      const int n = x.size();
      Eigen::Matrix<double, -1, 1> result(n);

      if (rate == 0) {  // bolus dose
        x0[cmt_ - 1] += amt;
        std::vector<double> pred = integrator_(F(), x0, t0, ii, theta, x_r, x_i)[0];
        for (int i = 0; i < n; ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else if (ii > 0) {  // multiple truncated infusions
        double dt = amt/rate;
        // consider dosing internal > ii
        int n = int(std::floor(dt / ii) + 0.1);

        std::vector<double> rate_vec(ncmt_, 0.0);

        double dt1 = dt - n * ii;
        rate_vec[cmt_ - 1] = (n + 1) * rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, double, double> f1(theta, rate_vec);
        x0 = integrator_(f1, to_array_1d(x), t0, dt1, f1.adaptor.adapted_param(), x_r, x_i)[0];

        double dt2 = (n + 1) * ii_(y_vec) - dt;
        rate_vec[cmt_ - 1] = n * rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, double, double> f2(theta, rate_vec);
        std::vector<double> pred = integrator_(f2, x0, t0, dt2, f2.adaptor.adapted_param(), x_r, x_i)[0];
        for (int i = 0; i < result.size(); ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else {  // constant infusion
        stan::math::check_less_or_equal(func, "AMT", amt, 0.0);
        std::vector<double> rate_vec(ncmt_, 0.0);
        rate_vec[cmt_ - 1] = rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, double, double> f(theta, rate_vec);
#ifdef TORSTEN_AS_FP
        stan::math::throw_domain_error(func, "algebra_solver_fp used for ", 1, "constant infusion");
#else
        result = to_vector(f(0, to_array_1d(x), f.adaptor.adapted_param(), x_r, x_i, 0));
#endif
      }

      return result;
    }

    /** 
     * With indepedent variable being <code>var</code>, this version
     * of the functor is used in KINSOL newton solver's Jacobian calculation
     * 
     * 
     * @return functor evaluation
     */
    inline Eigen::Matrix<stan::math::var, -1, 1>
    operator()(const Eigen::Matrix<stan::math::var, -1, 1>& x,
               const Eigen::VectorXd& y,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      using stan::math::value_of;
      using stan::math::var;

      static const char* func("Steady State Event");
      auto y_vec(stan::math::to_array_1d(y));
      std::vector<double> theta = value_of(theta_.theta_);
      double amt = value_of(amt_.amt_);
      double rate = value_of(rate_.rate_);
      double ii = value_of(ii_.ii_);


      double t0 = 0;            // FIXME: ODE explicitly depdends on time

      std::vector<var> x0 = stan::math::to_array_1d(x);
      const int n = x.size();
      Eigen::Matrix<var, -1, 1> result(n);

      if (rate == 0) {  // bolus dose
        x0[cmt_ - 1] += amt;
        std::vector<var> pred = integrator_(F(), x0, t0, ii, theta, x_r, x_i)[0];
        for (int i = 0; i < n; ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else if (ii > 0) {  // multiple truncated infusions
        double dt = amt/rate;
        // consider dosing internal > ii
        int n = int(std::floor(dt / ii) + 0.1);

        std::vector<double> rate_vec(ncmt_, 0.0);

        double dt1 = dt - n * ii;
        rate_vec[cmt_ - 1] = (n + 1) * rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, double, double> f1(theta, rate_vec);
        x0 = integrator_(f1, to_array_1d(x), t0, dt1, f1.adaptor.adapted_param(), x_r, x_i)[0];

        double dt2 = (n + 1) * ii_(y_vec) - dt;
        rate_vec[cmt_ - 1] = n * rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, double, double> f2(theta, rate_vec);
        std::vector<var> pred = integrator_(f2, x0, t0, dt2, f2.adaptor.adapted_param(), x_r, x_i)[0];
        for (int i = 0; i < result.size(); ++i) {
#ifdef TORSTEN_AS_FP
          result(i) = pred[i];
#else
          result(i) = x(i) - pred[i];
#endif
        }
      } else {  // constant infusion
        stan::math::check_less_or_equal(func, "AMT", amt, 0); 
        std::vector<double> rate_vec(ncmt_, 0.0);
        rate_vec[cmt_ - 1] = rate_(y_vec);
        const PMXOdeFunctorRateAdaptor<F, double, double> f(theta, rate_vec);
#ifdef TORSTEN_AS_FP
        stan::math::throw_domain_error(func, "algebra_solver_fp used for ", 1, "constant infusion");
#else
        result = to_vector(f(0, to_array_1d(x), f.adaptor.adapted_param(), x_r, x_i, 0));
#endif
      }

      return result;
    }
  };

  /**
   * ODE-based PKPD models.
   *
   * @tparam T_par PK parameters type
   * @tparam F ODE functor
   */
  template<typename T_par, typename F>
  class PKODEModel {
    static const double dt_min; /**< min step to move ode solution */
    const std::vector<double> x_r_dummy; /**< dummy data to point to*/
    const std::vector<int> x_i_dummy; /**< dummy data to point to */
    const std::vector<T_par> &par_; /**< parameters */
    const std::vector<double>& x_r_; /**< real data */
    const std::vector<int>& x_i_; /**< integer data */
    const F &f_;                /**< ODE functor */
    const int ncmt_;            /**< dim of ODE system */
  public:
    using par_type    = T_par;
    using f_type      = F;

    /**
     * Constructor
     *
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par, int ncmt, const F& f) :
      x_r_dummy(), x_i_dummy(),
      par_(par), x_r_(x_r_dummy), x_i_(x_i_dummy),f_(f), ncmt_(ncmt)
    {}

    /**
     * Constructor
     *
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i,
               int ncmt, const F& f) :
      par_(par), x_r_(x_r), x_i_(x_i),f_(f), ncmt_(ncmt)
    {}

    /**
     * Constructor
     *
     * @param par model parameters
     * @param f ODE functor
     * @param ncmt the ODE size.
     */
    PKODEModel(const std::vector<T_par> &par,
               const std::vector<double>& x_r,
               int ncmt, const F& f) :
      x_r_dummy(), x_i_dummy(),
      par_(par), x_r_(x_r), x_i_(x_i_dummy),f_(f), ncmt_(ncmt)
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
      x_r_(m.x_r_), x_i_(m.x_i_),
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
      if (is_var<T2>::value) res += rate.size();
      if (is_var<T3>::value) res += par.size();
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
        for (size_t j = 0; j < rate.size(); ++j) {
          res[i] = rate[j];
          i++;
        }
      }
      if (is_var<T3>::value) {
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
    template<typename integrator_type>
    Eigen::Matrix<double, Eigen::Dynamic, 1>
    integrate(double t1,
              const std::vector<double> &rate,
              const Eigen::Matrix<double, 1, Eigen::Dynamic>& y0,
              const double& dt,
              const integrator_type& integrator) const {
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
        const std::vector<double> pars{value_of(par_)};
        PMXOdeFunctorRateAdaptor<F, double, double> f(pars, rate);
        std::vector<std::vector<double> > res_v =
          integrator(f, y, t0, ts, f.adaptor.adapted_param(), x_r_, x_i_);
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
    template<typename Tt0, typename Tt1, typename T, typename T1, typename integrator_type>
    void solve(Eigen::Matrix<T, -1, 1>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const integrator_type& integrator) const {
      const double t0_d = stan::math::value_of(t0);
      std::vector<Tt1> ts(time_step(t0, t1));
      PMXOdeFunctorRateAdaptor<F, T_par, T1> f_rate(par_, rate);

      if (t1 - t0 > dt_min) {
        auto y_vec = stan::math::to_array_1d(y);
        std::vector<std::vector<T> > res_v =
          integrator(f_rate, y_vec, t0_d, ts,
                     f_rate.adaptor.adapted_param(), x_r_, x_i_);
        y = stan::math::to_vector(res_v[0]);
      }
    }

    /** 
     * Torsten's integrators can return 
     * results in form of data directly,
     * thanks to @c pmx_cvodes/arkode/odeint_integrator's observer implmenetation.
     * 
     * @param yd output data
     * @param y initial condition
     * @param t0 starting time
     * @param t1 end tme
     * @param rate infusion rate
     * @param integrator ODE integrator
     */
    template<typename T0, typename T, typename T1, typename integrator_type>
    void solve_d(Eigen::VectorXd& yd,
                 const PKRec<T>& y,
                 const T0& t0, const T0& t1,
                 const std::vector<T1>& rate,
                 const integrator_type& integrator) const {
      using stan::math::var;
      using stan::math::value_of;
      using stan::math::to_var;

      const double t0_d = value_of(t0);
      std::vector<T0> ts(time_step(t0, t1));
      PMXOdeFunctorRateAdaptor<F, T_par, T1> f_rate(par_, rate);

      if (t1 - t0 > dt_min) {
        auto y1d = stan::math::to_array_1d(y);
        yd = integrator.solve_d(f_rate, y1d, t0_d, ts,
                                f_rate.adaptor.adapted_param(), x_r_, x_i_).col(0);
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
    template<typename integrator_type, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_par, T_ii>, Eigen::Dynamic, 1> // NOLINT
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const integrator_type& integrator) const {
      using stan::math::value_of;
      using stan::math::algebra_solver_powell;
      using stan::math::algebra_solver_newton;
      using stan::math::algebra_solver_fp;

      typedef typename stan::return_type_t<T_amt, T_r, T_par, T_ii> scalar;

      double ii_dbl = value_of(ii);
      Eigen::Matrix<double, 1, -1> init_dbl(Eigen::Matrix<double, 1, -1>::Zero(ncmt_));
      std::vector<double> rate_vec(ncmt_, 0);
      // std::vector<int> x_i{cmt, ncmt_, int(par_.size()), int(integrator.max_num_step)};

      if (rate == 0) {                     // bolus dose
        init_dbl(cmt - 1) = value_of(amt); // bolus as initial condition
      } else {                             // infusion
        rate_vec[cmt - 1] = value_of(rate);
      }

      const double init_dt = (rate == 0.0 || ii > 0) ? ii_dbl : 24.0;
      PMXOdeFunctorSSAdaptor<integrator_type, T_par, T_amt, T_r, T_ii, F>
        fss(par_, amt, rate, ii, cmt, ncmt_, integrator);
      try {
#ifdef TORSTEN_AS_POWELL
      return algebra_solver_powell(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                                   fss.adapted_param(),
                                   x_r_, x_i_, 0,
                                   integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step);
#elif defined(TORSTEN_AS_FP)
      std::vector<double> u_scale(ncmt_, 1.0);
      std::vector<double> f_scale(ncmt_, 1.0);
      return algebra_solver_fp(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                               fss.adapted_param(),
                               x_r_, x_i_,
                               u_scale, f_scale, 0,
                               integrator.as_atol, integrator.as_max_num_step);
#else
      return pmx_algebra_solver_newton(fss, integrate(t0, rate_vec, init_dbl, init_dt, integrator),
                                       fss.adapted_param(),
                                       x_r_, x_i_, 0,
                                       integrator.as_rtol, integrator.as_atol, integrator.as_max_num_step);
#endif
      } catch (const std::exception& e) {
        const char *text =
          "Torsten failed to find steady state, due to "
          "either system's intrinsic lack of such a state "
          "or improper algebra solver controls. Details: ";
        throw std::runtime_error(std::string(text) + e.what());
      }
    }
  };

  template<typename T_par, typename F>
  const double PKODEModel<T_par, F>::dt_min = 1.e-12;
}

#endif
