#ifndef STAN_MATH_TORSTEN_COUPLED_MODEL_HPP
#define STAN_MATH_TORSTEN_COUPLED_MODEL_HPP

#include <stan/math/torsten/pmx_ode_integrator.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>

namespace torsten {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using torsten::PMXOdeIntegrator;

  /**
   * Packer in charge of packing & unpacking data when
   * adatpor is used for integration and algebra solutioin.
   * @tparam T_rate type of RATE
   * @tparam T_amt type of AMT
   * @tparam T_ii type of dosing interval II
   */
  template<typename T_rate, typename T_amt, typename T_ii>
  struct PMXOdeFunctorCouplingAdaptorPacker;

  /**
   * Packer in charge of packing & unpacking data when
   * adatpor is used for integration and algebra solutioin
   * when AMT & II are data.
   * @tparam T_rate type of RATE
   */
  template<typename T_rate>
  struct PMXOdeFunctorCouplingAdaptorPacker<T_rate, double, double> {

    /*
     * when solving coupled model with @c var rate, we
     * append rate and PK initial condition to
     * parameter vector. 
     * FIXME: spurious @c var parameters will be generated
     * if the original parameters are data.
     */
    template<typename T, typename T0, int R, int C>
    static std::vector<typename stan::return_type_t<T, T0, T_rate>>
    adapted_param(const std::vector<T> &par, const std::vector<T_rate> &rate,
                  const Eigen::Matrix<T0, R, C>& y0_pk) {
      std::vector<stan::math::var> theta;
      theta.reserve(par.size() + rate.size() + y0_pk.size());
      theta.insert(theta.end(), par.begin(), par.end());
      theta.insert(theta.end(), rate.begin(), rate.end());
      for (int i = 0; i < y0_pk.size(); ++i) {
        theta.push_back(y0_pk(i)); 
      }
      return theta;
    }

    /*
     * when solving coupled model with @c var rate, we
     * append rate and PK initial condition to
     * parameter vector. 
     * FIXME: spurious @c var parameters will be generated
     * if the original parameters are data.
     */
    template<typename T, typename T0, int R, int C>
    static std::vector<typename stan::return_type_t<T, T0, T_rate>>
    adapted_param(const std::vector<T> &par, int cmt, int ncmt, const T_rate& rate,
                  const Eigen::Matrix<T0, R, C>& y0_pk) {
      std::vector<stan::math::var> theta(par.size() + ncmt + y0_pk.size(), 0.0);
      std::copy(par.begin(), par.end(), theta.begin());
      theta[par.size() + cmt - 1] = rate;
      const int nPK = y0_pk.size();
      for (int i = 0; i < nPK; ++i) {
        *(theta.rbegin() + i) = y0_pk(nPK - i - 1);
        // theta[par.size() + ]
        // theta.push_back(y0_pk(i)); 
      }
      return theta;
    }

    /*
     * when solving coupled model with @c var rate, @c x_r
     * is filled with initial time
     */
    template<typename T0, int R, int C>
    static std::vector<double>
    adapted_x_r(const std::vector<T_rate> &rate,
                const Eigen::Matrix<T0, R, C>& y0_pk, double t0) {
      return {t0};
    }

    /*
     * when solving coupled model with @c var rate, @c x_r
     * is filled with initial time
     */
    template<typename T0, int R, int C>
    static std::vector<double>
    adapted_x_r(int cmt, int ncmt, const T_rate& rate,
                const Eigen::Matrix<T0, R, C>& y0_pk, double t0) {
      return {t0};
    }
  };

  /**
   * Packer in charge of packing & unpacking data when
   * adatpor is used for integration and algebra solutioin
   * when RATE, AMT & II are data.
   * @tparam T_rate type of RATE
   */
  template<>
  struct PMXOdeFunctorCouplingAdaptorPacker<double, double, double> {
    /*
     * when solving coupled model with @c var rate, we
     * append rate and PK initial condition to
     * parameter vector. 
     * FIXME: spurious @c var parameters will be generated
     * if the original parameters are data.
     */
    template<typename T, typename T0, int R, int C>
    static std::vector<typename stan::return_type_t<T, T0>>
    adapted_param(const std::vector<T> &par, const std::vector<double> &rate,
                  const Eigen::Matrix<T0, R, C>& y0_pk) {
      std::vector<typename stan::return_type_t<T, T0>> theta;
      theta.reserve(par.size() + y0_pk.size());
      theta.insert(theta.end(), par.begin(), par.end());
      for (int i = 0; i < y0_pk.size(); ++i) {
        theta.push_back(y0_pk(i)); 
      }
      return theta;
    }

    /*
     * when solving coupled model with @c var rate, we
     * append rate and PK initial condition to
     * parameter vector. 
     * FIXME: spurious @c var parameters will be generated
     * if the original parameters are data.
     */
    template<typename T, typename T0, int R, int C>
    static std::vector<typename stan::return_type_t<T, T0>>
    adapted_param(const std::vector<T> &par, int cmt, int ncmt, double rate,
                  const Eigen::Matrix<T0, R, C>& y0_pk) {
      return adapted_param(par, std::vector<double>(), y0_pk);
    }

    /*
     * when solving coupled model with @c var rate, @c x_r
     * is filled with initial time
     */
    template<typename T0, int R, int C>
    static std::vector<double>
    adapted_x_r(const std::vector<double> &rate,
                const Eigen::Matrix<T0, R, C>& y0_pk, double t0) {
      std::vector<double> res(rate);
      res.push_back(t0);
      return res;
    }

    /*
     * when solving coupled model with @c var rate, @c x_r
     * is filled with initial time
     */
    template<typename T0, int R, int C>
    static std::vector<double>
    adapted_x_r(int cmt, int ncmt, double rate,
                const Eigen::Matrix<T0, R, C>& y0_pk, double t0) {
      std::vector<double> res(ncmt + 1, 0.0);
      res[cmt - 1] = rate;
      res.back() = t0;
      return res;
    }
  };

  /**
   * Packer in charge of packing & unpacking data when
   * adatpor is used for integration and algebra solutioin.
   * @tparam T_rate type of RATE
   * @tparam T_amt type of AMT
   * @tparam T_ii type of dosing interval II
   */
  template<typename T_rate, typename T_amt, typename T_ii>
  struct PMXOdeFunctorCouplingSSAdaptorPacker;

  /**
   * Packer in charge of packing & unpacking data when
   * adatpor is used for integration and algebra solutioin
   * when RATE, AMT & II are data.
   */
  template<>
  struct PMXOdeFunctorCouplingSSAdaptorPacker<double, double, double> {
    /*
     * when solving coupled model with @c var rate, we
     * append rate and PK initial condition to
     * parameter vector. 
     * FIXME: spurious @c var parameters will be generated
     * if the original parameters are data.
     */
    template<typename T, typename T0, int R, int C>
    static Eigen::Matrix<typename stan::return_type_t<T, T0>, -1, 1>
    adapted_param(const std::vector<T> &par, int cmt, int ncmt, double rate,
                  const Eigen::Matrix<T0, R, C>& y0_pk) {
      Eigen::Matrix<typename stan::return_type_t<T, T0>, -1, 1> theta(par.size() + y0_pk.size());
      for (size_t i = 0; i < par.size(); ++i) {
        theta(i) = par[i];
      }
      for (int i = 0; i < y0_pk.size(); ++i) {
        theta(i + par.size()) = y0_pk(i);
      }
      return theta;
    }

    /*
     * when solving coupled model with @c var rate, @c x_r
     * is filled with initial time
     */
    template<typename T0, int R, int C>
    static std::vector<double>
    adapted_x_r(int cmt, int ncmt, double rate, double amt, double ii,
                const Eigen::Matrix<T0, R, C>& y0_pk, double t0) {
      std::vector<double> res(PMXOdeFunctorCouplingAdaptorPacker<double, double, double>::
                              adapted_x_r(cmt, ncmt, rate, y0_pk, t0));
      res.push_back(amt);
      return res;
    }
  };

  /**
   * In a coupled model's ODE functor, we first solve the PK
   * model using analytical solution, then pass it to the
   * numerical integrator to solve the ODE model. This
   * requires adapting the coupled model's functor to an ODE
   * functor that can be passed to ODE integrators.
   * The fuctor returns the derivative of the base ODE system when
   * @c rate is param.
   * The base PK component is calculated analytically and passed into
   * ODE functor as arguments.
   *  
   * @param t time
   * @param y indepedent PD variable
   * @param theta adapted coupled ODE params in the order of
   *    - original ODE param (PK param followed by PD param)
   *    - PK rate
   *    - PD rate
   *    - initial condition of PK ODE
   * @param x_r real data in the order of
   *    - dim PK states
   *    - initial time of PK ODE
   * @param x_i null int data
   */
  template <template<typename...> class T_m, typename F0, typename T_rate>
  struct PMXOdeFunctorCouplingAdaptor {
    PMXOdeFunctorCouplingAdaptor() {}

    /**
     *  default case: rate is a parameter, stored in theta.
     *  Theta contains in this order: ODE parameters, rates, and
     *  initial base PK states.
     */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename stan::return_type_t<T0, T1, T2, T3>>
    operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
      using std::vector;
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      typedef typename stan::return_type_t<T0, T1, T2, T3>
        scalar;
      typedef typename stan::return_type_t<T0, T2, T3>
        T_pk;  // return object of fTwoCpt  doesn't depend on T1
      using T_pkmodel = T_m<double>;

      size_t nPK = T_pkmodel::Ncmt;              // number of base PK states (for rates and inits)
      size_t nPD = y.size();                     // number of other states
      size_t nODEparms = theta.size() - 2 * nPK - nPD; // number of ODE parameters

      vector<T2> thetaPK(theta.begin(), theta.begin() + T_pkmodel::Npar);
      vector<T2> ratePK(theta.begin() + nODEparms, theta.begin() + nODEparms + nPK);

      // The last elements of theta contain the initial base PK states
      torsten::PKRec<T2> init_pk(nPK);
      for (int i = 0; i < nPK; ++i) {
        init_pk(nPK - i - 1) = *(theta.rbegin() + i);
      }

      // Last element of x_r contains the initial time
      T_m<T2> pkmodel(thetaPK);
      pkmodel.solve(init_pk, x_r.back(), t, ratePK);
      vector<scalar> dydt = F0()(t, y,
                                 to_array_1d(init_pk),
                                 theta, x_r, x_i, pstream_);

      for (size_t i = 0; i < dydt.size(); i++)
        dydt[i] += theta[nODEparms + nPK + i];

      return dydt;
    }
  };

  template <template<typename...> class T_m, typename F0>
  struct PMXOdeFunctorCouplingAdaptor<T_m, F0, double> {
    PMXOdeFunctorCouplingAdaptor() {}

    /**
     * Returns the derivative of the base ODE system when @c rate is data.
     * The base PK component is calculated analytically and passed into
     * ODE functor as arguments.
     *  
     * @param t time
     * @param y indepedent PD  variable
     * @param theta adapted coupled ODE params in the order of
     *    - original ODE param
     *    - initial condition of PK ODE
     * @param x_r real data in the order of
     *    - PK rate
     *    - PD rate
     *    - initial time of PK ODE
     * @param x_i null int data
     */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename stan::return_type_t<T0, T1, T2, T3>>
    operator()(const T0& t,
             const std::vector<T1>& y,
             const std::vector<T2>& theta,
             const std::vector<T3>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream_) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      using scalar = typename stan::return_type_t<T0, T1, T2, T3>;
      using T_pk = typename stan::return_type_t<T0, T2, T3>;
      using T_pkmodel = T_m<double>;

      // Get initial PK states stored at the end of @c theta
      int nPK = T_pkmodel::Ncmt;
      torsten::PKRec<T2> init_pk(nPK);
      for (int i = 0; i < nPK; ++i) {
        init_pk(nPK - i - 1) = *(theta.rbegin() + i);
      }

      // The last element of x_r contains the initial time,
      // and the beginning of theta are for PK params.
      T_m<T2> pkmodel(theta);
      std::vector<double> pk_rate(x_r.begin(), x_r.begin() + nPK);
      pkmodel.solve(init_pk, x_r.back(), t, pk_rate);
      // move PK RHS to current time then feed the solution to PD ODE
      std::vector<T_pk> y_pk = to_array_1d(init_pk);
      std::vector<scalar> dydt = F0()(t, y, y_pk, theta, x_r, x_i, pstream_);

      // x_r: {pk rate, pd rate, t0}
      for (size_t i = 0; i < dydt.size(); ++i) {
        dydt[i] += x_r[nPK + i];
      }

      return dydt;
    }
  };

  /**
   * A structure to store the algebraic system
   * which gets solved when computing the steady
   * state solution for coupled ODE models.
   *

   * @tparam T_amt @c amt type
   * @tparam T_rate @c rate type
   * @tparam F ODE RHS functor type
   * @tparam F2 type of the ODE that has analytical solution
   * @tparam T_integrator integrator type
   */
  template <template<typename...> class T_m, typename F, typename T_amt, typename T_rate, typename T_ii, typename T_integrator>
  struct PMXOdeFunctorCouplingSSAdaptor;

  /**
   * A structure to store the algebraic system
   * which gets solved when computing the steady
   * state solution.
   * 
   * specialization: both amt and rate are data
   */
  template <template<typename...> class T_m, typename F, typename T_integrator>
  struct PMXOdeFunctorCouplingSSAdaptor<T_m, F, double, double, double, T_integrator> {
    F f_;
    double ii_;
    int cmt_;  // dosing compartment
    const T_integrator integrator_;
    int nPK_;

    PMXOdeFunctorCouplingSSAdaptor() {}

    PMXOdeFunctorCouplingSSAdaptor(const F& f, double ii, int cmt,
                                   const T_integrator& integrator, int nPK)
      : f_(f), ii_(ii), cmt_(cmt), integrator_(integrator),
        nPK_(nPK) {}

    /**
     *  dd regime.
     *  dat contains the rates in each compartment followed
     *  by the adjusted amount (biovar * amt).
     */
    template <typename T0, typename T1>
    inline
    Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                  Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& dat,
               const std::vector<int>& x_i,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;

      typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
      typedef typename stan::return_type<T0, T1>::type T_deriv;

      double t0 = 0;

      std::vector<scalar> x0(x.size());
      for (size_t i = 0; i < x0.size(); i++) x0[i] = x(i);
      double amt = dat[dat.size() - 1];
      double rate = dat[cmt_ - 1];

      // real data, which gets passed to the integrator, shoud not have
      // amt in it, as the transient coupled solver needs to
      // have the last element to be the absolute time (in this case, 0).
      std::vector<double> x_r = dat;
      x_r.pop_back();

      Eigen::Matrix<scalar, -1, 1> result(x.size());
      std::vector<T1> y_vec(to_array_1d(y));
      std::vector<T0> x_vec(to_array_1d(x));

      if (rate == 0) {  // bolus dose
        if ((cmt_ - nPK_) >= 0) x0[cmt_ - nPK_ - 1] += amt;

        std::vector<scalar> pred = integrator_(f_, x0, t0, ii_, y_vec, x_r, x_i)[0];
        for (int i = 0; i < pred.size(); ++i) {
          pred[i] = x(i) - pred[i];          
        }
        result = stan::math::to_vector(pred);
      } else if (ii_ > 0) {  // multiple truncated infusions
        double delta = amt / rate;

        static const char* function("Steady State Event");
        torsten::check_mti(amt, delta, ii_, function);

        // PD solution of infusion, note that @c f_ is coupled adaptor functor
        x0 = integrator_(f_, x_vec, t0, delta, y_vec, x_r, x_i)[0];

        PKRec<T1> x_pk(nPK_);
        int nParms = y.size() - nPK_;
        for (int i = 0; i < nPK_; i++) x_pk(i) = y(nParms + i);

        std::vector<double> pk_rate(x_r.begin(), x_r.begin() + nPK_);
        T_m<T1>(y_vec).solve(x_pk, t0, delta, pk_rate);

        // no more infusion after amt/rate
        x_r[cmt_ - 1] = 0;

        std::vector<T1> y2(y.data(), y.data() + y.size());
        // Matrix<T1, -1, 1> y2(y);
        for (int i = 0; i < nPK_; ++i) y2[nParms + i] = x_pk(i);
        std::vector<scalar> pred = integrator_(f_, x0, t0, ii_ - delta, y2, x_r, x_i)[0];

        for (int i = 0; i < result.size(); i++)
          result(i) = x(i) - pred[i];
      } else {  // constant infusion
        std::vector<T_deriv> derivative = f_(0, x_vec, y_vec, x_r, x_i, 0);
        result = stan::math::to_vector(derivative);
      }

      return result;
    }
  };

  /**
   * A structure to store the algebraic system
   * which gets solved when computing the steady
   * state solution.
   * 
   * specialization: rate is data. When biovar is a parameter,
  *  amt will be a transformed parameter.
  *  The last element of y is contains amt.
  *  dat stores the rate.
   */
  template <template<typename...> class T_m, typename F, typename T_amt, typename T_integrator>
  struct PMXOdeFunctorCouplingSSAdaptor<T_m, F, T_amt, double, double, T_integrator> {
    F f_;
    double ii_;
    int cmt_;  // dosing compartment
    const T_integrator integrator_;
    int nPK_;

    PMXOdeFunctorCouplingSSAdaptor() {}

    PMXOdeFunctorCouplingSSAdaptor(const F& f, double ii, int cmt,
                                   const T_integrator& integrator, int nPK)
      : f_(f), ii_(ii), cmt_(cmt), integrator_(integrator),
        nPK_(nPK) {}

    template <typename T0, typename T1>
    inline
    Eigen::Matrix<typename boost::math::tools::promote_args<T0, T1>::type,
                  Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& dat,
               const std::vector<int>& dat_int,
               std::ostream* msgs) const {
      using stan::math::to_array_1d;
      using stan::math::to_vector;
      using std::vector;
      using stan::math::to_vector;
      using stan::math::invalid_argument;

      typedef typename boost::math::tools::promote_args<T0, T1>::type scalar;
      typedef typename stan::return_type<T0, T1>::type T_deriv;

      double t0 = 0;
      std::vector<double> ts(1);
      std::vector<double> rate_v(dat.size(), 0);
      for (size_t i = 0; i < rate_v.size(); i++) rate_v[i] = dat[i];

    std:vector<scalar> x0(x.size());
      for (size_t i = 0; i < x0.size(); i++) x0[i] = x(i);
      scalar amt = y(y.size() - 1);
      double rate = ((cmt_ - 1) >= 0) ? dat[cmt_ - 1] : 0;

      Eigen::Matrix<scalar, Eigen::Dynamic, 1> result(x.size());
      std::vector<scalar> parms(y.size() - 1);
      for (size_t i = 0; i < parms.size(); i++) parms[i] = y(i);

      if (rate == 0) {  // bolus dose
        if ((cmt_ - nPK_) >= 0) x0[cmt_ - nPK_ - 1] += amt;
        ts[0] = ii_;
        std::vector<scalar> pred = integrator_(f_, x0, t0, ts, parms,
                                          dat, dat_int)[0];

        for (int i = 0; i < result.size(); i++)
          result(i) = x(i) - pred[i];

      } else if (ii_ > 0) {  // multiple truncated infusions
        stan::math::invalid_argument("Steady State Event",
                         "Current version does not handle the case of",
                         "", " multiple truncated infusions ",
                         "(i.e ii > 0 and rate > 0) when F * amt is a parameter.");  // NOLINT

      } else {  // constant infusion
        if (amt != 0)
          stan::math::invalid_argument("Steady State Event",
                           "amt should be 0 when specifying a constant",
                           "", " infusion (i.e. rate > 0 and ii = 0).",
                           "");

        std::vector<T_deriv> derivative = f_(0, to_array_1d(x), parms,
                                        dat, dat_int, 0);
        result = to_vector(derivative);
      }

      return result;
    }
  };

  /**
   * Coupled model.
   *
   * @tparam T_m1 type of the 1st model in the coupling.
   * @tparam T_m2 type of the 2nd model in the coupling.
   */
  template <typename T_m1, typename T_m2>
  class PKCoupledModel;

  /**
   * Specialization of coupled model with 2nd model
   * being @c PKODEModel.
   *
   * @tparam T_m type of 1st model, choose among
   *             @c PMXOneCptmodel, @c PMXTwoCptmodel, @c PMXLinODEModel.
   * @tparam T_time type of time
   * @tparam T_rate type of dosing rate.
   * @tparam T_par type of parameter.
   * @tparam F type of ODE functor for @c PKODEModel.
   */
  template <template<typename...> class T_m,
            typename T_rate, typename T_par, typename F> // NOLINT
  class PKCoupledModel<T_m<T_par>,
                       torsten::PKODEModel<T_par,
                                           PMXOdeFunctorCouplingAdaptor<T_m, F, T_rate>> > { // NOLINT
    PMXOdeFunctorCouplingAdaptor<T_m, F, T_rate> f;
    int n_ode;

  public:
    using Fa = PMXOdeFunctorCouplingAdaptor<T_m, F, T_rate>;
    const T_m<T_par> pk_model;
    const torsten::PKODEModel<T_par, Fa> ode_model;

    using par_type    = T_par;
    using rate_type   = T_rate;

    /**
     * Coupled model constructor
     *
     * @param t0 initial time
     * @param y0 initial condition, with PK model's initial
     *           condition followed by ODE model's initial condition.
     * @param rate dosing rate
     * @param par model parameters
     * @param f ODE functor
     * @param n_ode the size of ode_model's ODE system
     */
    PKCoupledModel(const std::vector<T_par> & par,
                   const F& f0,
                   const int n_ode_) :
      f(),
      pk_model(par),
      ode_model(par, n_ode_, f),
      n_ode(n_ode_)
    {}
    
  private:

  /**
   * Steady state solver to
   * calculate the amount in each compartment at the end of a
   * steady-state dosing interval or during a steady-state
   * constant input (if ii = 0). The function is overloaded
   * to address the cases where amt or rate may be fixed or
   * random variables (yielding a total of 4 cases).
   * 
   * Case 1 (dd): amt and rate are fixed.
   *
   *	 @tparam T_time type of scalar for time
   *	 @tparam T_ii type of scalar for interdose interval
   *	 @tparam T_parameters type of scalar for ODE parameters
   *   @tparam T_biovar type of scalar for bio-availability
   *	 @tparam F type of ODE system function
   *	 @param[in] parameter model parameters at current event
   *	 @param[in] rate
   *	 @param[in] ii interdose interval
   *	 @param[in] cmt compartment in which the event occurs
   *	 @param[in] f functor for base ordinary differential equation
   *              that defines compartment model
   *   @return an eigen vector that contains predicted amount in each
   *           compartment at the current event.
   */
  template<typename T_ii, typename T_integrator>
  torsten::PKRec<typename stan::return_type_t<T_ii, T_par>>
  integrate(const double& amt,
            const double& rate,
            const T_ii& ii,
            const int& cmt,
            const T_integrator& integrator) const {
    typedef typename stan::return_type_t<T_ii, T_par> scalar;

    const double ii_dbl = stan::math::value_of(ii);
    const int nPK = pk_model.ncmt();
    const int nPD = n_ode;
    // const double t0 = stan::math::value_of(ode_model.t0());
    const double t0 = 0.0;      // FIXME: rm dummy

    Eigen::Matrix<T_par, -1, 1> predPK;
    if (cmt <= nPK) {  // check dosing occurs in a base state
      T_m<T_par> pkmodel(pk_model.par());
      predPK = pkmodel.solve(t0, amt, rate, ii_dbl, cmt);
    } else {
      predPK = Eigen::Matrix<scalar, -1, 1>::Zero(nPK);
    }

    using F_c = PMXOdeFunctorCouplingAdaptor<T_m, F, double>;

    // Tuning parameters for algebraic solver
    double rel_tol = 1e-10;  // default
    double f_tol = 1e-4;  // empirical
    long int max_num_steps = 1e3;  // default // NOLINT

    // construct algebraic system functor: note we adjust cmt
    // such that 1 corresponds to the first state we compute
    // numerically.
    PMXOdeFunctorCouplingSSAdaptor<T_m, F_c, double, double, T_ii, T_integrator>
      system(F_c(), ii_dbl, cmt, integrator, nPK);

    Eigen::Matrix<double, -1, 1> predPD_guess;
    Eigen::Matrix<scalar, 1, -1> predPD;

    std::vector<double> init_pd(nPD, 0.0);
    if (rate == 0) {  // bolus dose
      if (cmt > nPK) {
        init_pd[cmt - 1] = amt; 
      } else {
        predPK(cmt - 1) += amt;        
      }
    }

    using packer_t = PMXOdeFunctorCouplingAdaptorPacker<double, double, double>;
    using ss_packer_t = PMXOdeFunctorCouplingSSAdaptorPacker<double, double, double>;
    std::vector<T_par> pkpar(packer_t::adapted_param(ode_model.par(), cmt, nPK + nPD, rate, predPK));
    std::vector<double> x_r(packer_t::adapted_x_r(cmt, nPK + nPD, rate, predPK, t0));
    std::vector<int> x_i;

    // for time step of initail guess, in const infusion we
    // dose for 24 hours, otherwise use dosing interval.
    double init_dt = (rate > 0.0 && ii == 0.0) ? 24.0 : ii_dbl;
    predPD_guess = stan::math::to_vector(integrator(F_c(), init_pd, 0.0, init_dt,
                                                    packer_t::adapted_param(stan::math::value_of(ode_model.par()),
                                                                            cmt, nPK + nPD, rate,
                                                                            stan::math::value_of(predPK)),
                                                    packer_t::adapted_x_r(cmt, nPK + nPD, rate, predPK, t0),
                                                    x_i)[0]);
    predPD = stan::math::algebra_solver_powell(system, predPD_guess,
                                        ss_packer_t::adapted_param(ode_model.par(), cmt, nPK + nPD, rate, predPK),
                                        ss_packer_t::adapted_x_r(cmt, nPK + nPD, rate, amt, ii, predPK, t0),
                                        x_i,
                                        0, rel_tol, f_tol, max_num_steps);

    if (rate == 0.0) {
      if (cmt <= nPK) predPK(cmt - 1) -= amt;
    }

    torsten::PKRec<scalar> pred(nPK + nPD);
    for (int i = 0; i < nPK; i++) pred(i) = predPK(i);
    for (int i = 0; i < nPD; i++) pred(nPK + i) = predPD(i);

    return pred;
  }

  /**
   * Case 2 (vd): amt is random, rate is fixed.
   */
  template<typename T_ii, typename T_amt, typename T_integrator>
  torsten::PKRec<typename stan::return_type<T_ii, T_amt, T_par>::type>
  integrate(const T_amt& amt,
            const double& rate,
            const T_ii& ii,
            const int& cmt,
            const T_integrator& integrator) const {
    typedef typename stan::return_type_t<T_ii, T_amt, T_par> scalar;

    const double ii_dbl = stan::math::value_of(ii);
    const int nPK = pk_model.ncmt();
    const int nPD = n_ode;

    // Compute solution for base 1cpt PK
    Eigen::Matrix<scalar, -1, 1> predPK;
    std::vector<T_par> pkpar = ode_model.par();
    if (cmt <= nPK) {  // check dosing occurs in a base state
      const double t0 = 0.0;
      T_m<T_par> pkmodel(pk_model.par());
      predPK = pkmodel.solve(t0, amt, rate, ii_dbl, cmt);
      predPK(cmt - 1) += amt;
    } else {
      predPK = Eigen::Matrix<scalar, -1, 1>::Zero(nPK);
    }

    std::vector<scalar> theta2;
    theta2.insert(theta2.end(), pkpar.begin(), pkpar.end());
    theta2.insert(theta2.end(), predPK.data(), predPK.data() + predPK.size());
    theta2.push_back(amt);

    // Arguments for ODE integrator (and initial guess)
    std::vector<double> init_pd(nPD, 0.0);
    std::vector<double> x_r(nPK + nPD, 0);  // rate for the full system
    x_r.push_back(0);  // include initial time (at SS, t0 = 0)
    std::vector<int> x_i;

    // Tuning parameters for algebraic solver
    double rel_tol = 1e-10;  // default
    double f_tol = 1e-4;  // empirical
    long int max_num_steps = 1e3;  // default // NOLINT

    using F_c = PMXOdeFunctorCouplingAdaptor<T_m, F, double>;
    F_c f_coupled;

    PMXOdeFunctorCouplingSSAdaptor<T_m, F_c, T_amt, double, T_ii, T_integrator>
      system(F_c(), ii_dbl, cmt, integrator, nPK);

    Eigen::Matrix<double, -1, 1> predPD_guess;
    Eigen::Matrix<scalar, 1, -1> predPD;

    double init_dt = (rate > 0.0 && ii == 0.0) ? 24.0 : ii_dbl;
    if (rate == 0) {  // bolus dose
      predPD_guess = stan::math::to_vector(integrator(F_c(), init_pd,
                                        0.0, init_dt,
                                                      stan::math::value_of(theta2),
                                        x_r, x_i)[0]);

      predPD = algebra_solver_powell(system, predPD_guess,
                              to_vector(theta2),
                              x_r, x_i,
                              0, rel_tol, f_tol, max_num_steps);

      if (cmt <= nPK) predPK(cmt - 1) -= amt;
    } else if (ii > 0) {
      stan::math::invalid_argument("Steady State Event",
                       "Current version does not handle the case of",
                       "", " multiple truncated infusions ",
                       "(i.e ii > 0 and rate > 0) when F * amt is a parameter.");  // NOLINT
    } else {
      x_r[cmt - 1] = rate;

      predPD_guess = stan::math::to_vector(integrator(F_c(), init_pd,
                                                      0.0, init_dt,
                                                      stan::math::value_of(theta2),
                                                      x_r, x_i)[0]);

      predPD = algebra_solver_powell(system, predPD_guess,
                              stan::math::to_vector(theta2),
                              x_r, x_i,
                              0, rel_tol, f_tol, max_num_steps);
    }

    Eigen::Matrix<scalar, -1, 1> pred(nPK + nPD);
    for (int i = 0; i < nPK; i++) pred(i) = predPK(i);
    for (int i = 0; i < nPD; i++) pred(nPK + i) = predPD(i);

    return pred;
  }

  public:
    /*
     * solve the coupled model.
     */
    // template<torsten::PMXOdeIntegratorId It>
    // torsten::PKRec<scalar_type>
    // solve(const T_time& t_next,
    //       const torsten::PMXOdeIntegrator<It>& integrator) const {
    template<typename Tt0, typename Tt1, typename T, typename T1, PMXOdeIntegratorId It>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const PMXOdeIntegrator<It>& integrator) const {
      using std::vector;

      // pass fixed times to the integrator. FIX ME - see issue #30
      Tt1 t = t1;
      vector<double> t_dbl{stan::math::value_of(t)};
      double t0_dbl = stan::math::value_of(t0);

      PKRec<T> pred;
      if (t_dbl[0] > t0_dbl) {
        size_t nPK = pk_model.ncmt();
        // 
        // 
        PKRec<T> xPK(y.head(y.size() - n_ode));
        std::vector<T1> pk_rate(rate.begin(), rate.begin() + nPK);
        pk_model.solve(xPK, t0, t, pk_rate);
        PKRec<T> y_pk(y.head(y.size() - n_ode));
        PKRec<T> y_ode(y.segment(y_pk.size(), n_ode));

        // create vector with PD initial states
        vector<T> y0_PD(stan::math::to_array_1d(y_ode));
        PMXOdeFunctorCouplingAdaptor<T_m, F, T1> f_coupled;
        using packer_t = PMXOdeFunctorCouplingAdaptorPacker<T1, double, double>;
        vector<vector<T> >
          pred_V = integrator(f_coupled, y0_PD, t0_dbl, t_dbl,
                              packer_t::adapted_param(ode_model.par(), rate, y_pk),
                              packer_t::adapted_x_r(rate, y_pk, t0_dbl),
                              vector<int>());

        size_t nOde = pred_V[0].size();
        pred.resize(nPK + nOde);
        for (size_t i = 0; i < nPK; i++) pred(i) = xPK(i);
        for (size_t i = 0; i < nOde; i++) pred(nPK + i) = pred_V[0][i];
        y = pred;
      }
    }

    /* 
     * solve the coupled model, steady state. We delegate
     * the solution to @c integrate, in which the type of @c
     * amt will be used for template partial specification.
     */
    template<torsten::PMXOdeIntegratorId It, typename T_ii, typename T_amt>
    torsten::PKRec<typename stan::return_type_t<T_amt, T_ii, T_par>>
    solve(double t0, const T_amt& amt, const double& rate, const T_ii& ii, const int& cmt,
          const torsten::PMXOdeIntegrator<It>& integrator) const {
      return integrate(amt, rate, ii, cmt, integrator);
    }
  };

  template<typename T_rate, typename T_par, typename F> // NOLINT
  using PkOneCptOdeModel =
    PKCoupledModel< PMXOneCptModel<T_par>,
                    torsten::PKODEModel<T_par,
                               PMXOdeFunctorCouplingAdaptor<PMXOneCptModel, F, T_rate>> >; // NOLINT

  template<typename T_rate, typename T_par, typename F> // NOLINT
  using PkTwoCptOdeModel =
    PKCoupledModel< PMXTwoCptModel<T_par>,
                    torsten::PKODEModel<T_par,
                               PMXOdeFunctorCouplingAdaptor<PMXTwoCptModel, F, T_rate>> >; // NOLINT
}


#endif
