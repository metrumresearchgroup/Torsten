#ifndef STAN_MATH_TORSTEN_COUPLED_MODEL_HPP
#define STAN_MATH_TORSTEN_COUPLED_MODEL_HPP

#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/torsten/pmx_check.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>

namespace torsten {
  using Eigen::Matrix;
  using Eigen::Dynamic;

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
  template <template<typename...> class T_m, typename F>
  struct PMXOdeFunctorCouplingAdaptor {
    F const& f_;
    
    PMXOdeFunctorCouplingAdaptor(F const& f) : f_(f) {}

    /**
     */
    template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
    Eigen::Matrix<typename stan::return_type_t<T0, T1, T2, T3, T4, T5>, -1, 1>
    operator()(const T0& t,
               const Eigen::Matrix<T1, -1, 1>& y,
               std::ostream* pstream_,
               const std::vector<T2>& theta,
               const std::vector<T3>& rate,
               const T5& init_t,
               const Eigen::Matrix<T4, -1, 1>& init_pk,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      typedef typename stan::return_type_t<T0, T1, T2, T3, T4> scalar;
      using T_pkmodel = T_m<T2>;

      size_t nPK = T_pkmodel::Ncmt;
      size_t nPD = y.size();

      // theta & rate constains PK params followed by PD params
      std::vector<T2> thetaPK(theta.begin(), theta.begin() + T_pkmodel::Npar);
      std::vector<T3> ratePK(rate.begin(), rate.begin() + nPK);

      T_pkmodel pk(thetaPK);
      torsten::PKRec<typename stan::return_type_t<T0,T2,T3,T4,T5>> y_pk(nPK);
      for (auto i = 0; i < nPK; ++i) {
        y_pk[i] = init_pk[i];
      }
      pk.solve(y_pk, init_t, t, ratePK);

      Eigen::Matrix<scalar, -1, 1> dydt = f_(t, y, y_pk, theta, x_r, x_i, pstream_);
      for (auto i = 0; i < nPD; ++i) {
        dydt[i] += rate[nPK + i];
      }

      return dydt;
    }
  };

  /**
   * A structure to store the algebraic system
   * which gets solved when computing the steady
   * state solution for coupled ODE models.
   *
   */
  template <template<typename...> class T_m, typename F, typename T_integrator>
  struct PMXOdeFunctorCouplingSSAdaptor {
    F f_;
    const T_integrator integrator_;

    PMXOdeFunctorCouplingSSAdaptor(F const& f, T_integrator const& integrator)
    : f_(f), integrator_(integrator)
    {}
    
    /**
     *  dd regime.
     *  dat contains the rates in each compartment followed
     *  by the adjusted amount (biovar * amt).
     */
    template <typename T1, typename T2, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T1, T2, T_amt, T_r, T_ii>, -1, 1>
    operator()(const Eigen::Matrix<T1, -1, 1>& y,
               std::ostream* msgs,
               const std::vector<T2>& theta,
               const T_amt& amt, const T_r& rate_ss, const T_ii& ii, const int& cmt,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      using scalar = typename stan::return_type_t<T1, T2, T_amt, T_r, T_ii>;

      double init_t = 0;

      T_m<T2> pk_model(theta);
      auto init_pk = pk_model.solve(init_t, amt, rate_ss, ii, cmt);

      size_t nPK = T_m<double>::Ncmt;
      size_t nPD = y.size();
      Eigen::Matrix<scalar, -1, 1> residue(nPD);
      PMXOdeFunctorCouplingAdaptor<T_m, F> f1(f_);
      std::vector<T_r> rate(nPK + nPD, 0.0);

      if (rate_ss == 0) {  // bolus dose
        Eigen::Matrix<typename stan::return_type_t<T1, T_amt>, -1, 1> x0(nPD);
        for (size_t i = 0; i < nPD; i++) {
          x0[i] = y(i);
        }
        if (cmt > nPK) {
          x0[cmt - nPK - 1] += amt;
        } else {
          init_pk[cmt - 1] += amt;
        }
        Eigen::Matrix<scalar, -1, 1>  pred = integrator_(f1, x0, init_t, ii, theta, rate, init_t, init_pk, x_r, x_i);
        residue = y - pred;
      } else if (ii > 0) {  // multiple truncated infusions
        auto dt = amt / rate_ss;

        // FIXME: allow delta > ii
        torsten::check_mti(amt, dt, ii, "Steady State Event");

        // move PD forward, using PK SS solution as parameter
        rate[cmt - 1] = rate_ss;
        Eigen::Matrix<scalar, -1, 1> pred = integrator_(f1, y, init_t, dt, theta, rate, init_t, init_pk, x_r, x_i);

        // move PK forward, using PK SS solution as init condition
        if (cmt <= nPK) {
          std::vector<T_r> ratePK(rate.begin(), rate.begin() + nPK);
          pk_model.solve(init_pk, init_t, dt, ratePK);
        }

        // move PD forward with no more infusion
        rate[cmt - 1] = 0;
        pred = integrator_(f1, pred, init_t, ii - dt, theta, rate, init_t, init_pk, x_r, x_i);
        residue = y - pred;
      } else {  // constant infusion
        // rate[cmt - 1] = rate_ss;
        residue = f1(init_t, y, nullptr, theta, rate, init_t, init_pk, x_r, x_i);
      }

      return residue;
    }
  };

  /**
   * Coupled model.
   *
   * @tparam T_m1 type of the 1st model in the coupling.
   * @tparam T_m2 type of the 2nd model in the coupling.
   */
  template <typename T_m, typename F>
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
  template <template<typename...> class T_m, typename T_par, typename F>
  class PKCoupledModel<T_m<T_par>, F> {
    const std::vector<double> x_r_dummy; /**< dummy data to point to*/
    const std::vector<int> x_i_dummy; /**< dummy data to point to */
    const std::vector<T_par> &par_; /**< parameters */
    const std::vector<double>& x_r_; /**< real data */
    const std::vector<int>& x_i_; /**< integer data */
    F const& f_;
    int n_ode_;

  public:
    using par_type    = T_par;

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
    PKCoupledModel(const std::vector<T_par> & par, int n_ode, const F& f0) :
      par_(par),
      x_r_(x_r_dummy),
      x_i_(x_i_dummy),
      f_(f0),
      n_ode_(n_ode)
    {}
    
    PKCoupledModel(const std::vector<T_par> & par,
                   const std::vector<double>& x_r,
                   const std::vector<int>& x_i,
                   int n_ode, const F& f0) :
      par_(par),
      x_r_(x_r),
      x_i_(x_i),
      f_(f0),
      n_ode_(n_ode)
    {}

  public:
    /**
     * solve the transient coupled model.
     */
    template<typename Tt0, typename Tt1, typename T, typename T1, typename integrator_type>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const integrator_type& integrator) const {
      const double t0_d = stan::math::value_of(t0);
      // std::vector<Tt1> ts(time_step(t0, t1));
      // PMXOdeFunctorRateAdaptor<F> f_rate;
      const int nPK = T_m<T_par>::Ncmt;
      Eigen::Matrix<T, -1, 1> y_pk(y.head(nPK));
      Eigen::Matrix<T, -1, 1> y_pd(y.tail(n_ode_));
      const double dt_min = 1.e-12;
      if (t1 - t0 > dt_min) {
        PMXOdeFunctorCouplingAdaptor<T_m, F> f(f_);
        y_pd = integrator(f, y_pd, t0_d, t1, par_, rate, t0, y_pk, x_r_, x_i_);
        T_m<T_par> pk_model(par_);
        std::vector<T1> ratePK(rate.begin(), rate.begin() + nPK);
        pk_model.solve(y_pk, t0, t1, ratePK);
        y.head(T_m<T_par>::Ncmt) = y_pk;
        y.tail(n_ode_) = y_pd;
      }
    }

    /**
     * solve the coupled model, steady state. We delegate
     * the solution to @c integrate, in which the type of @c
     * amt will be used for template partial specification.
     */
    template<typename integrator_type, typename T_ii, typename T_r, typename T_amt>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_ii, T_par, T_r>, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const integrator_type& integrator) const {
      using stan::math::algebra_solver_newton;

      using scalar = typename stan::return_type_t<T_amt, T_r, T_par, T_ii>;

      size_t nPK = T_m<double>::Ncmt;
      size_t nPD = n_ode_;
      double ii_dbl = stan::math::value_of(ii);
      Eigen::Matrix<double, 1, -1> init_dbl(Eigen::Matrix<double, 1, -1>::Zero(nPK + nPD));
      std::vector<double> rate_vec(nPK + nPD, 0);

      if (rate == 0) {                     // bolus dose
        init_dbl(cmt - 1) = stan::math::value_of(amt); // bolus as initial condition
      } else {                             // infusion
        rate_vec[cmt - 1] = stan::math::value_of(rate);
      }

      const double init_dt = (rate == 0.0 || ii > 0) ? ii_dbl : 24.0;

      try {
        PMXOdeFunctorCouplingAdaptor<T_m, F> f(f_);
        PMXOdeFunctorCouplingSSAdaptor<T_m, F, integrator_type> fss(f_, integrator);
        Eigen::VectorXd y_pk_init(init_dbl.head(T_m<double>::Ncmt));
        Eigen::VectorXd y_pd_init(init_dbl.tail(n_ode_));
        std::vector<double> par_dbl(stan::math::value_of(par_));
        y_pd_init = integrator(f, y_pd_init, t0, t0 + init_dt, par_dbl, rate_vec, t0, y_pk_init, x_r_, x_i_);

        Eigen::Matrix<scalar, -1, 1> res(nPK + nPD);

        std::vector<double> scaling(nPK + nPD, 1.0);
        double rtol = integrator.as_rtol;
        double atol = integrator.as_atol;
        int max_num_step = integrator.as_max_num_step;
        res.tail(nPD) = pmx_algebra_solver_newton_tol(fss, y_pd_init, scaling, scaling,
                                                      rtol, atol, max_num_step,                                                                         
                                                      nullptr, par_, amt, rate, ii, cmt, x_r_, x_i_);

        T_m<T_par> pk_model(par_);
        res.head(nPK) = pk_model.solve(t0, amt, rate, ii, cmt);

        return res;
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
  using PkOneCptOdeModel = PKCoupledModel<PMXOneCptModel<T_par>, F>;

  template<typename T_par, typename F>
  using PkTwoCptOdeModel = PKCoupledModel<PMXTwoCptModel<T_par>, F>;
}


#endif
