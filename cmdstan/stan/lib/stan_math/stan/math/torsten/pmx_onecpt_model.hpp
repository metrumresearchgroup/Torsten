#ifndef STAN_MATH_TORSTEN_ONECPT_MODEL_HPP
#define STAN_MATH_TORSTEN_ONECPT_MODEL_HPP

#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/err/check_positive_finite.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/err/check_nonnegative.hpp>

namespace torsten {

  /**
   * standard one compartment with 1st order absorption PK ODE functor
   */
  struct PMXOneCptODE {
  /**
   * standard one compartment PK ODE RHS function
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
               const std::vector<T3>& rate,
               const std::vector<int>& dummy,
               std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

      scalar CL = parms.at(0), V1 = parms.at(1), ka = parms.at(2), k10 = CL / V1;
      std::vector<scalar> y(2, 0);

      y.at(0) = -ka * x.at(0);
      y.at(1) = ka * x.at(0) - k10 * x.at(1);

      return y;
    }

    /** 
     * Eigen::Matrix version
     * 
     */
    template <typename T0, typename T1, typename T2>
    inline
    Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
    operator()(const T0& t,
               const Eigen::Matrix<T1, -1, 1>& x,
               std::ostream* pstream__,
               const std::vector<T2>& parms,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      typedef typename stan::return_type<T0, T1, T2>::type scalar;

      T2 CL = parms.at(0), V1 = parms.at(1), ka = parms.at(2), k10 = CL / V1;
      Eigen::Matrix<scalar, -1, 1> y(2);

      y(0) = -ka * x(0);
      y(1) = ka * x(0) - k10 * x(1);

      return y;
    }
  };

  /**
   * One-compartment PK model. The static memebers provide
   * universal information, i.e. nb. of compartments,
   * nb. of parameters, and the RHS functor.
   *
   * @tparam T_time t type
   * @tparam T_par PK parameters type
   */
  template<typename T_par>
  class PMXOneCptModel {
    const T_par &CL_;
    const T_par &V2_;
    const T_par &ka_;
    const T_par k10_;
    const std::vector<T_par> alpha_;
    const std::vector<T_par> par_;

  public:
    static constexpr int Ncmt = 2;
    static constexpr int Npar = 3;
    static constexpr PMXOneCptODE f_ = PMXOneCptODE();

    using par_type    = T_par;

    /**
     * One-compartment PK model constructor
     *
     * @param t0 initial time
     * @param y0 initial condition
     * @param rate dosing rate
     * @param par model parameters
     * @param CL clearance
     * @param V2 central cpt vol
     * @param ka absorption
     */
    PMXOneCptModel(const T_par& CL,
                   const T_par& V2,
                   const T_par& ka) :
      CL_(CL),
      V2_(V2),
      ka_(ka),
      k10_(CL / V2),
      alpha_{k10_, ka_},
      par_{CL_, V2_, ka_}
    {
      const char* fun = "PMXOneCptModel";
      stan::math::check_positive_finite(fun, "CL", CL_);
      stan::math::check_positive_finite(fun, "V2", V2_);
      stan::math::check_nonnegative(fun, "ka", ka_);
      stan::math::check_finite(fun, "ka", ka_);
    }

    /**
     * One-compartment PK model constructor
     *
     * @param t0 initial time
     * @param y0 initial condition
     * @param rate dosing rate
     * @param par model parameters
     */
    PMXOneCptModel(const std::vector<T_par> & par) :
      PMXOneCptModel(par.at(0), par.at(1), par.at(2))
    {}

    /*
     * return @c vars that will be steady-state
     * solution. For SS solution @c rate_ or @ y0_ will not
     * be in the solution.
     */
    template<typename T_a, typename T_r, typename T_ii>
    std::vector<stan::math::var> vars(const T_a& a, const T_r& r, const T_ii& ii) {
      return torsten::dsolve::pk_vars(a, r, ii, par_);
    }

    /**
     * One-compartment PK model get methods
     */
    const std::vector<T_par>  & par()     const { return par_;   }
    const PMXOneCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }

    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(Eigen::Matrix<T, -1, 1>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::exp;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;
      Eigen::Matrix<T, -1, 1> pred = torsten::PKRec<T>::Zero(Ncmt);

      if (ka_ > 0.0) {
        Eigen::Matrix<T_par, -1, -1> p(Ncmt, Ncmt), p_inv(Ncmt, Ncmt);
        Eigen::Matrix<T_par, -1, 1> diag(Ncmt);
        p << -(ka_ - k10_)/ka_, 0, 1, 1;
        p_inv << -ka_ / (ka_ - k10_), 0, ka_ / (ka_ - k10_), 1;
        diag << -ka_, -k10_;
        LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p, diag, p_inv);
        PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
        linode_model.solve(y, t0, t1, rate, integ);
      } else {
        typename stan::return_type_t<T_par, Tt0, Tt1> exp1 = exp(-k10_ * dt);

        // contribution from cpt 1 bolus dose
        pred(1) += y(1) * exp1;

        // contribution from cpt 1 infusion dose
        pred(1) += rate[1] * (1 - exp1) / k10_;

        pred(0) += y(0) + rate[0] * dt;
        y = pred;
      }
    }

    /**
     * Solve two-cpt model: analytical solution for benchmarking & testing
     */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve_analytical(Eigen::Matrix<T, -1, 1>& y,
                          const Tt0& t0, const Tt1& t1,
                          const std::vector<T1>& rate,
                          const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::exp;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;
      Eigen::Matrix<T, -1, 1> pred = torsten::PKRec<T>::Zero(Ncmt);

      typename stan::return_type_t<T_par, Tt0, Tt1> exp1 = exp(-k10_ * dt);
      typename stan::return_type_t<T_par, Tt0, Tt1> exp2 = exp(-ka_ * dt);

      // contribution from cpt 1 bolus dose
      pred(1) += y(1) * exp1;

      // contribution from cpt 1 infusion dose
      pred(1) += rate[1] * (1 - exp1) / k10_;

      if (ka_ > 0.0) {
        // contribution from cpt 0 bolus dose
        pred(0) += y(0) * exp2;
        pred(1) += y(0) * ka_ / (ka_ - k10_) * (exp1 - exp2);

        // contribution from cpt 0 infusion dose
        pred(0) += rate[0] * (1 - exp2) / ka_;
        pred(1) += rate[0] * ka_ / (ka_ - k10_) * ((1 - exp1) / k10_ - (1 - exp2) / ka_);
      } else {
        // no absorption, GUT is accumulating dosages.
        pred(0) += y(0) + rate[0] * dt;
      }

      y = pred;
    }

    /**
     * Solve two-cpt model: analytical solution
     */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(Eigen::Matrix<T, -1, 1>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate) const {
      const dsolve::PMXAnalyiticalIntegrator integ;
      solve(y, t0, t1, rate, integ);
    }

    /**
     * Solve two-cpt model: analytical solution for bencharmking & testing
     */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve_analytical(Eigen::Matrix<T, -1, 1>& y,
                          const Tt0& t0, const Tt1& t1,
                          const std::vector<T1>& rate) const {
      const dsolve::PMXAnalyiticalIntegrator integ;
      solve_analytical(y, t0, t1, rate, integ);
    }

    /**
     * Solve one-cpt model: steady state solution
     *
     * @tparam T_time dosing interval time type
     * @tparam T_model ODE model type
     * @tparam T_amt dosing amount type
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, Eigen::Dynamic, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, int cmt) const { //NOLINT
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;

      const double inf = std::numeric_limits<double>::max();  // "infinity"

      stan::math::check_positive_finite("steady state one-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state one-cpt solver", "cmt", cmt, 3);
      stan::math::check_positive_finite("steady state one-cpt solver", "ka", ka_);

      using ss_scalar_type = typename stan::return_type_t<T_par, T_amt, T_r, T_ii>;
      PKRec<ss_scalar_type> pred = PKRec<ss_scalar_type>::Zero(Ncmt);

      Eigen::Matrix<T_par, -1, -1> p(Ncmt, Ncmt), p_inv(Ncmt, Ncmt);
      Eigen::Matrix<T_par, -1, 1> diag(Ncmt);
      p.setZero();
      p_inv.setZero();
      diag.setZero();
      p << -(ka_ - k10_)/ka_, 0, 1, 1;
      p_inv << -ka_ / (ka_ - k10_), 0, ka_ / (ka_ - k10_), 1;
      diag << -ka_, -k10_;
      LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p, diag, p_inv);
      PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
      pred = linode_model.solve(t0, amt, rate, ii, cmt);

      return pred;
    }

    /*
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<typename integrator_type, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_ii, T_par>, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const integrator_type& integrator) const {
      return solve(t0, amt, rate, ii, cmt);
    }

    /**
     * analytical solution for SS case, for testing 
     * 
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, Eigen::Dynamic, 1>
    solve_analytical(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, int cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;

      const double inf = std::numeric_limits<double>::max();  // "infinity"

      stan::math::check_positive_finite("steady state one-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state one-cpt solver", "cmt", cmt, 3);
      stan::math::check_positive_finite("steady state one-cpt solver", "ka", ka_);

      using ss_scalar_type = typename stan::return_type<T_par, T_amt, T_r, T_ii>::type;

      std::vector<ss_scalar_type> a(2, 0);
      Matrix<ss_scalar_type, -1, 1> pred = Matrix<ss_scalar_type, -1, 1>::Zero(2);
      if (rate == 0) {  // bolus dose
        switch (cmt) {
        case 1:
          pred(0) = amt / (exp(ka_ * ii) - 1.0);
          pred(1) = amt * ka_ / (ka_ - k10_) * (1.0 / (exp(k10_ * ii)-1.0) - 1.0 / (exp(ka_ * ii)-1.0));
          break;
        case 2:
          pred(1) = amt / (exp(k10_ * ii) - 1.0);
          break;
        }
      } else if (ii > 0) {  // multiple truncated infusions
        typename stan::return_type_t<T_amt, T_r> dt_infus = amt/rate;
        static const char* function("Steady State Event");
        torsten::check_mti(amt, stan::math::value_of(dt_infus), ii, function);

        // since we've checked ii > dt_infus
        // TODO: ii < dt_infus
        switch (cmt) {
        case 1:
          pred(0) = rate / ka_ * (1 - exp(-ka_ * dt_infus)) * exp(-ka_ * (ii - dt_infus)) / (1 - exp(-ka_ * ii));
          pred(1) = rate * ka_ / (k10_ * (ka_ - k10_)) * (1 - exp(-k10_ * dt_infus)) * exp(-k10_ * (ii - dt_infus)) / (1 - exp(-k10_ * ii))
            - rate / (ka_ - k10_) * (1 - exp(-ka_ * dt_infus)) * exp(-ka_ * (ii - dt_infus)) / (1 - exp(-ka_ * ii));
          break;
        case 2:
          pred(1) = rate / k10_ * (1 - exp(-k10_ * dt_infus)) * exp(-k10_ * (ii - dt_infus)) / (1 - exp(-k10_ * ii));
          break;
        }
      } else {  // constant infusion
        switch (cmt) {
        case 1:
          pred(0) = rate / ka_;
          pred(1) = rate * (ka_ / (k10_ * (ka_ - k10_)) - 1.0 / (ka_ - k10_));
          break;
        case 2:
          pred(1) = rate / k10_;
          break;
        }
      }
      return pred;
    }
  };



  template<typename T_par>
  constexpr int PMXOneCptModel<T_par>::Ncmt;

  template<typename T_par>
  constexpr int PMXOneCptModel<T_par>::Npar;

  template<typename T_par>
  constexpr PMXOneCptODE PMXOneCptModel<T_par>::f_;

}

#endif
