#ifndef STAN_MATH_TORSTEN_TWOCPT_MODEL_HPP
#define STAN_MATH_TORSTEN_TWOCPT_MODEL_HPP

#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/prim/err/check_positive_finite.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/err/check_nonnegative.hpp>

namespace torsten {
  /**
   * standard two compartment PK ODE functor.
   */
  struct PMXTwoCptODE {
  /**
   * standard two compartment PK ODE RHS function
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
               const std::vector<int>& dummy, std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

      scalar
        CL = parms.at(0),
        Q = parms.at(1),
        V1 = parms.at(2),
        V2 = parms.at(3),
        ka = parms.at(4),
        k10 = CL / V1,
        k12 = Q / V1,
        k21 = Q / V2;

      std::vector<scalar> y(3, 0);
      y.at(0) = -ka * x.at(0);
      y.at(1) = ka * x.at(0) - (k10 + k12) * x.at(1) + k21 * x.at(2);
      y.at(2) = k12 * x.at(1) - k21 * x.at(2);

      return y;
    }

    /** 
     * Eigen::Matrix vection
     * 
     */
    template <typename T0, typename T1, typename T2>
    inline
    Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
    operator()(const T0& t,
               const Eigen::Matrix<T1, -1, 1>& x,
               std::ostream* msg,
               const std::vector<T2>& parms,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      typedef typename stan::return_type_t<T0, T1, T2> scalar;

      T2 CL = parms.at(0);
      T2 Q = parms.at(1);
      T2 V1 = parms.at(2);
      T2 V2 = parms.at(3);
      T2 ka = parms.at(4);
      T2 k10 = CL / V1;
      T2 k12 = Q / V1;
      T2 k21 = Q / V2;

      Eigen::Matrix<scalar, -1, 1> y(3);
      y(0) = -ka * x(0);
      y(1) = ka * x(0) - (k10 + k12) * x(1) + k21 * x(2);
      y(2) = k12 * x(1) - k21 * x(2);

      return y;
    }
  };

  /**
   * two-compartment PK model. The static memebers provide
   * universal information, i.e. nb. of compartments,
   * nb. of parameters, and the RHS functor. Containing RHS
   * functor @c PMXTwoCptODE makes @c PMXTwoCptModel solvable
   * using general ODE solvers, which makes testing easier.
   *
   * @tparam T_time t type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   */
  template<typename T_par>
  class PMXTwoCptModel {
    const T_par &CL_;
    const T_par &Q_;
    const T_par &V2_;
    const T_par &V3_;
    const T_par &ka_;
    const T_par k10_;
    const T_par k12_;
    const T_par k21_;
    const T_par ksum_;
    std::vector<T_par> alpha_;
    const std::vector<T_par> par_;
    Eigen::Matrix<T_par, -1, -1> p_;
    Eigen::Matrix<T_par, -1, 1> diag_;
    Eigen::Matrix<T_par, -1, -1> p_inv_;

  public:
    static constexpr int Ncmt = 3;
    static constexpr int Npar = 5;
    static constexpr PMXTwoCptODE f_ = PMXTwoCptODE();

    using par_type    = T_par;

  /**
   * Two-compartment PK model constructor
   *
   * @param rate dosing rate
   * @param par model parameters
   * @param CL clearance
   * @param Q distributed amt
   * @param V2 central cpt vol
   * @param V3 peri cpt vol
   * @param ka absorption
   */
    PMXTwoCptModel(const T_par& CL,
                   const T_par& Q,
                   const T_par& V2,
                   const T_par& V3,
                   const T_par& ka) :
      CL_(CL),
      Q_(Q),
      V2_(V2),
      V3_(V3),
      ka_(ka),
      k10_(CL_ / V2_),
      k12_(Q_ / V2_),
      k21_(Q_ / V3_),
      ksum_(k10_ + k12_ + k21_),
      alpha_{0.5 * (ksum_ + stan::math::sqrt(ksum_ * ksum_ - 4 * k10_ * k21_)),
        0.5 * (ksum_ - stan::math::sqrt(ksum_ * ksum_ - 4 * k10_ * k21_)),
        ka_},
      par_{CL_, Q_, V2_, V3_, ka_},
      p_{Ncmt, Ncmt},
      diag_{Ncmt},
      p_inv_{Ncmt, Ncmt} {
      const char* fun = "PMXTwoCptModel";
      stan::math::check_positive_finite(fun, "CL", CL_);
      stan::math::check_positive_finite(fun, "Q", Q_);
      stan::math::check_positive_finite(fun, "V2", V2_);
      stan::math::check_positive_finite(fun, "V3", V3_);
      stan::math::check_nonnegative(fun, "ka", ka_);
      stan::math::check_finite(fun, "ka", ka_);

      T_par s = stan::math::sqrt(-4.0 * k10_ * k21_ + (k10_ + k12_ + k21_) * (k10_ + k12_ + k21_));
      T_par q = ka_ * ka_ - ka_ * k10_ - ka_ * k12_ - ka_ * k21_ + k10_ * k21_;
      T_par w = k10_ + k12_ - k21_;
      if (ka_ > 0) {
        p_ << q / (ka_ * k12_), 0, 0,
          -(ka_ - k21_)/k12_, -0.5 * (w + s) / k12_, -0.5 * (w - s) / k12_, 
          1, 1, 1;
        p_inv_ << ka_ * k12_/q, 0, 0,
          -ka_ * k12_ * ( 2.0 * ka_ - k10_ - k12_ - k21_ + s) / (2.0 * q * s), -k12_ / s, 0.5 * (s - w) / s,
          -ka_ * k12_ * (-2.0 * ka_ + k10_ + k12_ + k21_ + s) / (2.0 * q * s),  k12_ / s, 0.5 * (s + w) / s;
        diag_ << -ka_,
          -0.5 * (k10_ + k12_ + k21_ + s),
          -0.5 * (k10_ + k12_ + k21_ - s);
      } else {
        p_.resize(Ncmt-1, Ncmt-1);
        p_inv_.resize(Ncmt-1, Ncmt-1);
        diag_.resize(Ncmt-1);
        p_ << -0.5 * (w + s) / k12_, -0.5 * (w - s) / k12_, 1, 1;
        p_inv_ << -k12_/s,  0.5 * (s - w) / s, k12_/s, 0.5 * (s + w) / s;
        diag_ << -0.5 * (k12_ + k10_ + k21_ + s),
          -0.5 * (k12_ + k10_ + k21_ - s);
      }
    }

  /**
   * two-compartment PK model constructor
   *
   * @param t0 initial time
   * @param rate dosing rate
   * @param par model parameters
   */
    PMXTwoCptModel(const std::vector<T_par> & par) :
      PMXTwoCptModel(par[0], par[1], par[2], par[3], par[4])
    {}

    /**
     * two-compartment PK model get methods
     */
    const std::vector<T_par>  & par()     const { return par_;   }
    const PMXTwoCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }
    const int                 & npar ()   const { return Npar;   }


    Eigen::Matrix<T_par, -1, -1> to_linode_par() const {
      Eigen::Matrix<T_par, -1, -1> linode_par(Ncmt, Ncmt);
      linode_par << -ka_, 0.0, 0.0, ka_, -(k10_ + k12_), k21_, 0.0, k12_, -k21_;
      return linode_par;
    }

  /**
   * Solve two-cpt model: analytical solution
   */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::exp;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;

      if (ka_ > 0.0) {
        LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p_, diag_, p_inv_);
        PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
        linode_model.solve(y, t0, t1, rate, integ);
      } else {
        y(0) += rate[0] * dt;
        LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p_, diag_, p_inv_);
        PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
        PKRec<T> y2 = y.tail(Ncmt - 1);
        std::vector<T1> rate2(rate.begin() + 1, rate.end());
        linode_model.solve(y2, t0, t1, rate2, integ);
        y.tail(Ncmt - 1) = y2;
      }
    }

  /**
   * Solve two-cpt model: analytical solution for benchmarking & testing
   */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve_analytical(PKRec<T>& y,
                          const Tt0& t0, const Tt1& t1,
                          const std::vector<T1>& rate,
                          const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::exp;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;

      std::vector<T> a(Ncmt, 0);
      Eigen::Matrix<T, -1, 1> pred = torsten::PKRec<T>::Zero(Ncmt);

      // contribution from cpt 0
      {
        if (ka_ > 0.0) {
          const T_par a1 = ka_ * (k21_ - alpha_[0]) / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          const T_par a2 = ka_ * (k21_ - alpha_[1]) / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          const T_par a3 = -(a1 + a2);
          const T_par a4 = ka_ * k12_ / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          const T_par a5 = ka_ * k12_ / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          const T_par a6 = -(a4 + a5);

          // bolus
          pred(0) += y(0) * exp(-ka_ * dt);
          pred(1) += y(0) * (a1 * exp(-alpha_[0] * dt) + a2 * exp(-alpha_[1] * dt) + a3 * exp(-alpha_[2] * dt));
          pred(2) += y(0) * (a4 * exp(-alpha_[0] * dt) + a5 * exp(-alpha_[1] * dt) + a6 * exp(-alpha_[2] * dt));

          // infusion
          pred(0) += rate[0] * (1 - exp(-ka_ * dt)) / ka_;
          pred(1) += rate[0] * (a1 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a2 * (1 - exp(-alpha_[1] * dt)) / alpha_[1] + a3 * (1 - exp(-alpha_[2] * dt)) / alpha_[2]);
          pred(2) += rate[0] * (a4 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a5 * (1 - exp(-alpha_[1] * dt)) / alpha_[1] + a6 * (1 - exp(-alpha_[2] * dt)) / alpha_[2]);
        } else {
          pred(0) += y(0) + rate[0] * dt;
        }
      }

      // contribution from cpt 1
      {
        const T_par a1 = (k21_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
        const T_par a2 = (k21_ - alpha_[1]) / (alpha_[0] - alpha_[1]);

        // bolus
        pred(1) += y(1) * (a1 * exp(-alpha_[0] * dt) + a2 * exp(-alpha_[1] * dt));
        pred(2) += y(1) * k12_ / (alpha_[1] - alpha_[0]) * (exp(-alpha_[0] * dt) - exp(-alpha_[1] * dt));

        // infusion
        pred(1) += rate[1] * (a1 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a2 * (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
        pred(2) += rate[1] * k12_ / (alpha_[1] - alpha_[0]) * ((1 - exp(-alpha_[0] * dt)) / alpha_[0] - (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
      }

      // contribution from cpt 2
      {
        const T_par a1 = (k10_ + k12_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
        const T_par a2 = (k10_ + k12_ - alpha_[1]) / (alpha_[0] - alpha_[1]);

        // bolus
        pred(1) += y(2) * k21_ / (alpha_[1] - alpha_[0]) * (exp(-alpha_[0] * dt) - exp(-alpha_[1] * dt));
        pred(2) += y(2) * (a1 * exp(-alpha_[0] * dt) + a2 * exp(-alpha_[1] * dt));

        // infusion
        pred(1) += rate[2] * k21_ / (alpha_[1] - alpha_[0]) * ((1 - exp(-alpha_[0] * dt)) / alpha_[0] - (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
        pred(2) += rate[2] * (a1 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a2 * (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
      }
      y = pred;
    }

  /**
   * Solve two-cpt model: analytical solution
   */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate) const {
      const dsolve::PMXAnalyiticalIntegrator integ;
      solve(y, t0, t1, rate, integ);
    }

  /**
   * Solve two-cpt model: analytical solution used for benchmarking & testing
   */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve_analytical(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate) const {
      const dsolve::PMXAnalyiticalIntegrator integ;
      solve_analytical(y, t0, t1, rate, integ);
    }

  /**
   * Solve two-cpt steady state model. We have to consider
   * different scenarios: bolus/multiple truncated infusion/const infusion
   *
   * @tparam T_amt amt type
   * @param amt dosing amount
   * @param ii dosing interval
   * @param cmt dosing compartment
   */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;
      using stan::math::exp;
      using stan::math::matrix_exp;
      using stan::math::value_of;
      using stan::math::mdivide_left;
      using stan::math::multiply;

      using ss_scalar_type = typename stan::return_type<T_par, T_amt, T_r, T_ii>::type;

      stan::math::check_positive_finite("steady state two-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state two-cpt solver", "cmt", cmt, 4);
      stan::math::check_positive_finite("steady state two-cpt solver", "ka", ka_);

      LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p_, diag_, p_inv_);
      PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
      PKRec<ss_scalar_type> pred = linode_model.solve(t0, amt, rate, ii, cmt);

      return pred;
    }

    /*
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const dsolve::PMXAnalyiticalIntegrator& integrator) const {
      return solve(t0, amt, rate, ii, cmt);
    }

    /**
     * analytical solution used for testing
     * 
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, -1, 1>
    solve_analytical(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;
      using stan::math::exp;

      using ss_scalar_type = typename stan::return_type<T_par, T_amt, T_r, T_ii>::type;

      const double inf = std::numeric_limits<double>::max();

      stan::math::check_positive_finite("steady state two-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state two-cpt solver", "cmt", cmt, 4);
      stan::math::check_positive_finite("steady state two-cpt solver", "ka", ka_);

      std::vector<ss_scalar_type> a(3, 0);
      Matrix<ss_scalar_type, -1, 1> pred = Matrix<ss_scalar_type, 1, Dynamic>::Zero(3);

      if (rate == 0) {  // bolus dose
        switch (cmt) {
        case 1:
          pred(0) = amt / (exp(ka_ * ii) - 1.0);
          a[0] = ka_ * (k21_ - alpha_[0]) / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          a[1] = ka_ * (k21_ - alpha_[1]) / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          a[2] = -(a[0] + a[1]);
          pred(1) = amt * (a[0] / (exp(alpha_[0] * ii)-1.0) + a[1] / (exp(alpha_[1] * ii)-1.0) + a[2] / (exp(alpha_[2] * ii)-1.0));
          a[0] = ka_ * k12_ / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          a[1] = ka_ * k12_ / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          a[2] = -(a[0] + a[1]);
          pred(2) = amt * (a[0] / (exp(alpha_[0] * ii)-1.0) + a[1] / (exp(alpha_[1] * ii)-1.0) + a[2] / (exp(alpha_[2] * ii)-1.0));
          break;
        case 2:
          a[0] = (k21_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
          a[1] = (k21_ - alpha_[1]) / (alpha_[0] - alpha_[1]);
          pred(1) = amt * (a[0] / (exp(alpha_[0] * ii)-1.0) + a[1] / (exp(alpha_[1] * ii)-1.0));
          a[0] = k12_ / (alpha_[1] - alpha_[0]);
          a[1] = -a[0];
          pred(2) = amt * (a[0] / (exp(alpha_[0] * ii)-1.0) + a[1] / (exp(alpha_[1] * ii)-1.0));
          break;
        case 3:
          a[0] = k21_ / (alpha_[1] - alpha_[0]);
          a[1] = -a[0];
          pred(1) = amt * (a[0] / (exp(alpha_[0] * ii)-1.0) + a[1] / (exp(alpha_[1] * ii)-1.0));
          a[0] = (k10_ + k12_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
          a[1] = (k10_ + k12_ - alpha_[1]) / (alpha_[0] - alpha_[1]);
          pred(2) = amt * (a[0] / (exp(alpha_[0] * ii)-1.0) + a[1] / (exp(alpha_[1] * ii)-1.0));
          break;
        }
      } else if (ii > 0) {  // multiple truncated infusions
        typename stan::return_type_t<T_amt, T_r> dt_infus = amt/rate;
        static const char* function("Steady State Event");
        torsten::check_mti(amt, stan::math::value_of(dt_infus), ii, function);
        switch (cmt) {
        case 1:
          pred(0) = rate * trunc_infus_ss(alpha_[2], dt_infus, ii);
          a[0] = ka_ * (k21_ - alpha_[0]) / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          a[1] = ka_ * (k21_ - alpha_[1]) / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          a[2] = - (a[0] + a[1]);
          pred(1) = rate * (a[0] * trunc_infus_ss(alpha_[0], dt_infus, ii) +
                            a[1] * trunc_infus_ss(alpha_[1], dt_infus, ii) +
                            a[2] * trunc_infus_ss(alpha_[2], dt_infus, ii) );
          a[0] = ka_ * k12_ / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          a[1] = ka_ * k12_ / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          a[2] = -(a[0] + a[1]);
          pred(2) = rate * (a[0] * trunc_infus_ss(alpha_[0], dt_infus, ii) +
                            a[1] * trunc_infus_ss(alpha_[1], dt_infus, ii) +
                            a[2] * trunc_infus_ss(alpha_[2], dt_infus, ii) );
          break;
        case 2:
          a[0] = (k21_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
          a[1] = (k21_ - alpha_[1]) / (alpha_[0] - alpha_[1]);
          pred(1) = rate * (a[0] * trunc_infus_ss(alpha_[0], dt_infus, ii) +
                            a[1] * trunc_infus_ss(alpha_[1], dt_infus, ii) );
          a[0] = k12_ / (alpha_[1] - alpha_[0]);
          a[1] = -a[0];
          pred(2) = rate * (a[0] * trunc_infus_ss(alpha_[0], dt_infus, ii) +
                            a[1] * trunc_infus_ss(alpha_[1], dt_infus, ii) );
          break;
        case 3:
          a[0] = k21_ / (alpha_[1] - alpha_[0]);
          a[1] = -a[0];
          pred(1) = rate * (a[0] * trunc_infus_ss(alpha_[0], dt_infus, ii) +
                            a[1] * trunc_infus_ss(alpha_[1], dt_infus, ii) );
          a[0] = (k10_ + k12_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
          a[1] = (k10_ + k12_ - alpha_[1]) / (alpha_[0] - alpha_[1]);
          pred(2) = rate * (a[0] * trunc_infus_ss(alpha_[0], dt_infus, ii) +
                            a[1] * trunc_infus_ss(alpha_[1], dt_infus, ii) );
          break;
        }
      } else {  // constant infusion
        switch (cmt) {
        case 1:
          pred(0) = rate / alpha_[2];
          a[0] = ka_ * (k21_ - alpha_[0]) / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          a[1] = ka_ * (k21_ - alpha_[1]) / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          a[2] = -(a[0] + a[1]);
          pred(1) = rate * (a[0] / alpha_[0] + a[1] / alpha_[1] + a[2] / alpha_[2]);
          a[0] = ka_ * k12_ / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
          a[1] = ka_ * k12_ / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
          a[2] = -(a[0] + a[1]);
          pred(2) = rate * (a[0] / alpha_[0] + a[1] / alpha_[1] + a[2] / alpha_[2]);
          break;
        case 2:
          a[0] = (k21_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
          a[1] = (k21_ - alpha_[1]) / (alpha_[0] - alpha_[1]);
          pred(1) = rate * (a[0] / alpha_[0] + a[1] / alpha_[1] + a[2] / alpha_[2]);
          a[0] = k12_ / (alpha_[1] - alpha_[0]);
          a[1] = -a[0];
          pred(2) = rate * (a[0] / alpha_[0] + a[1] / alpha_[1] + a[2] / alpha_[2]);
          break;
        case 3:
          a[0] = k21_ / (alpha_[1] - alpha_[0]);
          a[1] = -a[0];
          pred(1) = rate * (a[0] / alpha_[0] + a[1] / alpha_[1] + a[2] / alpha_[2]);
          a[0] = (k10_ + k12_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
          a[1] = (k10_ + k12_ - alpha_[1]) / (alpha_[0] - alpha_[1]);
          pred(2) = rate * (a[0] / alpha_[0] + a[1] / alpha_[1] + a[2] / alpha_[2]);
        }
      }
      return pred;
    }

    template<typename T1, typename T2>
    inline typename stan::return_type_t<T1, T2, T_par>
    trunc_infus_ss(const T_par& p, const T1& dt, const T2& ii) const {
      return (1 - exp(-p * dt)) * exp(-p * (ii - dt)) / (1 - exp(-p * ii)) / p;
    }
  };

  template<typename T_par>
  constexpr int PMXTwoCptModel<T_par>::Ncmt;

  template<typename T_par>
  constexpr int PMXTwoCptModel<T_par>::Npar;

  template<typename T_par>
  constexpr PMXTwoCptODE PMXTwoCptModel<T_par>::f_;

}

#endif
