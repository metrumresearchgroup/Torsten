#ifndef STAN_MATH_TORSTEN_ONECPT_EFFCPT_MODEL_HPP
#define STAN_MATH_TORSTEN_ONECPT_EFFCPT_MODEL_HPP

#include <stan/math/torsten/pmx_linode_model.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/dsolve/pmx_ode_integrator.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/prim/err/check_positive_finite.hpp>
#include <stan/math/prim/err/check_less_or_equal.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/err/check_nonnegative.hpp>

namespace torsten {
  /**
   * one-compartment model coupled with an effective compartment PK model.
   */
  struct PMXOneCptEffCptODE {
  /**
   * The RHS function of the coupled linear ODE model
   * @tparam T0 t time
   * @tparam T1 initial condition type
   * @tparam T2 parameter type
   * @tparam T3 real data/rate type
   * @param t time
   * @param x initial condition
   * @param parms parameters
   * @param rate dosing rate
   * @param x_r real data
   * @param x_i integer data
   */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    std::vector<typename stan::return_type<T0, T1, T2, T3>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<T3>& x_r,
               const std::vector<int>& x_i, std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2, T3>::type scalar;

      scalar CL = parms.at(0);
      scalar V  = parms.at(1);
      scalar ka = parms.at(2);
      scalar ke = parms.at(3);
      scalar k10 = CL / V;

      std::vector<scalar> y(3, 0);
      y.at(0) = -ka * x.at(0);
      y.at(1) = ka * x.at(0) - k10 * x.at(1);
      y.at(2) = ke * x.at(1) - ke * x.at(2);

      return y;
    }

  /**
   * Eigen version where the statement variable is an Eigen vector.
   */
    template <typename T0, typename T1, typename T2, typename T3>
    inline
    Eigen::Matrix<typename stan::return_type_t<T0, T1, T2, T3>, -1, 1>
    operator()(const T0& t,
               const Eigen::Matrix<T1, -1, 1> & x,
               std::ostream* pstream__,
               const std::vector<T2>& parms,
               const std::vector<T3>& x_r,
               const std::vector<int>& x_i) const {
      typedef typename stan::return_type_t<T0, T1, T2, T3> scalar;

      T2 CL = parms.at(0);
      T2 V  = parms.at(1);
      T2 ka = parms.at(2);
      T2 ke = parms.at(3);
      T2 k10 = CL / V;

      Eigen::Matrix<scalar, -1, 1> y(3);
      y << -ka * x(0),
        ka * x(0) - k10 * x(1),
        ke * x(1) - ke * x(2);

      return y;
    }
  };

  /**
   * onecpt-effcpt coupled PK model. The static memebers provide
   * universal information: # of compartments,
   * # of parameters, and the RHS functor.
   *
   * @tparam T_par PK parameters type
   */
  template<typename T_par>
  class PMXOneCptEffCptModel {
    const T_par &CL_;
    const T_par &V_;
    const T_par &ka_;
    const T_par &ke_;
    const T_par k10_;
    const std::vector<T_par> par_;

  public:
    static constexpr int Ncmt = 3;
    static constexpr int Npar = 4;
    static constexpr PMXOneCptEffCptODE f_ = PMXOneCptEffCptODE();

    using par_type    = T_par;

  /**
   * Constructor
   *
   * @param par model parameters
   * @param CL clearance
   * @param V central cpt vol
   * @param ka absorption
   * @param ke effect cpt
   */
    PMXOneCptEffCptModel(const T_par& CL,
                         const T_par& V,
                         const T_par& ka,
                         const T_par& ke) :
      CL_(CL),
      V_(V),
      ka_(ka),
      ke_(ke),
      k10_(CL_ / V_),
      par_{CL_, V_, ka_, ke_}
    {
      const char* fun = "PMXOneCptEffCptModel";
      stan::math::check_positive_finite(fun, "CL", CL_);
      stan::math::check_positive_finite(fun, "V", V_);
      stan::math::check_nonnegative(fun, "ka", ka_);
      stan::math::check_positive_finite(fun, "ke", ke_);
      stan::math::check_finite(fun, "ka", ka_);
    }

  /**
   * constructor
   *
   * @param par model parameters
   */
    PMXOneCptEffCptModel(const std::vector<T_par> & par) :
      PMXOneCptEffCptModel(par[0], par[1], par[2], par[3])
    {}

    /**
     * Get methods
     */
    const std::vector<T_par>  & par()     const { return par_;   }
    const PMXOneCptEffCptODE  & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }
    const int                 & npar ()   const { return Npar;   }

    Eigen::Matrix<T_par, -1, -1> to_linode_par() const {
      Eigen::Matrix<T_par, -1, -1> linode_par(Ncmt, Ncmt);
      linode_par << -ka_, 0.0, 0.0, ka_, -k10_, 0.0, 0.0, ke_, -ke_;
      return linode_par;
    }

  /**
   * analytical solution from eigen decomposition
   */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::exp;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;

      std::vector<T> a(Ncmt, 0);
      Eigen::Matrix<T, -1, 1> pred = torsten::PKRec<T>::Zero(Ncmt);

      if (dt < 1.e-12) return;

      if (ka_ > 0.0) {
        Eigen::Matrix<T_par, -1, -1> p(Ncmt, Ncmt), p_inv(Ncmt, Ncmt);
        p << (ka_ - k10_) * (ka_ - ke_) / ( ka_ * ke_ ), 0, 0,
          - (ka_ - ke_) / ke_, - (k10_ - ke_) / ke_, 0,
          1, 1, 1;
        p_inv << ka_ * ke_/ (ka_ - ke_) / (ka_ - k10_), 0, 0,
          - ka_ * ke_/ (ka_ - k10_) / (k10_ - ke_), -ke_ / (k10_ - ke_), 0,
            ka_ * ke_/ (ka_ - ke_) / (k10_ - ke_), ke_ / (k10_ - ke_), 1;
        Eigen::Matrix<T_par, -1, 1> diag(Ncmt);
        diag << -ka_, -k10_, -ke_;
        LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p, diag, p_inv);
        PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
        linode_model.solve(y, t0, t1, rate, integ);
      } else {
        y(0) += rate[0] * dt;
        Eigen::Matrix<T_par, -1, -1> p(Ncmt-1, Ncmt-1), p_inv(Ncmt-1, Ncmt-1);
        Eigen::Matrix<T_par, -1, 1> diag(Ncmt-1);
        p << -(k10_ - ke_)/ke_, 0, 1, 1;
        p_inv << -ke_ / (k10_ - ke_), 0, ke_ / (k10_ - ke_), 1;
        diag << -k10_, -ke_;
        LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p, diag, p_inv);
        PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
        PKRec<T> y2 = y.tail(Ncmt - 1);
        std::vector<T1> rate2(rate.begin() + 1, rate.end());
        linode_model.solve(y2, t0, t1, rate2, integ);
        y.tail(Ncmt - 1) = y2;
      }
    }

  /**
   * Solve the coupled model: analytical solution
   */
    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate) const {
      const dsolve::PMXAnalyiticalIntegrator integ;
      solve(y, t0, t1, rate, integ);
    }

  /**
   * Steady state model. We have to consider
   * different scenarios: bolus/multiple truncated infusion/const infusion
   *
   * @tparam T_amt amt type
   * @tparam T_r rate type
   * @tparam T_ii dosing interval type
   * @param t0 start time for steady state events
   * @param amt dosing amount
   * @param ii dosing interval
   * @param rate infusion rate
   * @param cmt dosing compartment
   */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const {
      using ss_scalar_type = typename stan::return_type<T_par, T_amt, T_r, T_ii>::type;

      stan::math::check_positive_finite("steady state one-cpt/eff-cpt solver", "cmt", cmt);
      stan::math::check_less_or_equal("steady state one-cpt/eff-cpt solver", "cmt", cmt, Ncmt);
      stan::math::check_positive_finite("steady state one-cpt/eff-cpt solver", "ka", ka_);

      std::vector<ss_scalar_type> a(3, 0);
      PKRec<ss_scalar_type> pred = PKRec<ss_scalar_type>::Zero(Ncmt);

      Eigen::Matrix<T_par, -1, -1> p(Ncmt, Ncmt), p_inv(Ncmt, Ncmt);
      p << (ka_ - k10_) * (ka_ - ke_) / ( ka_ * ke_ ), 0, 0,
        - (ka_ - ke_) / ke_, - (k10_ - ke_) / ke_, 0,
        1, 1, 1;
      p_inv << ka_ * ke_/ (ka_ - ke_) / (ka_ - k10_), 0, 0,
        - ka_ * ke_/ (ka_ - k10_) / (k10_ - ke_), -ke_ / (k10_ - ke_), 0,
        - ka_ * ke_/ (ka_ - ke_) / (k10_ - ke_), ke_ / (k10_ - ke_), 0;
      Eigen::Matrix<T_par, -1, 1> diag(Ncmt);
      diag << -ka_, -k10_, -ke_;
      LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p, diag, p_inv);
      PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
      pred = linode_model.solve(t0, amt, rate, ii, cmt);

      return pred;
    }

    /**
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const dsolve::PMXAnalyiticalIntegrator& integrator) const {
      return solve(t0, amt, rate, ii, cmt);
    }
  };

  template<typename T_par>
  constexpr int PMXOneCptEffCptModel<T_par>::Ncmt;

  template<typename T_par>
  constexpr int PMXOneCptEffCptModel<T_par>::Npar;

  template<typename T_par>
  constexpr PMXOneCptEffCptODE PMXOneCptEffCptModel<T_par>::f_;
}

#endif
