#ifndef STAN_MATH_TORSTEN_TWOCPT_EFFCPT_MODEL_HPP
#define STAN_MATH_TORSTEN_TWOCPT_EFFCPT_MODEL_HPP

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
   * one-compartment model coupled with an effective compartment PK model.
   */
  struct PMXTwocptEffCptODE {
  /**
   * standard two compartment PK ODE RHS function
   * @tparam T0 t type
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

      scalar
        CL = parms.at(0),
        Q = parms.at(1),
        V1 = parms.at(2),
        V2 = parms.at(3),
        ka = parms.at(4),
        ke = parms.at(5),
        k10 = CL / V1,
        k12 = Q / V1,
        k21 = Q / V2;

      std::vector<scalar> y(4, 0);
      y.at(0) = -ka * x.at(0);
      y.at(1) = ka * x.at(0) - (k10 + k12) * x.at(1) + k21 * x.at(2);
      y.at(2) = k12 * x.at(1) - k21 * x.at(2);
      y.at(3) = ke * x.at(1) - ke * x.at(3);

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
      T2 Q  = parms.at(1);
      T2 V1 = parms.at(2);
      T2 V2 = parms.at(3);
      T2 ka = parms.at(4);
      T2 ke = parms.at(5);
      T2 k10 = CL / V1;
      T2 k12 = Q / V1;
      T2 k21 = Q / V2;

      Eigen::Matrix<scalar, -1, 1> y(4);
      y << -ka * x(0),
        ka * x(0) - (k10 + k12) * x(1) + k21 * x(2),
        k12 * x(1) - k21 * x(2),
        ke * x(1) - ke * x(3);

      return y;
    }
  };

  /**
   * twocpt-effcpt coupled model coupled with one effective compartment
   * model for PD effects. The static memebers provide
   * universal information: nb. of compartments,
   * nb. of parameters, and the RHS functor.
   *
   * @tparam T_par PK parameters type
   */
  template<typename T_par>
  class PMXTwocptEffCptModel {
    const T_par &CL_;
    const T_par &Q_;
    const T_par &V2_;
    const T_par &V3_;
    const T_par &ka_;
    const T_par &ke_;
    const T_par k10_;
    const T_par k12_;
    const T_par k21_;
    const std::vector<T_par> par_;
    Eigen::Matrix<T_par, -1, -1> p_;
    Eigen::Matrix<T_par, -1, 1> diag_;
    Eigen::Matrix<T_par, -1, -1> p_inv_;

  public:
    static constexpr int Ncmt = 4;
    static constexpr int Npar = 6;
    static constexpr PMXTwocptEffCptODE f_ = PMXTwocptEffCptODE();

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
    PMXTwocptEffCptModel(const T_par& CL,
                         const T_par& Q,
                         const T_par& V2,
                         const T_par& V3,
                         const T_par& ka,
                         const T_par& ke) :
      CL_(CL),
      Q_(Q),
      V2_(V2),
      V3_(V3),
      ka_(ka),
      ke_(ke),
      k10_(CL_ / V2_),
      k12_(Q_ / V2_),
      k21_(Q_ / V3_),
      par_{CL_, Q_, V2_, V3_, ka_, ke_},
      p_{Ncmt, Ncmt},
      diag_{Ncmt},
      p_inv_{Ncmt, Ncmt}
    {
      const char* fun = "PMXTwocptEffCptModel";
      stan::math::check_positive_finite(fun, "CL", CL_);
      stan::math::check_positive_finite(fun, "Q", Q_);
      stan::math::check_positive_finite(fun, "V2", V2_);
      stan::math::check_positive_finite(fun, "V3", V3_);
      stan::math::check_nonnegative(fun, "ka", ka_);
      stan::math::check_finite(fun, "ka", ka_);
      stan::math::check_positive_finite(fun, "ke", ke_);

      const T_par& a = ka_;
      const T_par& b = k10_;
      const T_par& c = k12_;
      const T_par& d = k21_;
      const T_par& e = ke_;
      const T_par s = stan::math::sqrt((b + c + d) * (b + c + d) - 4.0 * b * d);
      const T_par s1 = s + b + c + d;
      const T_par s2 = s - b - c - d;
      const T_par s3 = s + b + c - d;
      const T_par s4 = s - b - c + d;
      const T_par s5 = -s2;
      const T_par s6 = a * a - a * b - a * c - a * d + b * d;
      const T_par s7 = b * d - b * e - c * e - d * e + e * e;

      if (ka_ > 0.0) {
        p_ << (a - e) * s6 / (a * e * (a - d)), 0, 0, 0,
          1 - a/e, - 0.5 * (s1 - 2 * e)/e, - 0.5 * (s5 - 2 * e)/e, 0,
          c * (a - e) / e / (a - d), c * (s1 - 2 * e) / e / s3, -c * (s5 - 2 * e) / e / s4, 0,
          1, 1, 1 ,1;
        p_inv_ << 1.0/p_(0, 0), 0, 0, 0,
          - 0.5 * a * e * s3 * (s2 + 2 * a)/(s * s6 * (s1 - 2 * e)), - e * s3 / (s * (s1 - 2 * e)), 2 * d * e / (s * (s1 - 2 * e)), 0,
          0.5 * a * e * s4 * (s1 - 2 * a)/(s * s6 * (s5 - 2 * e)), - e * s4 / (s * (s5 - 2 * e)), -2 * d * e / (s * (s5 - 2 * e)), 0,
          a * e * (d - e) / ((a - e) * s7), e * (d - e) / s7, e * d / s7, 1;
        diag_ << -a, -0.5 * s1, 0.5 * s2, -e;
      } else {                  // ka = 0
        p_.resize(Ncmt-1, Ncmt-1);
        p_inv_.resize(Ncmt-1, Ncmt-1);
        diag_.resize(Ncmt-1);
        p_ << -0.5 * (s1 - 2 * e)/e, -0.5 * (s5 - 2 * e)/e, 0,
          0.25 * s4 * (s1 - 2 * e)/(e * d), 0.25 * s3 * (s2 + 2 * e)/(e * d), 0,
          1, 1, 1;
        T_par s8 = stan::math::sqrt(b * b + 2 * b * c - 2 * b * d + (c + d) * (c + d));
        p_inv_ <<
          -0.5 * e * (d - e) * ((d - s8) + b * e - b * d + c * d + c * e)/(s * s7), 0.5 * e * d * (s2 + 2 * e)/(s * s7), 0,
          -0.5 * e * (d - e) * ((d + s8) + b * e - b * d + c * d + c * e)/(s * s7), 0.5 * e * d * (s1 - 2 * e)/(s * s7), 0,
          e * (d -e) / s7, e * d/s7, 1;
        diag_ << -0.5 * s1, 0.5 * s2, -e;
      }
    }

  /**
   * constructor
   *
   * @param par model parameters
   */
    PMXTwocptEffCptModel(const std::vector<T_par> & par) :
      PMXTwocptEffCptModel(par[0], par[1], par[2], par[3], par[4], par[5])
    {}

    /**
     * Get methods
     */
    const std::vector<T_par>  & par()     const { return par_;   }
    const PMXTwocptEffCptODE  & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }
    const int                 & npar ()   const { return Npar;   }

    Eigen::Matrix<T_par, -1, -1> to_linode_par() const {
      Eigen::Matrix<T_par, -1, -1> linode_par(Ncmt, Ncmt);
      linode_par << -ka_, 0.0, 0.0, 0.0,
        ka_, -(k10_ + k12_), k21_, 0.0,
        0.0, k12_, -k21_, 0.0,
        0.0, ke_, 0.0, -ke_;
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
        LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p_, diag_, p_inv_);
        PMXLinOdeEigenDecompModel<T_par> linode_model(pdp);
        linode_model.solve(y, t0, t1, rate, integ);
      } else {                  // ka = 0
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
   * Solve two-cpt + effcpt model: analytical solution
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

      std::vector<ss_scalar_type> a(3, 0);
      PKRec<ss_scalar_type> pred = PKRec<ss_scalar_type>::Zero(Ncmt);

      LinOdeEigenDecomp<T_par> pdp = std::forward_as_tuple(p_, diag_, p_inv_);
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
  constexpr int PMXTwocptEffCptModel<T_par>::Ncmt;

  template<typename T_par>
  constexpr int PMXTwocptEffCptModel<T_par>::Npar;

  template<typename T_par>
  constexpr PMXTwocptEffCptODE PMXTwocptEffCptModel<T_par>::f_;
}

#endif
