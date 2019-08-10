#ifndef STAN_MATH_TORSTEN_TWOCPT_MODEL_HPP
#define STAN_MATH_TORSTEN_TWOCPT_MODEL_HPP

#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/pmx_ode_integrator.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/pk_nvars.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using torsten::PMXOdeIntegrator;

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
  };

  /**
   * two-compartment PK model. The static memebers provide
   * universal information, i.e. nb. of compartments,
   * nb. of parameters, and the RHS functor. Containing RHS
   * functor @c PMXTwoCptODE makes @c PMXTwoCptModel solvable
   * using general ODE solvers, which makes testing easier.
   *
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate dosing rate type
   * @tparam T_par PK parameters type
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  class PMXTwoCptModel {
    const T_time &t0_;
    const refactor::PKRec<T_init>& y0_;
    const std::vector<T_rate> &rate_;
    const T_par &CL_;
    const T_par &Q_;
    const T_par &V2_;
    const T_par &V3_;
    const T_par &ka_;
    const T_par k10_;
    const T_par k12_;
    const T_par k21_;
    const T_par ksum_;
    const std::vector<T_par> alpha_;
    const std::vector<T_par> par_;

  public:
    static constexpr int Ncmt = 3;
    static constexpr int Npar = 5;
    static constexpr PMXTwoCptODE f_ = PMXTwoCptODE();

    using scalar_type = typename stan::return_type<T_time, T_init, T_rate, T_par>::type;
    using init_type   = T_init;
    using time_type   = T_time;
    using par_type    = T_par;
    using rate_type   = T_rate;

  /**
   * Two-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param CL clearance
   * @param Q distributed amt
   * @param V2 central cpt vol
   * @param V3 peri cpt vol
   * @param ka absorption
   */
    PMXTwoCptModel(const T_time& t0,
                  const refactor::PKRec<T_init>& y0,
                  const std::vector<T_rate> &rate,
                  const T_par& CL,
                  const T_par& Q,
                  const T_par& V2,
                  const T_par& V3,
                  const T_par& ka) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      CL_(CL),
      Q_(Q),
      V2_(V2),
      V3_(V3),
      ka_(ka),
      k10_(CL_ / V2_),
      k12_(Q_ / V2_),
      k21_(Q_ / V3_),
      ksum_(k10_ + k12_ + k21_),
      alpha_{0.5 * (ksum_ + sqrt(ksum_ * ksum_ - 4 * k10_ * k21_)),
        0.5 * (ksum_ - sqrt(ksum_ * ksum_ - 4 * k10_ * k21_)),
        ka_},
      par_{CL_, Q_, V2_, V3_, ka_}
    {
      using stan::math::check_positive_finite;
      using stan::math::check_finite;
      const char* fun = "PMXTwoCptModel";
      check_positive_finite(fun, "CL", CL_);
      check_positive_finite(fun, "Q", Q_);
      check_positive_finite(fun, "V2", V2_);
      check_positive_finite(fun, "V3", V3_);
      check_positive_finite(fun, "ka", ka_);
    }

  /**
   * two-Compartment PK model constructor
   * FIXME need to remove parameter as this is for linode only.
   *
   * @tparam T_mp parameters class
   * @tparam Ts parameter types
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   * @param parameter ModelParameter type
   */
    template<template<typename...> class T_mp, typename... Ts>
    PMXTwoCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par,
                  const T_mp<Ts...> &parameter) :
      PMXTwoCptModel(t0, y0, rate, par[0], par[1], par[2], par[3], par[4])
    {}

  /**
   * two-compartment PK model constructor
   *
   * @param t0 initial time
   * @param y0 initial condition
   * @param rate dosing rate
   * @param par model parameters
   */
    PMXTwoCptModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par> & par) :
      PMXTwoCptModel(t0, y0, rate, par[0], par[1], par[2], par[3], par[4])
    {}


    /*
     * calculate number of @c vars for transient dosing.
     */
    static int nvars(int ncmt, int npar) {
      using stan::is_var;
      int n = 0;
      if (is_var<T_time>::value) n++; // t0
      if (is_var<T_init>::value) n += Ncmt; // y0 is fixed for twocpt model
      if (is_var<T_rate>::value) n += Ncmt; // rate is fixed for twocpt model
      if (is_var<T_par>::value) n += Npar; // par is fixed for twocpt model
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
      if (is_var<T_par>::value) n += Npar; // par is fixed for twocpt model
      return n;
    }

    /*
     * return the number @c var that will be the parameters
     * of the trasient dosing event's solution
     */
    template<typename T0>
    int nvars(const T0& t0) {
      return torsten::pk_nvars(t0, y0_, rate_, par_);
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
     * return @c vars that will be solution
     */
    template<typename T0>
    std::vector<stan::math::var> vars(const T0 t1) {
      return torsten::dsolve::pk_vars(t1, y0_, rate_, par_);
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

  /**
   * two-compartment PK model get methods
   */
    const T_time              & t0()      const { return t0_;    }
    const PKRec<T_init>    & y0()         const { return y0_;    }
    const std::vector<T_rate> & rate()    const { return rate_;  }
    const T_par               & CL()      const { return CL_;    }
    const T_par               & Q()       const { return Q_;     }
    const T_par               & V2()      const { return V2_;    }
    const T_par               & V3()      const { return V3_;    }
    const T_par               & ka()      const { return ka_;    }
    const T_par               & k10()     const { return k10_;   }
    const T_par               & k12()     const { return k12_;   }
    const T_par               & k21()     const { return k21_;   }
    const std::vector<T_par>  & par()     const { return par_;   }
    const std::vector<T_par>  & alpha()   const { return alpha_; }
    const PMXTwoCptODE         & f()       const { return f_;     }
    const int                 & ncmt ()   const { return Ncmt;   }
    const int                 & npar ()   const { return Npar;   }

  /**
   * Solve two-cpt model: analytical solution
   */
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> 
    solve(const T_time& t_next) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;

      T_time dt = t_next - t0_;

      std::vector<scalar_type> a(Ncmt, 0);
      Matrix<scalar_type, -1, 1> pred = PKRec<scalar_type>::Zero(Ncmt);

      // contribution from cpt 0
      {
        const T_par a1 = ka_ * (k21_ - alpha_[0]) / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
        const T_par a2 = ka_ * (k21_ - alpha_[1]) / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
        const T_par a3 = -(a1 + a2);
        const T_par a4 = ka_ * k12_ / ((ka_ - alpha_[0]) * (alpha_[1] - alpha_[0]));
        const T_par a5 = ka_ * k12_ / ((ka_ - alpha_[1]) * (alpha_[0] - alpha_[1]));
        const T_par a6 = -(a4 + a5);

        // bolus
        pred(0) += y0_[0] * exp(-ka_ * dt);
        pred(1) += y0_[0] * (a1 * exp(-alpha_[0] * dt) + a2 * exp(-alpha_[1] * dt) + a3 * exp(-alpha_[2] * dt));
        pred(2) += y0_[0] * (a4 * exp(-alpha_[0] * dt) + a5 * exp(-alpha_[1] * dt) + a6 * exp(-alpha_[2] * dt));

        // infusion
        pred(0) += rate_[0] * (1 - exp(-ka_ * dt)) / ka_;
        pred(1) += rate_[0] * (a1 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a2 * (1 - exp(-alpha_[1] * dt)) / alpha_[1] + a3 * (1 - exp(-alpha_[2] * dt)) / alpha_[2]);
        pred(2) += rate_[0] * (a4 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a5 * (1 - exp(-alpha_[1] * dt)) / alpha_[1] + a6 * (1 - exp(-alpha_[2] * dt)) / alpha_[2]);
      }

      // contribution from cpt 1
      {
        const T_par a1 = (k21_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
        const T_par a2 = (k21_ - alpha_[1]) / (alpha_[0] - alpha_[1]);

        // bolus
        pred(1) += y0_[1] * (a1 * exp(-alpha_[0] * dt) + a2 * exp(-alpha_[1] * dt));
        pred(2) += y0_[1] * k12_ / (alpha_[1] - alpha_[0]) * (exp(-alpha_[0] * dt) - exp(-alpha_[1] * dt));

        // infusion
        pred(1) += rate_[1] * (a1 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a2 * (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
        pred(2) += rate_[1] * k12_ / (alpha_[1] - alpha_[0]) * ((1 - exp(-alpha_[0] * dt)) / alpha_[0] - (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
      }

      // contribution from cpt 2
      {
        const T_par a1 = (k10_ + k12_ - alpha_[0]) / (alpha_[1] - alpha_[0]);
        const T_par a2 = (k10_ + k12_ - alpha_[1]) / (alpha_[0] - alpha_[1]);

        // bolus
        pred(1) += y0_[2] * k21_ / (alpha_[1] - alpha_[0]) * (exp(-alpha_[0] * dt) - exp(-alpha_[1] * dt));
        pred(2) += y0_[2] * (a1 * exp(-alpha_[0] * dt) + a2 * exp(-alpha_[1] * dt));

        // infusion
        pred(1) += rate_[2] * k21_ / (alpha_[1] - alpha_[0]) * ((1 - exp(-alpha_[0] * dt)) / alpha_[0] - (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
        pred(2) += rate_[2] * (a1 * (1 - exp(-alpha_[0] * dt)) / alpha_[0] + a2 * (1 - exp(-alpha_[1] * dt)) / alpha_[1]);
      }

      return pred;
    }

    /*
     * Solve the transient problem and return the result in
     * form of data, arranged as (solution value, grad1, grad2...)
     */
    Eigen::VectorXd solve_d(const T_time& t_next) const {
      return torsten::model_solve_d(*this, t_next);
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
    Eigen::Matrix<typename stan::return_type<T_par, T_amt, T_r, T_ii>::type, Eigen::Dynamic, 1>
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;

      using ss_scalar_type = typename stan::return_type<T_par, T_amt, T_r, T_ii>::type;

      const double inf = std::numeric_limits<double>::max();

      stan::math::check_positive_finite("steady state two-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state two-cpt solver", "cmt", cmt, 4);

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
        typename torsten::return_t<T_amt, T_r>::type dt_infus = amt/rate;
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

    /*
     * Solve the transient problem and return the result in
     * form of data, arranged as (solution value, grad1, grad2...)
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::VectorXd solve_d(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt) const {
      return torsten::model_solve_d(*this, amt, rate, ii, cmt);
    }

    /*
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<torsten::PMXOdeIntegratorId It, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const PMXOdeIntegrator<It>& integrator) const {
      return solve(amt, rate, ii, cmt);
    }

    template<typename T1, typename T2>
    inline typename torsten::return_t<T1, T2, T_par>::type
    trunc_infus_ss(const T_par& p, const T1& dt, const T2& ii) const {
      return (1 - exp(-p * dt)) * exp(-p * (ii - dt)) / (1 - exp(-p * ii)) / p;
    }
  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PMXTwoCptModel<T_time, T_init, T_rate, T_par>::Ncmt;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr int PMXTwoCptModel<T_time, T_init, T_rate, T_par>::Npar;

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PMXTwoCptODE PMXTwoCptModel<T_time, T_init, T_rate, T_par>::f_;

}

#endif
