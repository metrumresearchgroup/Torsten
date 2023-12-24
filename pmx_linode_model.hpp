#ifndef STAN_MATH_PMX_LINODE_MODEL_HPP
#define STAN_MATH_PMX_LINODE_MODEL_HPP

#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/rev/fun/multiply.hpp>
#include <stan/math/prim/fun/matrix_exp.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/rev/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/torsten/model_solve_d.hpp>
#include <stan/math/torsten/PKModel/functors/check_mti.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>

namespace torsten {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using Eigen::Map;

  struct PMXLinODE {
    /**
     * standard two compartment PK ODE RHS function
     * @tparam T0 t type
     * @tparam T1 initial condition type
     * @tparam T2 parameter type
     * @tparam T3 real data/rate type
     * @param t type
     * @param x initial condition type
     * @param parms parameters that form a n x n matrix, where
     * @c (n = x.size()). We assume the matrix is col-majored.
     * @param rate dosing rate
     * @param dummy dummy
     */
    template <typename T0, typename T1, typename T2>
    inline
    std::vector<typename stan::return_type<T0, T1, T2>::type>
    operator()(const T0& t,
               const std::vector<T1>& x,
               const std::vector<T2>& parms,
               const std::vector<double>& rate,
               const std::vector<int>& dummy,
               std::ostream* pstream__) const {
      typedef typename stan::return_type<T0, T1, T2>::type scalar;
      
      size_t n = x.size();
      std::vector<scalar> res(n);
      Matrix<scalar, Dynamic, Dynamic> m(n, n);
      Matrix<scalar, Dynamic, 1> v(n);
      for (size_t i = 0; i < n * n; ++i) m(i) = parms[i];
      for (size_t i = 0; i < n; ++i) v(i) = x[i];
      torsten::PMXLin<scalar> rv = stan::math::multiply(m, v);
      for (size_t i = 0; i < n; ++i) res[i] = rv(i);

      return res;
    }

    /** 
     * Eigen::Matrix version
     * 
     */
    template <typename T0, typename T1, typename T2>
    Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
    operator()(const T0& t,
               const Eigen::Matrix<T1, -1, 1>& x,
               std::ostream* pstream__,
               const std::vector<T2>& parms,
               const std::vector<double>& x_r,
               const std::vector<int>& x_i) const {
      typedef typename stan::return_type_t<T0, T1, T2> scalar;
      
      size_t n = x.size();
      Eigen::Matrix<scalar, -1, 1> res(n);
      Matrix<scalar, Dynamic, Dynamic> m(n, n);
      Matrix<scalar, Dynamic, 1> v(n);
      for (size_t i = 0; i < n * n; ++i) m(i) = parms[i];
      for (size_t i = 0; i < n; ++i) v(i) = x[i];
      torsten::PMXLin<scalar> rv = stan::math::multiply(m, v);
      for (size_t i = 0; i < n; ++i) res[i] = rv(i);

      return res;
    }
  };

  /**
   * Linear ODE model.
   *
   * @tparam T_time type of time
   * @tparam T_rate type of dosing rate
   * @tparam T_par  type of parameters.
   */
  template<typename T_par>
  class PMXLinODEModel {
  protected:
    const Eigen::Matrix<T_par, -1, -1> & par_;
    const int ncmt_;

  public:
    static constexpr PMXLinODE f_ = PMXLinODE();
    using par_type    = T_par;

    PMXLinODEModel(const Eigen::Matrix<T_par, -1, -1>& par,
                   int ncmt) :
      par_(par), ncmt_(ncmt)
    {}

    PMXLinODEModel(const Eigen::Matrix<T_par, -1, -1>& par) :
      par_(par), ncmt_(par.rows())
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

    const Eigen::Matrix<T_par,-1,-1>  & par ()  const { return par_; }
    const PMXLinODE            & f()     const { return f_; }

    const int ncmt () const {
      return ncmt_;
    }

    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::matrix_exp;
      using stan::math::mdivide_left;
      using stan::math::multiply;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;

      const int nCmt = par_.cols();
      Eigen::Matrix<T, -1, 1> y0t(nCmt);
      for (int i = 0; i < nCmt; ++i) y0t(i) = y(i);

      Eigen::Matrix<T, -1, 1> pred(nCmt);

      if (std::any_of(rate.begin(), rate.end(),
                      [](T1 r){return r != 0;})) {
        Eigen::Matrix<T, -1, 1> rate_vec(rate.size()), x(nCmt), x2(nCmt);
        for (size_t i = 0; i < rate.size(); ++i) rate_vec(i) = rate[i];
        x = mdivide_left(par_, rate_vec);
        x2 = x + y;
        Eigen::Matrix<T, -1, -1> dt_system = multiply(dt, par_);
        pred = matrix_exp(dt_system) * x2;
        pred -= x;
      } else {
        Eigen::Matrix<T, -1, -1> dt_system = multiply(dt, par_);
        pred = matrix_exp(dt_system) * y0t;
      }
      y = pred;
    }

    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate) const {
      const dsolve::PMXAnalyiticalIntegrator integ;
      solve(y, t0, t1, rate, integ);
    }

    /*
     * solve the linear ODE: steady state version
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_ii, T_par>, -1, 1> 
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii,
          const int& cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using stan::math::matrix_exp;
      using stan::math::mdivide_left;
      using stan::math::multiply;
      using stan::math::value_of;
      using std::vector;

      using T0 = typename stan::return_type_t<T_ii, T_par>;
      using scalar = typename stan::return_type_t<T_amt, T_r, T_ii, T_par>; //NOLINT

      int nCmt = ncmt();
      Matrix<T0, Dynamic, Dynamic> workMatrix;
      Matrix<T0, Dynamic, Dynamic> ii_system = multiply(ii, par_);
      Matrix<scalar, 1, Dynamic> pred(nCmt);
      pred.setZero();
      Matrix<scalar, Dynamic, 1> amounts(nCmt);
      amounts.setZero();

      if (rate == 0) {  // bolus dose
        amounts(cmt - 1) = amt;
        workMatrix = - matrix_exp(ii_system);
        for (int i = 0; i < nCmt; i++) workMatrix(i, i) += 1;
        amounts = mdivide_left(workMatrix, amounts);
        pred = multiply(matrix_exp(ii_system), amounts);

      } else if (ii > 0) {  // multiple truncated infusions
        /**
         * Each II has two sections due to overlapping infusion time,
         * e.g. when II < infusion time < 2 x II the first part of II
         * infusion has dosing rate <code>2 * rate</code>  while the rest of II has
         * dosing rate <code>rate</code>.
         *
         * Let A be linODE matrix, B the dosing rate at 1st section of
          * II, C the dosing rate at 2nd section of II, and y the
         * initial condition(steady state solution), then the ODE
         * solution at the end of 1st section is
         *
         * Exp(A * t1) * (y + A^{-1} * B) - A^{-1} * B
         *
         * where t1 is the length of 1st section(so that the length of
         * 2nd section t2 = II - t1). The final solution at the end of
         * 2nd section is obtained by using the above expression as new
         * y and apply t2 & C to it.
         *
         */

        typename stan::return_type_t<T_amt, T_r> dt = amt / rate;
        amounts(cmt - 1) = rate;

        int n = int(std::floor(value_of(dt) / value_of(ii)) + 0.1);
        typename stan::return_type_t<T_amt, T_r, T_ii> t1 = dt - n * ii;
        typename stan::return_type_t<T_amt, T_r, T_ii> t2 = (n + 1) * ii - dt;
        PKRec<T_r> rate1 = PKRec<T_r>::Zero(nCmt);
        PKRec<T_r> rate2 = PKRec<T_r>::Zero(nCmt);
        rate1(cmt - 1) = (n + 1) * rate;
        rate2(cmt - 1) = n * rate;
        PKRec<typename stan::return_type_t<T_r, T_par>> par_rate1 = mdivide_left(par_, rate1);
        PKRec<typename stan::return_type_t<T_r, T_par>> par_rate2 = mdivide_left(par_, rate2);
        Eigen::Matrix<scalar, -1, -1> exp_par_t1(matrix_exp(multiply(t1, par_)));
        Eigen::Matrix<scalar, -1, -1> exp_par_t2(matrix_exp(multiply(t2, par_)));
        Eigen::Matrix<scalar, -1, -1> lhs =
          multiply(exp_par_t1, exp_par_t2) - Eigen::Matrix<scalar, -1, -1>::Identity(nCmt, nCmt);
        PKRec<scalar> rhs = par_rate2 + multiply(exp_par_t2, par_rate1 - par_rate2) -
          multiply(exp_par_t2, multiply(exp_par_t1, par_rate1));
        pred = mdivide_left(lhs, rhs);
      } else {  // constant infusion
        amounts(cmt - 1) -= rate;
        pred = mdivide_left(par_, amounts);
      }
      return pred;
    }

    /*
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_ii, T_par>, -1, 1>
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const dsolve::PMXAnalyiticalIntegrator& integrator) const {
      return solve(t0, amt, rate, ii, cmt);
    }
  };

  template<typename T_par>
  constexpr PMXLinODE PMXLinODEModel<T_par>::f_;

  template<typename T>
  using LinOdeEigenDecomp = std::tuple<const Eigen::Matrix<T, -1, -1>&,
                                       const Eigen::Matrix<T, -1, 1>&,
                                       const Eigen::Matrix<T, -1, -1>& >;

  /**
   * linear ode with eigen decomposition
   * 
   */
  template<typename T_par>
  class PMXLinOdeEigenDecompModel {
    const Eigen::Matrix<T_par, -1, -1> & p_;
    const Eigen::Matrix<T_par, -1, 1>& diag_;
    const Eigen::Matrix<T_par, -1, -1> & p_inv_;
    const int ncmt_;

  public:
    PMXLinOdeEigenDecompModel(const Eigen::Matrix<T_par, -1, -1>& p,
                              const Eigen::Matrix<T_par, -1, 1>& diag,
                              const Eigen::Matrix<T_par, -1, -1>& p_inv) :
      p_(p), diag_(diag), p_inv_(p_inv), ncmt_(diag_.size())
    {}

    PMXLinOdeEigenDecompModel(const LinOdeEigenDecomp<T_par>& pdp) :
      PMXLinOdeEigenDecompModel(std::get<0>(pdp), std::get<1>(pdp), std::get<2>(pdp))
    {}

    template<typename Tt0, typename Tt1, typename T, typename T1>
    void solve(PKRec<T>& y,
               const Tt0& t0, const Tt1& t1,
               const std::vector<T1>& rate,
               const dsolve::PMXAnalyiticalIntegrator& integ) const {
      using stan::math::multiply;

      typename stan::return_type_t<Tt0, Tt1> dt = t1 - t0;

      Eigen::Matrix<T_par, -1, -1> work(ncmt_, ncmt_);
      work.setZero();
      for (int i = 0; i < ncmt_; ++i) {
        work(i, i) = 1.0 / diag_(i);
      }
      work = multiply(work, p_inv_);
      PKRec<T1> r = stan::math::to_vector(rate);
      PKRec<stan::return_type_t<T1, T_par> > p_inv_r = multiply(work, r);
      PKRec<T> pred = multiply(p_inv_, y) + p_inv_r;

      Eigen::Matrix<stan::return_type_t<T_par, Tt1, Tt0>, -1, -1> work1(ncmt_,ncmt_);
      work1.setZero();
      for (int i = 0; i < ncmt_; ++i) {
        work1(i, i) = stan::math::exp(diag_(i) * dt);
      }

      y = multiply(multiply(p_, work1), pred) - multiply(p_, p_inv_r);
    }

    /*
     * solve the linear ODE: steady state version
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<typename stan::return_type_t<T_amt, T_r, T_ii, T_par>, -1, 1> 
    solve(double t0, const T_amt& amt, const T_r& rate, const T_ii& ii,
          const int& cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      // using stan::math::matrix_exp;
      using stan::math::mdivide_left;
      using stan::math::multiply;
      using stan::math::value_of;
      using std::vector;

      using T0 = typename stan::return_type_t<T_ii, T_par>;
      using scalar = typename stan::return_type_t<T_amt, T_r, T_ii, T_par>; //NOLINT

      Matrix<scalar, -1, 1> pred(ncmt_);
      pred.setZero();
      Matrix<scalar, Dynamic, 1> amounts(ncmt_);
      amounts.setZero();

      if (rate == 0) {  // bolus dose
        /**
         * exp(At)(u + bolus) = u
         *
         * thus (I - exp(At))*u = exp(At)*bolus
         *
         */
        Matrix<T0, -1, -1> diag = Matrix<T0, -1, -1>::Zero(ncmt_, ncmt_);
        for (int i = 0; i < ncmt_; ++i) {
          diag(i, i) = stan::math::exp(ii * diag_(i));
        }
        PKRec<T_amt> bolus = PKRec<T_amt>::Zero(ncmt_);
        bolus(cmt - 1) = amt;
        pred = multiply(multiply(diag, p_inv_), bolus);
        for (int i = 0; i < ncmt_; ++i) {
          pred(i) /= (1.0 - diag(i, i));
        }
        pred = multiply(p_, pred);
      } else if (ii > 0) {  // multiple truncated infusions
        /**
         * with eigen-decomp A= P * diag * P_inv, change of variable: 
         * y = P_inv * u decouples original ODE based on u.
         */
        typename stan::return_type_t<T_amt, T_r> dt = amt / rate;
        int n = int(std::floor(value_of(dt) / value_of(ii)) + 0.1);
        typename stan::return_type_t<T_amt, T_r, T_ii> t1 = dt - n * ii;
        typename stan::return_type_t<T_amt, T_r, T_ii> t2 = (n + 1) * ii - dt;

        Eigen::Matrix<T_par, -1, -1> a(ncmt_, ncmt_);
        a.setZero();
        for (int i = 0; i < ncmt_; ++i) {
          a(i, i) = 1.0 / diag_(i);
        }
        // a = multiply(multiply(multiply(p_, a), p_inv_), p_inv_);
        a = multiply(a, p_inv_);
        PKRec<stan::return_type_t<T_par, T_r>> y1(ncmt_), y2(ncmt_);
        y1.setZero();
        y2.setZero();
        y1(cmt - 1) = (n + 1) * rate;
        y2(cmt - 1) = n * rate;
        y1 = multiply(a, y1);
        y2 = multiply(a, y2);
        Eigen::Matrix<scalar, -1, -1> work = Eigen::Matrix<scalar, -1, -1>::Zero(ncmt_, ncmt_);
        for (int i = 0; i < ncmt_; ++i) {
          work(i, i) = stan::math::exp(diag_(i) * t2);
        }
        pred = y2;
        y2 = y1 - y2;
        pred += multiply(work, y2);
        for (int i = 0; i < ncmt_; ++i) {
          work(i, i) = stan::math::exp(diag_(i) * (t1 + t2));
        }
        pred -= multiply(work, y1);
        work -= Eigen::Matrix<scalar, -1, -1>::Identity(ncmt_, ncmt_);
        for (int i = 0; i < ncmt_; ++i) {
          work(i, i) = 1.0 / work(i, i);
        }
        pred = multiply(p_, multiply(work, pred));
      } else {  // constant infusion
        amounts(cmt - 1) -= rate;
        PKRec<T_r> rvec = PKRec<T_r>::Zero(ncmt_);
        rvec(cmt - 1) -= rate;
        pred = multiply(p_inv_, rvec);
        Eigen::Matrix<T_par, -1, -1> work(ncmt_, ncmt_);
        work.setZero();
        for (int i = 0; i < ncmt_; ++i) {
          work(i, i) = 1.0 / diag_(i);
        }
        pred = multiply(p_, multiply(work, pred));
      }
      return pred;
    }
  };
}

#endif
