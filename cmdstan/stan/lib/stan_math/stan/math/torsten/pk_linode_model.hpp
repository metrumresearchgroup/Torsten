#ifndef PK_LINODE_MODEL_HPP
#define PK_LINODE_MODEL_HPP

#include <stan/math/torsten/torsten_def.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using Eigen::Map;

  struct PKLinODE {
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
      PKLin<scalar> rv = m * v;
      for (size_t i = 0; i < n; ++i) res[i] = rv(i);

      return res;
    }
  };

  /**
   * Linear ODE model.
   *
   * @tparam T_time type of time
   * @tparam T_init type of init condition
   * @tparam T_rate type of dosing rate
   * @tparam T_par  type of parameters.
   */
  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  class PKLinODEModel {
    const T_time &t0_;
    const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0_;
    const std::vector<T_rate> &rate_;
    const std::vector<T_par> & par_;

  public:
    static constexpr PKLinODE f_ = PKLinODE();
    using scalar_type = typename
      stan::return_type<T_time, T_init, T_rate, T_par>::type;

    /*
     * FIX ME: we need to get rid of @c ModelParameters in
     * @c Pred2
     */
    template<typename T0, typename T1, typename T2, typename T3>
    static std::vector<T1> 
    get_param(const torsten::ModelParameters<T0, T1, T2, T3>& p) {
      auto k = p.get_K();
      std::vector<T1> res(k.size());
      for (size_t i = 0; i < res.size(); ++i) {
        res[i] = k(i);
      }
      return res;
    }
  /**
   * Constructor
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
    // template<template<typename...> class T_mp, typename... Ts>
    // PKLinODEModel(const T_time& t0,
    //               const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
    //               const std::vector<T_rate> &rate,
    //               const std::vector<T_par> & par,
    //               const T_mp<Ts...> &parameter) :
    //   t0_(t0),
    //   y0_(y0),
    //   rate_(rate),
    // {}

    PKLinODEModel(const T_time& t0,
                  const Eigen::Matrix<T_init, 1, Eigen::Dynamic>& y0,
                  const std::vector<T_rate> &rate,
                  const std::vector<T_par>& par) :
      t0_(t0),
      y0_(y0),
      rate_(rate),
      par_(par)
    {}

    const T_time              & t0()    const { return t0_; }
    const PKRec<T_init>       & y0()    const { return y0_; }
    const std::vector<T_rate> & rate()  const { return rate_; }
    const std::vector<T_par>  & par ()  const { return par_; }
    const PKLinODE            & f()     const { return f_; }

    const int ncmt () const {
      return y0_.size();
    }
    const PKLin<T_par> coef () const {
      return Eigen::Map<const PKLin<T_par> >(par_.data(), ncmt(), ncmt());
    }

    /*
     * solve linear ODE model using matrix exponential function
     */
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& dt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using stan::math::value_of;
      using stan::math::matrix_exp;
      using stan::math::mdivide_left;
      using stan::math::multiply;
      using stan::math::scale_matrix_exp_multiply;

      auto system = coef();

      const int nCmt = system.cols();
      Matrix<scalar_type, Dynamic, 1> y0t(nCmt);
      for (int i = 0; i < nCmt; ++i) y0t(i) = y0_(i);

      if (std::any_of(rate_.begin(), rate_.end(),
                      [](T_rate r){return r != 0;})) {
        Matrix<scalar_type, Dynamic, 1> rate_vec(rate_.size()), x(nCmt), x2(nCmt);
        for (size_t i = 0; i < rate_.size(); i++) rate_vec(i) = rate_[i];
        x = mdivide_left(system, rate_vec);
        x2 = x + y0_.transpose();
        Matrix<scalar_type, Dynamic, Dynamic> dt_system = multiply(dt, system);
        Matrix<scalar_type, Dynamic, 1> pred = matrix_exp(dt_system) * x2;
        pred -= x;
        return pred.transpose();
      } else {
        // return scale_matrix_exp_multiply(value_of(dt), system, y0t);
        Matrix<scalar_type, Dynamic, Dynamic> dt_system = multiply(dt, system);
        Matrix<scalar_type, Dynamic, 1> pred = matrix_exp(dt_system) * y0t;
        return pred.transpose();
      }
    }

    /*
     * solve the linear ODE: steady state version
     */
    template<typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1> 
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii,
          const int& cmt) const {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using stan::math::matrix_exp;
      using stan::math::mdivide_left;
      using stan::math::multiply;
      using std::vector;

      typedef typename promote_args<T_ii, T_par>::type T0;
      typedef typename promote_args<T_amt, T_r, T_ii, T_par>::type scalar; //NOLINT

      Matrix<T_par, Dynamic, Dynamic> system = coef();
      int nCmt = ncmt();
      Matrix<T0, Dynamic, Dynamic> workMatrix;
      Matrix<T0, Dynamic, Dynamic> ii_system = multiply(ii, system);
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
        scalar delta = amt / rate;
        static const char* function("Steady State Event");
        torsten::check_mti(amt, delta, ii, function);

        amounts(cmt - 1) = rate;
        scalar t = delta;
        amounts = mdivide_left(system, amounts);
        Matrix<scalar, Dynamic, Dynamic> t_system = multiply(delta, system);
        pred = matrix_exp(t_system) * amounts;
        pred -= amounts;

        workMatrix = - matrix_exp(ii_system);
        for (int i = 0; i < nCmt; i++) workMatrix(i, i) += 1;

        Matrix<scalar, Dynamic, 1> pred_t = pred.transpose();
        pred_t = mdivide_left(workMatrix, pred_t);
        t = ii - t;
        t_system = multiply(t, system);
        pred_t = matrix_exp(t_system) * pred_t;
        pred = pred_t.transpose();

      } else {  // constant infusion
        amounts(cmt - 1) -= rate;
        pred = mdivide_left(system, amounts);
      }
      return pred;
    }

    /*
     * wrapper to fit @c PrepWrapper's call signature
     */
    template<PkOdeIntegratorId It>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_time& dt,
          const PkOdeIntegrator<It>& integrator) const {
      return solve(dt);
    }

    template<PkOdeIntegratorId It, typename T_amt, typename T_r, typename T_ii>
    Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>
    solve(const T_amt& amt, const T_r& rate, const T_ii& ii, const int& cmt,
          const PkOdeIntegrator<It>& integrator) const {
      return solve(amt, rate, ii, cmt);
    }
  };

  template<typename T_time, typename T_init, typename T_rate, typename T_par>
  constexpr PKLinODE PKLinODEModel<T_time, T_init, T_rate, T_par>::f_;

}

#endif
