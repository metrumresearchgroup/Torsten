#ifndef STAN_MATH_TORSTEN_MODEL_SOLVE_D_HPP
#define STAN_MATH_TORSTEN_MODEL_SOLVE_D_HPP

#include <stan/math/torsten/torsten_def.hpp>
#include <stan/math/torsten/val_and_grad_nested.hpp>
#include <stan/math/torsten/pmx_ode_integrator.hpp>

namespace torsten {
  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For transient solution, when the return type is already
   * data, we pass on to model's solver.
   */
  template<typename T_model,
           typename std::enable_if_t<!stan::is_var<typename T_model::scalar_type>::value >* = nullptr> //NOLINT
  Eigen::VectorXd model_solve_d(const T_model& pkmodel, const double& t_next) {
    return pkmodel.solve(t_next);
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For transient solution with prediction parameter in
   * solve call, such as an ODE integrator. When the return type is already
   * data, we pass on to model's solver.
   */
  template<typename T_model, typename T_pred,
           typename std::enable_if_t<!stan::is_var<typename T_model::scalar_type>::value >* = nullptr> //NOLINT
  Eigen::VectorXd model_solve_d(const T_model& pkmodel, const double& t_next, const T_pred& pred_par) {
    return pkmodel.solve(t_next, pred_par);
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For transient solution. When the solution is @c var
   * type, we take gradient using autodiff.
   */
  template<typename T_model>
  Eigen::VectorXd model_solve_d(const T_model& pkmodel,
                                typename T_model::time_type const& t_next) {
    using std::vector;
    using Eigen::VectorXd;
    using Eigen::Matrix;
    using stan::math::value_of;
    using stan::math::var;

    using T_time = typename T_model::time_type; 
    using T_init = typename T_model::init_type;
    using T_rate = typename T_model::rate_type;
    using T_par = typename T_model::par_type;

    const Matrix<T_init, 1, -1>& y0 = pkmodel.y0();
    const vector<T_rate>& rate = pkmodel.rate();
    const vector<T_par>& par = pkmodel.par();      

    VectorXd res_d;

    stan::math::start_nested();

    try {
      Matrix<T_init, 1, -1> y0_new(y0.size());
      vector<T_rate> rate_new(rate.size());
      vector<T_par> par_new(par.size());

      for (int i = 0; i < y0_new.size(); ++i) {y0_new(i) = value_of(y0(i));}
      for (size_t i = 0; i < rate_new.size(); ++i) {rate_new[i] = value_of(rate[i]);}
      for (size_t i = 0; i < par_new.size(); ++i) {par_new[i] = value_of(par[i]);}

      T_time t0 = value_of(pkmodel.t0());
      T_time t1 = value_of(t_next);
      T_model pkmodel_new(t0, y0_new, rate_new, par_new);

      auto res = pkmodel_new.solve(t1);
      vector<var> var_new(pkmodel_new.vars(t1));
      res_d = val_and_grad_nested(res, var_new);
    } catch (const std::exception& e) {
      stan::math::recover_memory_nested();
      throw;
    }
    stan::math::recover_memory_nested();

    return res_d;
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For transient solution of models that requires an
   * integrator to solve. When the solution is @c var
   * type, we take gradient using autodiff.
   */
  template<typename T_model, PMXOdeIntegratorId It,
           typename std::enable_if_t<stan::is_var<typename T_model::scalar_type>::value >* = nullptr> //NOLINT
  Eigen::VectorXd model_solve_d(const T_model& pkmodel,
                                typename T_model::time_type const& t_next,
                                PMXOdeIntegrator<It> const& integrator) {
    using std::vector;
    using Eigen::VectorXd;
    using Eigen::Matrix;
    using stan::math::value_of;
    using stan::math::var;

    using T_time = typename T_model::time_type; 
    using T_init = typename T_model::init_type;
    using T_rate = typename T_model::rate_type;
    using T_par = typename T_model::par_type;

    const Matrix<T_init, 1, -1>& y0 = pkmodel.y0();
    const vector<T_rate>& rate = pkmodel.rate();
    const vector<T_par>& par = pkmodel.par();      

    VectorXd res_d;

    stan::math::start_nested();

    try {
      Matrix<T_init, 1, -1> y0_new(y0.size());
      vector<T_rate> rate_new(rate.size());
      vector<T_par> par_new(par.size());

      for (int i = 0; i < y0_new.size(); ++i) {y0_new(i) = value_of(y0(i));}
      for (size_t i = 0; i < rate_new.size(); ++i) {rate_new[i] = value_of(rate[i]);}
      for (size_t i = 0; i < par_new.size(); ++i) {par_new[i] = value_of(par[i]);}

      T_time t0 = value_of(pkmodel.t0());
      T_time t1 = value_of(t_next);
      T_model pkmodel_new(t0, y0_new, rate_new, par_new, pkmodel.f());

      auto res = pkmodel_new.solve(t1, integrator);
      vector<var> var_new(pkmodel_new.vars(t1));
      res_d = val_and_grad_nested(res, var_new);
    } catch (const std::exception& e) {
      stan::math::recover_memory_nested();
      throw;
    }
    stan::math::recover_memory_nested();

    return res_d;
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For steady-state solution. when the return type is already
   * data, we pass on to model's solver.
   */
  template<typename T_model,
           typename std::enable_if_t<!stan::is_var<typename T_model::par_type>::value >* = nullptr>
    Eigen::VectorXd model_solve_d(const T_model& pkmodel, const double& amt, const double& rate, const double& ii, const int& cmt) { // NOLINT
      return pkmodel.solve(amt, rate, ii, cmt);
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For steady-state solution with prediction parameter in
   * solve call, such as an ODE integrator. when the return type is already
   * data, we pass on to model's solver.
   */
  template<typename T_model, typename T_pred,
           typename std::enable_if_t<!stan::is_var<typename T_model::par_type>::value >* = nullptr>
  Eigen::VectorXd model_solve_d(const T_model& pkmodel, const double& amt, const double& rate, const double& ii, const int& cmt, // NOLINT
                                const T_pred& pred_par) {
    return pkmodel.solve(amt, rate, ii, cmt, pred_par);
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For steady-state solution. When the solution is @c var
   * type, we take gradient using autodiff.
   */
  template<typename T_model, typename T_amt, typename T_rate, typename T_ii>
  Eigen::VectorXd model_solve_d(const T_model& pkmodel, const T_amt& amt, const T_rate& r, const T_ii& ii, const int& cmt) { // NOLINT
    using std::vector;
    using Eigen::VectorXd;
    using Eigen::Matrix;
    using stan::math::value_of;
    using stan::math::var;

    VectorXd res_d;

    using T_par = typename T_model::par_type;
    const std::vector<T_par>& par(pkmodel.par());

    stan::math::start_nested();
    try {
      std::vector<T_par> par_new(par.size());
      for (size_t i = 0; i < par.size(); ++i) par_new[i] = value_of(par[i]);

      T_model pkmodel_new(pkmodel.t0(), pkmodel.y0(), pkmodel.rate(), par_new);

      T_amt amt_new = value_of(amt);
      T_rate r_new = value_of(r);
      T_ii ii_new = value_of(ii);
      auto res = pkmodel_new.solve(amt_new, r_new, ii_new, cmt);
      vector<var> var_new(pkmodel_new.vars(amt_new, r_new, ii_new));
      res_d = val_and_grad_nested(res, var_new);
    } catch (const std::exception& e) {
      stan::math::recover_memory_nested();
      throw;
    }
    stan::math::recover_memory_nested();

    return res_d;
  }

  /*
   * Solve PK models and return the results in form of data,
   * arrange as
   *
   *    sol value y1, dy1/dp1, dy1/dp2..., sol value y2, dy2/dp2, dy2/dp2...
   *
   * For steady-state solution in a model that requires an
   * integrator to solve. When the solution is @c var
   * type, we take gradient using autodiff.
   */
  template<typename T_model, typename T_amt, typename T_ii, PMXOdeIntegratorId It,
           typename std::enable_if_t<stan::is_var<typename stan::return_type<typename T_model::par_type, T_amt, T_ii>::type>::value >* = nullptr> // NOLINT
  Eigen::VectorXd model_solve_d(const T_model& pkmodel, const T_amt& amt, const double& r, const T_ii& ii, const int& cmt, // NOLINT
                                PMXOdeIntegrator<It> const& integrator) {
      using std::vector;
      using Eigen::VectorXd;
      using Eigen::Matrix;
      using stan::math::value_of;
      using stan::math::var;

      VectorXd res_d;

      using T_par = typename T_model::par_type;
      const std::vector<T_par>& par(pkmodel.par());

      stan::math::start_nested();
      try {
        std::vector<T_par> par_new(par.size());
        for (size_t i = 0; i < par.size(); ++i) par_new[i] = value_of(par[i]);

        T_model pkmodel_new(pkmodel.t0(), pkmodel.y0(), pkmodel.rate(), par_new, pkmodel.f());

        T_amt amt_new = value_of(amt);
        double r_new = r;
        T_ii ii_new = value_of(ii);
        auto res = pkmodel_new.solve(amt_new, r_new, ii_new, cmt, integrator);
        vector<var> var_new(pkmodel_new.vars(amt_new, r_new, ii_new));
        res_d = val_and_grad_nested(res, var_new);
      } catch (const std::exception& e) {
        stan::math::recover_memory_nested();
        throw;
      }
      stan::math::recover_memory_nested();

      return res_d;
    }
}
#endif
