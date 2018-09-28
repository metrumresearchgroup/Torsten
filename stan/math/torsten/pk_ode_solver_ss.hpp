#ifndef STAN_MATH_TORSTEN_REFACTOR_ODE_SOLVER_SS_HPP
#define STAN_MATH_TORSTEN_REFACTOR_ODE_SOLVER_SS_HPP

#include <stan/math/torsten/PKModel/functors/general_functor.hpp>
#include <stan/math/torsten/pk_ss_system.hpp>

namespace refactor {

  using boost::math::tools::promote_args;
  using namespace torsten;

  /**
   * Steady state ODE solver
   *
   */
  class PKODEModelSolverSS {
    torsten::TorstenIntegrator integrator_;
  public:
  /**
   * Steady state ODE solver constructor
   *
   */
    PKODEModelSolverSS(const double& rel_tol,
                       const double& abs_tol,
                       const long int& max_num_steps,
                       std::ostream* msgs,
                       const int& integratorType) :
      integrator_(rel_tol, abs_tol, max_num_steps, msgs, integratorType)
    {}

  /**
   * Solve steady state ODE-based PKPD model.
   *
   * @tparam T_ii dosing interval type
   * @tparam T_model ODE model type
   * @tparam Ts_par type of parameters
   * @param pkmodel Linear ODE model
   * @param amt dosing amount
   * @param rate dosing rate
   * @param ii dosing interval
   * @param cmt compartment where the dosing occurs
   * @return col vector of the steady state ODE solution
   */
    template<typename T_ii,
             template <class, class... > class T_model, class... Ts_par>
    Eigen::Matrix<typename promote_args<T_ii,
                                        typename
                                        T_model<Ts_par...>::par_type>::type, 
                  Eigen::Dynamic, 1>
    solve(const T_model<Ts_par...> &pkmodel,
          const double& amt,
          const double& rate,
          const T_ii& ii,
          const int& cmt) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using Eigen::VectorXd;
    using std::vector;
    using stan::math::algebra_solver;
    using stan::math::to_vector;

    typedef typename promote_args<T_ii,
                                  typename
                                  T_model<Ts_par...>::par_type>::type scalar;

    using F = typename T_model<Ts_par...>::f_type;
    const F& f = pkmodel.rhs_fun();
    const int& ncmt = pkmodel.ncmt();

    Matrix<scalar, Dynamic, 1> pred;

    // Arguments for ODE integrator (and initial guess)
    double ii_dbl = unpromote(ii);
    Matrix<double, 1, Dynamic> init_dbl(ncmt);
    for (int i = 0; i < ncmt; i++) init_dbl(i) = 0;
    vector<double> x_r(ncmt, 0);
    vector<int> x_i(0);

    // Arguments for algebraic solver
    Matrix<double, 1, Dynamic> y;
    double rel_tol = 1e-10;  // default
    double f_tol = 1e-4;  // empirical
    long int max_num_steps = 1e3;  // default // NOLINT

    // construct algebraic function
    general_functor<F> f1(f);
    using SS_functor = ode_rate_dbl_functor<general_functor<F>>;
    SteadyStateSys_dd<SS_functor, void>
      system(SS_functor(f1), 
             ii_dbl, cmt, integrator_);

    refactor::PKODEModelSolver sol(integrator_);

    if (rate == 0) {  // bolus dose
      // compute initial guess
      init_dbl(cmt - 1) = amt;
      y = sol.solve(general_functor<F>(f),
                    stan::math::value_of(pkmodel.t0()),
                    init_dbl,
                    ii_dbl,
                    stan::math::value_of(pkmodel.par()),
                    x_r);
      x_r.push_back(amt);
      pred = algebra_solver(system, y,
                            to_vector(pkmodel.par()),
                            x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
      // DEV - what tuning parameters should we use for the algebra solver?
      // DEV - update initial guess or tuning parameters if result not good?
    }  else if (ii > 0) {  // multiple truncated infusions
      x_r[cmt - 1] = rate;
      // compute initial guess
      y = sol.solve(general_functor<F>(f),
                    stan::math::value_of(pkmodel.t0()),
                    init_dbl,
                    ii_dbl,
                    stan::math::value_of(pkmodel.par()),
                    x_r);
      x_r.push_back(amt);
      pred = algebra_solver(system, y,
                            to_vector(pkmodel.par()),
                            x_r, x_i,
                            0, rel_tol, 1e-3, max_num_steps);  // FIX ME
                                                               // use ftol
    } else {  // constant infusion
      x_r[cmt - 1] = rate;
      y = sol.solve(general_functor<F>(f),
                    stan::math::value_of(pkmodel.t0()),
                    init_dbl,
                    100.0,
                    stan::math::value_of(pkmodel.par()),
                    x_r);

      x_r.push_back(amt);
      pred = algebra_solver(system, y,
                            to_vector(pkmodel.par()),
                            x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    }

    return pred;
  }

    template<// typename F,
           typename T_amt,
           typename T_ii,
           template <class, class... > class T_model, class... Ts_par>
  Eigen::Matrix<typename
  boost::math::tools::promote_args<T_ii,
                                   T_amt,
                                   typename T_model<Ts_par...>::par_type>
                ::type, Eigen::Dynamic, 1>
    solve(T_model<Ts_par...> pkmodel,
          const T_amt& amt,
          const double& rate,
          const T_ii& ii,
          const int& cmt) const
{
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using Eigen::VectorXd;
    using std::vector;
    using stan::math::algebra_solver;
    using stan::math::to_vector;
    using stan::math::invalid_argument;

    typedef typename promote_args<T_ii,
                                  T_amt,
                                  typename
                                  T_model<Ts_par...>::time_type,
                                  typename
                                  T_model<Ts_par...>::par_type>::type scalar;

    using par_type = typename T_model<Ts_par...>::par_type;
    using F = typename T_model<Ts_par...>::f_type;
    const F& f = pkmodel.rhs_fun();
    const int& ncmt = pkmodel.ncmt();
    const std::vector<par_type>& pars = pkmodel.par();

    Matrix<scalar, Dynamic, 1> pred;

    // Arguments for the ODE integrator
    double ii_dbl = unpromote(ii);
    Matrix<double, 1, Dynamic> init_dbl(ncmt);
    for (int i = 0; i < ncmt; i++) init_dbl(i) = 0;
    vector<double> x_r(ncmt, 0);
    vector<int> x_i(0);

    // Arguments for algebraic solver
    Matrix<double, 1, Dynamic> y;
    double rel_tol = 1e-10;  // default
    double f_tol = 5e-4;  // empirical (note: differs from other function)
    long int max_num_steps = 1e4;  // default  // NOLINT

    // construct algebraic function
    general_functor<F> f1(f);
    using SS_functor = ode_rate_dbl_functor<general_functor<F>>;
    SteadyStateSys_vd<SS_functor>
      system(SS_functor(f1), ii_dbl, cmt, integrator_);

    refactor::PKODEModelSolver sol(integrator_);

    int npar = pars.size();
    Matrix<scalar, Dynamic, 1> parms(npar + 1);
    for (int i = 0; i < npar; i++) parms(i) = pars[i];
    parms(npar) = amt;

    if (rate == 0) {  // bolus dose
      // compute initial guess
      init_dbl(cmt - 1) = unpromote(amt);
      y = sol.solve(general_functor<F>(f),
                    stan::math::value_of(pkmodel.t0()),
                    init_dbl,
                    ii_dbl,
                    stan::math::value_of(pkmodel.par()),
                    x_r);

      pred = algebra_solver(system, y, parms, x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    }  else if (ii > 0) {  // multiple truncated infusions
      // compute initial guess
      x_r[cmt - 1] = rate;
      y = sol.solve(general_functor<F>(f),
                    stan::math::value_of(pkmodel.t0()),
                    init_dbl,
                    ii_dbl,
                    stan::math::value_of(pkmodel.par()),
                    x_r);

      pred = algebra_solver(system, y, parms, x_r, x_i,
                            0, rel_tol, 1e-3, max_num_steps);  // use ftol
    } else {  // constant infusion
      x_r[cmt - 1] = rate;
      y = sol.solve(general_functor<F>(f),
                    stan::math::value_of(pkmodel.t0()),
                    init_dbl,
                    100.0,
                    stan::math::value_of(pkmodel.par()),
                    x_r);

      pred = algebra_solver(system, y, parms, x_r, x_i,
                            0, rel_tol, f_tol, max_num_steps);
    }
    return pred;
  }
  };

}



#endif
