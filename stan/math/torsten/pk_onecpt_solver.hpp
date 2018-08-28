#ifndef STAN_MATH_TORSTEN_ONECPT_SOLVER_HPP
#define STAN_MATH_TORSTEN_ONECPT_SOLVER_HPP

#include <stan/math/torsten/pk_onecpt_model.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>

namespace refactor {

  using boost::math::tools::promote_args;

  /**
   * standard one compartment PK ODE solver based on
   * analytical solution.
   */
  class PKOneCptModelSolver {
  public:
    static constexpr int Ncmt =
      PKOneCptModel<double, double, double, double>::Ncmt;
    static constexpr int Npar =
      PKOneCptModel<double, double, double, double>::Npar;

    PKOneCptModelSolver() {}

  /**
   * standard one compartment PK model attached to this solver.
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate rate type
   * @tparam T_par parameter type
   */
    // template<typename T_time, typename T_init, typename T_rate, typename T_par>
    // using default_model = PKOneCptModel<T_time, T_init, T_rate, T_par>;

  /**
   * Solve one-cpt model using analytical solution.
   *
   * @tparam T_time time type
   * @tparam T_model ODE model type
   */
    template<typename T_time, typename T_model>
    static
    Eigen::Matrix<torsten::scalar_t<T_model>, Eigen::Dynamic, 1>
    solve(const T_model &pkmodel, const T_time& dt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using scalar_type = torsten::scalar_t<T_model>;

      auto init = pkmodel.y0()   ;
      auto rate = pkmodel.rate() ;
      auto ka   = pkmodel.ka()   ;
      auto alpha= pkmodel.alpha();

      std::vector<scalar_type> a(Ncmt, 0);
      Matrix<scalar_type, 1, Dynamic> pred = PKRec<scalar_type>::Zero(Ncmt);

      if ((init[0] != 0) || (rate[0] != 0)) {
        pred(0, 0) = init[0] * exp(-ka * dt) + rate[0] * (1 - exp(-ka * dt)) / ka;
        a[0] = ka / (ka - alpha[0]);
        a[1] = -a[0];
        pred(0, 1) += torsten::PolyExp(dt, init[0], 0, 0, 0, false, a, alpha, 2) +
          torsten::PolyExp(dt, 0, rate[0], dt, 0, false, a, alpha, 2);
      }

      if ((init[1] != 0) || (rate[1] != 0)) {
        a[0] = 1;
        pred(0, 1) += torsten::PolyExp(dt, init[1], 0, 0, 0, false, a, alpha, 1) +
          torsten::PolyExp(dt, 0, rate[1], dt, 0, false, a, alpha, 1);
      }
      return pred;
    }

  /**
   * Solve one-cpt model: steady state solution
   *
   * @tparam T_time dosing interval time type
   * @tparam T_model ODE model type
   * @tparam T_amt dosing amount type
   */
    template<typename T_time, typename T_model, typename T_amt>
    static
    Eigen::Matrix<torsten::scalar_t<T_model>, Eigen::Dynamic, 1>
    solve(const T_model& pkmodel,
          const T_amt& amt,
          const T_time& ii,
          const int& cmt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;
      using scalar_type = scalar_t<T_model>;

      const double inf = std::numeric_limits<double>::max();  // "infinity"

      stan::math::check_positive("steady state one-cpt solver", "cmt", cmt);
      stan::math::check_less("steady state one-cpt solver", "cmt", cmt, 3);

      const auto rate = pkmodel.rate().at(cmt - 1);
      auto ka   = pkmodel.ka()   ;
      auto alpha= pkmodel.alpha();

      std::vector<scalar_type> a(2, 0);
      Matrix<scalar_type, 1, Dynamic> pred = Matrix<scalar_type, 1, Dynamic>::Zero(2);
      if (stan::math::value_of(rate) == 0) {  // bolus dose
        if (cmt == 1) {
          a[0] = 0;
          a[1] = 1;
          pred(0) = torsten::PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 2);
          a[0] = ka / (ka - alpha[0]);
          a[1] = -a[0];
          pred(1) = torsten::PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 2);
        } else {  // cmt=2
          a[0] = 1;
          pred(1) = torsten::PolyExp(ii, amt, 0, 0, ii, true, a, alpha, 1);
        }
      } else if (ii > 0) {  // multiple truncated infusions
        double delta = torsten::unpromote(amt / rate);
        static const char* function("Steady State Event");
        torsten::check_mti(amt, delta, ii, function);

        if (cmt == 1) {
          a[0] = 0;
          a[1] = 1;
          pred(0) = torsten::PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 2);
          a[0] = ka / (ka - alpha[0]);
          a[1] = -a[0];
          pred(1) = torsten::PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 2);
        } else {  // cmt = 2
          a[0] = 1;
          pred(1) = torsten::PolyExp(ii, 0, rate, amt / rate, ii, true, a, alpha, 1);
        }
      } else {  // constant infusion
        if (cmt == 1) {
          a[0] = 0;
          a[1] = 1;
          pred(0) = torsten::PolyExp(0, 0, rate, inf, 0, true, a, alpha, 2);
          a[0] = ka / (ka - alpha[0]);
          a[1] = -a[0];
          pred(1) = torsten::PolyExp(0, 0, rate, inf, 0, true, a, alpha, 2);
        } else {  // cmt = 2
          a[0] = 1;
          pred(1) = torsten::PolyExp(0, 0, rate, inf, 0, true, a, alpha, 1);
        }
      }
      return pred;
    }
  };

}

#endif
