#ifndef STAN_MATH_TORSTEN_ONECPT_SOLVER_HPP
#define STAN_MATH_TORSTEN_ONECPT_SOLVER_HPP

#include <stan/math/torsten/pk_onecpt_model.hpp>

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
   * Solve one-cpt model.
   *
   * @tparam T_time time type
   * @tparam T_model ODE model type
   */
    template<typename T_time, typename T_model>
    static
    Eigen::Matrix<typename T_model::scalar_type, Eigen::Dynamic, 1>
    solve(const T_model &pkmodel, const T_time& dt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using scalar_type = typename T_model::scalar_type;

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
  };

}

#endif
