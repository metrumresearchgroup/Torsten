#ifndef PK_TWOCPT_SOLVER_HPP
#define PK_TWOCPT_SOLVER_HPP

#include <stan/math/torsten/pk_twocpt_model.hpp>

namespace refactor {

  using boost::math::tools::promote_args;

  /**
   * standard two compartment PK ODE solver based on
   * analytical solution.
   */
  class PKTwoCptModelSolver {
  public:
    static constexpr int Ncmt =
      PKTwoCptModel<double, double, double, double>::Ncmt;
    static constexpr int Npar =
      PKTwoCptModel<double, double, double, double>::Npar;

    /**
     * Constructor
     */
    PKTwoCptModelSolver() {}

  /**
   * standard two compartment PK model attached to this solver.
   * @tparam T_time t type
   * @tparam T_init initial condition type
   * @tparam T_rate rate type
   * @tparam T_par parameter type
   */
    template<typename T_time, typename T_init, typename T_rate, typename T_par>
    using default_model = PKTwoCptModel<T_time, T_init, T_rate, T_par>;

  /**
   * Solve two-cpt model.
   *
   * @tparam T_time time type
   * @tparam T_model ODE model type
   * @tparam Ts_par type of parameters
   */
    template<typename T_time, typename T_model>
    static
    Eigen::Matrix<torsten::scalar_t<T_model>, Eigen::Dynamic, 1> 
    solve(const T_model& pkmodel, const T_time& dt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using std::vector;
      using scalar_type = torsten::scalar_t<T_model>;

      auto init = pkmodel.y0()   ;
      auto rate = pkmodel.rate() ;
      auto ka   = pkmodel.ka()   ;
      auto k10  = pkmodel.k10()  ;
      auto k12  = pkmodel.k12()  ;
      auto k21  = pkmodel.k21()  ;
      auto alpha= pkmodel.alpha();

      std::vector<scalar_type> a(Ncmt, 0);
      Matrix<scalar_type, 1, Dynamic> pred = PKRec<scalar_type>::Zero(3);

      if ((init[0] != 0) || (rate[0] != 0))  {
        pred(0, 0) = init[0] * exp(-ka * dt) + rate[0] *
          (1 - exp(-ka * dt)) / ka;
        a[0] = ka * (k21 - alpha[0]) / ((ka - alpha[0]) * 
                                        (alpha[1] - alpha[0]));
        a[1] = ka * (k21 - alpha[1]) / ((ka - alpha[1]) *
                                        (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 1) += torsten::PolyExp(dt, init[0], 0, 0, 0, false, a, alpha, 3)
          + torsten::PolyExp(dt, 0, rate[0], dt, 0, false, a, alpha, 3);
        a[0] = ka * k12 / ((ka - alpha[0]) * (alpha[1] - alpha[0]));
        a[1] = ka * k12 / ((ka - alpha[1]) * (alpha[0] - alpha[1]));
        a[2] = -(a[0] + a[1]);
        pred(0, 2) += torsten::PolyExp(dt, init[0], 0, 0, 0, false, a, alpha, 3)
          + torsten::PolyExp(dt, 0, rate[0], dt, 0, false, a, alpha, 3);
      }

      if ((init[1] != 0) || (rate[1] != 0)) {
        a[0] = (k21 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k21 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 1) += torsten::PolyExp(dt, init[1], 0, 0, 0, false, a, alpha, 2)
          + torsten::PolyExp(dt, 0, rate[1], dt, 0, false, a, alpha, 2);
        a[0] = k12 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 2) += torsten::PolyExp(dt, init[1], 0, 0, 0, false, a, alpha, 2)
          + torsten::PolyExp(dt, 0, rate[1], dt, 0, false, a, alpha, 2);
      }

      if ((init[2] != 0) || (rate[2] != 0)) {
        a[0] = k21 / (alpha[1] - alpha[0]);
        a[1] = -a[0];
        pred(0, 1) += torsten::PolyExp(dt, init[2], 0, 0, 0, false, a, alpha, 2)
          + torsten::PolyExp(dt, 0, rate[2], dt, 0, false, a, alpha, 2);
        a[0] = (k10 + k12 - alpha[0]) / (alpha[1] - alpha[0]);
        a[1] = (k10 + k12 - alpha[1]) / (alpha[0] - alpha[1]);
        pred(0, 2) += torsten::PolyExp(dt, init[2], 0, 0, 0, false, a, alpha, 2)
          + torsten::PolyExp(dt, 0, rate[2], dt, 0, false, a, alpha, 2);
      }

      return pred;
    }
  };

}

#endif
