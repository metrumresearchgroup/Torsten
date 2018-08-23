#ifndef PK_LINODE_SOLVER_HPP
#define PK_LINODE_SOLVER_HPP

namespace refactor {

  using boost::math::tools::promote_args;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  /**
   * Linear ODE solver based on matrix exponential.
   *
   */
  class PKLinODEModelSolver {
  public:
    PKLinODEModelSolver() {}

  /**
   * Solve Linear ODE model, or any model that can provide
   * RHS matrix (thus can be seen as linear ODE model)
   *
   * @tparam T_time  type of time
   * @tparam T_model type of model
   * @param pkmodel  linear ODE model
   * @param dt  time span
   */
    template<typename T_time, typename T_model>
    static
    Eigen::Matrix<torsten::scalar_t<T_model>, Eigen::Dynamic, 1> 
    solve(const T_model &pkmodel, const T_time& dt) {
      using Eigen::Matrix;
      using Eigen::Dynamic;
      using stan::math::value_of;
      using stan::math::matrix_exp;
      using stan::math::mdivide_left;
      using stan::math::multiply;
      using stan::math::scale_matrix_exp_multiply;
      using scalar_type = torsten::scalar_t<T_model>;

      auto init = pkmodel.y0()   ;
      auto rate = pkmodel.rate() ;
      auto system = pkmodel.coef();
      const int nCmt = system.cols();
      Matrix<scalar_type, Dynamic, 1> y0t(nCmt);
      for (int i = 0; i < nCmt; ++i) y0t(i) = init(i);

      if (std::any_of(rate.begin(), rate.end(),
                      [](torsten::rate_t<T_model> r){return r != 0;})) {
        Matrix<scalar_type, Dynamic, 1> rate_vec(rate.size()), x(nCmt), x2(nCmt);
        for (size_t i = 0; i < rate.size(); i++) rate_vec(i) = rate[i];
        x = mdivide_left(system, rate_vec);
        x2 = x + init.transpose();
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

      // if (dt == 0) { return init;
      // } else {


      //   bool rate_zeros = true;
      //   for (size_t i = 0; i < rate.size(); i++)
      //     if (rate[i] != 0) rate_zeros = false;

      //   // trick to promote dt, and dt_system
      //   scalar_type dt_s = dt;

      //   if (rate_zeros) {
      //     Matrix<scalar_type, Dynamic, Dynamic> dt_system = multiply(dt_s, system);
      //     Matrix<scalar_type, Dynamic, 1> pred = matrix_exp(dt_system)
      //       * init.transpose();
      //     return pred.transpose();
      //   } else {
      //     int nCmt = system.cols();
      //     Matrix<scalar_type, Dynamic, 1> rate_vec(rate.size()), x(nCmt), x2(nCmt);
      //     for (size_t i = 0; i < rate.size(); i++) rate_vec(i) = rate[i];
      //     x = mdivide_left(system, rate_vec);
      //     x2 = x + init.transpose();
      //     Matrix<scalar_type, Dynamic, Dynamic> dt_system = multiply(dt_s, system);
      //     Matrix<scalar_type, Dynamic, 1> pred = matrix_exp(dt_system) * x2;
      //     pred -= x;
      //     return pred.transpose();
      //   }
      // }
    }
  };
}

#endif
