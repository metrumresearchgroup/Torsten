#ifndef STAN_MATH_TORSTEN_MPI_PRECOMPUTED_GRADIENTS_HPP
#define STAN_MATH_TORSTEN_MPI_PRECOMPUTED_GRADIENTS_HPP

#include<stan/math/rev/core.hpp>
#include<stan/math/prim/functor/apply.hpp>
#include<stan/math/prim/fun/to_array_1d.hpp>
#include<stan/math/prim/meta/is_vector.hpp>
#include<stan/math/torsten/dsolve/ode_tuple_functor.hpp>
#include <nvector/nvector_serial.h>

namespace torsten {

  template<typename... T_par>
  inline stan::math::vari**
  varis_from_ode_pars(const std::tuple<const T_par&...>& arg_tuple) {
    using stan::math::vari;
    const size_t m = torsten::dsolve::count_vars_in_tuple(arg_tuple);
    vari** varis = stan::math::ChainableStack::instance_->memalloc_.alloc_array<vari*>(m);
    stan::math::apply([&](auto&&... args) {stan::math::save_varis(varis, args...);}, arg_tuple);
    return varis;
  }

  inline stan::math::vari**
  varis_from_ode_pars(const std::vector<stan::math::var>& theta) {
    using stan::math::vari;
    vari** varis = stan::math::ChainableStack::instance_->memalloc_.alloc_array<vari*>(theta.size());
    save_varis(varis, theta);
    return varis;
  }

  /**
   * 
   * @tparam Ty @c vector type, can be either std::vector or eigen::vector
   */
  template<typename Ty = std::vector<double>,
    typename Ty_out = std::vector<stan::math::var>,
    std::enable_if_t<stan::is_vector<Ty>::value>* = nullptr>
  inline Ty_out
  precomputed_gradients(const Ty& sol, const std::vector<stan::math::var>& theta) {
    const size_t n = sol.size()/(1 + theta.size());
    Ty_out y_res(n);
    std::vector<double> g(theta.size());

    stan::math::vari** varis = varis_from_ode_pars(theta);

    for (size_t j = 0; j < n; j++) {
      double* g = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(theta.size());
      for (size_t k = 0; k < theta.size(); k++) *(g + k) = sol[n + n * k + j];
      y_res[j] = new stan::math::precomputed_gradients_vari(sol[j], theta.size(), varis, g);
    }

    return y_res;
  }

  /** 
   * Generate <code>var vector</code> that has specified ODE solution
   * w.r.t. a tuple of args from variadic args.
   * 
   * @param sol data-only solution from ODE solvers
   * @param arg_tuple tuple from variadic args.
   * 
   * @return a <code>var vector</code> as ODE soution.
   */
  template<typename Ty, typename... T_par>
  Eigen::Matrix<stan::math::var, -1, 1>  
  precomputed_gradients(const Ty& sol, const std::tuple<const T_par&...>& arg_tuple) {
    using stan::math::ChainableStack;

    const size_t m = torsten::dsolve::count_vars_in_tuple(arg_tuple);
    const size_t n = sol.size()/(1 + m);
    Eigen::Matrix<stan::math::var, -1, 1>  y_res(n);

    stan::math::vari** varis = varis_from_ode_pars(arg_tuple);

    for (size_t j = 0; j < n; j++) {
      double* g = ChainableStack::instance_->memalloc_.alloc_array<double>(m);
      for (size_t k = 0; k < m; k++) *(g + k) = sol[n + n * k + j];
      y_res[j] = new stan::math::precomputed_gradients_vari(sol[j], m, varis, g);
    }

    return y_res;
  }

  /** 
   * Generate <code>var array1d</code> that has specified CVODES ODE solution
   * w.r.t. a tuple of args from variadic args.
   * 
   * @param y value solution from CVODES
   * @param ys sensitiy solution from CVODES
   * @param theta parameter
   * 
   * @return a <code>var array</code> of ODE solution
   */
  inline std::vector<stan::math::var>
  precomputed_gradients(const N_Vector& y, const N_Vector* ys,
                        const std::vector<stan::math::var>& theta) {
    const int n = NV_LENGTH_S(y);
    const int ns = theta.size();
    std::vector<double> g(ns);
    std::vector<stan::math::var> y_res(n);
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < ns; ++j) {
        g[j] = NV_Ith_S(ys[j], k);
      }
      y_res[k] = precomputed_gradients(NV_Ith_S(y, k), theta, g);
    }
    return y_res;
  }

  /** 
   * Generate <code>var vector</code> that has specified CVODES ODE solution
   * w.r.t. a tuple of args from variadic args.
   * 
   * @param y value solution from CVODES
   * @param ys sensitiy solution from CVODES
   * @param arg_tuple tuple of variadic parameters
   * 
   * @return a <code>var vector</code> of ODE solution
   */
  template<typename... T_par>
  inline Eigen::Matrix<stan::math::var, -1, 1>
  precomputed_gradients(const N_Vector& y, const N_Vector* ys,
                        const std::tuple<const T_par&...>& arg_tuple) {
    using stan::math::ChainableStack;
    const int n = NV_LENGTH_S(y);
    const size_t ns = torsten::dsolve::count_vars_in_tuple(arg_tuple);
    Eigen::Matrix<stan::math::var, -1, 1> y_res(n);

    stan::math::vari** varis = varis_from_ode_pars(arg_tuple);
    for (size_t k = 0; k < n; ++k) {
      double* g = ChainableStack::instance_->memalloc_.alloc_array<double>(ns);
      for (size_t j = 0; j < ns; ++j) {
        *(g + j) = NV_Ith_S(ys[j], k);
      }
      y_res[k] = new stan::math::precomputed_gradients_vari(NV_Ith_S(y, k), ns, varis, g);
    }
    return y_res;
  }

  template<typename Ty_out>
  inline Ty_out
  precomputed_gradients(const Eigen::MatrixXd& sol, int i_col, const std::vector<stan::math::var>& theta) {
    const size_t n = sol.rows()/(1 + theta.size());
    Ty_out y_res(n);
    std::vector<double> g(theta.size());
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < theta.size(); k++) g[k] = sol(n + n * k + j, i_col);
      y_res[j] = precomputed_gradients(sol(j, i_col), theta, g);
    }
    return y_res;
  }

  inline stan::math::matrix_v
  precomputed_gradients (const Eigen::MatrixXd& sol, const std::vector<stan::math::var>& theta) {
    const size_t n = sol.rows()/(1 + theta.size());
    stan::math::matrix_v y_res(n, sol.cols());
    for (int i = 0; i < sol.cols(); i++) {
      y_res.col(i) = precomputed_gradients<stan::math::vector_v>(sol, i, theta);
    }
    return y_res;
  }

  inline stan::math::matrix_v
  precomputed_gradients (const Eigen::MatrixXd& sol, const Eigen::Matrix<stan::math::var,-1,1>& theta) {
    const size_t n = sol.rows()/(1 + theta.size());
    stan::math::matrix_v y_res(n, sol.cols());
    std::vector<stan::math::var> theta_ = stan::math::to_array_1d(theta);
    for (int i = 0; i < sol.cols(); i++) {
      y_res.col(i) = precomputed_gradients<stan::math::vector_v>(sol, i, theta_);
    }
    return y_res;
  }

  namespace mpi {
    /**
     * Generate a Eigen::Vector with @c var entries that have given
     * value and gradients. The value and gradients are provided
     * through @c VectorXd consisting of multiple
     * records in the format (value, grad1, grad2, grad3...).
     * @param d input vector data with consisting of
     * multiple records, each record of format (value, grad1, grad2...)
     * @return @c var matrix with given value and gradients
     */
    inline stan::math::vector_v
    precomputed_gradients(const Eigen::VectorXd& sol, const std::vector<stan::math::var>& theta) {
      const size_t n = sol.size()/(1 + theta.size());
      stan::math::vector_v y_res(n);
      std::vector<double> g(theta.size());
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < theta.size(); k++) g[k] = sol(n + n * k + j);
        y_res(j) = precomputed_gradients(sol(j), theta, g);
      }
      return y_res;
    }
  }
}
#endif

    
