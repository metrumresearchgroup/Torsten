#ifndef STAN_MATH_TORSTEN_DSOLVE_PMX_ODE_INTEGRATOR_HPP
#define STAN_MATH_TORSTEN_DSOLVE_PMX_ODE_INTEGRATOR_HPP

#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/functor/jacobian.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/err/check_positive_finite.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_integrator.hpp>
#include <stan/math/torsten/dsolve/pmx_odeint_integrator.hpp>
#include <ostream>
#include <vector>

namespace torsten {
  namespace dsolve {
    struct PMXOdeIntegratorBase {
      const double rtol;
      const double atol;
      const long int max_num_step;
      const double as_rtol;
      const double as_atol;
      const long int as_max_num_step;
      std::ostream* msgs;      

      const char* caller = "pmx_integrate_ode";

      PMXOdeIntegratorBase(double rtol0,
                       double atol0,
                       long int max_num_step0,
                       double as_rtol0,
                       double as_atol0,
                       long int as_max_num_step0,
                       std::ostream* msgs0) :
        rtol(rtol0), atol(atol0), max_num_step(max_num_step0),
        as_rtol(as_rtol0), as_atol(as_atol0), as_max_num_step(as_max_num_step0),
        msgs(msgs0)
      {
        check_controls();        
      }
      
      PMXOdeIntegratorBase(double rtol0, double atol0, long int max_num_step0,
                       std::ostream* msgs0) :
        rtol(rtol0), atol(atol0), max_num_step(max_num_step0),
        as_rtol(1e-6), as_atol(1e-6), as_max_num_step(100),
        msgs(msgs0)
      {
        check_controls();        
      }

      PMXOdeIntegratorBase() :
        rtol(1e-10), atol(1e-10), max_num_step(1e8),
        as_rtol(1e-6), as_atol(1e-6), as_max_num_step(100),
        msgs(0)
      {
        check_controls();
      }

      inline void check_controls() {
        using stan::math::check_positive_finite;
        check_positive_finite(caller, "relative tolerance", rtol);
        check_positive_finite(caller, "absolute tolerance", atol);
        check_positive_finite(caller, "maximum number of steps", max_num_step);
        check_positive_finite(caller, "algebra solver relative tolerance", as_rtol);
        check_positive_finite(caller, "algebra solver absolute tolerance", as_atol);
        check_positive_finite(caller, "algebra solver maximum number of steps", as_max_num_step);
      }
    };

    template<template<typename...> class ode_type, typename integrator_t>
    struct PMXOdeIntegrator;

    template<typename integrator_t>
    struct PMXOdeIntegrator<PMXOdeSystem, integrator_t> : public PMXOdeIntegratorBase {
      PMXOdeIntegrator(double rtol0,
                       double atol0,
                       long int max_num_step0,
                       double as_rtol0,
                       double as_atol0,
                       long int as_max_num_step0,
                       std::ostream* msgs0) :
        PMXOdeIntegratorBase(rtol0, atol0, max_num_step0, as_rtol0,
                             as_atol0, as_max_num_step0, msgs0)
      {}
      
      PMXOdeIntegrator(double rtol0, double atol0, long int max_num_step0,
                       std::ostream* msgs0) :
        PMXOdeIntegratorBase(rtol0, atol0, max_num_step0,
                             1e-6, 1e-6, 100, msgs0)
      {}

      PMXOdeIntegrator() :
        PMXOdeIntegratorBase(1e-10, 1e-10, 1e8,
                             1e-6, 1e-6, 100, 0)
      {}

      /**
       * Integrator that takes "old" std::vector args
       * 
       * @return A 2D array of solutions.
       */
      template <typename F, typename Tt, typename T_initial, typename T_param>
      std::vector<std::vector<typename stan::return_type_t<Tt, T_initial, T_param>> >
      operator()(const F& f,
                 const std::vector<T_initial>& y0,
                 double t0,
                 const std::vector<Tt>& ts,
                 const std::vector<T_param>& theta,
                 const std::vector<double>& x_r,
                 const std::vector<int>& x_i) const {
        using Ode = PMXOdeSystem<F, Tt, T_initial, T_param>;
        Ode ode{f, t0, ts, y0, theta, x_r, x_i, msgs};
        integrator_t solver(rtol, atol, max_num_step);
        dsolve::OdeObserver<Ode> observer(ode);
        solver.integrate(ode, observer);
        return observer.y;
      }

      template <typename F, typename Tt, typename T_initial, typename T_param>
      std::vector<std::vector<typename stan::return_type_t<Tt, T_initial, T_param>> >
      operator()(const F& f,
                 const std::vector<T_initial>& y0,
                 double t0,
                 const Tt& t1,
                 const std::vector<T_param>& theta,
                 const std::vector<double>& x_r,
                 const std::vector<int>& x_i) const {
        std::vector<Tt> ts{t1};
        return (*this)(f, y0, t0, ts, theta, x_r, x_i);
      }

      template <typename F, typename Tt, typename T_initial, typename T_param>
      Eigen::MatrixXd
      solve_d(const F& f,
              const std::vector<T_initial>& y0,
              double t0,
              const std::vector<Tt>& t1,
              const std::vector<T_param>& theta,
              const std::vector<double>& x_r,
              const std::vector<int>& x_i) const {
        std::vector<Tt> ts{t1};
        using Ode = PMXOdeSystem<F, Tt, T_initial, T_param>;
        Ode ode{f, t0, ts, y0, theta, x_r, x_i, msgs};
        integrator_t solver(rtol, atol, max_num_step);
        dsolve::OdeDataObserver<Ode> observer(ode);
        solver.integrate(ode, observer);
        return observer.y;
      }
    };

    /** 
     * Specialization for ODE systems that take parameter pack
     * 
     * @param rtol0 
     * @param atol0 
     * @param max_num_step0 
     * @param as_rtol0 
     * @param as_atol0 
     * @param as_max_num_step0 
     * @param msgs0 
     * 
     */
    template<typename integrator_t>
    struct PMXOdeIntegrator<PMXVariadicOdeSystem, integrator_t> : public PMXOdeIntegratorBase {
      PMXOdeIntegrator(double rtol0,
                       double atol0,
                       long int max_num_step0,
                       double as_rtol0,
                       double as_atol0,
                       long int as_max_num_step0,
                       std::ostream* msgs0) :
        PMXOdeIntegratorBase(rtol0, atol0, max_num_step0, as_rtol0,
                             as_atol0, as_max_num_step0, msgs0)
      {}
      
      PMXOdeIntegrator(double rtol0, double atol0, long int max_num_step0,
                       std::ostream* msgs0) :
        PMXOdeIntegratorBase(rtol0, atol0, max_num_step0,
                             1e-6, 1e-6, 100, msgs0)
      {}

      PMXOdeIntegrator() :
        PMXOdeIntegratorBase(1e-10, 1e-10, 1e8,
                             1e-6, 1e-6, 100, 0)
      {}

      template <typename F, typename Tt, typename T_initial, typename... T_par>
      std::vector<Eigen::Matrix<typename stan::return_type_t<Tt, T_initial, T_par...>, -1, 1>>
      operator()(const F& f,
                 const Eigen::Matrix<T_initial, -1, 1>& y0,
                 double t0,
                 const std::vector<Tt>& ts,
                 const T_par&... args) const {
        using Ode = PMXVariadicOdeSystem<F, Tt, T_initial, T_par...>;
        Ode ode{f, t0, ts, y0, msgs, args...};
        integrator_t solver(rtol, atol, max_num_step);
        dsolve::OdeObserver<Ode> observer(ode);
        solver.integrate(ode, observer);
        return observer.y;
      }

      template <typename F, typename Tt, typename T_initial, typename... T_par>
      Eigen::Matrix<typename stan::return_type_t<Tt, T_initial, T_par...>, -1, 1>
      operator()(const F& f,
                 const Eigen::Matrix<T_initial, -1, 1>& y0,
                 double t0,
                 const Tt& t1,
                 const T_par&... args) const {
        const double dt_min = 1.e-10; // too small would cause CVODES -27 error
        if (stan::math::value_of(t1) - stan::math::value_of(t0) > dt_min) {
          std::vector<Tt> ts{t1};
          return (*this)(f, y0, t0, ts, args...)[0];
        } else {                // single-step forward Euler
          return y0 + (t1 - t0) * f(t0, y0, nullptr, args...);
        }
      }

      template <typename F, typename Tt, typename T_initial, typename... T_par>
      Eigen::MatrixXd
      solve_d(const F& f,
              const Eigen::Matrix<T_initial, -1, 1>& y0,
              double t0,
              const std::vector<Tt>& ts,
              const T_par&... args) const {
        using Ode = PMXVariadicOdeSystem<F, Tt, T_initial, T_par...>;
        Ode ode{f, t0, ts, y0, msgs, args...};
        integrator_t solver(rtol, atol, max_num_step);
        dsolve::OdeDataObserver<Ode> observer(ode);
        solver.integrate(ode, observer);
        return observer.y;
      }
    };

    /**
     * Dummy integrator type for analytical solutions such as one/two cpt.
     * 
     */
    struct PMXAnalyiticalIntegrator {};
  }

  template<typename integrator_type>
  struct has_data_only_output {
    static const bool value = false;
  };

  template<typename scheme_t>
  struct has_data_only_output<dsolve::PMXOdeIntegrator<dsolve::PMXVariadicOdeSystem, dsolve::PMXOdeintIntegrator<scheme_t>>> {
    static const bool value = true;
  };

  template<int lmm_type, int ism_type>
  struct
  has_data_only_output<dsolve::PMXOdeIntegrator<dsolve::PMXVariadicOdeSystem, dsolve::PMXCvodesIntegrator<lmm_type, ism_type>>> {
    static const bool value = true;
  };
}

#endif
