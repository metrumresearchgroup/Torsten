#ifndef STAN_MATH_TORSTEN_DSOLVE_ODE_SYSTEM_HPP
#define STAN_MATH_TORSTEN_DSOLVE_ODE_SYSTEM_HPP

#include <stan/math/torsten/dsolve/pk_vars.hpp>

namespace torsten {
  namespace dsolve {

    /**
     * General ODE system that contains informtion on residual
     * equation functor, sensitivity residual equation functor,
     * as well as initial conditions. This is a base type that
     * is intended to contain common values used by forward
     * sensitivity system.
     *
     * @tparam Tt scalar type of time steps
     * @tparam Ty scalar type of initial unknown values
     * @tparam Tp scalar type of parameters
     */
    template <typename Tt, typename Ty, typename Tp>
    struct PMXOdeSystem {
      virtual inline const std::vector<Tt>& ts() const = 0;
      virtual inline const std::vector<Ty>& y0() const = 0;
      virtual inline const std::vector<Tp>& theta() const = 0;
      virtual inline const size_t fwd_system_size() const = 0;
      const std::vector<stan::math::var> vars() const {
        return pk_vars(y0(), theta(), ts());
      }
    };

    template <>
    struct PMXOdeSystem<stan::math::var, double, double> {
      virtual inline const std::vector<stan::math::var>& ts() const = 0;
      virtual inline const std::vector<double>& y0() const = 0;
      virtual inline const std::vector<double>& theta() const = 0;
      virtual inline const size_t fwd_system_size() const = 0;
      const std::vector<stan::math::var>& vars() const {
        return ts();
      }
    };

    template <>
    struct PMXOdeSystem<double, stan::math::var, double> {
      virtual inline const std::vector<double>& ts() const = 0;
      virtual inline const std::vector<stan::math::var>& y0() const = 0;
      virtual inline const std::vector<double>& theta() const = 0;
      virtual inline const size_t fwd_system_size() const = 0;
      const std::vector<stan::math::var>& vars() const {
        return y0();
      }
    };

    template <>
    struct PMXOdeSystem<double, double, stan::math::var> {
      virtual inline const std::vector<double>& ts() const = 0;
      virtual inline const std::vector<double>& y0() const = 0;
      virtual inline const std::vector<stan::math::var>& theta() const = 0;
      virtual inline const size_t fwd_system_size() const = 0;
      const std::vector<stan::math::var>& vars() const {
        return theta();
      }
    };
  }
}

#endif
