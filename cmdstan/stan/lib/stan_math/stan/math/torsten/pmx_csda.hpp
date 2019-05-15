#ifndef TORSTEN_COMPLEX_STEP_DERIVATIVE_HPP
#define TORSTEN_COMPLEX_STEP_DERIVATIVE_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/core/precomp_v_vari.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/var.hpp>

#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>

namespace torsten {

/**
 * Return a var that has value of given functor F and derivative
 * of df/d(theta), using complex step derivative
 * approximation. "f" does not have to support "var"
 * type, as its signature should be
 * (complex, std::vector<double>, std::vector<int>, stream*) : complex
 *
 * @tparam F type of functor F
 * @param[in] f functor for the complex number evaluation,
 * must support @c std::complex<double> as arg.
 * @param[in] theta parameter where f and df/d(theta) is requested.
 * @param[in] x_r continuous data vector for the ODE.
 * @param[in] x_i integer data vector for the ODE.
 * @param[in] h complex step size
 * @param[out] msgs the print stream for warning messages.
 * @return a var with value f(theta.val()) and derivative at theta.
 */
template <typename F>
stan::math::var pk_csda(const F& f,
                        const stan::math::var& theta,
                        const std::vector<double>& x_r,
                        const std::vector<int>& x_i,
                        const double h,
                        std::ostream* msgs) {
  using stan::math::var;
  using std::complex;
  const double theta_d = theta.val();
  const complex<double> res = f(complex<double>(theta_d, h), x_r, x_i, msgs);
  const double fx = std::real(res);
  const double g = std::imag(res) / h;
  return var(new stan::math::precomp_v_vari(fx, theta.vi_, g));
}

/**
 * CSDA, calculate directional derivative of a vector function
 *
 */
template <typename F>
std::vector<double> pk_csda(const F& f,
                            const std::vector<double>& y,
                            const std::vector<double>& dy,
                            const std::vector<double>& x_r,
                            const std::vector<int>& x_i,
                            std::ostream* msgs) {
  using stan::math::var;
  using std::complex;
  using cplx = complex<double>;
  const double h = 1.E-20;
  std::vector<cplx> cplx_y(y.size());
  std::transform(y.begin(), y.end(), dy.begin(), cplx_y.begin(),
                 [&h](const double &r, const double &i) -> cplx {
                    return cplx(r, h * i); }
                 );
  const std::vector<cplx> res = f(cplx_y, x_r, x_i, msgs);
  std::vector<double> g(y.size());
  std::transform(res.begin(), res.end(), g.begin(),
                 [&h](cplx x) -> double { return std::imag(x)/h; });
  return g;
}

/**
 * CSDA, calculate directional derivative of a vector
 * function, without system input @c x_r, x_i, msgs
 */
template <typename F>
std::vector<double> pk_csda(const F& f,
                            const std::vector<double>& y,
                            const std::vector<double>& dy) {
  using std::complex;
  using cplx = complex<double>;
  const double h = 1.E-20;
  std::vector<cplx> cplx_y(y.size());
  std::transform(y.begin(), y.end(), dy.begin(), cplx_y.begin(),
                 [&h](const double &r, const double &i) -> cplx {
                    return cplx(r, h * i); }
                 );
  const std::vector<cplx> res = f(cplx_y);
  std::vector<double> g(y.size());
  std::transform(res.begin(), res.end(), g.begin(),
                 [&h](cplx x) -> double { return std::imag(x)/h; });
  return g;
}

/**
 * CSDA, default h version, with h = 1.E-20
 *
 * @tparam F type of functor F
 * @param[in] f functor for the complex number evaluation,
 * must support @c std::complex<double> as arg.
 * @param[in] theta parameter where f and df/d(theta) is requested.
 * @param[in] x_r continuous data vector for the ODE.
 * @param[in] x_i integer data vector for the ODE.
 * @param[out] msgs the print stream for warning messages.
 * @return a var with value f(theta.val()) and derivative at theta.
 */
template <typename F>
stan::math::var pk_csda(const F& f,
                        const stan::math::var& theta,
                        const std::vector<double>& x_r,
                        const std::vector<int>& x_i,
                        std::ostream* msgs) {
  return pk_csda(f, theta, x_r, x_i, 1.E-20, msgs);
}

}  // namespace torsten

#endif
