#include <gtest/gtest.h>
#include <stan/math/torsten/test/unit/test_util.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/meta/is_std_ode.hpp>
#include <stan/math/torsten/meta/is_eigen_ode.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/dsolve/ode_tuple_functor.hpp>
#include <stan/math/torsten/ode_rhs_ostream_adatpor.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>

TEST(Torsten, ode_signature) {
  static_assert(torsten::is_std_ode<harm_osc_ode_fun>::value, "Wrong ODE signature");
  static_assert(!torsten::is_std_ode<harm_osc_ode_fun_eigen>::value, "Wrong ODE signature");
  static_assert(torsten::is_eigen_ode<harm_osc_ode_fun_eigen,
                std::vector<double>, std::vector<double>, std::vector<int>>::value,
                "Wrong ODE signature");
  static_assert(torsten::is_eigen_ode<harm_osc_ode_fun_eigen,
                std::vector<stan::math::var>, std::vector<double>, std::vector<int>>::value,
                "Wrong ODE signature");
  static_assert(!torsten::is_eigen_ode<harm_osc_ode_fun,
                std::vector<double>, std::vector<double>, std::vector<int>>::value,
                "Wrong ODE signature");
  torsten::std_ode<harm_osc_ode_fun> ode1;
  torsten::eigen_ode<harm_osc_ode_fun_eigen,
                     std::vector<double>, std::vector<double>, std::vector<int>> ode2;
}

TEST(Torsten, ode_variadic_and_tuple_functor) {
  using torsten::dsolve::TupleOdeFunc;
  using torsten::dsolve::VariadicOdeFunc;

  Eigen::Matrix<stan::math::var, -1, 1> y_var(2);
  y_var << 1.0, 0.5;

  std::vector<stan::math::var> theta_var(1);
  theta_var[0] = 0.15;

  std::vector<double> x_r;
  std::vector<int> x_i;
  double t = 1.0;

  harm_osc_ode_fun_eigen f0;
  auto y0 = f0(t, y_var, nullptr, theta_var, x_r, x_i);

  VariadicOdeFunc<harm_osc_ode_fun_eigen, double, stan::math::var>
    f1(f0, t, y_var, nullptr);
  auto y1 = f1(theta_var, x_r, x_i);

  auto theta_tuple = std::make_tuple(theta_var, x_r, x_i);

  TupleOdeFunc<harm_osc_ode_fun_eigen> f(f0);
  auto y2 = f(1.0, y_var, nullptr, theta_tuple);

  torsten::test::test_grad(theta_var, y0, y1, 1.e-10, 1.e-10);
  torsten::test::test_grad(theta_var, y0, y2, 1.e-10, 1.e-10);
  std::vector<stan::math::var> y(stan::math::to_array_1d(y_var));
  torsten::test::test_grad(y, y0, y1, 1.e-10, 1.e-10);
  torsten::test::test_grad(y, y0, y2, 1.e-10, 1.e-10);
}

struct Func
{
  template <typename T0, typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type_t<T0, T1, T2>, -1, 1>
  operator()(const T0& t,
             const Eigen::Matrix<T1, -1, 1>& x,
             const std::vector<T2>& parms,
             const std::vector<double>& x_r,
             const std::vector<int>& x_i,
             std::ostream* pstream__) const {
    typedef typename stan::return_type<T0, T1, T2>::type scalar;
    return Eigen::Matrix<scalar, -1, 1>(2);
  }
};

TEST(Torsten, ode_ostream_arg_placement_adaptor) {
  static_assert(torsten::pmx_ode_ostream_not_adapted<Func>::value, "Wrong ODE ostream argument location");
  static_assert(!torsten::pmx_ode_ostream_adapted<Func>::value, "Wrong ODE ostream argument location");
  static_assert(torsten::pmx_ode_ostream_adapted<torsten::PMXOneCptODE>::value, "Wrong ODE ostream argument location");
  static_assert(!torsten::pmx_ode_ostream_not_adapted<torsten::PMXOneCptODE>::value, "Wrong ODE ostream argument location");
}
