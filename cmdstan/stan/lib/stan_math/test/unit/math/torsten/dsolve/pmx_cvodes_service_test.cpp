#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/torsten/dsolve/sundials_check.hpp>
#include <stan/math/torsten/dsolve/pmx_integrate_ode_adams.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_integrator.hpp>
#include <stan/math/torsten/dsolve/cvodes_service.hpp>
#include <stan/math/torsten/dsolve/cvodes_rhs.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <ostream>
#include <vector>


struct Sho {
  harm_osc_ode_fun f;
  std::vector<double> y;
  std::vector<double> fy;
  const std::vector<double>& p;
  const std::vector<double>& xr;
  const std::vector<int>& xi;
  std::ostream* m;    
  Sho(const std::vector<double>& theta,
      const std::vector<double>& x_r,
      const std::vector<int>& x_i,
      std::ostream* msgs) :
    y(2), fy(2), p(theta), xr(x_r), xi(x_i), m(msgs)
  {}

  void eval_rhs(double t, N_Vector& nvy, N_Vector& nvydot) {
    for (size_t i = 0; i < 2; ++i) y[i] = NV_Ith_S(nvy, i);
    fy = f(t, y, p, xr, xi, m);
    for (size_t i = 0; i < 2; ++i) NV_Ith_S(nvydot, i) = fy[i];
  }

  static CVRhsFn rhs() {
    return torsten::dsolve::cvodes_rhs<Sho>();
  }

  static constexpr int lmm_type       = CV_ADAMS;
  static constexpr bool is_var_y0     = false;
  static constexpr bool is_var_par    = false;
  static constexpr bool need_fwd_sens = false;
};

TEST_F(TorstenOdeTest_sho, service) {
  using torsten::dsolve::PMXCvodesIntegrator;
  using torsten::dsolve::PMXOdeService;

  PMXOdeService<Sho> serv(2, 1);

  size_t n = 2;
  N_Vector& y = serv.nv_y;
  double t1 = t0;
  auto yy = stan::math::integrate_ode_adams(f, y0, t0,
                                            ts, theta, x_r, x_i);
  Sho ode(theta, x_r, x_i, msgs);
  void* user_data = static_cast<void*>(&ode);
  for (int i = 0; i < 10; ++i) {
    void* mem = serv.mem;
    for (size_t i = 0; i < n; ++i) NV_Ith_S(y, i) = y0[i];
    CHECK_SUNDIALS_CALL(CVodeReInit(mem, 0.0, y));
    CHECK_SUNDIALS_CALL(CVodeSStolerances(mem, rtol, atol));
    CHECK_SUNDIALS_CALL(CVodeSetUserData(mem, user_data));
    CHECK_SUNDIALS_CALL(CVodeSetMaxNumSteps(mem, max_num_steps));
    CHECK_SUNDIALS_CALL(CVode(mem, ts[0], y, &t1, CV_NORMAL));

    for (size_t i = 0; i < n; ++i) EXPECT_FLOAT_EQ(NV_Ith_S(y, i), yy[0][i]);
  }
}
