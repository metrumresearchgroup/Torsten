#ifdef TORSTEN_BRAID
#include <stan/math/torsten/mpi/cvodes_braid.hpp>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <memory>

using torsten::dsolve::PMXCvodesFwdSystem;
using torsten::dsolve::PMXCvodesIntegrator;
using torsten::dsolve::PMXOdeService;
using torsten::PMXCvodesSensMethod;
using torsten::mpi::CVBraidVec;
using torsten::mpi::CVBraidApp;
using stan::math::var;

TEST_F(TorstenOdeTest_neutropenia, braid_app) {
  int rank, size;
  MPI_Init(NULL, NULL);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int nt = 10;
  double h = 0.1;
  ts.resize(nt);
  for (int i = 0; i < nt; ++i) ts[i] = (i + 1) * h;
  auto y = torsten::pmx_integrate_ode_bdf(f, y0, t0, ts, theta, x_r, x_i);

  PMXCvodesIntegrator solver(rtol, atol, 1e8);
  using Ode = PMXCvodesFwdSystem<F, double, double, double, CV_BDF, torsten::AD>;
  PMXOdeService<typename Ode::Ode> s1(y0.size(), theta.size());
  Ode ode{s1, f, t0, ts, y0, theta, x_r, x_i, msgs};
  auto mem       = ode.mem();

  CVBraidApp app(mem, ode.to_user_data(), ode.nv_y(), MPI_COMM_WORLD, t0, ts.back(), nt);
  BraidUtil util;
  BraidCore core(comm, &app);
  core.SetMaxLevels(2);

  core.Drive();

  std::vector<std::vector<double> > res;
  if (app.rank == 0) {
    res = app.get_result();
    torsten::test::test_val(res, y, 1.e-3, 1.e-4);
  }

  MPI_Finalize();
}
#endif
