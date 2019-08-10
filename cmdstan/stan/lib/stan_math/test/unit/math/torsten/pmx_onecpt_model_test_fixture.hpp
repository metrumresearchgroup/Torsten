#ifndef TEST_UNIT_TORSTEN_PK_ONECPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_ONECPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_ode_model.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct TorstenOneCptModelTest : public testing::Test {
  double t0;
  std::vector<double> ts;
  Eigen::Matrix<double, 1, Eigen::Dynamic> y0;
  std::vector<double> rate;
  double CL;
  double V2;
  double ka;
  const std::vector<double> x_r;
  const std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
  TorstenOneCptModelTest() :
    t0(0.0),
    ts{0.1, 0.5, 1.0},
    y0(2),
    rate(2, 0.0),
    CL(50.0),
    V2(80.0),
    ka(1.2),
    msgs{nullptr} {
      y0 << 0.0, 0.0;
      SetUp();
  }
};

struct TorstenOneCptTest : public testing::Test {
  int nCmt;
  int nt;
  std::vector<std::vector<double> > pMatrix;
  std::vector<std::vector<double> > biovar;
  std::vector<std::vector<double> > tlag;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }

  TorstenOneCptTest() :
    nCmt(2),
    nt(10),
    pMatrix(1, {10, 80, 1.2}),
    biovar(1, {1, 1}),
    tlag(1, {0, 0}),
    time(nt),
    amt(10, 0),
    rate(10, 0),
    cmt(10, 2),
    evid(10, 0),
    ii(10, 0),
    addl(10, 0),
    ss(10, 0)
  {
    time[0] = 0.0;
    for(int i = 0; i < nt - 1; i++) time[i] = i * 0.25;
    time[nt - 1] = 4.0;

    amt[0] = 1000.0;

    cmt[0] = 1;

    evid[0] = 1;

    ii[0] = 12;

    addl[0] = 14;

    SetUp();
  }

  void resize(int n) {
    nt = n;
    time.resize(nt);
    amt .resize(nt);
    rate.resize(nt);
    cmt .resize(nt);
    evid.resize(nt);
    ii  .resize(nt);
    addl.resize(nt);
    ss  .resize(nt);
  }
};

#endif
