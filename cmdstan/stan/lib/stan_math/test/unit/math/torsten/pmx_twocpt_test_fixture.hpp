#ifndef TEST_UNIT_TORSTEN_PK_TWOCPT_MODEL_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_TWOCPT_MODEL_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

class TorstenTwoCptTest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
public:
  TorstenTwoCptTest() :
    nCmt(3),
    nt(10),
    time(nt, 0.0),
    amt(nt, 0),
    rate(nt, 0),
    cmt(nt, 3),
    evid(nt, 0),
    ii(nt, 0),
    addl(nt, 0),
    ss(nt, 0),
    pMatrix{ {5, 8, 20, 70, 1.2 } },
    biovar{ { 1, 1, 1 } },
    tlag{ { 0, 0, 0 } }
  {
    SetUp();

    time[0] = 0;
    for(int i = 1; i < 9; i++) time[i] = time[i - 1] + 0.25;
    time[9] = 4.0;
    amt[0] = 1000;
    cmt[0] = 1;    
    evid[0] = 1;
    ii[0] = 12;
    addl[0] = 14;
  }

  const int nCmt;
  int nt;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  std::vector<std::vector<double> > pMatrix;  // CL, VC, Ka
  std::vector<std::vector<double> > biovar;
  std::vector<std::vector<double> > tlag;

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
