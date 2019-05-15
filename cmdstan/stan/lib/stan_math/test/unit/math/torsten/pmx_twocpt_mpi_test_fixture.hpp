#ifndef TEST_UNIT_TORSTEN_PK_TWOCPT_MODEL_MPI_TEST_FIXTURE
#define TEST_UNIT_TORSTEN_PK_TWOCPT_MODEL_MPI_TEST_FIXTURE

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

class TorstenPopulationPMXTwoCptTest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
    torsten::mpi::Envionment::init();
  }
public:
  TorstenPopulationPMXTwoCptTest() :
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
    tlag{ { 0, 0, 0 } },
    np(10),
    time_m   (np * nt),
    amt_m    (np * nt),
    rate_m   (np * nt),
    cmt_m    (np * nt),
    evid_m   (np * nt),
    ii_m     (np * nt),
    addl_m   (np * nt),
    ss_m     (np * nt),
    pMatrix_m(np * pMatrix.size()),
    biovar_m (np * biovar.size()),
    tlag_m   (np * tlag.size()),
    len(np),
    len_pMatrix(np),
    len_biovar(np),
    len_tlag(np)
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

    // population data
    for (int i = 0; i < np; ++i) {
      for (int j = 0; j < nt; ++j) {      
        time_m[i * nt + j] = time[j];
        amt_m [i * nt + j] = amt[j]; 
        rate_m[i * nt + j] = rate[j];
        cmt_m [i * nt + j] = cmt[j]; 
        evid_m[i * nt + j] = evid[j];
        ii_m  [i * nt + j] = ii[j];  
        addl_m[i * nt + j] = addl[j];
        ss_m  [i * nt + j] = ss[j];  
      }
    }

    // population param
    for (int i = 0; i < np; ++i) {
      pMatrix_m[i] = pMatrix[0];
      biovar_m[i] = biovar[0];
      tlag_m[i] = tlag[0];
    }

    // population lenght
    for (int i = 0; i < np; ++i) {
      len[i] = nt;
      len_pMatrix[i] = 1;
      len_biovar[i] = 1;
      len_tlag[i] = 1;
    }
  }

  const int nt;
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
  const int np;
  std::vector<double> time_m   ;
  std::vector<double> amt_m    ;
  std::vector<double> rate_m   ;
  std::vector<int   > cmt_m    ;
  std::vector<int   > evid_m   ;
  std::vector<double> ii_m     ;
  std::vector<int   > addl_m   ;
  std::vector<int   > ss_m     ;
  std::vector<std::vector<double> > pMatrix_m;
  std::vector<std::vector<double> > biovar_m ;
  std::vector<std::vector<double> > tlag_m   ;
  std::vector<int> len;
  std::vector<int> len_pMatrix;
  std::vector<int> len_biovar;
  std::vector<int> len_tlag;
};

#endif
