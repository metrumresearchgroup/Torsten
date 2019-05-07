#ifndef STAN_MATH_TORSTEN_PK_NEUT_MPI_TEST_FIXTURE_HPP
#define STAN_MATH_TORSTEN_PK_NEUT_MPI_TEST_FIXTURE_HPP

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <test/unit/math/torsten/pmx_ode_test_fixture.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_fwd_system.hpp>
#include <stan/math/torsten/dsolve/pmx_cvodes_integrator.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/lorenz.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

class TorstenPopulationNeutropeniaTest : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
    torsten::mpi::Envionment::init();
  }
public:
  TorstenPopulationNeutropeniaTest() :
    f(),
    nCmt(8),
    nt(27),
    time{ 0.000,  0.000,  0.083,  0.167,  0.250,  0.500,  0.750,  1.000,  1.500,  2.000,
      3.000,  4.000,  6.000,  8.000, 12.000, 12.083, 12.167, 12.250, 12.500, 12.750,
      13.000, 13.500, 14.000, 15.000, 16.000, 18.000, 20.000},
    amt(nt, 0),
    rate(nt, 0),
    cmt(nt, 2),
    evid(nt, 0),
    ii(nt, 0),
    addl(nt, 0),
    ss(nt, 0),
    theta(1, {10.74924, 16.83211, 37.33329, 98.06352, 2.001674, 119.6034, 4.787322, 0.201094, 0.000278} ),
    biovar(1, std::vector<double>(nCmt, 1.0)),
    tlag(1, std::vector<double>(nCmt, 0.0)),
    np(8),
    time_m   (np * nt),
    amt_m    (np * nt),
    rate_m   (np * nt),
    cmt_m    (np * nt),
    evid_m   (np * nt),
    ii_m     (np * nt),
    addl_m   (np * nt),
    ss_m     (np * nt),
    theta_m  (np * theta.size()),
    biovar_m (np * biovar.size()),
    tlag_m   (np * tlag.size()),
    len(np)
  {
    SetUp();

    amt[1] = 80000;
    cmt[1] = 1;
    evid[1] = 1;
    ii[1] = 12;
    addl[1] = 1;

    setup_population();
  }

  void setup_population() {
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
      theta_m[i] = theta[0];
      biovar_m[i] = biovar[0];
      tlag_m[i] = tlag[0];
    }

    // population length
    for (int i = 0; i < np; ++i) {
      len[i] = nt;
    }    
  }

  /*
   * setup population given length of each individual,
   * so that individual data that are greater than
   * given length is modified by @c std::vector::resize()
   */
  void setup_population(const std::vector<int>& length) {
    time_m.clear();
    amt_m.clear();
    rate_m.clear();
    cmt_m.clear(); 
    evid_m.clear();
    ii_m.clear();
    addl_m.clear();
    ss_m.clear();
    theta_m.clear();
    biovar_m.clear();
    tlag_m.clear();

    // population data
    for (size_t i = 0; i < length.size(); ++i) {
      assert(length[i] <= nt);
      for (int j = 0; j < length[i]; ++j) {
        time_m.push_back(time[j]);
        amt_m .push_back(amt[j]);
        rate_m.push_back(rate[j]);
        cmt_m .push_back(cmt[j]);
        evid_m.push_back(evid[j]);
        ii_m  .push_back(ii[j]);
        addl_m.push_back(addl[j]);
        ss_m  .push_back(ss[j]);
      }
    }

    // population param
    for (size_t i = 0; i < length.size(); ++i) {
      for (int j = 0; j < length[i]; ++j) {
        theta_m.push_back(theta[0]);
        biovar_m.push_back(biovar[0]);
        tlag_m.push_back(tlag[0]);
      }
    }

    // population length
    len.resize(length.size());
    for (size_t i = 0; i < length.size(); ++i) {
      len[i] = length[i];
    }

    np = length.size();
  }

  const TwoCptNeutModelODE f;
  const int nCmt;
  const int nt;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> rate;
  std::vector<int> cmt;
  std::vector<int> evid;
  std::vector<double> ii;
  std::vector<int> addl;
  std::vector<int> ss;
  std::vector<std::vector<double> > theta;
  std::vector<std::vector<double> > biovar;
  std::vector<std::vector<double> > tlag;
  int np;
  std::vector<double> time_m   ;
  std::vector<double> amt_m    ;
  std::vector<double> rate_m   ;
  std::vector<int   > cmt_m    ;
  std::vector<int   > evid_m   ;
  std::vector<double> ii_m     ;
  std::vector<int   > addl_m   ;
  std::vector<int   > ss_m     ;
  std::vector<std::vector<double> > theta_m;
  std::vector<std::vector<double> > biovar_m ;
  std::vector<std::vector<double> > tlag_m   ;
  std::vector<int> len;
};

#endif
