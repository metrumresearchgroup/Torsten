#include <stan/math/rev/core.hpp>
#include <stan/math/torsten/pmx_onecpt_model.hpp>
#include <stan/math/torsten/pmx_twocpt_model.hpp>
#include <stan/math/torsten/ev_manager.hpp>
#include <stan/math/torsten/test/unit/pmx_onecpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_twocpt_test_fixture.hpp>
#include <stan/math/torsten/test/unit/pmx_twocpt_mpi_test_fixture.hpp>
#include <stan/math/torsten/test/unit/test_macros.hpp>
#include <Eigen/Dense>
#include <Eigen/src/Core/NumTraits.h>
#include <gtest/gtest.h>
#include <vector>

using stan::math::var;
using torsten::EventsManager;
using torsten::NONMENEventsRecord;
using torsten::NonEventParameters;

TEST_F(TorstenOneCptTest, lag_time) {
  // 2 subjects
  nt = 10;
  resize(nt);

  // two subjects
  std::vector<int> len{6, 4};
  evid = {0, 1, 1, 1, 0, 1, 0, 1, 1, 1};
  cmt = std::vector<int>(nt, 1);
  cmt[3] = 2;
  ii = std::vector<double>(nt, 0.0);
  addl = std::vector<int>(nt, 0);
  time = {0, 12, 24, 24, 48, 36, 0, 6, 18, 30};
  ss = std::vector<int>(nt, 0);
  rate = std::vector<double>(nt, 0.0);
  amt = std::vector<double>{0, 1100.0, 1100.0, 900.0, 0.0, 1100.0,
    0, 800.0, 800.0, 800.0};

  tlag.resize(nt);
  for (auto& lag : tlag) {
    lag = std::vector<double>{0.0, 0.0};
  }

  tlag[2] = std::vector<double>{6.0, 6.0};
  tlag[8] = std::vector<double>{6.0, 6.0};

  const NONMENEventsRecord<double, double, double, double>
    events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss);

  using EM = EventsManager<NONMENEventsRecord<double, double, double, double>,
                           NonEventParameters<double, double, std::vector, std::tuple<double, double> >>;

  // subject 1
  EM em1(0, events_rec, pMatrix, biovar, tlag);
  auto ev1 = em1.events();

  // subject 2
  EM em2(1, events_rec, pMatrix, biovar, tlag);
  auto ev2 = em2.events();

  EXPECT_EQ(ev1.size(), len[0] + 1);
  EXPECT_EQ(ev2.size(), len[1] + 1);

  EXPECT_FLOAT_EQ(ev1.time(2), 24.0);
  EXPECT_FLOAT_EQ(ev1.time(3), 24.0);
  EXPECT_FLOAT_EQ(ev1.time(4), 30.0);
  EXPECT_EQ(ev1.evid(2), 9);
  EXPECT_EQ(ev1.evid(4), 1);
  EXPECT_EQ(ev1.keep(4), 0);

  EXPECT_FLOAT_EQ(ev2.time(2), 18.0);
  EXPECT_FLOAT_EQ(ev2.time(3), 24.0);
  EXPECT_EQ(ev2.evid(2), 9);
  EXPECT_EQ(ev2.evid(3), 1);
  EXPECT_EQ(ev2.keep(3), 0);

  // std::cout << "taki test: " << events_rec.has_positive_param(0, tlag) << "\n";
  int nev1 = EM::num_events(0, events_rec, pMatrix, biovar, tlag);
  int nev2 = EM::num_events(1, events_rec, pMatrix, biovar, tlag);
  EXPECT_EQ(nev1, len[0] + 1);
  EXPECT_EQ(nev2, len[1] + 1);
}

