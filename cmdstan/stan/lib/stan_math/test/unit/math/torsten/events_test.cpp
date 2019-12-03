#include <test/unit/math/torsten/pmx_onecpt_test_fixture.hpp>
#include <test/unit/math/torsten/pmx_twocpt_test_fixture.hpp>
#include <test/unit/math/torsten/pmx_twocpt_mpi_test_fixture.hpp>
#include <test/unit/math/torsten/test_util.hpp>
#include <Eigen/Dense>
#include <Eigen/src/Core/NumTraits.h>
#include <gtest/gtest.h>
#include <vector>

using std::vector;
using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::var;
using refactor::PKRec;
using torsten::EventsManager;
using torsten::NONMENEventsRecord;

TEST_F(TorstenOneCptTest, lag_time) {
  using stan::math::var;

  nt = 2;
  time.resize(nt);
  amt.resize(nt);
  rate.resize(nt);
  cmt.resize(nt);
  evid.resize(nt);
  ii.resize(nt);
  addl.resize(nt);
  ss.resize(nt);
  evid[0] = 1;
  cmt[0] = 1;
  ii[0] = 0;
  addl[0] = 0;
  time[0] = 0.0;
  tlag[0][0] = 1.5;

  time[1] = 2.5;

  const NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>
    events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EventsManager<NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double> >
    em(events_rec);
  auto ev = em.events();

  EXPECT_EQ(ev.num_state_times(), nt + 1);  

  EXPECT_FLOAT_EQ(ev.time(0), 0.0);
  EXPECT_FLOAT_EQ(ev.time(1), 1.5);
  EXPECT_FLOAT_EQ(ev.time(2), 2.5);
  EXPECT_EQ(ev.evid(0), 2);
  EXPECT_EQ(ev.evid(1), 1);
  EXPECT_EQ(ev.evid(2), 0);

  std::vector<std::vector<var> > tlag_v(1, stan::math::to_var(tlag[0]));
  const NONMENEventsRecord<double, double, double, double, std::vector<double>, double, var>
    events_rec_v(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag_v);
  EventsManager<NONMENEventsRecord<double, double, double, double, std::vector<double>, double, var> >
    em1(events_rec_v);
  auto ev1 = em1.events();
  stan::math::set_zero_all_adjoints();
  std::vector<double> g;
  ev1.time(1).grad(tlag_v[0], g);
  EXPECT_FLOAT_EQ(g[0], 1.0);
  EXPECT_FLOAT_EQ(g[1], 0.0);
}

TEST_F(TorstenTwoCptTest, events_addl) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 5;
  addl[3] = 3;

  {
    ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
    EM em(events_rec);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0]);
    EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );
  }

  {
    ii[3] = 4.0;
    ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
    EM em(events_rec);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0]);
    EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );
  }

  amt[3] = 400.0;
  evid[3] = 1;
  ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EM em(events_rec);
  auto ev = em.events();
  EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0] + addl[3]);
  EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );

  EXPECT_EQ(ev.time(0 ), 0    );  EXPECT_EQ(ev.amt(0 ), 1000);  EXPECT_EQ(ev.evid(0 ), 1);
  EXPECT_EQ(ev.time(1 ), 0.25 );  EXPECT_EQ(ev.amt(1 ), 0   );  EXPECT_EQ(ev.evid(1 ), 0);
  EXPECT_EQ(ev.time(2 ), 0.5  );  EXPECT_EQ(ev.amt(2 ), 0   );  EXPECT_EQ(ev.evid(2 ), 0);
  EXPECT_EQ(ev.time(3 ), 0.75 );  EXPECT_EQ(ev.amt(3 ), 400 );  EXPECT_EQ(ev.evid(3 ), 1);
  EXPECT_EQ(ev.time(4 ), 1    );  EXPECT_EQ(ev.amt(4 ), 0   );  EXPECT_EQ(ev.evid(4 ), 0);
  EXPECT_EQ(ev.time(5 ), 1.25 );  EXPECT_EQ(ev.amt(5 ), 0   );  EXPECT_EQ(ev.evid(5 ), 0);
  EXPECT_EQ(ev.time(6 ), 1.5  );  EXPECT_EQ(ev.amt(6 ), 0   );  EXPECT_EQ(ev.evid(6 ), 0);
  EXPECT_EQ(ev.time(7 ), 1.75 );  EXPECT_EQ(ev.amt(7 ), 0   );  EXPECT_EQ(ev.evid(7 ), 0);
  EXPECT_EQ(ev.time(8 ), 2    );  EXPECT_EQ(ev.amt(8 ), 0   );  EXPECT_EQ(ev.evid(8 ), 0);
  EXPECT_EQ(ev.time(9 ), 4    );  EXPECT_EQ(ev.amt(9 ), 0   );  EXPECT_EQ(ev.evid(9 ), 0);
  EXPECT_EQ(ev.time(10), 4.75 );  EXPECT_EQ(ev.amt(10), 400 );  EXPECT_EQ(ev.evid(10), 1);
  EXPECT_EQ(ev.time(11), 8.75 );  EXPECT_EQ(ev.amt(11), 400 );  EXPECT_EQ(ev.evid(11), 1);
  EXPECT_EQ(ev.time(12), 12   );  EXPECT_EQ(ev.amt(12), 1000);  EXPECT_EQ(ev.evid(12), 1);
  EXPECT_EQ(ev.time(13), 12.75);  EXPECT_EQ(ev.amt(13), 400 );  EXPECT_EQ(ev.evid(13), 1);
  EXPECT_EQ(ev.time(14), 24   );  EXPECT_EQ(ev.amt(14), 1000);  EXPECT_EQ(ev.evid(14), 1);
  EXPECT_EQ(ev.time(15), 36   );  EXPECT_EQ(ev.amt(15), 1000);  EXPECT_EQ(ev.evid(15), 1);
  EXPECT_EQ(ev.time(16), 48   );  EXPECT_EQ(ev.amt(16), 1000);  EXPECT_EQ(ev.evid(16), 1);
  EXPECT_EQ(ev.time(17), 60   );  EXPECT_EQ(ev.amt(17), 1000);  EXPECT_EQ(ev.evid(17), 1);
}

TEST_F(TorstenTwoCptTest, events_addl_singled_ragged_array) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 5;
  addl[3] = 3;

  {
    ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
    EM em(0, events_rec);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0]);
    EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );
  }

  {
    ii[3] = 4.0;
    ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
    EM em(0, events_rec);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0]);
    EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );
  }

  amt[3] = 400.0;
  evid[3] = 1;
  ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EM em(0, events_rec);
  auto ev = em.events();
  EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0] + addl[3]);
  EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );

  EXPECT_EQ(ev.time(0 ), 0    );  EXPECT_EQ(ev.amt(0 ), 1000);  EXPECT_EQ(ev.evid(0 ), 1);
  EXPECT_EQ(ev.time(1 ), 0.25 );  EXPECT_EQ(ev.amt(1 ), 0   );  EXPECT_EQ(ev.evid(1 ), 0);
  EXPECT_EQ(ev.time(2 ), 0.5  );  EXPECT_EQ(ev.amt(2 ), 0   );  EXPECT_EQ(ev.evid(2 ), 0);
  EXPECT_EQ(ev.time(3 ), 0.75 );  EXPECT_EQ(ev.amt(3 ), 400 );  EXPECT_EQ(ev.evid(3 ), 1);
  EXPECT_EQ(ev.time(4 ), 1    );  EXPECT_EQ(ev.amt(4 ), 0   );  EXPECT_EQ(ev.evid(4 ), 0);
  EXPECT_EQ(ev.time(5 ), 1.25 );  EXPECT_EQ(ev.amt(5 ), 0   );  EXPECT_EQ(ev.evid(5 ), 0);
  EXPECT_EQ(ev.time(6 ), 1.5  );  EXPECT_EQ(ev.amt(6 ), 0   );  EXPECT_EQ(ev.evid(6 ), 0);
  EXPECT_EQ(ev.time(7 ), 1.75 );  EXPECT_EQ(ev.amt(7 ), 0   );  EXPECT_EQ(ev.evid(7 ), 0);
  EXPECT_EQ(ev.time(8 ), 2    );  EXPECT_EQ(ev.amt(8 ), 0   );  EXPECT_EQ(ev.evid(8 ), 0);
  EXPECT_EQ(ev.time(9 ), 4    );  EXPECT_EQ(ev.amt(9 ), 0   );  EXPECT_EQ(ev.evid(9 ), 0);
  EXPECT_EQ(ev.time(10), 4.75 );  EXPECT_EQ(ev.amt(10), 400 );  EXPECT_EQ(ev.evid(10), 1);
  EXPECT_EQ(ev.time(11), 8.75 );  EXPECT_EQ(ev.amt(11), 400 );  EXPECT_EQ(ev.evid(11), 1);
  EXPECT_EQ(ev.time(12), 12   );  EXPECT_EQ(ev.amt(12), 1000);  EXPECT_EQ(ev.evid(12), 1);
  EXPECT_EQ(ev.time(13), 12.75);  EXPECT_EQ(ev.amt(13), 400 );  EXPECT_EQ(ev.evid(13), 1);
  EXPECT_EQ(ev.time(14), 24   );  EXPECT_EQ(ev.amt(14), 1000);  EXPECT_EQ(ev.evid(14), 1);
  EXPECT_EQ(ev.time(15), 36   );  EXPECT_EQ(ev.amt(15), 1000);  EXPECT_EQ(ev.evid(15), 1);
  EXPECT_EQ(ev.time(16), 48   );  EXPECT_EQ(ev.amt(16), 1000);  EXPECT_EQ(ev.evid(16), 1);
  EXPECT_EQ(ev.time(17), 60   );  EXPECT_EQ(ev.amt(17), 1000);  EXPECT_EQ(ev.evid(17), 1);
}

TEST_F(TorstenTwoCptTest, events_addl_multiple_identical_ragged_array) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 5;
  addl[3] = 3;
  ii[3] = 4.0;
  amt[3] = 400.0;
  evid[3] = 1;

  size_t n = time.size();
  time.resize (2 * n) ; std::copy_n(time.begin() , n, time.begin() + n);
  amt.resize  (2 * n) ; std::copy_n(amt.begin()  , n, amt.begin()  + n);
  rate.resize (2 * n) ; std::copy_n(rate.begin() , n, rate.begin() + n);
  ii.resize   (2 * n) ; std::copy_n(ii.begin()   , n, ii.begin()   + n);
  evid.resize (2 * n) ; std::copy_n(evid.begin() , n, evid.begin() + n);
  cmt.resize  (2 * n) ; std::copy_n(cmt.begin()  , n, cmt.begin()  + n);
  addl.resize (2 * n) ; std::copy_n(addl.begin() , n, addl.begin() + n);
  ss.resize   (2 * n) ; std::copy_n(ss.begin()   , n, ss.begin()   + n);

  size_t n1 = pMatrix.size();
  size_t n2 = biovar.size();
  size_t n3 = tlag.size();
  pMatrix.resize(2 * n1) ; std::copy_n(pMatrix.begin(), n1, pMatrix.begin() + n1);
  biovar.resize (2 * n2) ; std::copy_n(biovar.begin(),  n2, biovar.begin()  + n2);
  tlag.resize   (2 * n3) ; std::copy_n(tlag.begin(),    n3, tlag.begin()    + n3);

  std::vector<int> len(2, n);
  ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  {
    EM em(0, events_rec, 0, n1, 0, n2, 0, n3);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), EM::num_events(0, events_rec) );
    EXPECT_EQ(ev.time(0 ), 0    );  EXPECT_EQ(ev.amt(0 ), 1000);  EXPECT_EQ(ev.evid(0 ), 1);
    EXPECT_EQ(ev.time(1 ), 0.25 );  EXPECT_EQ(ev.amt(1 ), 0   );  EXPECT_EQ(ev.evid(1 ), 0);
    EXPECT_EQ(ev.time(2 ), 0.5  );  EXPECT_EQ(ev.amt(2 ), 0   );  EXPECT_EQ(ev.evid(2 ), 0);
    EXPECT_EQ(ev.time(3 ), 0.75 );  EXPECT_EQ(ev.amt(3 ), 400 );  EXPECT_EQ(ev.evid(3 ), 1);
    EXPECT_EQ(ev.time(4 ), 1    );  EXPECT_EQ(ev.amt(4 ), 0   );  EXPECT_EQ(ev.evid(4 ), 0);
    EXPECT_EQ(ev.time(5 ), 1.25 );  EXPECT_EQ(ev.amt(5 ), 0   );  EXPECT_EQ(ev.evid(5 ), 0);
    EXPECT_EQ(ev.time(6 ), 1.5  );  EXPECT_EQ(ev.amt(6 ), 0   );  EXPECT_EQ(ev.evid(6 ), 0);
    EXPECT_EQ(ev.time(7 ), 1.75 );  EXPECT_EQ(ev.amt(7 ), 0   );  EXPECT_EQ(ev.evid(7 ), 0);
    EXPECT_EQ(ev.time(8 ), 2    );  EXPECT_EQ(ev.amt(8 ), 0   );  EXPECT_EQ(ev.evid(8 ), 0);
    EXPECT_EQ(ev.time(9 ), 4    );  EXPECT_EQ(ev.amt(9 ), 0   );  EXPECT_EQ(ev.evid(9 ), 0);
    EXPECT_EQ(ev.time(10), 4.75 );  EXPECT_EQ(ev.amt(10), 400 );  EXPECT_EQ(ev.evid(10), 1);
    EXPECT_EQ(ev.time(11), 8.75 );  EXPECT_EQ(ev.amt(11), 400 );  EXPECT_EQ(ev.evid(11), 1);
    EXPECT_EQ(ev.time(12), 12   );  EXPECT_EQ(ev.amt(12), 1000);  EXPECT_EQ(ev.evid(12), 1);
    EXPECT_EQ(ev.time(13), 12.75);  EXPECT_EQ(ev.amt(13), 400 );  EXPECT_EQ(ev.evid(13), 1);
    EXPECT_EQ(ev.time(14), 24   );  EXPECT_EQ(ev.amt(14), 1000);  EXPECT_EQ(ev.evid(14), 1);
    EXPECT_EQ(ev.time(15), 36   );  EXPECT_EQ(ev.amt(15), 1000);  EXPECT_EQ(ev.evid(15), 1);
    EXPECT_EQ(ev.time(16), 48   );  EXPECT_EQ(ev.amt(16), 1000);  EXPECT_EQ(ev.evid(16), 1);
    EXPECT_EQ(ev.time(17), 60   );  EXPECT_EQ(ev.amt(17), 1000);  EXPECT_EQ(ev.evid(17), 1);    
  }

  {
    EM em(1, events_rec, n1, n1, n2, n2, n3, n3);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), EM::num_events(1, events_rec));
    EXPECT_EQ(ev.time(0 ), 0    );  EXPECT_EQ(ev.amt(0 ), 1000);  EXPECT_EQ(ev.evid(0 ), 1);
    EXPECT_EQ(ev.time(1 ), 0.25 );  EXPECT_EQ(ev.amt(1 ), 0   );  EXPECT_EQ(ev.evid(1 ), 0);
    EXPECT_EQ(ev.time(2 ), 0.5  );  EXPECT_EQ(ev.amt(2 ), 0   );  EXPECT_EQ(ev.evid(2 ), 0);
    EXPECT_EQ(ev.time(3 ), 0.75 );  EXPECT_EQ(ev.amt(3 ), 400 );  EXPECT_EQ(ev.evid(3 ), 1);
    EXPECT_EQ(ev.time(4 ), 1    );  EXPECT_EQ(ev.amt(4 ), 0   );  EXPECT_EQ(ev.evid(4 ), 0);
    EXPECT_EQ(ev.time(5 ), 1.25 );  EXPECT_EQ(ev.amt(5 ), 0   );  EXPECT_EQ(ev.evid(5 ), 0);
    EXPECT_EQ(ev.time(6 ), 1.5  );  EXPECT_EQ(ev.amt(6 ), 0   );  EXPECT_EQ(ev.evid(6 ), 0);
    EXPECT_EQ(ev.time(7 ), 1.75 );  EXPECT_EQ(ev.amt(7 ), 0   );  EXPECT_EQ(ev.evid(7 ), 0);
    EXPECT_EQ(ev.time(8 ), 2    );  EXPECT_EQ(ev.amt(8 ), 0   );  EXPECT_EQ(ev.evid(8 ), 0);
    EXPECT_EQ(ev.time(9 ), 4    );  EXPECT_EQ(ev.amt(9 ), 0   );  EXPECT_EQ(ev.evid(9 ), 0);
    EXPECT_EQ(ev.time(10), 4.75 );  EXPECT_EQ(ev.amt(10), 400 );  EXPECT_EQ(ev.evid(10), 1);
    EXPECT_EQ(ev.time(11), 8.75 );  EXPECT_EQ(ev.amt(11), 400 );  EXPECT_EQ(ev.evid(11), 1);
    EXPECT_EQ(ev.time(12), 12   );  EXPECT_EQ(ev.amt(12), 1000);  EXPECT_EQ(ev.evid(12), 1);
    EXPECT_EQ(ev.time(13), 12.75);  EXPECT_EQ(ev.amt(13), 400 );  EXPECT_EQ(ev.evid(13), 1);
    EXPECT_EQ(ev.time(14), 24   );  EXPECT_EQ(ev.amt(14), 1000);  EXPECT_EQ(ev.evid(14), 1);
    EXPECT_EQ(ev.time(15), 36   );  EXPECT_EQ(ev.amt(15), 1000);  EXPECT_EQ(ev.evid(15), 1);
    EXPECT_EQ(ev.time(16), 48   );  EXPECT_EQ(ev.amt(16), 1000);  EXPECT_EQ(ev.evid(16), 1);
    EXPECT_EQ(ev.time(17), 60   );  EXPECT_EQ(ev.amt(17), 1000);  EXPECT_EQ(ev.evid(17), 1);    
  }
}

TEST_F(TorstenTwoCptTest, events_addl_rate) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 5;
  addl[3] = 2;
  amt[3] = 1200.0;
  ii[3] = 5;
  cmt[3] = 2;
  rate[3] = 400;
  evid[3] = 1;

  ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EM em(events_rec);
  auto ev = em.events();

  /* each IV dose has an end event.*/
  EXPECT_EQ(ev.num_state_times(), evid.size() + addl[0] + addl[3] * 2 + 1);
  EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );

  EXPECT_FLOAT_EQ(ev.time(0 ), 0    );   EXPECT_EQ(ev.evid(0 ), 1);
  EXPECT_FLOAT_EQ(ev.time(1 ), 0.25 );   EXPECT_EQ(ev.evid(1 ), 0);
  EXPECT_FLOAT_EQ(ev.time(2 ), 0.5  );   EXPECT_EQ(ev.evid(2 ), 0);
  EXPECT_FLOAT_EQ(ev.time(3 ), 0.75 );   EXPECT_EQ(ev.evid(3 ), 1);
  EXPECT_FLOAT_EQ(ev.time(4 ), 1    );   EXPECT_EQ(ev.evid(4 ), 0);
  EXPECT_FLOAT_EQ(ev.time(5 ), 1.25 );   EXPECT_EQ(ev.evid(5 ), 0);
  EXPECT_FLOAT_EQ(ev.time(6 ), 1.5  );   EXPECT_EQ(ev.evid(6 ), 0);
  EXPECT_FLOAT_EQ(ev.time(7 ), 1.75 );   EXPECT_EQ(ev.evid(7 ), 0);
  EXPECT_FLOAT_EQ(ev.time(8 ), 2    );   EXPECT_EQ(ev.evid(8 ), 0);
  EXPECT_FLOAT_EQ(ev.time(9 ), 3.75 );   EXPECT_EQ(ev.evid(9 ), 2);
  EXPECT_FLOAT_EQ(ev.time(10), 4    );   EXPECT_EQ(ev.evid(10), 0);
  EXPECT_FLOAT_EQ(ev.time(11), 5.75 );   EXPECT_EQ(ev.evid(11), 1);
  EXPECT_FLOAT_EQ(ev.time(12), 8.75 );   EXPECT_EQ(ev.evid(12), 2);
  EXPECT_FLOAT_EQ(ev.time(13), 10.75);   EXPECT_EQ(ev.evid(13), 1);
  EXPECT_FLOAT_EQ(ev.time(14), 12   );   EXPECT_EQ(ev.evid(14), 1);
  EXPECT_FLOAT_EQ(ev.time(15), 13.75);   EXPECT_EQ(ev.evid(15), 2);
  EXPECT_FLOAT_EQ(ev.time(16), 24   );   EXPECT_EQ(ev.evid(16), 1);
  EXPECT_FLOAT_EQ(ev.time(17), 36   );   EXPECT_EQ(ev.evid(17), 1);
  EXPECT_FLOAT_EQ(ev.time(18), 48   );   EXPECT_EQ(ev.evid(18), 1);
  EXPECT_FLOAT_EQ(ev.time(19), 60   );   EXPECT_EQ(ev.evid(19), 1);

  EXPECT_FLOAT_EQ(ev.fractioned_rates(0 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(1 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(2 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(3 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(4 )[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(5 )[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(6 )[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(7 )[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(8 )[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(9 )[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(10)[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(11)[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(12)[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(13)[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(14)[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(15)[1], 400);
  EXPECT_FLOAT_EQ(ev.fractioned_rates(16)[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(17)[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(18)[1], 0  );
  EXPECT_FLOAT_EQ(ev.fractioned_rates(19)[1], 0  );
}

TEST_F(TorstenTwoCptTest, events_addl_rate_multiple_identical_ragged_array) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 5;
  addl[3] = 2;
  amt[3] = 1200.0;
  ii[3] = 5;
  cmt[3] = 2;
  rate[3] = 400;
  evid[3] = 1;

  size_t n = time.size();
  time.resize (2 * n) ; std::copy_n(time.begin() , n, time.begin() + n);
  amt.resize  (2 * n) ; std::copy_n(amt.begin()  , n, amt.begin()  + n);
  rate.resize (2 * n) ; std::copy_n(rate.begin() , n, rate.begin() + n);
  ii.resize   (2 * n) ; std::copy_n(ii.begin()   , n, ii.begin()   + n);
  evid.resize (2 * n) ; std::copy_n(evid.begin() , n, evid.begin() + n);
  cmt.resize  (2 * n) ; std::copy_n(cmt.begin()  , n, cmt.begin()  + n);
  addl.resize (2 * n) ; std::copy_n(addl.begin() , n, addl.begin() + n);
  ss.resize   (2 * n) ; std::copy_n(ss.begin()   , n, ss.begin()   + n);

  size_t n1 = pMatrix.size();
  size_t n2 = biovar.size();
  size_t n3 = tlag.size();
  pMatrix.resize(2 * n1) ; std::copy_n(pMatrix.begin(), n1, pMatrix.begin() + n1);
  biovar.resize (2 * n2) ; std::copy_n(biovar.begin(),  n2, biovar.begin()  + n2);
  tlag.resize   (2 * n3) ; std::copy_n(tlag.begin(),    n3, tlag.begin()    + n3);

  std::vector<int> len(2, n);
  ER events_rec(nCmt, len, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  {
    EM em(0, events_rec, 0, n1, 0, n2, 0, n3);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), EM::num_events(0, events_rec) );
    EXPECT_FLOAT_EQ(ev.time(0 ), 0    );   EXPECT_EQ(ev.evid(0 ), 1);
    EXPECT_FLOAT_EQ(ev.time(1 ), 0.25 );   EXPECT_EQ(ev.evid(1 ), 0);
    EXPECT_FLOAT_EQ(ev.time(2 ), 0.5  );   EXPECT_EQ(ev.evid(2 ), 0);
    EXPECT_FLOAT_EQ(ev.time(3 ), 0.75 );   EXPECT_EQ(ev.evid(3 ), 1);
    EXPECT_FLOAT_EQ(ev.time(4 ), 1    );   EXPECT_EQ(ev.evid(4 ), 0);
    EXPECT_FLOAT_EQ(ev.time(5 ), 1.25 );   EXPECT_EQ(ev.evid(5 ), 0);
    EXPECT_FLOAT_EQ(ev.time(6 ), 1.5  );   EXPECT_EQ(ev.evid(6 ), 0);
    EXPECT_FLOAT_EQ(ev.time(7 ), 1.75 );   EXPECT_EQ(ev.evid(7 ), 0);
    EXPECT_FLOAT_EQ(ev.time(8 ), 2    );   EXPECT_EQ(ev.evid(8 ), 0);
    EXPECT_FLOAT_EQ(ev.time(9 ), 3.75 );   EXPECT_EQ(ev.evid(9 ), 2);
    EXPECT_FLOAT_EQ(ev.time(10), 4    );   EXPECT_EQ(ev.evid(10), 0);
    EXPECT_FLOAT_EQ(ev.time(11), 5.75 );   EXPECT_EQ(ev.evid(11), 1);
    EXPECT_FLOAT_EQ(ev.time(12), 8.75 );   EXPECT_EQ(ev.evid(12), 2);
    EXPECT_FLOAT_EQ(ev.time(13), 10.75);   EXPECT_EQ(ev.evid(13), 1);
    EXPECT_FLOAT_EQ(ev.time(14), 12   );   EXPECT_EQ(ev.evid(14), 1);
    EXPECT_FLOAT_EQ(ev.time(15), 13.75);   EXPECT_EQ(ev.evid(15), 2);
    EXPECT_FLOAT_EQ(ev.time(16), 24   );   EXPECT_EQ(ev.evid(16), 1);
    EXPECT_FLOAT_EQ(ev.time(17), 36   );   EXPECT_EQ(ev.evid(17), 1);
    EXPECT_FLOAT_EQ(ev.time(18), 48   );   EXPECT_EQ(ev.evid(18), 1);
    EXPECT_FLOAT_EQ(ev.time(19), 60   );   EXPECT_EQ(ev.evid(19), 1);

    EXPECT_FLOAT_EQ(ev.fractioned_rates(0 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(1 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(2 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(3 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(4 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(5 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(6 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(7 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(8 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(9 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(10)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(11)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(12)[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(13)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(14)[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(15)[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(16)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(17)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(18)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(19)[1], 0  );
  }

{
    EM em(1, events_rec, n1, n1, n2, n2, n3, n3);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), EM::num_events(1, events_rec) );
    EXPECT_FLOAT_EQ(ev.time(0 ), 0    );   EXPECT_EQ(ev.evid(0 ), 1);
    EXPECT_FLOAT_EQ(ev.time(1 ), 0.25 );   EXPECT_EQ(ev.evid(1 ), 0);
    EXPECT_FLOAT_EQ(ev.time(2 ), 0.5  );   EXPECT_EQ(ev.evid(2 ), 0);
    EXPECT_FLOAT_EQ(ev.time(3 ), 0.75 );   EXPECT_EQ(ev.evid(3 ), 1);
    EXPECT_FLOAT_EQ(ev.time(4 ), 1    );   EXPECT_EQ(ev.evid(4 ), 0);
    EXPECT_FLOAT_EQ(ev.time(5 ), 1.25 );   EXPECT_EQ(ev.evid(5 ), 0);
    EXPECT_FLOAT_EQ(ev.time(6 ), 1.5  );   EXPECT_EQ(ev.evid(6 ), 0);
    EXPECT_FLOAT_EQ(ev.time(7 ), 1.75 );   EXPECT_EQ(ev.evid(7 ), 0);
    EXPECT_FLOAT_EQ(ev.time(8 ), 2    );   EXPECT_EQ(ev.evid(8 ), 0);
    EXPECT_FLOAT_EQ(ev.time(9 ), 3.75 );   EXPECT_EQ(ev.evid(9 ), 2);
    EXPECT_FLOAT_EQ(ev.time(10), 4    );   EXPECT_EQ(ev.evid(10), 0);
    EXPECT_FLOAT_EQ(ev.time(11), 5.75 );   EXPECT_EQ(ev.evid(11), 1);
    EXPECT_FLOAT_EQ(ev.time(12), 8.75 );   EXPECT_EQ(ev.evid(12), 2);
    EXPECT_FLOAT_EQ(ev.time(13), 10.75);   EXPECT_EQ(ev.evid(13), 1);
    EXPECT_FLOAT_EQ(ev.time(14), 12   );   EXPECT_EQ(ev.evid(14), 1);
    EXPECT_FLOAT_EQ(ev.time(15), 13.75);   EXPECT_EQ(ev.evid(15), 2);
    EXPECT_FLOAT_EQ(ev.time(16), 24   );   EXPECT_EQ(ev.evid(16), 1);
    EXPECT_FLOAT_EQ(ev.time(17), 36   );   EXPECT_EQ(ev.evid(17), 1);
    EXPECT_FLOAT_EQ(ev.time(18), 48   );   EXPECT_EQ(ev.evid(18), 1);
    EXPECT_FLOAT_EQ(ev.time(19), 60   );   EXPECT_EQ(ev.evid(19), 1);

    EXPECT_FLOAT_EQ(ev.fractioned_rates(0 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(1 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(2 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(3 )[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(4 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(5 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(6 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(7 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(8 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(9 )[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(10)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(11)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(12)[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(13)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(14)[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(15)[1], 400);
    EXPECT_FLOAT_EQ(ev.fractioned_rates(16)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(17)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(18)[1], 0  );
    EXPECT_FLOAT_EQ(ev.fractioned_rates(19)[1], 0  );
  }
}

TEST_F(TorstenTwoCptTest, events_addl_const_tlag) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 1;
  addl[3] = 2;
  amt[3] = 1200.0;
  ii[3] = 3;
  cmt[3] = 2;
  evid[3] = 4;

  // lag for cmt 1, this shouldn't affect the dose on cmt 2.
  tlag[0][0] = 0.25;
  {
    ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
    EM em(events_rec);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), time.size() + addl[0] * 2 + 1 + addl[3]);
    EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );

    EXPECT_FLOAT_EQ(ev.time(0 ), 0    ); EXPECT_FLOAT_EQ(ev.amt(0 ), 1000.);
    EXPECT_FLOAT_EQ(ev.time(1 ), 0.25 ); EXPECT_FLOAT_EQ(ev.amt(1 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(2 ), 0.25 ); EXPECT_FLOAT_EQ(ev.amt(2 ), 1000.);
    EXPECT_FLOAT_EQ(ev.time(3 ), 0.5  ); EXPECT_FLOAT_EQ(ev.amt(3 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(4 ), 0.75 ); EXPECT_FLOAT_EQ(ev.amt(4 ), 1200.);
    EXPECT_FLOAT_EQ(ev.time(5 ), 1    ); EXPECT_FLOAT_EQ(ev.amt(5 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(6 ), 1.25 ); EXPECT_FLOAT_EQ(ev.amt(6 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(7 ), 1.5  ); EXPECT_FLOAT_EQ(ev.amt(7 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(8 ), 1.75 ); EXPECT_FLOAT_EQ(ev.amt(8 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(9 ), 2    ); EXPECT_FLOAT_EQ(ev.amt(9 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(10), 3.75 ); EXPECT_FLOAT_EQ(ev.amt(10), 1200.);
    EXPECT_FLOAT_EQ(ev.time(11), 4    ); EXPECT_FLOAT_EQ(ev.amt(11), 0.   );
    EXPECT_FLOAT_EQ(ev.time(12), 6.75 ); EXPECT_FLOAT_EQ(ev.amt(12), 1200.);
    EXPECT_FLOAT_EQ(ev.time(13), 12   ); EXPECT_FLOAT_EQ(ev.amt(13), 1000.);
    EXPECT_FLOAT_EQ(ev.time(14), 12.25); EXPECT_FLOAT_EQ(ev.amt(14), 1000.);

    EXPECT_EQ(ev.cmt(0 ), 1);     EXPECT_EQ(ev.evid(0 ), 2);
    EXPECT_EQ(ev.cmt(1 ), 3);     EXPECT_EQ(ev.evid(1 ), 0);
    EXPECT_EQ(ev.cmt(2 ), 1);     EXPECT_EQ(ev.evid(2 ), 1);
    EXPECT_EQ(ev.cmt(3 ), 3);     EXPECT_EQ(ev.evid(3 ), 0);
    EXPECT_EQ(ev.cmt(4 ), 2);     EXPECT_EQ(ev.evid(4 ), 4);
    EXPECT_EQ(ev.cmt(5 ), 3);     EXPECT_EQ(ev.evid(5 ), 0);
    EXPECT_EQ(ev.cmt(6 ), 3);     EXPECT_EQ(ev.evid(6 ), 0);
    EXPECT_EQ(ev.cmt(7 ), 3);     EXPECT_EQ(ev.evid(7 ), 0);
    EXPECT_EQ(ev.cmt(8 ), 3);     EXPECT_EQ(ev.evid(8 ), 0);
    EXPECT_EQ(ev.cmt(9 ), 3);     EXPECT_EQ(ev.evid(9 ), 0);
    EXPECT_EQ(ev.cmt(10), 2);     EXPECT_EQ(ev.evid(10), 4);
    EXPECT_EQ(ev.cmt(11), 3);     EXPECT_EQ(ev.evid(11), 0);
    EXPECT_EQ(ev.cmt(12), 2);     EXPECT_EQ(ev.evid(12), 4);
    EXPECT_EQ(ev.cmt(13), 1);     EXPECT_EQ(ev.evid(13), 2);
    EXPECT_EQ(ev.cmt(14), 1);     EXPECT_EQ(ev.evid(14), 1);
  }

  // lag for cmt 2, this shouldn't affect the dose on cmt 1.
  tlag[0][0] = 0.0;
  tlag[0][1] = 0.2;
  {
    ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
    EM em(events_rec);
    auto ev = em.events();
    EXPECT_EQ(ev.num_state_times(), time.size() + addl[3] * 2 + 1 + addl[0]);
    EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );
    
    EXPECT_FLOAT_EQ(ev.time(0 ), 0   );   EXPECT_FLOAT_EQ(ev.amt(0 ), 1000.);
    EXPECT_FLOAT_EQ(ev.time(1 ), 0.25);   EXPECT_FLOAT_EQ(ev.amt(1 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(2 ), 0.5 );   EXPECT_FLOAT_EQ(ev.amt(2 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(3 ), 0.75);   EXPECT_FLOAT_EQ(ev.amt(3 ), 1200.);
    EXPECT_FLOAT_EQ(ev.time(4 ), 0.95);   EXPECT_FLOAT_EQ(ev.amt(4 ), 1200.);
    EXPECT_FLOAT_EQ(ev.time(5 ), 1   );   EXPECT_FLOAT_EQ(ev.amt(5 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(6 ), 1.25);   EXPECT_FLOAT_EQ(ev.amt(6 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(7 ), 1.5 );   EXPECT_FLOAT_EQ(ev.amt(7 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(8 ), 1.75);   EXPECT_FLOAT_EQ(ev.amt(8 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(9 ), 2   );   EXPECT_FLOAT_EQ(ev.amt(9 ), 0.   );
    EXPECT_FLOAT_EQ(ev.time(10), 3.75);   EXPECT_FLOAT_EQ(ev.amt(10), 1200.);
    EXPECT_FLOAT_EQ(ev.time(11), 3.95);   EXPECT_FLOAT_EQ(ev.amt(11), 1200.);
    EXPECT_FLOAT_EQ(ev.time(12), 4   );   EXPECT_FLOAT_EQ(ev.amt(12), 0.   );
    EXPECT_FLOAT_EQ(ev.time(13), 6.75);   EXPECT_FLOAT_EQ(ev.amt(13), 1200.);
    EXPECT_FLOAT_EQ(ev.time(14), 6.95);   EXPECT_FLOAT_EQ(ev.amt(14), 1200.);
    EXPECT_FLOAT_EQ(ev.time(15), 12  );   EXPECT_FLOAT_EQ(ev.amt(15), 1000.);

    EXPECT_EQ(ev.cmt(0 ), 1);     EXPECT_EQ(ev.evid(0 ), 1);
    EXPECT_EQ(ev.cmt(1 ), 3);     EXPECT_EQ(ev.evid(1 ), 0);
    EXPECT_EQ(ev.cmt(2 ), 3);     EXPECT_EQ(ev.evid(2 ), 0);
    EXPECT_EQ(ev.cmt(3 ), 2);     EXPECT_EQ(ev.evid(3 ), 2);
    EXPECT_EQ(ev.cmt(4 ), 2);     EXPECT_EQ(ev.evid(4 ), 4);
    EXPECT_EQ(ev.cmt(5 ), 3);     EXPECT_EQ(ev.evid(5 ), 0);
    EXPECT_EQ(ev.cmt(6 ), 3);     EXPECT_EQ(ev.evid(6 ), 0);
    EXPECT_EQ(ev.cmt(7 ), 3);     EXPECT_EQ(ev.evid(7 ), 0);
    EXPECT_EQ(ev.cmt(8 ), 3);     EXPECT_EQ(ev.evid(8 ), 0);
    EXPECT_EQ(ev.cmt(9 ), 3);     EXPECT_EQ(ev.evid(9 ), 0);
    EXPECT_EQ(ev.cmt(10), 2);     EXPECT_EQ(ev.evid(10), 2);
    EXPECT_EQ(ev.cmt(11), 2);     EXPECT_EQ(ev.evid(11), 4);
    EXPECT_EQ(ev.cmt(12), 3);     EXPECT_EQ(ev.evid(12), 0);
    EXPECT_EQ(ev.cmt(13), 2);     EXPECT_EQ(ev.evid(13), 2);
    EXPECT_EQ(ev.cmt(14), 2);     EXPECT_EQ(ev.evid(14), 4);
    EXPECT_EQ(ev.cmt(15), 1);     EXPECT_EQ(ev.evid(15), 1);
  }
}

TEST_F(TorstenTwoCptTest, events_addl_rate_const_tlag) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 1;
  addl[3] = 2;
  amt[3] = 1200.0;
  ii[3] = 5;
  cmt[3] = 2;
  rate[3] = 400;
  evid[3] = 1;

  tlag[0][1] = 0.25;
  ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EM em(events_rec);
  auto ev = em.events();
  EXPECT_EQ(ev.num_state_times(), time.size() + addl[3] * 3 + 2 + addl[0]);
  EXPECT_EQ(ev.num_state_times(), EM::num_events(events_rec) );

  EXPECT_FLOAT_EQ(ev.time(0 ), 0    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(0 )[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(1 ), 0.25 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(1 )[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(2 ), 0.5  ); EXPECT_FLOAT_EQ(ev.fractioned_rates(2 )[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(3 ), 0.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(3 )[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(4 ), 1    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(4 )[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(5 ), 1    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(5 )[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(6 ), 1.25 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(6 )[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(7 ), 1.5  ); EXPECT_FLOAT_EQ(ev.fractioned_rates(7 )[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(8 ), 1.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(8 )[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(9 ), 2    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(9 )[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(10), 4    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(10)[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(11), 4    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(11)[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(12), 5.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(12)[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(13), 6    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(13)[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(14), 9    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(14)[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(15), 10.75); EXPECT_FLOAT_EQ(ev.fractioned_rates(15)[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(16), 11   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(16)[1], 0.  );
  EXPECT_FLOAT_EQ(ev.time(17), 12   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(17)[1], 400.);
  EXPECT_FLOAT_EQ(ev.time(18), 14   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(18)[1], 400.);

  EXPECT_EQ(ev.evid(0 ), 1);  EXPECT_EQ(ev.cmt(0 ), 1);
  EXPECT_EQ(ev.evid(1 ), 0);  EXPECT_EQ(ev.cmt(1 ), 3);
  EXPECT_EQ(ev.evid(2 ), 0);  EXPECT_EQ(ev.cmt(2 ), 3);
  EXPECT_EQ(ev.evid(3 ), 2);  EXPECT_EQ(ev.cmt(3 ), 2);
  EXPECT_EQ(ev.evid(4 ), 0);  EXPECT_EQ(ev.cmt(4 ), 3);
  EXPECT_EQ(ev.evid(5 ), 1);  EXPECT_EQ(ev.cmt(5 ), 2);
  EXPECT_EQ(ev.evid(6 ), 0);  EXPECT_EQ(ev.cmt(6 ), 3);
  EXPECT_EQ(ev.evid(7 ), 0);  EXPECT_EQ(ev.cmt(7 ), 3);
  EXPECT_EQ(ev.evid(8 ), 0);  EXPECT_EQ(ev.cmt(8 ), 3);
  EXPECT_EQ(ev.evid(9 ), 0);  EXPECT_EQ(ev.cmt(9 ), 3);
  EXPECT_EQ(ev.evid(10), 0);  EXPECT_EQ(ev.cmt(10), 3);
  EXPECT_EQ(ev.evid(11), 2);  EXPECT_EQ(ev.cmt(11), 2);
  EXPECT_EQ(ev.evid(12), 2);  EXPECT_EQ(ev.cmt(12), 2);
  EXPECT_EQ(ev.evid(13), 1);  EXPECT_EQ(ev.cmt(13), 2);
  EXPECT_EQ(ev.evid(14), 2);  EXPECT_EQ(ev.cmt(14), 2);
  EXPECT_EQ(ev.evid(15), 2);  EXPECT_EQ(ev.cmt(15), 2);
  EXPECT_EQ(ev.evid(16), 1);  EXPECT_EQ(ev.cmt(16), 2);
  EXPECT_EQ(ev.evid(17), 1);  EXPECT_EQ(ev.cmt(17), 1);
  EXPECT_EQ(ev.evid(18), 2);  EXPECT_EQ(ev.cmt(18), 2);
}

TEST_F(TorstenPopulationPMXTwoCptTest, events_addl_rate_const_tlag) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  int id = 1;
  addl_m[id * nt + 0] = 1;
  addl_m[id * nt + 3] = 2;
  amt_m[id * nt + 3] = 1200.0;
  ii_m[id * nt + 3] = 5;
  cmt_m[id * nt + 3] = 2;
  rate_m[id * nt + 3] = 400;
  evid_m[id * nt + 3] = 1;
  tlag_m[id * 1][1] = 0.25;

  id = 5;
  addl_m[id * nt + 3] = 1;
  amt_m[id * nt + 3] = 1000.0;
  ii_m[id * nt + 3] = 5;
  cmt_m[id * nt + 3] = 2;
  rate_m[id * nt + 3] = 100;
  evid_m[id * nt + 3] = 1;
  tlag_m[id * 1][1] = 0.12;

  std::vector<int> len(np, time.size());
  ER all_events_rec(nCmt, len, time_m, amt_m, rate_m, ii_m, evid_m, cmt_m, addl_m, ss_m, pMatrix_m, biovar_m, tlag_m);

  int ibegin = 0, ibegin_pMatrix = 0, ibegin_biovar = 0, ibegin_tlag = 0;
  int isize = 0, isize_pMatrix = 0, isize_biovar = 0, isize_tlag = 0;
  for (int id = 0; id < np; ++id) {
    ibegin += isize;
    ibegin_pMatrix += isize_pMatrix;
    ibegin_biovar += isize_biovar;
    ibegin_tlag += isize_tlag;
    isize = len[id];
    isize_pMatrix = len_pMatrix[id];
    isize_biovar = len_biovar[id];
    isize_tlag = len_tlag[id];
    
    std::vector<double> time_i(time_m.begin() + ibegin, time_m.begin() + ibegin + isize);
    std::vector<double> amt_i(amt_m.begin() + ibegin, amt_m.begin() + ibegin + isize);
    std::vector<double> rate_i(rate_m.begin() + ibegin, rate_m.begin() + ibegin + isize);
    std::vector<double> ii_i(ii_m.begin() + ibegin, ii_m.begin() + ibegin + isize);
    std::vector<int> evid_i(evid_m.begin() + ibegin, evid_m.begin() + ibegin + isize);
    std::vector<int> cmt_i(cmt_m.begin() + ibegin, cmt_m.begin() + ibegin + isize);
    std::vector<int> addl_i(addl_m.begin() + ibegin, addl_m.begin() + ibegin + isize);
    std::vector<int> ss_i(ss_m.begin() + ibegin, ss_m.begin() + ibegin + isize);
    std::vector<std::vector<double> > pMatrix_i(pMatrix_m.begin() + ibegin_pMatrix, pMatrix_m.begin() + ibegin_pMatrix + isize_pMatrix);
    std::vector<std::vector<double> > biovar_i(biovar_m.begin() + ibegin_biovar, biovar_m.begin() + ibegin_biovar + isize_biovar);
    std::vector<std::vector<double> > tlag_i(tlag_m.begin() + ibegin_tlag, tlag_m.begin() + ibegin_tlag + isize_tlag);

    ER events_rec(nCmt, time_i, amt_i, rate_i, ii_i, evid_i, cmt_i, addl_i, ss_i, pMatrix_i, biovar_i, tlag_i);
    EM em(events_rec);
    EM em_i(id, all_events_rec, ibegin_pMatrix, isize_pMatrix, ibegin_biovar, isize_biovar, ibegin_tlag, isize_tlag);
    EXPECT_EQ(em.events().num_state_times(), em_i.events().num_state_times());
    EXPECT_EQ(em.events().num_state_times(), EM::num_events(id, all_events_rec) );
    for (size_t j = 0; j < em.events().num_state_times(); ++j) {
      EXPECT_FLOAT_EQ(em.events().time(j), em_i.events().time(j));
      EXPECT_FLOAT_EQ(em.events().amt(j), em_i.events().amt(j));
      EXPECT_FLOAT_EQ(em.events().fractioned_rates(j)[1], em_i.events().fractioned_rates(j)[1]);
      EXPECT_EQ(em.events().evid(j), em_i.events().evid(j));
    }
  }
}

TEST_F(TorstenTwoCptTest, events_addl_rate_tlag) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 0;
  addl[3] = 2;
  amt[3] = 1200.0;
  ii[3] = 5;
  cmt[3] = 2;
  rate[3] = 400;
  evid[3] = 1;

  tlag.resize(time.size());
  for (auto& l : tlag) l.resize(nCmt);
  tlag[3][1] = 0.25;

  ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EM em(events_rec);
  auto ev = em.events();

  // Last Observation Carried Forward for the parameters for
  // the addl events. Since the last param has no lags, the
  // addl events have no lags.
  EXPECT_EQ(ev.num_state_times(), time.size() + addl[3] * 2 + 2);

  EXPECT_FLOAT_EQ(ev.time(0 ), 0    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(0 )[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(1 ), 0.25 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(1 )[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(2 ), 0.5  ); EXPECT_FLOAT_EQ(ev.fractioned_rates(2 )[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(3 ), 0.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(3 )[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(4 ), 1    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(4 )[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(5 ), 1    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(5 )[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(6 ), 1.25 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(6 )[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(7 ), 1.5  ); EXPECT_FLOAT_EQ(ev.fractioned_rates(7 )[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(8 ), 1.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(8 )[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(9 ), 2    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(9 )[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(10), 4    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(10)[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(11), 4    ); EXPECT_FLOAT_EQ(ev.fractioned_rates(11)[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(12), 5.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(12)[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(13), 8.75 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(13)[1],  400.0);
  EXPECT_FLOAT_EQ(ev.time(14), 10.75); EXPECT_FLOAT_EQ(ev.fractioned_rates(14)[1],    0.0);
  EXPECT_FLOAT_EQ(ev.time(15), 13.75); EXPECT_FLOAT_EQ(ev.fractioned_rates(15)[1],  400.0);
}

TEST_F(TorstenTwoCptTest, events_addl_rate_tlag_LOCF) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  addl[0] = 0;
  addl[3] = 2;
  amt[3] = 1200.0;
  ii[3] = 2.1;
  cmt[3] = 2;
  rate[3] = 400;
  evid[3] = 1;
  time[time.size() - 2] = 3.0;

  tlag.resize(time.size());
  for (auto& l : tlag) l.resize(nCmt);
  tlag[time.size() - 1][1] = 0.25;

  ER events_rec(nCmt, time, amt, rate, ii, evid, cmt, addl, ss, pMatrix, biovar, tlag);
  EM em(events_rec);
  auto ev = em.events();

  // Last Observation Carried Forward for the parameters for
  // the addl events. Since the last param has a lag, the
  // addl events that are after it have a lag. Since the
  // first of the two addl dosing event is at t = 2.85, and
  // it's a time that are not in original time points, the
  // parameter for this event is the one of the subsequent
  // event, which has no lag. Therefore the 1st addl event
  // has no lag, but the 2nd has.
  EXPECT_EQ(ev.num_state_times(), time.size() + 6);

  EXPECT_FLOAT_EQ(ev.time(0 ), 0   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(0 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.time(1 ), 0.25); EXPECT_FLOAT_EQ(ev.fractioned_rates(1 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.time(2 ), 0.5 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(2 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.time(3 ), 0.75); EXPECT_FLOAT_EQ(ev.fractioned_rates(3 )[1], 0  );
  EXPECT_FLOAT_EQ(ev.time(4 ), 1   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(4 )[1], 400);
  EXPECT_FLOAT_EQ(ev.time(5 ), 1.25); EXPECT_FLOAT_EQ(ev.fractioned_rates(5 )[1], 400);
  EXPECT_FLOAT_EQ(ev.time(6 ), 1.5 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(6 )[1], 400);
  EXPECT_FLOAT_EQ(ev.time(7 ), 1.75); EXPECT_FLOAT_EQ(ev.fractioned_rates(7 )[1], 400);
  EXPECT_FLOAT_EQ(ev.time(8 ), 2.85); EXPECT_FLOAT_EQ(ev.fractioned_rates(8 )[1], 400);
  EXPECT_FLOAT_EQ(ev.time(9 ), 3   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(9 )[1], 800);
  EXPECT_FLOAT_EQ(ev.time(10), 3.75); EXPECT_FLOAT_EQ(ev.fractioned_rates(10)[1], 800);
  EXPECT_FLOAT_EQ(ev.time(11), 4   ); EXPECT_FLOAT_EQ(ev.fractioned_rates(11)[1], 400);
  EXPECT_FLOAT_EQ(ev.time(12), 4.95); EXPECT_FLOAT_EQ(ev.fractioned_rates(12)[1], 400);
  EXPECT_FLOAT_EQ(ev.time(13), 5.2 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(13)[1], 400);
  EXPECT_FLOAT_EQ(ev.time(14), 5.85); EXPECT_FLOAT_EQ(ev.fractioned_rates(14)[1], 800);
  EXPECT_FLOAT_EQ(ev.time(15), 8.2 ); EXPECT_FLOAT_EQ(ev.fractioned_rates(15)[1], 400);

  EXPECT_EQ(ev.evid(0 ), 1);
  EXPECT_EQ(ev.evid(1 ), 0);
  EXPECT_EQ(ev.evid(2 ), 0);
  EXPECT_EQ(ev.evid(3 ), 1);
  EXPECT_EQ(ev.evid(4 ), 0);
  EXPECT_EQ(ev.evid(5 ), 0);
  EXPECT_EQ(ev.evid(6 ), 0);
  EXPECT_EQ(ev.evid(7 ), 0);
  EXPECT_EQ(ev.evid(8 ), 1);
  EXPECT_EQ(ev.evid(9 ), 0);
  EXPECT_EQ(ev.evid(10), 2);
  EXPECT_EQ(ev.evid(11), 0);
  EXPECT_EQ(ev.evid(12), 2);
  EXPECT_EQ(ev.evid(13), 1);
  EXPECT_EQ(ev.evid(14), 2);
  EXPECT_EQ(ev.evid(15), 2);
}

TEST_F(TorstenPopulationPMXTwoCptTest, events_addl_rate_tlag_LOCF) {
  using ER = NONMENEventsRecord<double, double, double, double, std::vector<double>, double, double>;
  using EM = EventsManager<ER>;

  int nCmt = refactor::PMXTwoCptModel<double, double, double, double>::Ncmt;
  int id = 1;
  addl_m[id * nt + 0] = 0;
  addl_m[id * nt + 3] = 2;
  amt_m[id * nt + 3] = 1200.0;
  ii_m[id * nt + 3] = 2.1;
  cmt_m[id * nt + 3] = 2;
  rate_m[id * nt + 3] = 400;
  evid_m[id * nt + 3] = 1;
  time_m[id * nt + time.size() - 2] = 3.0;

  id  = 5;
  addl_m[id * nt + 0] = 0;
  addl_m[id * nt + 3] = 2;
  amt_m[id * nt + 3] = 1000.0;
  ii_m[id * nt + 3] = 2.1;
  cmt_m[id * nt + 3] = 2;
  rate_m[id * nt + 3] = 400;
  evid_m[id * nt + 3] = 1;
  time_m[id * nt + time.size() - 2] = 3.0;

  // tlag now contains both constant data and time-dependent
  // data(for id 1 & 5)
  tlag_m.clear();
  for (int id = 0; id < np; ++id) {
    if (id == 1) {
      for (int j = 0; j < nt; ++j) {
        std::vector<double> tlag_i(nCmt, 0.0);
        if (j == len[id] - 1) tlag_i[1] = 0.25;
        tlag_m.push_back(tlag_i);
      }
      len_tlag[id] = nt;
    } else if (id == 5) {
      for (int j = 0; j < nt; ++j) {
        std::vector<double> tlag_i(nCmt, 0.0);
        if (j == len[id] - 1) tlag_i[1] = 0.13;
        tlag_m.push_back(tlag_i);
      }
      len_tlag[id] = nt;
    } else {
      std::vector<double> tlag_i(nCmt, 0.0);
      tlag_m.push_back(tlag_i);
      len_tlag[id] = 1;
    }
  }
  
  std::vector<int> len(np, time.size());
  ER all_events_rec(nCmt, len, time_m, amt_m, rate_m, ii_m, evid_m, cmt_m, addl_m, ss_m, pMatrix_m, biovar_m, tlag_m);

  int ibegin = 0, ibegin_pMatrix = 0, ibegin_biovar = 0, ibegin_tlag = 0;
  int isize = 0, isize_pMatrix = 0, isize_biovar = 0, isize_tlag = 0;
  for (int id = 0; id < np; ++id) {
    ibegin += isize;
    ibegin_pMatrix += isize_pMatrix;
    ibegin_biovar += isize_biovar;
    ibegin_tlag += isize_tlag;
    isize = len[id];
    isize_pMatrix = len_pMatrix[id];
    isize_biovar = len_biovar[id];
    isize_tlag = len_tlag[id];

    std::vector<double> time_i(time_m.begin() + ibegin, time_m.begin() + ibegin + isize);
    std::vector<double> amt_i(amt_m.begin() + ibegin, amt_m.begin() + ibegin + isize);
    std::vector<double> rate_i(rate_m.begin() + ibegin, rate_m.begin() + ibegin + isize);
    std::vector<double> ii_i(ii_m.begin() + ibegin, ii_m.begin() + ibegin + isize);
    std::vector<int> evid_i(evid_m.begin() + ibegin, evid_m.begin() + ibegin + isize);
    std::vector<int> cmt_i(cmt_m.begin() + ibegin, cmt_m.begin() + ibegin + isize);
    std::vector<int> addl_i(addl_m.begin() + ibegin, addl_m.begin() + ibegin + isize);
    std::vector<int> ss_i(ss_m.begin() + ibegin, ss_m.begin() + ibegin + isize);
    std::vector<std::vector<double> > pMatrix_i(pMatrix_m.begin() + ibegin_pMatrix, pMatrix_m.begin() + ibegin_pMatrix + isize_pMatrix);
    std::vector<std::vector<double> > biovar_i(biovar_m.begin() + ibegin_biovar, biovar_m.begin() + ibegin_biovar + isize_biovar);
    std::vector<std::vector<double> > tlag_i(tlag_m.begin() + ibegin_tlag, tlag_m.begin() + ibegin_tlag + isize_tlag);

    ER events_rec(nCmt, time_i, amt_i, rate_i, ii_i, evid_i, cmt_i, addl_i, ss_i, pMatrix_i, biovar_i, tlag_i);
    EM em(events_rec);
    EM em_i(id, all_events_rec, ibegin_pMatrix, isize_pMatrix, ibegin_biovar, isize_biovar, ibegin_tlag, isize_tlag);
    EXPECT_EQ(em.events().num_state_times(), em_i.events().num_state_times());
    for (size_t j = 0; j < em.events().num_state_times(); ++j) {
      EXPECT_FLOAT_EQ(em.events().time(j), em_i.events().time(j));
      EXPECT_FLOAT_EQ(em.events().fractioned_rates(j)[1], em_i.events().fractioned_rates(j)[1]);
      EXPECT_EQ(em.events().evid(j), em_i.events().evid(j));
    }
  }
}
