#ifndef STAN_MATH_TORSTEN_EVENTS_MANAGER_HPP
#define STAN_MATH_TORSTEN_EVENTS_MANAGER_HPP

#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/ev_history.hpp>
#include <stan/math/torsten/ev_record.hpp>
#include <stan/math/torsten/event.hpp>
#include <boost/math/tools/promotion.hpp>
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace torsten {
  template <typename T_event_record>
  struct EventsManager;

  template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
            template<class...> class theta_container>
  struct EventsManager<NONMENEventsRecord<T0, T1, T2, T3, T4, T5, T6,
                                          theta_container> > {
    using ER = NONMENEventsRecord<T0, T1, T2, T3, T4, T5, T6, theta_container>;
    using T_scalar = typename ER::T_scalar;
    using T_time   = typename ER::T_time;
    using T_rate   = typename ER::T_rate;
    using T_amt    = typename ER::T_amt;
    using T_par    = typename ER::T_par;
    using T_par_rate = typename ER::T_par_rate;
    using T_par_ii   = typename ER::T_par_ii;

    EventHistory<T0, T1, T2, T3, T4, T5, T6, theta_container> event_his;

    const int nKeep;
    const int ncmt;

    static int nCmt(const ER& rec) {
      return rec.ncmt;
    }

    /*
     * the index in the result/input where subject @c id begins.
     */
    static int begin(int id, const ER& rec) {
      return rec.begin_.at(id);
    }

    EventsManager(const ER& rec) : EventsManager(0, rec) {}

    EventsManager(int id, const ER& rec) :
      EventsManager(id, rec,
                    rec.begin_param(id), rec.len_param(id),
                    rec.begin_biovar(id), rec.len_biovar(id),
                    rec.begin_tlag(id), rec.len_tlag(id))
    {}

    /*
     * For population models, we need generate events using
     * ragged arrays.
     */
    EventsManager(int id, const ER& rec,
                  int ibegin_theta, int isize_theta,
                  int ibegin_biovar, int isize_biovar,
                  int ibegin_tlag, int isize_tlag) :
      event_his(rec.ncmt, rec.begin_[id], rec.len_[id], rec.time_, rec.amt_, rec.rate_, rec.ii_, rec.evid_, rec.cmt_, rec.addl_, rec.ss_,
                ibegin_theta, isize_theta, rec.pMatrix_,
                ibegin_biovar, isize_biovar, rec.biovar_,
                ibegin_tlag, isize_tlag, rec.tlag_),
      nKeep(event_his.num_event_times),
      ncmt(rec.ncmt)
    {}

    const EventHistory<T0, T1, T2, T3, T4, T5, T6, theta_container>& events() const {
      return event_his;
    }

    /*
     * number of events for a sinlge individual
     */
    static int num_events(const ER& rec) {
      return num_events(0, rec);
    }

    Event<T_time, T3, T_amt, T_rate, T2> event(int i) const {
      int id;
      switch (event_his.evid(i)) {
      case 2:                   // "other" type given "-cmt" indicates turn-off/reset
        if (event_his.cmt(i) < 0) {
            id = 5;
        } else {
          id = 0;
        }
        break;
      case 3:                   // reset
        id = 1;
        break;
      case 4:                   // reset + dosing
        if (event_his.is_ss_dosing(i)) {
          // since it's reset, SS reset is irrelevant 
          id = 4;
        } else {
          id = 2;
        }
        break;
      case 8:                   // mrgsolve: "evid=9" overwrite cmt
        if (event_his.cmt(i) > 0) {
          id = 6;
        }
        break;
      default:
        if (event_his.is_ss_dosing(i)) {
          if (event_his.ss(i) == 2) {
            id = 3;
          } else {
            id = 4;
          }
        } else {
          id = 0;
        }
      }

      T_time t0, t1;
      if (event_his.is_ss_dosing(i)) {
        // t0 = event_his.time(i);
        // t1 = event_his.ii(i);
        t0 = i == 0 ? event_his.time(0) : event_his.time(i-1);
        t1 = event_his.time(i);
      } else {
        t0 = i == 0 ? event_his.time(0) : event_his.time(i-1);
        t1 = event_his.time(i);
      }
      PKRec<T_amt> amt = PKRec<T_amt>::Zero(ncmt);
      if (event_his.is_bolus_dosing(i) || event_his.is_ss_dosing(i)) {
        amt(event_his.cmt(i) - 1) = event_his.fractioned_amt(i);
      }
      std::vector<T_rate> rate(event_his.fractioned_rates(i));
      return {id, t0, t1, event_his.ii(i),
              amt, rate, event_his.rate(i), event_his.cmt(i)};
    }

    /*
     * number of events for a sinlge individual when given a
     * population and individual id.
     */
    static int num_events(int id, const ER& rec) {
      using stan::math::value_of;

      int res;
      bool has_lag = rec.has_lag(id);

      if (!has_lag) {
        int n = rec.len_[id];
        for (int i = rec.begin_[id]; i < rec.begin_[id] + rec.len_[id]; ++i) {
          if (rec.evid_[i] == 1 || rec.evid_[i] == 4) {      // is dosing event
            if (rec.addl_[i] > 0 && rec.ii_[i] > 0) {        // has addl doses
              if (rec.rate_[i] > 0 && rec.amt_[i] > 0) {
                n++;                               // end event for original IV dose
                n += 2 * rec.addl_[i];                  // end event for addl IV dose
              } else {
                n += rec.addl_[i];
              }
            } else if (rec.rate_[i] > 0 && rec.amt_[i] > 0) {
              n++;                                 // end event for IV dose
            }
          }
        }
        res = n;
      } else if (rec.len_tlag(id) == 1) {
        int n = rec.len_[id];
        std::vector<std::tuple<double, int>> dose;
        dose.reserve(rec.len_[id]);
        for (int i = rec.begin_[id]; i < rec.begin_[id] + rec.len_[id]; ++i) {
          if (rec.evid_[i] == 1 || rec.evid_[i] == 4) {      // is dosing event
            if (rec.tlag_[rec.begin_tlag(id)][rec.cmt_[i] - 1] > 0.0) {       // tlag dose
              n++;
            }
            if (rec.addl_[i] > 0 && rec.ii_[i] > 0) {        // has addl doses
              if (rec.rate_[i] > 0 && rec.amt_[i] > 0) {
                n++;                               // end ev for IV dose
                n += 2 * rec.addl_[i];                  // end ev for addl IV dose
              } else {
                n += rec.addl_[i];
              }
              if (rec.tlag_[rec.begin_tlag(id)][rec.cmt_[i] - 1] > 0.0) {     // tlag dose
                n += rec.addl_[i];
              }
            } else if (rec.rate_[i] > 0 && rec.amt_[i] > 0) {
              n++;                                 // end event for IV dose
            }
          }
        }
        res = n;
      } else {
        // FIXME
        res = EventsManager(id, rec).events().num_state_times();
      }

      return res;
    }
  };

}

#endif
