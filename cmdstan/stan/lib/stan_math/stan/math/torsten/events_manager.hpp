#ifndef STAN_MATH_TORSTEN_EVENTS_MANAGER_HPP
#define STAN_MATH_TORSTEN_EVENTS_MANAGER_HPP

#include <stan/math/torsten/PKModel/PKModel.hpp>
#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <stan/math/torsten/events_record.hpp>
#include <boost/math/tools/promotion.hpp>
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace torsten {
  template <typename T_event_record>
  struct EventsManager;

  template <typename T0, typename T1, typename T2, typename T3, typename T4_container, typename T5, typename T6>
    struct EventsManager<NONMENEventsRecord<T0, T1, T2, T3, T4_container, T5, T6> > {
    using ER = NONMENEventsRecord<T0, T1, T2, T3, T4_container, T5, T6>;
    using T_scalar = typename ER::T_scalar;
    using T_time   = typename ER::T_time;
    using T_rate   = typename ER::T_rate;
    using T_amt    = typename ER::T_amt;
    using T_par    = typename ER::T_par;
    using T_par_rate = typename ER::T_par_rate;
    using T_par_ii   = typename ER::T_par_ii;
    using T4 = typename stan::math::value_type<T4_container>::type;

    EventHistory<T0, T1, T2, T3, T4_container, T5, T6> event_his;

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

    const EventHistory<T0, T1, T2, T3, T4_container, T5, T6>& events() const {
      return event_his;
    }

    /*
     * number of events for a sinlge individual
     */
    static int num_events(const ER& rec) {
      return num_events(0, rec);
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
