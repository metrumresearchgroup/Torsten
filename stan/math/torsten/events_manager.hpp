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

  template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    struct EventsManager<NONMENEventsRecord<T0, T1, T2, T3, T4, T5, T6> > {
    using ER = NONMENEventsRecord<T0, T1, T2, T3, T4, T5, T6>;
    using T_scalar = typename ER::T_scalar;
    using T_time   = typename ER::T_time;
    using T_rate   = typename ER::T_rate;
    using T_amt    = typename ER::T_amt;
    using T_par    = typename ER::T_par;
    using T_par_rate = typename ER::T_par_rate;
    using T_par_ii   = typename ER::T_par_ii;

    EventHistory<T_time, T1, T2, T3> event_his;
    std::vector<std::vector<T_rate> > rate_v;
    std::vector<T_amt> amt_v;
    std::vector<std::vector<T_par> > par_v;

    const int nKeep;
    const int ncmt;

    static int nCmt(const ER& rec) {
      return rec.ncmt;
    }

    static int population_size(const ER& rec) {
      return rec.len_.size();
    }

    static int solution_size(const ER& rec) {
      return solution_size(0, rec);
    }

    static int solution_size(int id, const ER& rec) {
      return rec.len_[id];
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
      event_his(rec.begin_[id], rec.len_[id], rec.time_, rec.amt_, rec.rate_, rec.ii_, rec.evid_, rec.cmt_, rec.addl_, rec.ss_),
      nKeep(event_his.size()),
      ncmt(rec.ncmt)
    {
      event_his.Sort();

      ModelParameterHistory<T_time, T4, T5, T6>
        param_his(rec.begin_[id], rec.len_[id], rec.time_,
                  ibegin_theta, isize_theta, rec.pMatrix_, rec.systems_,
                  ibegin_biovar, isize_biovar, rec.biovar_,
                  ibegin_tlag, isize_tlag, rec.tlag_);
      param_his.Sort();

      event_his.AddlDoseEvents();
      param_his.CompleteParameterHistory(event_his);

      event_his.AddLagTimes(param_his, ncmt);
      RateHistory<T_time, T2> rate_history(event_his, ncmt);
      param_his.CompleteParameterHistory(event_his);

      int iRate = 0;
      for (size_t i = 0; i < event_his.size(); i++) {

        // Use index iRate instead of i to find rate at matching time, given there
        // is one rate per time, not per event.
        if (rate_history.time(iRate) != event_his.time(i)) iRate++;
        std::vector<T_rate> rate_i(ncmt);
        for (int j = 0; j < ncmt; ++j) {
          rate_i[j] = rate_history.rate(iRate, j) * param_his.GetValueBio(i, j);
        }
        rate_v.push_back(rate_i);

        amt_v.push_back(param_his.GetValueBio(i, event_his.cmt(i) - 1) * event_his.amt(i));
      }

      par_v.resize(event_his.size());
      for (size_t i = 0; i < event_his.size(); ++i) {
        auto p = param_his.GetModelParameters(i);
        par_v[i] = p.get_RealParameters(param_his.has_matrix_param);
      }
    }

    const EventHistory<T_time, T1, T2, T3>& events() const {
      return event_his;
    }

    const std::vector<std::vector<T_rate> >& rates() const {
      return rate_v;
    }

    const std::vector<T_amt>& amts() const {
      return amt_v;
    }

    const std::vector<std::vector<T_par> >& pars() const {
      return par_v;
    }

    /*
     * check the exisitence of SS dosing events
     */
    static bool has_ss_dosing(int id, const ER& rec) {
      return rec.has_ss_dosing(id);
    }

    /*
     * number of parameter for a given model is supposed to
     * be constant across the population.
     */
    static int parameter_size(const ER& rec) {
      return rec.pMatrix_.empty() ? rec.systems_[0].size() : rec.pMatrix_[0].size();
    }

    static int num_events(int id, const ER& rec) {
      if (rec.pMatrix_.empty()) {
        // TODO
        return 0;
      } else {
        return num_events(rec.begin_[id], rec.len_[id],
                       rec.time_, rec.amt_, rec.rate_, rec.ii_, rec.evid_, rec.cmt_, rec.addl_, rec.ss_,
                       rec.begin_param(id), rec.len_param(id), rec.pMatrix_, 
                       rec.begin_biovar(id), rec.len_biovar(id), rec.biovar_, 
                       rec.begin_tlag(id), rec.len_tlag(id), rec.tlag_);        
      }
    }

    /*
     * calculate the total nb. of events without generating
     * events history.
     */
    template <typename T0_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_, typename T6_>
    static int num_events(const std::vector<T0_>& time,
                       const std::vector<T1_>& amt,
                       const std::vector<T2_>& rate,
                       const std::vector<T3_>& ii,
                       const std::vector<int>& evid,
                       const std::vector<int>& cmt,
                       const std::vector<int>& addl,
                       const std::vector<int>& ss,
                       const std::vector<std::vector<T4_> >& pMatrix,
                       const std::vector<std::vector<T5_> >& biovar,
                       const std::vector<std::vector<T6_> >& tlag) {
      return num_events(0, time.size(), time, amt, rate, ii, evid, cmt, addl, ss,
                     0, pMatrix.size(), pMatrix,
                     0, biovar.size(), biovar,
                     0, tlag.size(), tlag);
    }

    template <typename T0_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_, typename T6_>
    static int num_events(int ibegin, int isize,
                       const std::vector<T0_>& time,
                       const std::vector<T1_>& amt,
                       const std::vector<T2_>& rate,
                       const std::vector<T3_>& ii,
                       const std::vector<int>& evid,
                       const std::vector<int>& cmt,
                       const std::vector<int>& addl,
                       const std::vector<int>& ss,
                       int ibegin_pMatrix, int isize_pMatrix,
                       const std::vector<std::vector<T4_> >& pMatrix,
                       int ibegin_biovar, int isize_biovar,
                       const std::vector<std::vector<T5_> >& biovar,
                       int ibegin_tlag, int isize_tlag,
                       const std::vector<std::vector<T6_> >& tlag) {
      using stan::math::value_of;

      int res;
      bool has_lag = std::any_of(tlag.begin() + ibegin_tlag, tlag.begin() + ibegin_tlag + isize_tlag,
                                 [](const std::vector<T6_>& v) {
                                   return std::any_of(v.begin(), v.end(), [](const T6_& x) { return std::abs(value_of(x)) > 1.E-10; });
                                 });

      if (!has_lag) {
        int n = isize;
        for (int i = ibegin; i < ibegin + isize; ++i) {
          if (evid[i] == 1 || evid[i] == 4) {      // is dosing event
            if (addl[i] > 0 && ii[i] > 0) {        // has addl doses
              if (rate[i] > 0 && amt[i] > 0) {
                n++;                               // end event for original IV dose
                n += 2 * addl[i];                  // end event for addl IV dose
              } else {
                n += addl[i];
              }
            } else if (rate[i] > 0 && amt[i] > 0) {
              n++;                                 // end event for IV dose
            }
          }
        }
        res = n;
      } else if (isize_tlag == 1) {
        int n = isize;
        std::vector<std::tuple<double, int>> dose;
        dose.reserve(isize);
        for (int i = ibegin; i < ibegin + isize; ++i) {
          if (evid[i] == 1 || evid[i] == 4) {      // is dosing event
            if (tlag[ibegin_tlag][cmt[i] - 1] > 0.0) {       // tlag dose
              n++;
            }
            if (addl[i] > 0 && ii[i] > 0) {        // has addl doses
              if (rate[i] > 0 && amt[i] > 0) {
                n++;                               // end ev for IV dose
                n += 2 * addl[i];                  // end ev for addl IV dose
              } else {
                n += addl[i];
              }
              if (tlag[ibegin_tlag][cmt[i] - 1] > 0.0) {     // tlag dose
                n += addl[i];
              }
            } else if (rate[i] > 0 && amt[i] > 0) {
              n++;                                 // end event for IV dose
            }
          }
        }
        res = n;
      } else {
        // int n = time.size();
        // std::vector<std::tuple<double, int>> addl_dose;
        // for (size_t i = 0; i < time.size(); i++) {
        //   if (evid[i] == 1 || evid[i] == 4) {      // is dosing event
        //     if (tlag[0][cmt[i] - 1] > 0.0) {       // tlag dose
        //       n++;
        //     }
        //     if (addl[i] > 0 && ii[i] > 0) {        // has addl doses
        //       if (rate[i] > 0 && amt[i] > 0) {
        //         n++;                               // end ev for IV dose
        //         n += 2 * addl[i];                  // end ev for addl IV dose
        //       } else {
        //         n += addl[i];
        //       }
        //       for (int j = 0; j < addl[i]; ++j) {
        //         addl_dose.push_back(std::make_tuple(value_of(time[i]) + (1+j) * value_of(ii[i]), cmt[i]));
        //       }
        //     } else if (rate[i] > 0 && amt[i] > 0) {
        //       n++;                                 // end event for IV dose
        //     }
        //   }
        // }
        // std::sort(addl_dose.begin(), addl_dose.end(),
        //           [](std::tuple<double, int>& a, std::tuple<double, int>& b)
        //           {
        //             return std::get<0>(a) < std::get<0>(b);
        //           });
      }

      return res;
    }
  };

}

#endif
