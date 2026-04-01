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
  template <typename T_event_record, typename T_params>
  struct EventsManager;

  /** 
   * Manage NMTRAN events by containing an <code>EventHistory</code> member
   * The class is effectively the interface of
   * <code>EventHistory</code>.
   * 
   * @tparam NONMENEventsRecord<...> currently only support NMTRAN records
   * @tparam NonEventParameters<...> non-event parameters such as ODE params.
   * 
   */
  template <typename T0, typename T1, typename T2, typename T3, typename T4,
            template<class...> class theta_container, typename... tuple_pars_t, typename... Ts>
  struct EventsManager<NONMENEventsRecord<T0, T1, T2, T3>,
                       NonEventParameters<T0, T4, theta_container, std::tuple<tuple_pars_t...>, Ts...> > {
    using param_t = NonEventParameters<T0, T4, theta_container, std::tuple<tuple_pars_t...>, Ts...>;
    using ER = NONMENEventsRecord<T0, T1, T2, T3>;
    using T_scalar = typename stan::return_type_t<typename ER::T_scalar, T4, tuple_pars_t..., Ts...>;
    using T_time   = typename stan::return_type_t<typename ER::T_time, typename param_t::lag_t>;
    using T_rate   = typename stan::return_type_t<typename ER::T_rate, typename param_t::biovar_t>;
    using T_amt    = typename stan::return_type_t<typename ER::T_amt, typename param_t::biovar_t>;
    using T_par    = T4;
    using T_par_rate = T2;
    using T_par_ii   = T3;

    param_t params;             /**< non-event parameters */
    EventHistory<T0, T1, T2, T3, typename param_t::lag_t> event_his;

    int nKeep;                  /**< number of original events */
    int ncmt;                   /**< number of compartments of the PMX model */

    /** 
     * number of compartments in the NMTRAN record
     * 
     * @param rec NMTRAN record
     * 
     * @return number of compartments
     */
    static int nCmt(const ER& rec) {
      return rec.ncmt;
    }

    /** 
     * the index in the result/input where subject @c id begins.
     * 
     * @param id subject id
     * @param rec NMTRAN record
     * 
     * @return the beginning index of the subject in the record
     */
    static int begin(int id, const ER& rec) {
      return rec.begin_.at(id);
    }

    /**
     * For population models, we need generate events using
     * ragged arrays.
     */
    template<typename... Tss>
    EventsManager(const ER& rec,
                  const std::vector<theta_container<T4>>& theta,
                  const std::vector<std::vector<Tss>>&... non_event_params) :
      EventsManager(0, rec, theta, non_event_params...)
    {}

    /**
     * For population models, we need generate events using
     * ragged arrays.
     */
    template<typename... Tss>
    EventsManager(int id, const ER& rec,
                  const std::vector<theta_container<T4>>& theta,
                  const std::vector<std::vector<Tss>>&... non_event_params) :
      params(id, rec, theta, non_event_params...),
      event_his(rec.ncmt, rec.begin_[id], rec.len_[id], rec.time_, rec.amt_, rec.rate_, rec.ii_, rec.evid_, rec.cmt_, rec.addl_, rec.ss_)
    {
      ncmt = rec.ncmt;

      attach_event_parameters();
      insert_lag_dose();
      event_his.generate_rates(ncmt);
      attach_event_parameters();

      nKeep = event_his.num_event_times;
    }

    template<typename... Tss>
    EventsManager(int id, const ER& rec,
                  int ibegin_theta, int isize_theta,
                  int ibegin_biovar, int isize_biovar,
                  int ibegin_tlag, int isize_tlag,
                  const std::vector<theta_container<T4>>& theta,
                  const std::vector<std::vector<Tss>>&... non_event_params) :
      params(id, rec, ibegin_theta, isize_theta, ibegin_biovar, isize_biovar, ibegin_tlag, isize_tlag, theta, non_event_params...),
      event_his(rec.ncmt, rec.begin_[id], rec.len_[id], rec.time_, rec.amt_, rec.rate_, rec.ii_, rec.evid_, rec.cmt_, rec.addl_, rec.ss_)
    {
      ncmt = rec.ncmt;

      attach_event_parameters();
      insert_lag_dose();
      event_his.generate_rates(ncmt);
      attach_event_parameters();

      nKeep = event_his.num_event_times;
    }

    /**
     * Implement absorption lag times by modifying the times of the dosing events.
     * Two cases: parameters are either constant or vary with each event.
     * Function sorts events at the end of the procedure.
     * The old event is set with a special EVID = 9 and it introduces no action.
     */
    void insert_lag_dose() {
      // reverse loop so we don't process same lagged events twice
      int nEvent = event_his.size();
      int iEvent = nEvent - 1;
      while (iEvent >= 0) {
        if (event_his.is_dosing(iEvent)) {
          if (params.lag_time(iEvent, event_his.cmt(iEvent) - 1) > ER::lag_time_min) {
            event_his.insert_event(iEvent);
            event_his.gen_time.back() += params.lag_time(iEvent, event_his.cmt(iEvent) - 1);
            event_his.idx[iEvent][2] = 9;
          }
        }
        iEvent--;
      }
      event_his.sort_state_time();
    }

    /** 
     * Get the @c EventHistory member
     * 
     * 
     * @return @c EventHistory
     */
    const EventHistory<T0, T1, T2, T3, typename param_t::lag_t>& events() const {
      return event_his;
    }

    /** 
     * When new events are created we need to attach its corresponding
     * parameters. 
     * 
     */
    void attach_event_parameters() {
      int nEvent = event_his.size();
      assert(nEvent > 0);
      int len_Parameters = params.size();  // numbers of events for which parameters are determined
      assert(len_Parameters > 0);

      if (!params.is_ordered()) params.sort();
      params.pars.resize(nEvent);

      int iEvent = 0;
      for (int i = 0; i < len_Parameters - 1; ++i) {
        while (event_his.isnew(iEvent)) iEvent++;  // skip new events
        assert(params.get_par_time(i) == event_his.time(iEvent));  // compare time of "old' events to time of parameters.
        iEvent++;
      }

      if (len_Parameters == 1)  {
        for (int i = 0; i < nEvent; ++i) {
          params.set_par_time(i, stan::math::value_of(event_his.time(i)));
          params.set_par_array(i, params.get_par_array(0));
          event_his.idx[i][3] = 0;
        }
      } else {  // parameters are event dependent.
        std::vector<double> times(nEvent, 0);
        for (int i = 0; i < nEvent; ++i) times[i] = params.pars[i].first;
        iEvent = 0;

        using par_t = typename param_t::par_t;
        par_t newParameter;
        int j = 0;
        typename std::vector<par_t>::const_iterator lower = params.pars.begin();
        typename std::vector<par_t>::const_iterator it_param_end = params.pars.begin() + len_Parameters;
        for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
          if (event_his.isnew(iEvent)) {
            // Find the index corresponding to the time of the new event in the
            // times vector.
            const double t = stan::math::value_of(event_his.time(iEvent));
            lower = std::lower_bound(lower, it_param_end, t,
                                     [](const par_t& t1, const double& t2) {return t1.first < t2;});
            newParameter = lower == (it_param_end) ? params.pars[len_Parameters-1] : *lower;
            newParameter.first = t;
            params.pars[len_Parameters + j] = newParameter;
            event_his.idx[iEvent][3] = 0; /**< item is "new" no more, set "isnew" false */
            j++;
          }
        }
      }
      params.sort();
    }

    /** 
     * Get model param for certain subject
     * 
     * @param i subject id
     * 
     * @return model param
     */
    inline const theta_container<T4>& theta(int i) const {
      return params.theta(i);
    }

    template<size_t Is>
    inline auto& get_model_array_1d_param(int i) const {
      return params.template get_model_array_1d_param<Is>(i);
    }

    /** 
     * Get dosing rate adjusted with bioavailability for a subject
     * 
     * @param i subject id
     * 
     * @return dosing rate
     */
    inline std::vector<T_rate> fractioned_rates(int i) const {
      const int n = event_his.rates[0].second.size();
      const std::vector<T2>& r = event_his.rates[event_his.rate_index[i]].second;
      std::vector<T_rate> res(r.size());
      for (size_t j = 0; j < r.size(); ++j) {
        res[j] = r[j] * params.bioavailability(i, j);
      }
      return res;
    }

    /** 
     * Get dosing amount adjusted with bioavailability for a subject
     * 
     * @param i subject id
     * 
     * @return dosing amount
     */
    inline T_amt fractioned_amt(int i) const {
      // return bioavailability(i, cmt(i) - 1) * amt(i);
      return params.bioavailability(i, event_his.cmt(i) - 1) * event_his.amt(i);
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
        amt(event_his.cmt(i) - 1) = fractioned_amt(i);
      }
      std::vector<T_rate> rate(fractioned_rates(i));
      return {id, t0, t1, event_his.ii(i),
              amt, rate, event_his.rate(i), event_his.cmt(i)};
    }

    /*
     * number of events for a sinlge individual
     */
    template<typename... Tss>
    static int num_events(const ER& rec,
                          const std::vector<theta_container<T4>>& theta,
                          const std::vector<std::vector<Tss>>&... non_event_params) {
      return num_events(0, rec, theta, non_event_params...);
    }

    static int num_events(int id, const ER& rec,
                          const std::vector<theta_container<T4>>& theta) {
      int res;
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
      return res;
    }
    
    template<typename T5>
    static int num_events(int id, const ER& rec,
                          const std::vector<theta_container<T4>>& theta,
                          const std::vector<std::vector<T5> >& biovar) {
      return num_events(id, rec, theta);
    }

    template<typename T5, typename T6, typename... Tss>
    static int num_events(int id, const ER& rec,
                          const std::vector<theta_container<T4>>& theta,
                          const std::vector<std::vector<T5> >& biovar,
                          const std::vector<std::vector<T6> >& tlag,
                          const std::vector<std::vector<Tss>>&... non_event_params) {
      int res;
      bool has_lag = rec.has_positive_param(id, tlag);

      if (!has_lag) {
        res = num_events(id, rec, theta);
      } else if (rec.len_param(id, tlag) == 1) {
        int n = rec.len_[id];
        std::vector<std::tuple<double, int>> dose;
        dose.reserve(rec.len_[id]);
        for (int i = rec.begin_[id]; i < rec.begin_[id] + rec.len_[id]; ++i) {
          if (rec.evid_[i] == 1 || rec.evid_[i] == 4) {      // is dosing event
            if (tlag[rec.begin_param(id, tlag)][rec.cmt_[i] - 1] > 0.0) {       // tlag dose
              n++;
            }
            if (rec.addl_[i] > 0 && rec.ii_[i] > 0) {        // has addl doses
              if (rec.rate_[i] > 0 && rec.amt_[i] > 0) {
                n++;                               // end ev for IV dose
                n += 2 * rec.addl_[i];                  // end ev for addl IV dose
              } else {
                n += rec.addl_[i];
              }
              if (tlag[rec.begin_param(id, tlag)][rec.cmt_[i] - 1] > 0.0) {     // tlag dose
                n += rec.addl_[i];
              }
            } else if (rec.rate_[i] > 0 && rec.amt_[i] > 0) {
              n++;                                 // end event for IV dose
            }
          }
        }
        res = n;
      } else {
        // FIXME not to use brute force
        res = EventsManager(id, rec, theta, biovar, tlag, non_event_params...).events().size();
      }

      return res;
    }
  };

}

#endif
