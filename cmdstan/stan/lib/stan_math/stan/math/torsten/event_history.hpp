#ifndef STAN_MATH_TORSTEN_EVENT_HISTORY_HPP
#define STAN_MATH_TORSTEN_EVENT_HISTORY_HPP

#include <iomanip>
#include <stan/math/torsten/return_type.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/torsten/PKModel/functions.hpp>
#include <stan/math/torsten/pk_nsys.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <vector>

namespace torsten {

  /**
   * The EventHistory class defines objects that contain a vector of Events,
   * along with a series of functions that operate on them.
   */
  template<typename T0, typename T1, typename T2, typename T3, typename T4_container, typename T5, typename T6>
  struct EventHistory {
    using T4 = typename stan::math::value_type<T4_container>::type;
    using T_scalar = typename torsten::return_t<T0, T1, T2, T3, T4, T5, T6>::type;
    using T_time = typename torsten::return_t<T0, T1, T6, T2>::type;
    using T_rate = typename torsten::return_t<T2, T5>::type;
    using T_amt = typename torsten::return_t<T1, T5>::type;
    using Param = std::pair<double, std::array<int, 3> >;
    using rate_t = std::pair<double, std::vector<T2> >;

    const std::vector<T0>& time_;
    const std::vector<T1>& amt_;
    const std::vector<T2>& rate_;
    const std::vector<T3>& ii_;
    const std::vector<int>& evid_;
    const std::vector<int>& cmt_;
    const std::vector<int>& addl_;
    const std::vector<int>& ss_;

    std::vector<Param> param_index;
    const std::vector<T4_container>& theta_;
    const std::vector<std::vector<T5> >& biovar_;
    const std::vector<std::vector<T6> >& tlag_;

    // internally generated events
    const size_t num_event_times;
    std::vector<T_time> gen_time;
    std::vector<T1> gen_amt;
    std::vector<T2> gen_rate;
    std::vector<T3> gen_ii;
    std::vector<int> gen_cmt;
    std::vector<int> gen_addl;
    std::vector<int> gen_ss;

    using IDVec = std::array<int, 7>;
    // 0: original(0)/generated(1)
    // 1: index in original/generated arrays
    // 2: evid
    // 3: is new?(0/1)
    // 4: theta id
    // 5: biovar id
    // 6: tlag id
    std::vector<IDVec> index;

    // rate at distinct time
    std::vector<rate_t> rates;
    std::vector<int> rate_index;

    inline bool keep(const IDVec& id)  const { return id[0] == 0; }
    inline bool isnew(const IDVec& id) const { return id[3] == 1; }
    inline int evid (const IDVec& id) const { return id[2] ; }

    inline bool keep(int i)  const { return keep(index[i]); }
    inline bool isnew(int i) const { return isnew(index[i]); }
    inline int evid (int i) const { return evid(index[i]); }

    /*
     * for a population with data in ragged array form, we
     * form the events history using the population data and
     * the location of the individual in the ragged arrays.
     * In this constructor we assume @c p_ii.size() > 1 and
     * @c p_ss.size() > 1.
     */
    EventHistory(int ncmt, int ibegin, int isize,
                 const std::vector<T0>& p_time, const std::vector<T1>& p_amt,
                 const std::vector<T2>& p_rate, const std::vector<T3>& p_ii,
                 const std::vector<int>& p_evid, const std::vector<int>& p_cmt,
                 const std::vector<int>& p_addl, const std::vector<int>& p_ss,
                 int ibegin_theta, int isize_theta,
                 const std::vector<T4_container>& theta,
                 int ibegin_biovar, int isize_biovar,
                 const std::vector<std::vector<T5> >& biovar,
                 int ibegin_tlag, int isize_tlag,
                 const std::vector<std::vector<T6> >& tlag) :
      time_(p_time),
      amt_(p_amt),
      rate_(p_rate),
      ii_(p_ii),
      evid_(p_evid),
      cmt_(p_cmt),
      addl_(p_addl),
      ss_(p_ss),
      param_index(isize),
      theta_(theta),
      biovar_(biovar),
      tlag_(tlag),
      num_event_times(isize),
      index(isize, {0, 0, 0, 0, ibegin_theta, ibegin_biovar, ibegin_tlag})
    {
      const int iend = ibegin + isize;
      using stan::math::check_greater_or_equal;
      static const char* caller = "EventHistory::EventHistory";
      check_greater_or_equal(caller, "isize", isize , 1);
      check_greater_or_equal(caller, "time size", p_time.size() , size_t(iend));
      check_greater_or_equal(caller, "amt size", p_amt.size()   , size_t(iend));
      check_greater_or_equal(caller, "rate size", p_rate.size() , size_t(iend));
      check_greater_or_equal(caller, "ii size", p_ii.size()     , size_t(iend));
      check_greater_or_equal(caller, "evid size", p_evid.size() , size_t(iend));
      check_greater_or_equal(caller, "cmt size", p_cmt.size()   , size_t(iend));
      check_greater_or_equal(caller, "addl size", p_addl.size() , size_t(iend));
      check_greater_or_equal(caller, "ss size", p_ss.size()     , size_t(iend));
      for (size_t i = 0; i < isize; ++i) {
        index[i][1] = ibegin + i;
        index[i][2] = evid_[ibegin + i];
        index[i][4] = isize_theta   > 1 ? ibegin_theta  + i : ibegin_theta;
        index[i][5] = isize_biovar  > 1 ? ibegin_biovar + i : ibegin_biovar;
        index[i][6] = isize_tlag    > 1 ? ibegin_tlag   + i : ibegin_tlag;
      }
      insert_addl_dose();
      sort_state_time();

      for (int i = 0; i < isize; ++i) {
        int j = isize_theta   > 1 ? ibegin_theta  + i : ibegin_theta;
        int k = isize_biovar  > 1 ? ibegin_biovar + i : ibegin_biovar;
        int l = isize_tlag    > 1 ? ibegin_tlag   + i : ibegin_tlag;
        param_index[i] = std::make_pair<double, std::array<int, 3> >(stan::math::value_of(time_[ibegin + i]), {j, k, l });
      }
      param_sort();

      attach_event_parameters();
      insert_lag_dose();
      generate_rates(ncmt);
      attach_event_parameters();
    }

    EventHistory(int ncmt,
                 const std::vector<T0>& p_time, const std::vector<T1>& p_amt,
                 const std::vector<T2>& p_rate, const std::vector<T3>& p_ii,
                 const std::vector<int>& p_evid, const std::vector<int>& p_cmt,
                 const std::vector<int>& p_addl, const std::vector<int>& p_ss,
                 const std::vector<T4_container>& theta,
                 const std::vector<std::vector<T5> >& biovar,
                 const std::vector<std::vector<T6> >& tlag)
    : EventHistory(ncmt,
                   0, p_time.size(), p_time, p_amt, p_rate, p_ii, p_evid, p_cmt, p_addl, p_ss,
                   0, theta.size(), theta,
                   0, biovar.size(), biovar,
                   0, tlag.size(), tlag)
    {}

    void attach_event_parameters() {
      int nEvent = num_state_times();
      assert(nEvent > 0);
      int len_Parameters = param_index.size();  // numbers of events for which parameters are determined
      assert(len_Parameters > 0);

      if (!param_check()) param_sort();
      param_index.resize(nEvent);

      int iEvent = 0;
      for (int i = 0; i < len_Parameters - 1; i++) {
        while (isnew(iEvent)) iEvent++;  // skip new events
        assert(std::get<0>(param_index[i]) == time(iEvent));  // compare time of "old' events to time of parameters.
        iEvent++;
      }

      if (len_Parameters == 1)  {
        for (int i = 0; i < nEvent; i++) {
          param_index[i] = std::make_pair<double, std::array<int, 3> >(stan::math::value_of(time(i)) , std::array<int,3>(std::get<1>(param_index[0])));
          index[i][3] = 0;
        }
      } else {  // parameters are event dependent.
        std::vector<double> times(nEvent, 0);
        for (int i = 0; i < nEvent; i++) times[i] = param_index[i].first;
        iEvent = 0;

        Param newParameter;
        int j = 0;
        std::vector<Param>::const_iterator lower = param_index.begin();
        std::vector<Param>::const_iterator it_param_end = param_index.begin() + len_Parameters;
        for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
          if (isnew(iEvent)) {
            // Find the index corresponding to the time of the new event in the
            // times vector.
            const double t = stan::math::value_of(time(iEvent));
            lower = std::lower_bound(lower, it_param_end, t,
                                     [](const Param& t1, const double& t2) {return t1.first < t2;});
            newParameter = lower == (it_param_end) ? param_index[len_Parameters-1] : *lower;
            newParameter.first = t;
            param_index[len_Parameters + j] = newParameter;
            index[iEvent][3] = 0;
            j++;
          }
        }
      }
      param_sort();
    }

    bool is_reset(int i) const {
      return evid(i) == 3 || evid(i) == 4;
    }

    bool is_reset_event(const IDVec& id) {
      return evid(id) == 3 || evid(id) == 4;
    }

    bool is_dosing(int i) const {
      return evid(i) == 1 || evid(i) == 4;
    }

    /*
     * if an event is steady-state dosing event.
     */
    bool is_ss_dosing(int i) const {
      return (is_dosing(i) && (ss(i) == 1 || ss(i) == 2)) || ss(i) == 3;
    }

    static bool is_dosing(const std::vector<int>& event_id, int i) {
      return event_id[i] == 1 || event_id[i] == 4;
    }

    bool is_bolus_dosing(int i) const {
      const double eps = 1.0E-12;
      return is_dosing(i) && rate(i) < eps;
    }

    /*
     * use current event #i as template to @c push_back to
     * another event.
     */
    void insert_event(int i) {
      index.push_back({1, int(gen_time.size()), index[i][2], 1, index[i][4], index[i][5], index[i][6]});
      gen_time. push_back(time (i));
      gen_amt.  push_back(amt  (i));
      gen_rate. push_back(rate (i));
      gen_ii.   push_back(ii   (i));
      gen_cmt.  push_back(cmt  (i));
      gen_addl. push_back(addl (i));
      gen_ss.   push_back(ss   (i));
    }

    /**
     * Add events to EventHistory object, corresponding to additional dosing,
     * administered at specified inter-dose interval. This information is stored
     * in the addl and ii members of the EventHistory object.
     *
     * Events is sorted at the end of the procedure.
     */
    void insert_addl_dose() {
      for (int i = 0; i < num_state_times(); i++) {
        if (is_dosing(i) && ((addl(i) > 0) && (ii(i) > 0))) {
          for (int j = 1; j <= addl(i); j++) {
            insert_event(i);
            gen_time.back() += j * ii(i);
            gen_ii.back() = 0;
            gen_addl.back() = 0;
            gen_ss.back() = 0;
          }
        }
      }
    }

    /*
     * sort PMX events and nonevents times
     */
    void sort_state_time() { std::stable_sort(index.begin(), index.end(),
                                              [this](const IDVec &a, const IDVec &b) {
                                                using stan::math::value_of;
                                                double ta = keep(a) ? value_of(time_[a[1]]) : value_of(gen_time[a[1]]);
                                                double tb = keep(b) ? value_of(time_[b[1]]) : value_of(gen_time[b[1]]);
                                                return ta < tb;
                                              });
    }

    /*
     * return if an event is a "reset" event(evid=3 or 4)
     */
    void param_sort() {
      std::sort(param_index.begin(), param_index.end(),
                [](const Param& a, const Param& b)
                { return a.first < b.first; });
    }

    bool param_check() {
      // check that elements are in chronological order.
      int i = param_index.size() - 1;
      bool ordered = true;

      while (i > 0 && ordered) {
        ordered = (param_index[i].first >= param_index[i-1].first);
        i--;
      }
      return ordered;
    }

    void generate_rates(int nCmt) {
      using std::vector;
      using stan::math::value_of;

      const int n = num_state_times();
      for (size_t i = 0; i < n; ++i) {
        if ((is_dosing(i)) && (rate(i) > 0 && amt(i) > 0)) {
          insert_event(i);
          index.back()[2] = 2;    // reset evid
          gen_time. back() += amt(i)/rate(i);
          gen_amt.  back() = 0;
          gen_rate. back() = 0;
          gen_ii.   back() = 0;
          gen_addl. back() = 0;
          gen_ss.   back() = 0;
        }
      }
      sort_state_time();

      rate_t newRate{0.0, std::vector<T2>(nCmt, 0.0)};
      // unique_times is sorted
      std::vector<int> ut(unique_times());
      for (auto i : ut) {
        newRate.first = value_of(time(i));
        rates.push_back(newRate);
      }
      // check sorted?
      std::sort(rates.begin(), rates.end(), [](rate_t const &a, rate_t const &b) {
          return a.first < b.first;
        });

      for (size_t i = 0; i < num_state_times(); ++i) {
        if ((is_dosing(i)) && (rate(i) > 0 && amt(i) > 0)) {
          double t0 = value_of(time(i));
          double t1 = t0 + value_of(amt(i)/rate(i));
          for (auto&& r : rates) {
            if (r.first > t0 && r.first <= t1) {
              r.second[cmt(i) - 1] += rate(i);
            }
          }
        }
      }

      /*
       * rate index points to the rates for each state time,
       * since there is one rates vector per time, not per event.
       */
      rate_index.resize(index.size());
      int iRate = 0;
      for (size_t i = 0; i < index.size(); ++i) {
        if (rates[iRate].first != time(i)) iRate++;
        rate_index[i] = iRate;
      }
    }

    // Access functions
    inline T_time time (const IDVec& id) const { return keep(id) ? time_[id[1]] : gen_time[id[1]] ; }
    inline const T1& amt      (const IDVec& id) const { return keep(id) ? amt_[id[1]] : gen_amt[id[1]] ; }
    inline const T2& rate     (const IDVec& id) const { return keep(id) ? rate_[id[1]] : gen_rate[id[1]] ; }
    inline const T3& ii       (const IDVec& id) const { return keep(id) ? ii_[id[1]] : gen_ii[id[1]] ; }
    inline int cmt     (const IDVec& id) const { return keep(id) ? cmt_[id[1]] : gen_cmt[id[1]] ; }
    inline int addl    (const IDVec& id) const { return keep(id) ? addl_[id[1]] : gen_addl[id[1]] ; }
    inline int ss      (const IDVec& id) const { return keep(id) ? ss_[id[1]] : gen_ss[id[1]] ; }

    inline T_time time (int i) const { return time(index[i]); }
    inline const T1& amt      (int i) const { return amt (index[i]); }
    inline const T2& rate     (int i) const { return rate(index[i]); }
    inline const T3& ii       (int i) const { return ii  (index[i]); }
    inline int cmt     (int i) const { return cmt (index[i]); }
    inline int addl    (int i) const { return addl(index[i]); }
    inline int ss      (int i) const { return ss  (index[i]); }

    inline size_t num_state_times() const { return index.size(); }

    inline std::vector<T_rate> fractioned_rates(int i) const {
      const int n = rates[0].second.size();
      const std::vector<T2>& r = rates[rate_index[i]].second;
      std::vector<T_rate> res(r.size());
      for (size_t j = 0; j < r.size(); ++j) {
        res[j] = r[j] * bioavailability(i, j);
      }
      return res;
    }

    inline T_amt fractioned_amt(int i) const {
      return bioavailability(i, cmt(i) - 1) * amt(i);
    }

    std::vector<int> unique_times() {
      std::vector<int> t(index.size());
      std::iota(t.begin(), t.end(), 0);
      auto last = std::unique(t.begin(), t.end(),
                              [this](const int& i, const int& j) {return time(i) == time(j);});
      t.erase(last, t.end());
      return t;
    }

    /**
     * Implement absorption lag times by modifying the times of the dosing events.
     * Two cases: parameters are either constant or vary with each event.
     * Function sorts events at the end of the procedure.
     *
     * @tparam T_parameters type of scalar model parameters
     * @return - modified events that account for absorption lag times
     */
    void insert_lag_dose() {
      // reverse loop so we don't process same lagged events twice
      int nEvent = num_state_times();
      int iEvent = nEvent - 1;
      while (iEvent >= 0) {
        if (is_dosing(iEvent)) {
          if (GetValueTlag(iEvent, cmt(iEvent) - 1) != 0) {
            insert_event(iEvent);
            gen_time.back() += GetValueTlag(iEvent, cmt(iEvent) - 1);

            // Events[iEvent].evid = 2;  // Check
            index[iEvent][2] = 2;
            // The above statement changes events so that CleanEvents does
            // not return an object identical to the original. - CHECK
          }
        }
        iEvent--;
      }
      sort_state_time();
    }

    inline const T4_container& model_param(int i) const {
      return theta_[param_index[i].second[0]];
    }

    inline const T4& GetValue(int iEvent, int iParameter) const {
      return theta_[std::get<1>(param_index[iEvent])[0]][iParameter];
    }

    inline const T5& bioavailability(int iEvent, int iParameter) const {
      return biovar_[std::get<1>(param_index[iEvent])[1]][iParameter];
    }

    inline const T6& GetValueTlag(int iEvent, int iParameter) const {
      return tlag_[std::get<1>(param_index[iEvent])[2]][iParameter];
    }

    /*
     * Overloading the << Operator
     */
    friend std::ostream& operator<<(std::ostream& os, const EventHistory& ev) {
      const int w = 6;
      os << "\n";
      os << std::setw(w) << "time" <<
        std::setw(w) << "amt" <<
        std::setw(w) << "rate" <<
        std::setw(w) << "ii" <<
        std::setw(w) << "evid" <<
        std::setw(w) << "cmt" <<
        std::setw(w) << "addl" <<
        std::setw(w) << "ss" <<
        std::setw(w) << "keep" <<
        std::setw(w) << "isnew" << "\n";
      for (size_t i = 0; i < ev.size(); ++i) {
        os <<
          std::setw(w)   << ev.time(i) << " " <<
          std::setw(w-1) << ev.amt(i) << " " <<
          std::setw(w-1) << ev.rate(i) << " " <<
          std::setw(w-1) << ev.ii(i) << " " <<
          std::setw(w-1) << ev.evid(i) << " " <<
          std::setw(w-1) << ev.cmt(i) << " " <<
          std::setw(w-1) << ev.addl(i) << " " <<
          std::setw(w-1) << ev.ss(i) << " " <<
          std::setw(w-1) << ev.keep(i) << " " <<
          std::setw(w-1) << ev.isnew(i) << "\n";
      }
      return os;
    }
  };

}    // torsten namespace
#endif
