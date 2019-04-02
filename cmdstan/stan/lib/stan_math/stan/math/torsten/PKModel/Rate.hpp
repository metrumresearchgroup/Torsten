#ifndef STAN_MATH_TORSTEN_PKMODEL_RATE_HPP
#define STAN_MATH_TORSTEN_PKMODEL_RATE_HPP

#include <stan/math/torsten/PKModel/Event.hpp>
#include <stan/math/torsten/PKModel/functions.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace torsten {

/**
 * The RateHistory class defines objects that contain a vector of rates,
 * along with a series of functions that operate on them.
 */
template <typename T_time, typename T_rate>
struct RateHistory {

  template<typename T1, typename T2>
  struct Rate {
    T1 time;
    std::vector<T2> rate;
    Rate() {}
    Rate(T1 p_time, std::vector<T2>& p_rate) : time(p_time), rate(p_rate) {}
  };

  std::vector<Rate<T_time, T_rate> > Rates;

  /*
   * generate rates using event history
   */
  template <typename T_amt, typename T_ii>
  RateHistory(torsten::EventHistory<T_time, T_amt, T_rate, T_ii>& events, int nCmt) {
    using std::vector;

    if (!events.Check()) events.Sort();

    vector<T_rate> rate_init(nCmt, 0);
    Rate<T_time, T_rate> newRate(0, rate_init);
    for (size_t j = 0; j < events.size(); j++)
      if (j == 0 || events.time(j) != events.time(j - 1)) {
        newRate.time = events.time(j);
        Rates.push_back(newRate);
      }

    if (!std::is_sorted(Rates.begin(), Rates.end(), by_time())) sort();

    // Create time vector for rates
    vector<T_time> RateTimes(Rates.size(), 0);
    for (size_t j = 0; j < Rates.size(); j++) RateTimes[j] = Rates[j].time;

    // Create time vector for events
    vector<T_time> EventTimes(events.size(), 0);
    for (size_t j = 0; j < events.size(); j++) EventTimes[j] = events.time(j);

    size_t i = 0, k, l;
    T_time endTime;
    torsten::Event<T_time, T_amt, T_rate, T_ii> newEvent;
    while (i < events.size()) {
      if ((events.is_dosing(i)) && (events.rate(i) > 0 && events.amt(i) > 0)) {
          endTime = events.time(i) + events.amt(i)/events.rate(i);
          newEvent = newEvent(endTime, 0, 0, 0, 2, events.cmt(i), 0, 0, false, true);
          events.InsertEvent(newEvent);
          if (!events.Check()) events.Sort();
          EventTimes.push_back(endTime);
          std::sort(EventTimes.begin(), EventTimes.end());

          // Only create a new Rate if endTime does not correspond to a time
          // that is already in RateHistory. - CHECK
          if (!find_time(RateTimes, endTime)) {
            newRate.time = endTime;
            Rates.push_back(newRate);
            // InsertRate(newRate);
            if (!std::is_sorted(Rates.begin(), Rates.end(), by_time())) sort();
            RateTimes.push_back(endTime);
            std::sort(RateTimes.begin(), RateTimes.end());
          }

          // Find indexes at which time of event and endtime occur.
          l = SearchReal(RateTimes, events.size(), events.time(i));
          k = SearchReal(RateTimes, events.size(), endTime);

          // Compute Rates for each element between the two times
          for (size_t iRate = l ; iRate < k; iRate++)
            Rates[iRate].rate[events.cmt(i) - 1] += events.rate(i);
        }
        i++;
    }

    // Sort events and rates
    if (!std::is_sorted(Rates.begin(), Rates.end(), by_time())) sort();
    if (!events.Check()) events.Sort();
  }

  T_time time(int i) { return Rates[i].time; }

  T_rate rate(int i, int j) { return Rates[i].rate[j]; }

  struct by_time {
    bool operator()(Rate<T_time, T_rate> const &a, Rate<T_time, T_rate>
      const &b) {
      return a.time < b.time;
    }
  };

  void sort() { std::sort(Rates.begin(), Rates.end(), by_time()); }
};

}

#endif
