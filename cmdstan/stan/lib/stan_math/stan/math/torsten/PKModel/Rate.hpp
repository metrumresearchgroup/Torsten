#ifndef STAN_MATH_TORSTEN_PKMODEL_RATE_HPP
#define STAN_MATH_TORSTEN_PKMODEL_RATE_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/functions.hpp>
#include <algorithm>
#include <vector>

namespace torsten {

// forward declaration
template<typename T_time, typename T_rate> class RateHistory;

/**
 * The Rate class defines objects that contain the rate in each compartment
 * at each time of the event schedule (but not nescessarily at each event,
 * since two events may happen at the same time).
 */
template<typename T_time, typename T_rate>
class Rate {
private:
  T_time time;
public:
  std::vector<T_rate> rate;  // rate for each compartment

  Rate() {
    std::vector<T_rate> v(1, 0);
    time = 0;
    rate = v;
  }

  Rate(T_time p_time, std::vector<T_rate> p_rate) {
    time = p_time;
    rate = p_rate;
  }

  // access functions
  T_time get_time() const { return time; }
  std::vector<T_rate> get_rate() const { return rate; }

  // Overload = operator
  // Allows us to construct a rate of var from a rate of double
  template <typename T0, typename T1>
  void copy(const Rate<T0, T1>& rate1) {
    time = rate1.get_time();
    rate.resize(rate1.get_rate().size());
    for (size_t i = 0; i < rate.size(); i++)
      rate[i] = rate1.get_rate()[i];
  }

  void Print() {
    std::cout << time << " ";
    for (int i = 0; i < rate.size(); i++) std::cout << rate[i] << " ";
    std::cout << std::endl;
  }

  friend class RateHistory<T_time, T_rate>;
  template <typename T_amt, typename T_ii>
  friend void MakeRates(torsten::EventHistory<T_time, T_amt, T_rate, T_ii>&,
    RateHistory<T_time, T_rate>&);

  template <typename T_0, typename T_1, typename T_2, typename T_3,
    typename T_4, typename T_5, typename T_6, typename F_1, typename F_2>
  friend
  Eigen::Matrix<typename boost::math::tools::promote_args<T_0, T_1, T_2, T_3,
    typename boost::math::tools::promote_args<T_4, T_5, T_6>::type
    >::type, Eigen::Dynamic, Eigen::Dynamic>
  Pred(const std::vector<T_0>& time,
       const std::vector<T_1>& amt,
       const std::vector<T_2>& rate,
       const std::vector<T_3>& ii,
       const std::vector<int>& evid,
       const std::vector<int>& cmt,
       const std::vector<int>& addl,
       const std::vector<int>& ss,
       const std::vector<std::vector<T_4> >& pMatrix,
       const std::vector<std::vector<T_5> >& biovar,
       const std::vector<std::vector<T_6> >& tlag,
       const int& nCmt,
       const std::vector<Eigen::Matrix<T_4, Eigen::Dynamic, Eigen::Dynamic> >&
         system,
       const F_1& f1,
       const F_2& fss);
};

/**
 * The RateHistory class defines objects that contain a vector of rates,
 * along with a series of functions that operate on them.
 */
template <typename T_time, typename T_rate>
class RateHistory {
private:
  std::vector<Rate<T_time, T_rate> > Rates;

public:
  RateHistory() {
    Rate<T_time, T_rate> initRate;
    Rates.resize(1);
    Rates[0] = initRate;
  }

  template <typename T0, typename T1>
  RateHistory(std::vector<T0> p_time, std::vector<std::vector<T1> > p_rate) {
    int nRate = p_rate.size();
    Rates.resize(nRate);
    for (int i = 0; i < nRate; i++) Rates[i] = Rate<T_time, T_rate>(p_time[i],
      p_rate[i]);
  }

  explicit RateHistory(int nEvent) {
    Rate<T_time, T_rate> initRate;
    Rates.resize(nEvent);
    for (int i = 0; i < nEvent; i++) Rates[i] = initRate;
  }

  T_time get_time(int i) { return Rates[i].time; }
  std::vector<T_rate> get_rate(int i) { return Rates[i].rate; }

  bool Check() {
    int i = Rates.size() - 1;
    bool ordered = true;

    while ((i > 0) && (ordered)) {
      ordered = (Rates[i].time >= Rates[i - 1].time);
      i--;
    }
    return ordered;
  }

  Rate<T_time, T_rate> GetRate(int i) {
    Rate<T_time, T_rate> newRate(Rates[i].time, Rates[i].rate);
    return newRate;
  }

  void InsertRate(Rate<T_time, T_rate> p_Rate) { Rates.push_back(p_Rate); }

  void RemoveRate(int i) {
    assert(i >= 0);
    Rates.erase(Rates.begin() + i);
  }

  int Size() { return Rates.size(); }

  void Print(int j) {
    std::cout << Rates[j].time << " ";
    for (int i = 0; i < Rates[j].rate.size(); i++)
      std::cout << Rates[j].rate[i] << " ";
    std::cout << std::endl;
  }

  struct by_time {
    bool operator()(Rate<T_time, T_rate> const &a, Rate<T_time, T_rate>
      const &b) {
      return a.time < b.time;
    }
  };

  void Sort() { std::sort(Rates.begin(), Rates.end(), by_time()); }

  template <typename T_amt, typename T_ii>
  void MakeRates(torsten::EventHistory<T_time, T_amt, T_rate, T_ii>& events, int nCmt) {
    using std::vector;

    if (!events.Check()) events.Sort();

    vector<T_rate> rate_init(nCmt, 0);
    Rate<T_time, T_rate> newRate(0, rate_init);
    for (int j = 0; j < events.get_size(); j++)
      if (j == 0 || events.get_time(j) != events.get_time(j - 1)) {
        newRate.time = events.get_time(j);
        InsertRate(newRate);
      }

    RemoveRate(0);  // remove rate created by default constructor.

    if (!Check()) Sort();

    // Create time vector for rates
    vector<T_time> RateTimes(Size(), 0);
    for (int j = 0; j < Size(); j++) RateTimes[j] = Rates[j].time;

    // Create time vector for events
    vector<T_time> EventTimes(events.get_size(), 0);
    for (int j = 0; j < events.get_size(); j++)
      EventTimes[j] = events.get_time(j);

    int i = 0, k, l;
    T_time endTime;
    torsten::Event<T_time, T_amt, T_rate, T_ii> newEvent;
    while (i < events.get_size()) {
      if ((events.get_evid(i) == 1 || events.get_evid(i) == 4)
        && (events.get_rate(i) > 0 && events.get_amt(i) > 0)) {
          endTime = events.get_time(i) + events.get_amt(i)/events.get_rate(i);
          newEvent = newEvent(endTime, 0, 0, 0, 2, events.get_cmt(i), 0, 0,
            false, true);
          events.InsertEvent(newEvent);
          if (!events.Check()) events.Sort();
          EventTimes.push_back(endTime);
          std::sort(EventTimes.begin(), EventTimes.end());

          // Only create a new Rate if endTime does not correspond to a time
          // that is already in RateHistory. - CHECK
          if (!find_time(RateTimes, endTime)) {
            newRate.time = endTime;
            InsertRate(newRate);
            if (!Check()) Sort();
            RateTimes.push_back(endTime);
            std::sort(RateTimes.begin(), RateTimes.end());
          }

          // Find indexes at which time of event and endtime occur.
          l = SearchReal(RateTimes, events.get_size(), events.get_time(i));
          k = SearchReal(RateTimes, events.get_size(), endTime);

          // Compute Rates for each element between the two times
          for (int iRate = l ; iRate < k; iRate++)
            Rates[iRate].rate[events.get_cmt(i) - 1] += events.get_rate(i);
        }
        i++;
    }

    // Sort events and rates
    if (!Check()) Sort();
    if (!events.Check()) events.Sort();
  }

  // declare friends
  friend class Rate<T_time, T_rate>;
  template <typename T_amt, typename T_ii>
  friend void MakeRates(torsten::EventHistory<T_time, T_amt, T_rate, T_ii>&,
    RateHistory<T_time, T_amt>&);
};

}

#endif
