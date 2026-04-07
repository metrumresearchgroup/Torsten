/**
 * @file   ev_record.hpp
 * @author Yi Zhang <yiz@Yis-MacBook-Pro.local>
 * @date   Wed Aug 26 23:43:13 2020
 * 
 * @brief  raw event record
 * 
 * 
 */
#ifndef STAN_MATH_TORSTEN_NONMEN_EVENTS_RECORD_HPP
#define STAN_MATH_TORSTEN_NONMEN_EVENTS_RECORD_HPP

#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <Eigen/Dense>
#include <numeric>
#include <string>
#include <vector>

namespace torsten {

  /**
   * Raw events record in the form of NONMEN, except that
   * for the population input @c len consisting of data record
   * length for each individual.
   * @tparam T0 type of scalar for time of events. 
   * @tparam T1 type of scalar for amount at each event.
   * @tparam T2 type of scalar for rate at each event.
   * @tparam T3 type of scalar for inter-dose inteveral at each event.
   * @tparam theta_container type of container for parameter. <code>linode</code> model uses
   *         matrix as parameter, other than <code>std::vector</code> used by others.
   */
  template <typename T0, typename T1, typename T2, typename T3>
  struct NONMENEventsRecord {
    using T_scalar = typename stan::return_type_t<T0, T1, T2, T3>;
    using T_time = typename stan::return_type_t<T0, T1, T2, T3>;
    using T_rate = T2;
    using T_amt = T1;
    using T_par_rate = T2;
    using T_par_ii = T3;

  private:
    /// <code>len</code> needs actual variable to ref to.
    const std::vector<int> len_1_;

  public:
    static const double lag_time_min; /**< minimum lag time allowed */

    /// nb. of compartments
    const int ncmt;
    /// begin id for each subject in concat vector <code>time_</code>, etc
    std::vector<int> begin_;
    /// data length for each subject
    const std::vector<int>& len_;
    /// nb. of events in the entire population
    const int total_num_event_times;
    /// event time
    const std::vector<T0>& time_;
    /// dosing amount
    const std::vector<T1>& amt_;
    /// dosing rate
    const std::vector<T2>& rate_;
    /// dosing interval
    const std::vector<T3>& ii_;
    /// dosing event id\n
    ///      (0) observation\n
    ///      (1) dosing\n
    ///      (2) other\n
    ///      (3) reset\n
    ///      (4) reset AND dosing
    const std::vector<int>& evid_;
    /// event compartment
    const std::vector<int>& cmt_;
    /// nb. of additional doses
    const std::vector<int>& addl_;
    /// steady-states flag
    const std::vector<int>& ss_;

    /**
     * Constructor using population data with parameter give as
     * matrix, such as in linear ODE models
     * @param[in] n number of compartments in model
     * @param[in] len record length for each individual.
     * @param[in] time times of events  
     * @param[in] amt amount at each event
     * @param[in] rate rate at each event
     * @param[in] ii inter-dose interval at each event
     * @param[in] evid event identity: \n
     *                    (0) observation\n 
     *                    (1) dosing\n
     *                    (2) other \n
     *                    (3) reset \n
     *                    (4) reset AND dosing 
     * @param[in] cmt compartment number at each event 
     * @param[in] addl additional dosing at each event 
     * @param[in] ss steady state approximation at each event (0: no, 1: yes)
     */
    template <typename T0_, typename T1_, typename T2_, typename T3_>
    NONMENEventsRecord(int n,
                       const std::vector<int>& len,
                       const std::vector<T0_>& time,
                       const std::vector<T1_>& amt,
                       const std::vector<T2_>& rate,
                       const std::vector<T3_>& ii,
                       const std::vector<int>& evid,
                       const std::vector<int>& cmt,
                       const std::vector<int>& addl,
                       const std::vector<int>& ss) :
      len_1_(),
      ncmt(n),
      begin_(len.size()),
      len_(len),
      total_num_event_times(std::accumulate(len_.begin(), len_.end(), 0)),
      time_   (time  ),
      amt_    (amt   ),
      rate_   (rate   ),
      ii_     (ii     ),
      evid_   (evid   ),
      cmt_    (cmt    ),
      addl_   (addl   ),
      ss_     (ss     )
    {
      begin_[0] = 0;
      std::partial_sum(len.begin(), len.end() - 1, begin_.begin() + 1);
    }

    /**
     * Constructor using individual data with parameter give as
     * matrix, such as in linear ODE models
     * @param[in] n number of compartments in model
     * @param[in] time times of events  
     * @param[in] amt amount at each event
     * @param[in] rate rate at each event
     * @param[in] ii inter-dose interval at each event
     * @param[in] evid event identity: \n
     *                    (0) observation \n
     *                    (1) dosing\n
     *                    (2) other \n
     *                    (3) reset \n
     *                    (4) reset AND dosing 
     * @param[in] cmt compartment number at each event 
     * @param[in] addl additional dosing at each event 
     * @param[in] ss steady state approximation at each event (0: no, 1: yes)
     */
    template <typename T0_, typename T1_, typename T2_, typename T3_>
    NONMENEventsRecord(int n,
                       const std::vector<T0_>& time,
                       const std::vector<T1_>& amt,
                       const std::vector<T2_>& rate,
                       const std::vector<T3_>& ii,
                       const std::vector<int>& evid,
                       const std::vector<int>& cmt,
                       const std::vector<int>& addl,
                       const std::vector<int>& ss) :
      len_1_(1, time.size()),
      ncmt(n),
      begin_{0},
      len_(len_1_),
      total_num_event_times(std::accumulate(len_.begin(), len_.end(), 0)),
      time_   (time   ),
      amt_    (amt    ),
      rate_   (rate   ),
      ii_     (ii     ),
      evid_   (evid   ),
      cmt_    (cmt    ),
      addl_   (addl   ),
      ss_     (ss     )
    {}
  
    /** 
     * begin of parameters for a subject @c id
     * in @c pMatrix. It is assumed that all the paramter are
     * either constant or time dependent.
     *
     * @param id subject id
     * 
     * @return begin index in @c pMatrix for the subject
     */
    template<typename T>
    int begin_param(int id, const std::vector<T>& param) const {
      return param.size() == len_.size() ? id : begin_[id];
    }

    /**
     * length of parameters for a subject @c id
     * in @c pMatrix. It is assumed that all the paramter are
     * either constant or time dependent.
     *
     * @param id subject id
     * 
     * @return len in @c pMatrix for the subject
     */
    template<typename T>
    int len_param(int id, const std::vector<T>& param) const {
      return param.size() == len_.size() ? 1 : len_[id]; 
    }

    /** 
     * check the exisitence of steady state dosing events
     * for the 1st subject, used for single-individual data.
     * 
     * @return if the subject has any SS dosing event.
     */
    inline bool has_ss_dosing() const {
      return has_ss_dosing(0);
    }

    /** 
     * check the exisitence of steady state dosing events
     * 
     * @param id subject id
     * 
     * @return if the subject has any SS dosing event.
     */
    inline bool has_ss_dosing(int id) const {
      bool res = false;
      int begin = begin_[id];
      int end = size_t(id + 1) == len_.size() ? time_.size() : begin_[id + 1];
      for (int i = begin; i < end; ++i) {
        if ((evid_[i] == 1 || evid_[i] == 4) && ss_[i] != 0) {
          res = true;
          break;
        }
      }
      return res;
    }

    /** 
     * check the exisitence of time-lagged events
     * 
     * @param id subject id
     * 
     * @return if the subject has any time-lagged dosing event.
     */
    template<typename T>
    inline bool has_positive_param(int id, const std::vector<std::vector<T>>& params) const {
      using stan::math::value_of;
      return std::any_of(params.begin() + begin_param(id, params),
                         params.begin() + begin_param(id, params) + len_param(id, params),
                         [](const std::vector<T>& v) {
                           return std::any_of(v.begin(), v.end(), [](const T& x) { return std::abs(value_of(x)) > lag_time_min; });
                         });
    }

    /** 
     * check the exisitence of time-lagged dosing events
     * for the 1st subject, used for single-individual data.
     * 
     * @return if the subject has any time-lagged dosing event.
     */
    template<typename T>
    inline bool has_positive_param(const std::vector<std::vector<T>>& params) const {
      return has_positive_param(0, params);
    }

    /** 
     * @param id subject id
     * 
     * @return nb. of events for the subject
     */
    inline int num_event_times(int id) const {
      return len_.at(id);
    }

    /** 
     * @return nb. of events for subject 0.
     */
    inline int num_event_times() const {
      return len_.at(0);
    }

    /** 
     * @return population size.
     */
    inline int num_subjects() const {
      return len_.size();
    }
  };

  template <typename T0, typename T1, typename T2, typename T3>
  const double NONMENEventsRecord<T0, T1, T2, T3>::lag_time_min = 1.e-12;

}

#endif
