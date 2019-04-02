#ifndef STAN_MATH_TORSTEN_NONMEN_EVENTS_RECORD_HPP
#define STAN_MATH_TORSTEN_NONMEN_EVENTS_RECORD_HPP

#include <stan/math/torsten/dsolve/pk_vars.hpp>
#include <boost/math/tools/promotion.hpp>
#include <Eigen/Dense>
#include <numeric>
#include <string>
#include <vector>

namespace torsten {

  /*
   * Raw events record in the form of NONMEN, except that
   * for the population input @c len consisting of data record
   * length for each individual.
   * @tparam T0 type of scalar for time of events. 
   * @tparam T1 type of scalar for amount at each event.
   * @tparam T2 type of scalar for rate at each event.
   * @tparam T3 type of scalar for inter-dose inteveral at each event.
   * @tparam T4 type of scalars for the model parameters.
   * @tparam T5 type of scalars for the bio-variability parameters.
   * @tparam T6 type of scalars for the model tlag parameters.
   */
template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
struct NONMENEventsRecord {
  using T_scalar = typename stan::return_type<T0, T1, T2, T3, typename stan::return_type<T4, T5, T6>::type >::type;
  using T_time = typename stan::return_type<T0, T1, T6, T2>::type;
  using T_rate = typename stan::return_type<T2, T5>::type;
  using T_amt = typename stan::return_type<T1, T5>::type;
  using T_par = T4;
  using T_par_rate = T2;
  using T_par_ii = T3;

private:
  static const std::vector<Eigen::Matrix<T4, -1, -1>> dummy_param_matrix;
  static const std::vector<std::vector<T4> > dummy_param_vector;
  const std::vector<int> len_1_;  

public:
  const int nev;
  const int ncmt;
  std::vector<int> begin_;
  const std::vector<int>& len_;
  const std::vector<T0>& time_;
  const std::vector<T1>& amt_;
  const std::vector<T2>& rate_;
  const std::vector<T3>& ii_;
  const std::vector<int>& evid_;
  const std::vector<int>& cmt_;
  const std::vector<int>& addl_;
  const std::vector<int>& ss_;
  const std::vector<std::vector<T4> >& pMatrix_;
  const std::vector<Eigen::Matrix<T4, -1, -1>>& systems_;
  const std::vector<std::vector<T5> >& biovar_;
  const std::vector<std::vector<T6> >& tlag_;

  /*
   * Constructor using population data with parameter give as
   * matrix, such as in linear ODE models
   * @param[in] n number of compartments in model
   * @param[in] len record length for each individual.
   * @param[in] time times of events  
   * @param[in] amt amount at each event
   * @param[in] rate rate at each event
   * @param[in] ii inter-dose interval at each event
   * @param[in] evid event identity: 
   *                    (0) observation 
   *                    (1) dosing
   *                    (2) other 
   *                    (3) reset 
   *                    (4) reset AND dosing 
   * @param[in] cmt compartment number at each event 
   * @param[in] addl additional dosing at each event 
   * @param[in] ss steady state approximation at each event (0: no, 1: yes)
   * @param[in] pMatrix parameters at each event
   * @param[in] biovar bioavailability
   * @param[in] tlag lag time
   */
  template <typename T0_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_, typename T6_>
  NONMENEventsRecord(int n,
                     const std::vector<int>& len,
                     const std::vector<T0_>& time,
                     const std::vector<T1_>& amt,
                     const std::vector<T2_>& rate,
                     const std::vector<T3_>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<Eigen::Matrix<T4_, -1, -1>>& systems,
                     const std::vector<std::vector<T5_> >& biovar,
                     const std::vector<std::vector<T6_> >& tlag) :
    len_1_(),
    nev(time.size()),
    ncmt(n),
    begin_(len.size()),
    len_(len),
    time_   (time  ),
    amt_    (amt   ),
    rate_   (rate   ),
    ii_     (ii     ),
    evid_   (evid   ),
    cmt_    (cmt    ),
    addl_   (addl   ),
    ss_     (ss     ),
    pMatrix_(dummy_param_vector),
    systems_(systems),
    biovar_ (biovar ),
    tlag_   (tlag   )
  {
    begin_[0] = 0;
    std::partial_sum(len.begin(), len.end() - 1, begin_.begin() + 1);
  }

  /*
   * Constructor using individual data with parameter give as
   * matrix, such as in linear ODE models
   * @param[in] n number of compartments in model
   * @param[in] len record length for each individual.
   * @param[in] time times of events  
   * @param[in] amt amount at each event
   * @param[in] rate rate at each event
   * @param[in] ii inter-dose interval at each event
   * @param[in] evid event identity: 
   *                    (0) observation 
   *                    (1) dosing
   *                    (2) other 
   *                    (3) reset 
   *                    (4) reset AND dosing 
   * @param[in] cmt compartment number at each event 
   * @param[in] addl additional dosing at each event 
   * @param[in] ss steady state approximation at each event (0: no, 1: yes)
   * @param[in] pMatrix parameters at each event
   * @param[in] biovar bioavailability
   * @param[in] tlag lag time
   */
  template <typename T0_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_, typename T6_>
  NONMENEventsRecord(int n,
                     const std::vector<T0_>& time,
                     const std::vector<T1_>& amt,
                     const std::vector<T2_>& rate,
                     const std::vector<T3_>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<Eigen::Matrix<T4_, -1, -1>>& systems,
                     const std::vector<std::vector<T5_> >& biovar,
                     const std::vector<std::vector<T6_> >& tlag) :
    len_1_(1, time.size()),
    nev(time.size()),
    ncmt(n),
    begin_{0},
    len_(len_1_),
    time_   (time   ),
    amt_    (amt    ),
    rate_   (rate   ),
    ii_     (ii     ),
    evid_   (evid   ),
    cmt_    (cmt    ),
    addl_   (addl   ),
    ss_     (ss     ),
    pMatrix_(dummy_param_vector),
    systems_(systems),
    biovar_ (biovar ),
    tlag_   (tlag   )
  {}

  /*
   * Constructor using population data with parameter give as vectors.
   * @param[in] n number of compartments in model
   * @param[in] len record length for each individual.
   * @param[in] time times of events  
   * @param[in] amt amount at each event
   * @param[in] rate rate at each event
   * @param[in] ii inter-dose interval at each event
   * @param[in] evid event identity: 
   *                    (0) observation 
   *                    (1) dosing
   *                    (2) other 
   *                    (3) reset 
   *                    (4) reset AND dosing 
   * @param[in] cmt compartment number at each event 
   * @param[in] addl additional dosing at each event 
   * @param[in] ss steady state approximation at each event (0: no, 1: yes)
   * @param[in] pMatrix parameters at each event
   * @param[in] biovar bioavailability
   * @param[in] tlag lag time
   */
  template <typename T0_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_, typename T6_>
  NONMENEventsRecord(int n,
                     const std::vector<int>& len,
                     const std::vector<T0_>& time,
                     const std::vector<T1_>& amt,
                     const std::vector<T2_>& rate,
                     const std::vector<T3_>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<std::vector<T4_> >& pMatrix,
                     const std::vector<std::vector<T5_> >& biovar,
                     const std::vector<std::vector<T6_> >& tlag) :
    len_1_(),
    nev(time.size()),
    ncmt(n),
    begin_(len.size()),
    len_(len),
    time_   (time   ),
    amt_    (amt    ),
    rate_   (rate   ),
    ii_     (ii     ),
    evid_   (evid   ),
    cmt_    (cmt    ),
    addl_   (addl   ),
    ss_     (ss     ),
    pMatrix_(pMatrix),
    systems_(dummy_param_matrix),
    biovar_ (biovar ),
    tlag_   (tlag   )
  {
    begin_[0] = 0;
    std::partial_sum(len.begin(), len.end() - 1, begin_.begin() + 1);    
  }

  /*
   * Constructor using individual data with parameter give as vectors.
   * @param[in] n number of compartments in model
   * @param[in] len record length for each individual.
   * @param[in] time times of events  
   * @param[in] amt amount at each event
   * @param[in] rate rate at each event
   * @param[in] ii inter-dose interval at each event
   * @param[in] evid event identity: 
   *                    (0) observation 
   *                    (1) dosing
   *                    (2) other 
   *                    (3) reset 
   *                    (4) reset AND dosing 
   * @param[in] cmt compartment number at each event 
   * @param[in] addl additional dosing at each event 
   * @param[in] ss steady state approximation at each event (0: no, 1: yes)
   * @param[in] pMatrix parameters at each event
   * @param[in] biovar bioavailability
   * @param[in] tlag lag time
   */
  template <typename T0_, typename T1_, typename T2_, typename T3_, typename T4_, typename T5_, typename T6_>
  NONMENEventsRecord(int n,
                     const std::vector<T0_>& time,
                     const std::vector<T1_>& amt,
                     const std::vector<T2_>& rate,
                     const std::vector<T3_>& ii,
                     const std::vector<int>& evid,
                     const std::vector<int>& cmt,
                     const std::vector<int>& addl,
                     const std::vector<int>& ss,
                     const std::vector<std::vector<T4_> >& pMatrix,
                     const std::vector<std::vector<T5_> >& biovar,
                     const std::vector<std::vector<T6_> >& tlag) :
    len_1_(1, time.size()),
    nev(time.size()),
    ncmt(n),
    begin_{0},
    len_(len_1_),
    time_   (time   ),
    amt_    (amt    ),
    rate_   (rate   ),
    ii_     (ii     ),
    evid_   (evid   ),
    cmt_    (cmt    ),
    addl_   (addl   ),
    ss_     (ss     ),
    pMatrix_(pMatrix),
    systems_(dummy_param_matrix),
    biovar_ (biovar ),
    tlag_   (tlag   )
  {}

  /*
   * begin of the parameters for individual @c id
   * in @c pMatrix. It is assumed that all the paramter are
   * either constant or time dependent.
   */
  int begin_param(int id) const {
    if(pMatrix_.empty()) {
      return systems_.size() == len_.size() ? id : begin_[id];
    } else {
      return pMatrix_.size() == len_.size() ? id : begin_[id];
    }    
  }

  /*
   * length of the parameters for individual @c id
   * in @c pMatrix. It is assumed that all the paramter are
   * either constant or time dependent.
   */
  int len_param(int id) const {
    if(pMatrix_.empty()) {
      return systems_.size() == len_.size() ? 1 : len_[id]; 
    } else {
      return pMatrix_.size() == len_.size() ? 1 : len_[id]; 
    }
  }

  int begin_biovar(int id) const {
    return biovar_.size()  == len_.size() ? id : begin_[id];
  }

  int len_biovar(int id) const {
    return biovar_.size()  == len_.size() ? 1 : len_[id];
  }

  int begin_tlag(int id) const {
    return tlag_.size()  == len_.size() ? id : begin_[id];
  }

  int len_tlag(int id) const {
    return tlag_.size()  == len_.size() ? 1 : len_[id];
  }

  static bool has_ss_dosing() {
    return has_ss_dosing(0);
  }
    /*
     * check the exisitence of steady state dosing events
     */
  bool has_ss_dosing(int id) const {
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

};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
const std::vector<Eigen::Matrix<T4, -1, -1>>
NONMENEventsRecord<T0, T1, T2, T3, T4, T5, T6>::dummy_param_matrix = {};

template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
const std::vector<std::vector<T4> >
NONMENEventsRecord<T0, T1, T2, T3, T4, T5, T6>::dummy_param_vector = {};

}

#endif
