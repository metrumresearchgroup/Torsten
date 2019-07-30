#ifndef STAN_MATH_TORSTEN_PKMODEL_MODELPARAMETERS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_MODELPARAMETERS_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/torsten/event_history.hpp>
#include <stan/math/torsten/PKModel/ExtractVector.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace torsten {

/**
 * The ModelParameters class defines objects that contain the parameters of
 * a model at a given time.
 */
template<typename T_time, typename T_parameters, typename T_biovar, typename T_tlag>
struct ModelParameters {
  double time_;
  std::vector<T_parameters> theta_;
  int nrow, ncol;
  std::vector<T_biovar> biovar_;
  std::vector<T_tlag> tlag_;

  ModelParameters() {}

  ModelParameters(const T_time& time,
                  const Eigen::Matrix<T_parameters, Eigen::Dynamic, Eigen::Dynamic>& K,
                  const std::vector<T_biovar>& biovar,
                  const std::vector<T_tlag>& tlag)
    : time_(stan::math::value_of(time)), theta_(K.size()), nrow(K.rows()), ncol(K.cols()), biovar_(biovar), tlag_(tlag)
  {
    Eigen::Matrix<T_parameters, -1, -1>::Map(theta_.data(), nrow, ncol) = K;
  }

  ModelParameters(const T_time& time,
                  const std::vector<T_parameters>& theta,
                  const std::vector<T_biovar>& biovar,
                  const std::vector<T_tlag>& tlag)
    : time_(stan::math::value_of(time)), theta_(theta), biovar_(biovar), tlag_(tlag) {}

  /**
   * Adds parameters. Useful for the mixed solver, where
   * we want to augment the parameters with the intial PK
   * states when calling the numerical integrator.
   */
  template <typename T>
  ModelParameters<T_time, T, T_biovar, T_tlag>
  augment(const std::vector<T>& thetaAdd) const {
    std::vector<T> theta(theta_.size());
    for (size_t i = 0; i < theta.size(); i++) theta[i] = theta_[i];
    for (size_t i = 0; i < thetaAdd.size(); i++) theta.push_back(thetaAdd[i]);
    return
      ModelParameters<T_time, T, T_biovar, T_tlag>
        (time_, theta, biovar_, tlag_);
  }

  template <typename T>
  ModelParameters<T_time, T, T_biovar, T_tlag>
  augment(const Eigen::Matrix<T, Eigen::Dynamic, 1>& thetaAdd)
  const {
    return augment(stan::math::to_array_1d(thetaAdd));
  }

  /**
   * Edit time stored in parameter object.
   */
  void time(double time) {
    time_ = time;
  }

  int CountParameters() const {
    return theta_.size();
  }

  // access functions   // FIX ME - name should be get_theta.
  double get_time() const { return time_; }
  std::vector<T_parameters> get_RealParameters(bool return_matrix) const {
    if (return_matrix) {
      auto k = get_K();
      std::vector<T_parameters> par(k.size());
      for (size_t j = 0; j < par.size(); ++j) par[j] = k(j);
      return par;
    } else {
      return theta_;
    }
  }
  std::vector<T_biovar> get_biovar() const {
    return biovar_;
  }
  std::vector<T_tlag> get_tlag() const {
    return tlag_;
  }
  Eigen::Matrix<T_parameters, Eigen::Dynamic, Eigen::Dynamic> get_K() const {
    Eigen::Matrix<T_parameters, -1, -1> res(nrow, ncol);
    res = Eigen::Matrix<T_parameters, -1, -1>::Map(theta_.data(), nrow, ncol);
    return res;
  }
};

  template<typename T>
  struct is_matrix : std::false_type {};
  
  template<typename T, int R, int C>
  struct is_matrix<Eigen::Matrix<T, R, C> > : std::true_type {};
  
/**
 * The ModelParameterHistory class defines objects that contain a vector
 * of ModelParameters, along with a series of functions that operate on
 * them.
 */
template<typename T_time, typename T4_container, typename T5, typename T6>
struct ModelParameterHistory {
  using T4 = typename stan::math::value_type<T4_container>::type;
  using Param = std::pair<double, std::array<int, 3> >;

  static const bool has_matrix_param;

  std::vector<Param> index;
  const std::vector<T4_container>& theta_;
  const std::vector<std::vector<T5> >& biovar_;
  const std::vector<std::vector<T6> >& tlag_;

  template<typename T0>
   ModelParameterHistory(const std::vector<T0>& time,
                        const std::vector<T4_container>& theta,
                        const std::vector<std::vector<T5> >& biovar,
                        const std::vector<std::vector<T6> >& tlag) :
     ModelParameterHistory(0, time.size(), time,
                           0, theta.size(), theta,
                           0, biovar.size(), biovar,
                           0, tlag.size(), tlag)
  {}

  /*
   * For population data in form of ragged array, we need to
   * generate individual parameter history given the entire
   * population data and the location of the
   * inidividual. However, @c theta, @c biovar and @c tlag
   * could have different lengths, so for each variable we
   * need an index that points to the range that belongs to the individual.
   * Note that if all three variables are of size 1, their
   * time is set to be the first entry of the @c time vector
   */
  template<typename T0>
  ModelParameterHistory(int ibegin, int isize,
                        const std::vector<T0>& time,
                        int ibegin_theta, int isize_theta,
                        const std::vector<T4_container>& theta,
                        int ibegin_biovar, int isize_biovar,
                        const std::vector<std::vector<T5> >& biovar,
                        int ibegin_tlag, int isize_tlag,
                        const std::vector<std::vector<T6> >& tlag) :
    index(isize),
    theta_(theta),
    biovar_(biovar),
    tlag_(tlag)
  {
    for (int i = 0; i < isize; ++i) {
      int j = isize_theta   > 1 ? ibegin_theta  + i : ibegin_theta;
      int k = isize_biovar  > 1 ? ibegin_biovar + i : ibegin_biovar;
      int l = isize_tlag    > 1 ? ibegin_tlag   + i : ibegin_tlag;
      index[i] = std::make_pair<double, std::array<int, 3> >(stan::math::value_of(time[ibegin + i]), {j, k, l });
    }
    Sort();
  }

  const T4_container& model_param(int i) const {
    return theta_[index[i].second[0]];
  }

  ModelParameters<T_time, T4, T5, T6> GetModelParameters(int i) const {
    return ModelParameters<T_time, T4, T5, T6>(std::get<0>(index[i]), theta_[std::get<1>(index[i])[0]], biovar_[std::get<1>(index[i])[1]], tlag_[std::get<1>(index[i])[2]]);
  }

  /**
   * MPV.size gives us the number of events.
   * MPV[i].RealParameters.size gives us the number of
   * ODE parameters for the ith event.
   * 
   * FIX ME - rename this GetValueTheta
   */
  inline const T4& GetValue(int iEvent, int iParameter) const {
    return theta_[std::get<1>(index[iEvent])[0]][iParameter];
  }

  inline const T5& GetValueBio(int iEvent, int iParameter) const {
    return biovar_[std::get<1>(index[iEvent])[1]][iParameter];
  }

  inline const T6& GetValueTlag(int iEvent, int iParameter) const {
    return tlag_[std::get<1>(index[iEvent])[2]][iParameter];
  }

  inline int get_size() const {
    return index.size();
  }

  void Sort() {
    std::sort(index.begin(), index.end(),
              [](const Param& a, const Param& b)
              { return std::get<0>(a) < std::get<0>(b); });
  }

  bool Check() {
  // check that elements are in chronological order.
    int i = index.size() - 1;
    bool ordered = true;

    while (i > 0 && ordered) {
      ordered = (std::get<0>(index[i]) >= std::get<0>(index[i-1]));
      i--;
    }
    return ordered;
  }

  /**
   * COMPLETE MODEL PARAMETERS
   *
   * Completes parameters so that it contains model parameters for each event 
   * in events. If parameters contains only one set of parameters (case where
   * the parameters are constant), this set is replicated for each event in
   * events. Otherwise a new parameter vector is added for each new event 
   * (isnew = true). This parameter vector is identical to the parameter vector
   * at the subsequent event. If the new event occurs at a time posterior to
   * the time of the last event, than the new vector parameter equals the
   * parameter vector of the last event. This amounts to doing an LOCF
   * (Last Observation Carried Forward):
   * Three cases:
   * (a) The time of the new event is higher than the time of the last
   *     parameter vector in parameters (k = len_parameters).
   *     Create a parameter vector at the the time of the new event,
   *     with the parameters of the last parameter vector.
   *     (Last Observation Carried Forward)
   * (b) The time of the new event matches the time of a parameter vector
   *     in parameters. This parameter vector gets replicated.
   * (c) (a) is not verified and no parameter vector occurs at the time
   *     of the new event. A new parameter vector is created at the time
   *     of the new event, and its parameters are equal to the parameters
   *     of the subsequent parameter vector in parameters.
   *
   * Since both @c index and events @c index are sorted
   * in time, we always move paramtter pointer to events
   *
   * Events and Parameters are sorted at the end of the procedure.
   *
   * @param[in] parameters at each event
   * @param[in] events elements (following NONMEM convention) at each event
   * @return - modified parameters and events.
   */
  template<typename T0, typename T_p1, typename T_p2, typename T_p3>
  void CompleteParameterHistory(torsten::EventHistory<T0, T_p1, T_p2, T_p3, T4_container, T5, T6>& events) {
    int nEvent = events.size();
    assert(nEvent > 0);
    int len_Parameters = index.size();  // numbers of events for which parameters are determined
    assert(len_Parameters > 0);

    if (!Check()) Sort();
    if (!events.Check()) events.Sort();
    index.resize(nEvent);

    int iEvent = 0;
    for (int i = 0; i < len_Parameters - 1; i++) {
      while (events.isnew(iEvent)) iEvent++;  // skip new events
      assert(std::get<0>(index[i]) == events.time(iEvent));  // compare time of "old' events to time of parameters.
      iEvent++;
    }

    if (len_Parameters == 1)  {
      for (int i = 0; i < nEvent; i++) {
        index[i] = std::make_pair<double, std::array<int, 3> >(stan::math::value_of(events.time(i)) , std::array<int,3>(std::get<1>(index[0])));
        events.index[i][3] = 0;
      }
    } else {  // parameters are event dependent.
      std::vector<double> times(nEvent, 0);
      for (int i = 0; i < nEvent; i++) times[i] = index[i].first;
      iEvent = 0;

      Param newParameter;
      int j = 0;
      std::vector<Param>::const_iterator lower = index.begin();
      std::vector<Param>::const_iterator it_param_end = index.begin() + len_Parameters;
      for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (events.isnew(iEvent)) {
          // Find the index corresponding to the time of the new event in the
          // times vector.
          const double t = stan::math::value_of(events.time(iEvent));
          lower = std::lower_bound(lower, it_param_end, t,
                                   [](const Param& t1, const double& t2) {return t1.first < t2;});
          newParameter = lower == (it_param_end) ? index[len_Parameters-1] : *lower;
          newParameter.first = t;
          index[len_Parameters + j] = newParameter;
          events.index[iEvent][3] = 0;
          j++;
        }
      }
    }
    Sort();
  }
};

template<typename T_time, typename T4_container, typename T5, typename T6>
const bool ModelParameterHistory<T_time, T4_container, T5, T6>::has_matrix_param = is_matrix<T4_container>::value;

}

#endif
