#ifndef STAN_MATH_TORSTEN_PKMODEL_MODELPARAMETERS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_MODELPARAMETERS_HPP

#include <stan/math/prim/mat/fun/to_array_1d.hpp>
#include <stan/math/torsten/PKModel/Event.hpp>
#include <stan/math/torsten/PKModel/ExtractVector.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace torsten {

template<typename T_time,
         typename T_parameters,
         typename T_biovar,
         typename T_tlag>
  struct ModelParameterHistory;
/**
 * The ModelParameters class defines objects that contain the parameters of
 * a model at a given time.
 */
template<typename T_time,
         typename T_parameters,
         typename T_biovar,
         typename T_tlag>
struct ModelParameters {
  T_time time_;
  std::vector<T_parameters> theta_;
  int nrow, ncol;
  std::vector<T_biovar> biovar_;
  std::vector<T_tlag> tlag_;

  ModelParameters() {}

  ModelParameters(const T_time& time,
                  const std::vector<T_biovar>& biovar,
                  const std::vector<T_tlag>& tlag,
                  const Eigen::Matrix<T_parameters, Eigen::Dynamic, Eigen::Dynamic>& K)
    : time_(time), theta_(K.size()), nrow(K.rows()), ncol(K.cols()), biovar_(biovar), tlag_(tlag)
  {
    Eigen::Matrix<T_parameters, -1, -1>::Map(theta_.data(), nrow, ncol) = K;
  }

  ModelParameters(const T_time& time,
                  const std::vector<T_parameters>& theta,
                  const std::vector<T_biovar>& biovar,
                  const std::vector<T_tlag>& tlag)
    : time_(time), theta_(theta), biovar_(biovar), tlag_(tlag) {}

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
  T_time get_time() const { return time_; }
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

/**
 * The ModelParameterHistory class defines objects that contain a vector
 * of ModelParameters, along with a series of functions that operate on
 * them.
 */
template<typename T_time,
         typename T_parameters,
         typename T_biovar,
         typename T_tlag>
struct ModelParameterHistory{
  const bool has_matrix_param;
  std::vector<ModelParameters<T_time, T_parameters, T_biovar, T_tlag> > MPV_;

  template<typename T0, typename T1, typename T2, typename T3>
  ModelParameterHistory(std::vector<T0> time,
                        std::vector<std::vector<T1> > theta,
                        std::vector<std::vector<T2> > biovar,
                        std::vector<std::vector<T3> > tlag) :
    has_matrix_param(false),
    MPV_(std::max(theta.size(), std::max(biovar.size(), tlag.size())))
  {
    using std::max;
    int nParameters = max(theta.size(), max(biovar.size(), tlag.size()));
    int j, k, l;
    // FIX ME - is this the most efficient way of storing data?
    for (int i = 0; i < nParameters; i++) {
      (theta.size() == 1) ? j = 0 : j = i;
      (biovar.size() == 1) ? k = 0 : k = i;
      (tlag.size() == 1) ? l = 0 : l = i;
       MPV_[i] = ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
         (time[i], theta[j], biovar[k], tlag[l]);
    }
  }

  template<typename T0, typename T1, typename T2, typename T3>
  ModelParameterHistory(std::vector<T0> time,
                        std::vector< Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> > K,
                        std::vector<std::vector<T2> > biovar,
                        std::vector<std::vector<T3> > tlag) :
    has_matrix_param(true),
    MPV_(std::max(K.size(), std::max(biovar.size(), tlag.size())))
  {
    using std::max;
    int nParameters = max(K.size(), max(biovar.size(), tlag.size()));
    int k, l, m;
    // FIX ME - is this the most efficient way of storing data?
    for (int i = 0; i < nParameters; i++) {
      (biovar.size() == 1) ? k = 0 : k = i;
      (tlag.size() == 1) ? l = 0 : l = i;
      (K.size() == 1) ? m = 0 : m = i;
       MPV_[i] = ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
         (time[i], biovar[k], tlag[l], K[m]);
    }
  }

  template<typename T0, typename T1, typename T2, typename T3>
  ModelParameterHistory(std::vector<T0> time,
                        std::vector<std::vector<T1> > theta,
                        std::vector< Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> > K,
                        std::vector<std::vector<T2> > biovar,
                        std::vector<std::vector<T3> > tlag) :
    has_matrix_param(theta.empty()),
    MPV_(has_matrix_param ? std::max(K.size(), std::max(biovar.size(), tlag.size()))
         : std::max(theta.size(), std::max(biovar.size(), tlag.size())))
  {
    using std::max;
    if (has_matrix_param) {
      int nParameters = max(K.size(), max(biovar.size(), tlag.size()));
      int k, l, m;
      // FIX ME - is this the most efficient way of storing data?
      for (int i = 0; i < nParameters; i++) {
        (biovar.size() == 1) ? k = 0 : k = i;
        (tlag.size() == 1) ? l = 0 : l = i;
        (K.size() == 1) ? m = 0 : m = i;
        MPV_[i] = ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
          (time[i], biovar[k], tlag[l], K[m]);
      }      
    } else {
      int nParameters = max(theta.size(), max(biovar.size(), tlag.size()));
      int j, k, l;
      // FIX ME - is this the most efficient way of storing data?
      for (int i = 0; i < nParameters; i++) {
        (theta.size() == 1) ? j = 0 : j = i;
        (biovar.size() == 1) ? k = 0 : k = i;
        (tlag.size() == 1) ? l = 0 : l = i;
        MPV_[i] = ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
          (time[i], theta[j], biovar[k], tlag[l]);
      }
    }
  }

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
  template<typename T0, typename T1, typename T2, typename T3>
  ModelParameterHistory(int ibegin, int isize,
                        std::vector<T0> time,
                        int ibegin_theta, int isize_theta,
                        std::vector<std::vector<T1> > theta,
                        std::vector< Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> > K,
                        int ibegin_biovar, int isize_biovar,
                        std::vector<std::vector<T2> > biovar,
                        int ibegin_tlag, int isize_tlag,
                        std::vector<std::vector<T3> > tlag) :
    has_matrix_param(theta.empty()),
    MPV_(std::max(isize_theta, std::max(isize_biovar, isize_tlag)))
  {
    using std::max;
    int n1 = isize_theta;
    int n2 = isize_biovar;
    int n3 = isize_tlag;
    int nParameters = max(n1, max(n2, n3));
    static const char* caller = "ModelParameterHistory::ModelParameterHistory";
    stan::math::check_greater_or_equal(caller, "isize", isize, nParameters);
    stan::math::check_greater_or_equal(caller, "time size", time.size(), size_t(ibegin + nParameters - 1));
    int j, k, l;
    for (int i = 0; i < nParameters; i++) {
      (n1 == 1) ? j = ibegin_theta  : j = ibegin_theta  + i;
      (n2 == 1) ? k = ibegin_biovar : k = ibegin_biovar + i;
      (n3 == 1) ? l = ibegin_tlag   : l = ibegin_tlag   + i;
      if (has_matrix_param) {
        MPV_[i] = ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
          (time[i + ibegin], biovar[k], tlag[l], K[j]);
      } else {
        MPV_[i] = ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
          (time[i + ibegin], theta[j], biovar[k], tlag[l]);
      }
    }
  }

  ModelParameters<T_time, T_parameters, T_biovar, T_tlag>
    GetModelParameters(int i) const {
      return MPV_[i];
  }

  /**
   * MPV.size gives us the number of events.
   * MPV[i].RealParameters.size gives us the number of
   * ODE parameters for the ith event.
   * 
   * FIX ME - rename this GetValueTheta
   */
  T_parameters GetValue(int iEvent, int iParameter) {
    assert((iEvent >= 0) && ((size_t) iEvent < MPV_.size()));
    assert((iParameter >= 0) && ((size_t) iParameter
      < MPV_[iEvent].theta_.size()));
    return MPV_[iEvent].theta_[iParameter];
  }

  T_biovar GetValueBio(int iEvent, int iParameter) const {
    assert(iEvent >= 0 && (size_t) iEvent < MPV_.size());
    assert(iParameter >= 0 && (size_t) iParameter
             < MPV_[iEvent].biovar_.size());
    return MPV_[iEvent].biovar_[iParameter];
  }

  T_tlag GetValueTlag(int iEvent, int iParameter) {
    assert(iEvent >= 0 && (size_t) iEvent < MPV_.size());
    assert(iParameter >= 0 && (size_t) iParameter
             < MPV_[iEvent].tlag_.size());
    return MPV_[iEvent].tlag_[iParameter];
  }

  void InsertModelParameters(ModelParameters<T_time, T_parameters,
    T_biovar, T_tlag> M) {
    MPV_.push_back(M);
  }

  int get_size() {
    return MPV_.size();
  }

  struct by_time {
    bool operator()(ModelParameters<T_time, T_parameters, T_biovar,
                                    T_tlag> const &a,
                    ModelParameters<T_time, T_parameters, T_biovar,
                                    T_tlag> const &b) {
      return a.time_ < b.time_;
    }
  };

  void Sort() {
    std::sort(MPV_.begin(), MPV_.end(), by_time());
  }

  bool Check() {
  // check that elements are in chronological order.
    int i = MPV_.size() - 1;
    bool ordered = true;

    while (i > 0 && ordered) {
      ordered = (MPV_[i].time_ >= MPV_[i-1].time_);
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
   * (Last Observation Carried Forward).
   *
   * Events and Parameters are sorted at the end of the procedure.
   *
   * @param[in] parameters at each event
   * @param[in] events elements (following NONMEM convention) at each event
   * @return - modified parameters and events.
   */
  template<typename T0, typename T1, typename T2, typename T3>
  void CompleteParameterHistory(torsten::EventHistory<T0, T1, T2, T3>& events) {
    int nEvent = events.size();
    assert(nEvent > 0);
    int len_Parameters = MPV_.size();  // numbers of events for which parameters
                                       // are determined
    assert(len_Parameters > 0);

    if (!Check()) Sort();
    if (!events.Check()) events.Sort();
    MPV_.resize(nEvent);

    int iEvent = 0;
    for (int i = 0; i < len_Parameters - 1; i++) {
      while (events.isnew(iEvent)) iEvent++;  // skip new events
      assert(MPV_[i].time_ == events.time(iEvent));  // compare time of
                                                         // "old' events to
                                                         // time of
                                                         // parameters.
      iEvent++;
    }

    if (len_Parameters == 1)  {
      for (int i = 0; i < nEvent; i++) {
        // FIX ME - inefficient data storage
        MPV_[i].theta_ = MPV_[0].theta_;
        MPV_[i].biovar_ = MPV_[0].biovar_;
        MPV_[i].tlag_ = MPV_[0].tlag_;
        MPV_[i].nrow = MPV_[0].nrow;
        MPV_[i].ncol = MPV_[0].ncol;
        MPV_[i].time_ = events.time(i);
        events.Events[i].isnew = false;
      }
    } else {  // parameters are event dependent.
      std::vector<T_time> times(nEvent, 0);
      for (int i = 0; i < nEvent; i++) times[i] = MPV_[i].time_;
      iEvent = 0;

      int k, j = 0;
      ModelParameters<T_time, T_parameters, T_biovar, T_tlag> newParameter;

      for (int i = 0; i < nEvent; i++) {
        while (events.isnew(iEvent)) {
          /* Three cases:
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
           */
          // Find the index corresponding to the time of the new event in the
          // times vector.
          k = SearchReal(times, len_Parameters - 1, events.time(iEvent));

          if ((k == len_Parameters) ||
            (events.time(iEvent) == MPV_[k - 1].time_))
            newParameter = GetModelParameters(k - 1);
          else
            newParameter = GetModelParameters(k);

          newParameter.time_ = events.time(iEvent);
          MPV_[len_Parameters + j] = newParameter;
          events.Events[iEvent].isnew = false;
          if (iEvent < nEvent - 1) iEvent++;
          j++;
        }

        if (iEvent < nEvent - 1) iEvent++;
      }
    }
    if (!Check()) Sort();
  }
};

} 

#endif
