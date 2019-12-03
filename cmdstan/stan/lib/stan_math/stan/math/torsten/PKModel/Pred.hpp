#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/torsten/PKModel/ModelParameters.hpp>
#include <stan/math/rev/mat/fun/multiply.hpp>
#include <boost/math/tools/promotion.hpp>
#include <Eigen/Dense>
#include <vector>

namespace torsten{

/**
 * The Event class defines objects that contain the elements of an event,
 * following NONMEM conventions:
 *    time
 *    amt: amount
 *    rate
 *    ii: interdose interval
 *    evid: event identity
 *      (0) observation
 *      (1) dosing
 *      (2) other
 *      (3) reset
 *      (4) reset and dosing
 *    cmt: compartment in which the event occurs
 *    addl: additional doses
 *    ss: steady state approximation (0: no, 1: yes)
 *    keep: if TRUE, save the predicted amount at this event
 *          in the final output of the pred function.
 *    isnew: if TRUE, event was created when pred augmented
 *           the input data set
 */
template <typename T_time, typename T_amt, typename T_rate, typename T_ii>
struct Event{

  T_time time;
  T_amt amt;
  T_rate rate;
  T_ii ii;
  int evid, cmt, addl, ss;
  bool keep, isnew;

  Event() : time(0), amt(0), rate(0), ii(0), cmt(0), addl(0), ss(0), keep(false), isnew(false)
  {}

  Event(T_time p_time, T_amt p_amt, T_rate p_rate, T_ii p_ii, int p_evid,
        int p_cmt, int p_addl, int p_ss, bool p_keep, bool p_isnew) :
    time  (p_time ),
    amt   (p_amt  ),
    rate  (p_rate ),
    ii    (p_ii   ),
    evid  (p_evid ),
    cmt   (p_cmt  ),
    addl  (p_addl ),
    ss    (p_ss   ),
    keep  (p_keep ),
    isnew (p_isnew)
  {}

  /**
   * The function operator is handy when we need to define the same event
   * multiple times, as we might in a FOR loop.
   */
  Event operator()(T_time p_time, T_amt p_amt, T_rate p_rate, T_ii p_ii,
                   int p_evid, int p_cmt, int p_addl, int p_ss, bool p_keep,
                   bool p_isnew) {
    Event newEvent;
    newEvent.time = p_time;
    newEvent.amt = p_amt;
    newEvent.rate = p_rate;
    newEvent.ii = p_ii;
    newEvent.evid = p_evid;
    newEvent.cmt = p_cmt;
    newEvent.addl = p_addl;
    newEvent.ss = p_ss;
    newEvent.keep = p_keep;
    newEvent.isnew = p_isnew;
    return newEvent;
  }

  // Access functions
  T_time get_time() { return time; }
  T_amt get_amt() { return amt; }
  T_rate get_rate() { return rate; }
  T_ii get_ii() { return ii; }
  int get_evid() { return evid; }
  int get_cmt() { return cmt; }
  int get_addl() { return addl; }
  int get_ss() { return ss; }
  bool get_keep() { return keep; }
  bool get_isnew() { return isnew; }

  // declare friends
};

/**
 * Every Torsten function calls Pred.
 *
 * Predicts the amount in each compartment for each event,
 * given the event schedule and the parameters of the model.
 *
 * Proceeds in two steps. First, computes all the events that
 * are not included in the original data set but during which
 * amounts in the system get updated. Secondly, predicts
 * the amounts in each compartment sequentially by going
 * through the augmented schedule of events. The returned pred
 * Matrix only contains the amounts in the event originally
 * specified by the users.
 *
 * This function is valid for all models. What changes from one
 * model to the other are the Pred1 and PredSS functions, which
 * calculate the amount at an individual event.
 *
 * @tparam T_time type of scalar for time
 * @tparam T_amt type of scalar for amount
 * @tparam T_rate type of scalar for rate
 * @tparam T_ii type of scalar for interdose interval
 * @tparam T_parameters type of scalar for the ODE parameters
 * @tparam T_biovar type of scalar for bio-variability parameters
 * @tparam T_tlag type of scalar for lag times parameters
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
 * @param[in] cmt compartment number at each event (starts at 1)
 * @param[in] addl additional dosing at each event
 * @param[in] ss steady state approximation at each event
 * (0: no, 1: yes)
 * @param[in] pMatrix parameters at each event
 * @param[in] addParm additional parameters at each event
 * @parem[in] model basic info for ODE model and evolution operators
 * @param[in] SystemODE matrix describing linear ODE system that
 * defines compartment model. Used for matrix exponential solutions.
 * Included because it may get updated in modelParameters.
 * @return a matrix with predicted amount in each compartment
 * at each event.
 */
template<typename T_time,
        typename T_amt,
        typename T_rate,
        typename T_ii,
        typename T_parameters,
        typename T_biovar,
        typename T_tlag,
        typename F_one,
        typename F_SS>
Eigen::Matrix<typename boost::math::tools::promote_args<T_time, T_amt, T_rate,
  T_ii, typename boost::math::tools::promote_args<T_parameters, T_biovar,
  T_tlag>::type >::type, Eigen::Dynamic, Eigen::Dynamic>
Pred(const std::vector<T_time>& time,
     const std::vector<T_amt>& amt,
     const std::vector<T_rate>& rate,
     const std::vector<T_ii>& ii,
     const std::vector<int>& evid,
     const std::vector<int>& cmt,
     const std::vector<int>& addl,
     const std::vector<int>& ss,
     const std::vector<std::vector<T_parameters> >& pMatrix,
     const std::vector<std::vector<T_biovar> >& biovar,
     const std::vector<std::vector<T_tlag> >& tlag,
     const int& nCmt,
     const std::vector<Eigen::Matrix<T_parameters,
       Eigen::Dynamic, Eigen::Dynamic> >& system,
     const F_one& Pred1,
     const F_SS& PredSS) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using boost::math::tools::promote_args;
  using std::vector;
  using::stan::math::multiply;

  typedef typename promote_args<T_time, T_amt, T_rate, T_ii,
    typename promote_args<T_parameters, T_biovar, T_tlag>::type >::type scalar;
  typedef typename promote_args<T_time, T_amt, T_tlag, T_rate>::type T_tau;
  typedef typename promote_args<T_rate, T_biovar>::type T_rate2;

  // BOOK-KEEPING: UPDATE DATA SETS
  EventHistory<T_time, T_amt, T_rate, T_ii, std::vector<T_parameters>, T_biovar, T_tlag>
    events(nCmt, time, amt, rate, ii, evid, cmt, addl, ss);

  ModelParameterHistory<T_tau, std::vector<T_parameters>, T_biovar, T_tlag>
    parameters(time, pMatrix, biovar, tlag, system);

  events.Sort();
  parameters.Sort();
  int nKeep = events.size();

  events.insert_addl_dose();
  parameters.CompleteParameterHistory(events);

  events.insert_lag_dose();
  // RateHistory<T_time, T_amt, T_rate, T_ii, std::vector<T_parameters>, T_biovar, T_tlag> rates(events, nCmt);
  parameters.CompleteParameterHistory(events);

  Matrix<scalar, 1, Dynamic> zeros = Matrix<scalar, 1, Dynamic>::Zero(nCmt);
  Matrix<scalar, 1, Dynamic> init = zeros;

  // COMPUTE PREDICTIONS
  Matrix<scalar, Dynamic, Dynamic>
    pred = Matrix<scalar, Dynamic, Dynamic>::Zero(nKeep, nCmt);

  scalar Scalar = 1;  // trick to promote variables to scalar

  T_tau dt, tprev = events.get_time(0);
  Matrix<scalar, Dynamic, 1> pred1;
  Event<T_tau, T_amt, T_rate, T_ii> event;
  ModelParameters<T_tau, T_parameters, T_biovar, T_tlag> parameter;
  int iRate = 0, ikeep = 0;

  for (int i = 0; i < events.size(); i++) {
    event = events.GetEvent(i);

    // Use index iRate instead of i to find rate at matching time, given there
    // is one rate per time, not per event.
    // if (rates.get_time(iRate) != events.get_time(i)) iRate++;
    std::vector<T_rate2> rate2(nCmt);
    for (int j = 0; j < nCmt; ++j) {
      // rate2[j] = rates.Rates[iRate].rate[j] * parameters.GetValueBio(i, j);
    }

    parameter = parameters.GetModelParameters(i);

    if ((event.get_evid() == 3) || (event.get_evid() == 4)) {  // reset events
      dt = 0;
      init = zeros;
    } else {
      T_tau t_new = event.get_time();
      dt = t_new - tprev;
      pred1 = Pred1(dt, parameter, init, rate2.get_rate());
      init = pred1;
    }

    if (((event.get_evid() == 1 || event.get_evid() == 4)
      && (event.get_ss() == 1 || event.get_ss() == 2)) ||
      event.get_ss() == 3) {  // steady state event
      pred1 = multiply(PredSS(parameter,
                              parameters.GetValueBio(i, event.get_cmt() - 1)
                                * event.get_amt(),
                              event.get_rate(), event.get_ii(),
                              event.get_cmt()),
                       Scalar);

      // the object PredSS returns doesn't always have a scalar type. For
      // instance, PredSS does not depend on tlag, but pred does. So if
      // tlag were a var, the code must promote PredSS to match the type
      // of pred1. This is done by multiplying predSS by a Scalar.

      if (event.get_ss() == 2) init += pred1;  // steady state without reset
      else
        init = pred1;  // steady state with reset (ss = 1)
    }

    if (((event.get_evid() == 1) || (event.get_evid() == 4)) &&
      (event.get_rate() == 0)) {  // bolus dose
      init(0, event.get_cmt() - 1)
        += parameters.GetValueBio(i, event.get_cmt() - 1) * event.get_amt();
    }

    if (event.get_keep()) {
      pred.row(ikeep) = init;
      ikeep++;
    }
  tprev = event.get_time();
  }

  return pred;
}

}

#endif
