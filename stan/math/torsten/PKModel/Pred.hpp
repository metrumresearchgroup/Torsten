#ifndef STAN_MATH_TORSTEN_PKMODEL_PRED_HPP
#define STAN_MATH_TORSTEN_PKMODEL_PRED_HPP

#include <Eigen/Dense>
#include <vector>

/**
 * Every Torsten function calls Pred.
 *
 * Predicts the amount in each compartment for each event,
 * given the event schedule and the parameters of the model.
 *
 * Proceeds in two steps. First, computes all the events that
 * are not included in the original data set but during which
 * the amounts in the system are updated. Secondly, predicts
 * the amounts in each compartment sequentially by going
 * through the augmented schedule of events. The returned pred
 * Matrix only contains the amounts in the event originally
 * specified by the users.
 *
 * This function is valid for all models. What changes from one
 * model to the other are the Pred1 and PredSS functions, which
 * calculate the amount at an individual event.
 *
 * @tparam T_parameters type of scalar for the model parameters
 * @tparam T_time type of scalar for time
 * @tparam T_amt type of scalar for amount
 * @tparam T_rate type of scalar for rate
 * @tparam T_ii type of scalar for interdose interval
 * @tparam F type of ODE system function
 * @param[in] pMatrix parameters at each event
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
 * @param[in] ss steady state approximation at each event
 * (0: no, 1: yes)
 * @parem[in] model basic structural information on compartment
 * model
 * @param[in] f functor for base ordinary differential equation
 * that defines compartment model. Used for ODE integrators
 * (optional).
 * @param[in] SystemODE matrix describing linear ODE system that
 * defines compartment model. Used for matrix exponential solutions
 * (optional).
 * @return a matrix with predicted amount in each compartment
 * at each event.
 */
template <typename T_parameters, typename T_time, typename T_amt,
  typename T_rate, typename T_ii, typename F, typename T_system>
Eigen::Matrix<typename boost::math::tools::promote_args<T_parameters, T_time,
  T_amt, T_rate, typename boost::math::tools::promote_args<T_ii, T_system>::
  type >::type, Eigen::Dynamic, Eigen::Dynamic>
Pred(const std::vector<std::vector<T_parameters> >& pMatrix,
     const std::vector<T_time>& time,
     const std::vector<T_amt>& amt,
     const std::vector<T_rate>& rate,
     const std::vector<T_ii>& ii,
     const std::vector<int>& evid,
     const std::vector<int>& cmt,
     const std::vector<int>& addl,
     const std::vector<int>& ss,
     PKModel model,
     const F& f,
     const std::vector<Eigen::Matrix<T_system,
       Eigen::Dynamic, Eigen::Dynamic> >& system) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using boost::math::tools::promote_args;
  // using namespace stan::math;
  using std::vector;

  typedef typename promote_args<T_parameters, T_time, T_amt, T_rate,
    typename promote_args<T_ii, T_system>::type >::type scalar;

  /////////////////////////////////////////////////////////////////////////////
  // BOOK-KEEPING: UPDATE DATA SETS

  int nParameter = model.GetNParameter(),
    nCmt = model.GetNCmt(),
    F1Index = model.GetF1Index(),
    tlag1Index = model.GetTLagIndex();

  EventHistory<scalar, scalar, scalar, scalar>
    events(time, amt, rate, ii, evid, cmt, addl, ss);
  ModelParameterHistory<scalar, T_parameters, T_system>
    parameters(time, pMatrix, system);
  RateHistory<scalar, scalar> rates;

  events.Sort();
  parameters.Sort();

  int nKeep = events.get_size();
  int np = pMatrix[0].size() * pMatrix.size();

  assert((np > 0) && (np % nParameter == 0));
  events.AddlDoseEvents();

  parameters.CompleteParameterHistory(events);

  vector<int> tlagIndexes;
  tlagIndexes.assign(nCmt, 0);
  vector<int> tlagCmts;
  tlagCmts.assign(nCmt, 0);
  for (int i = 0; i < nCmt; i++) {
    tlagIndexes[i] = tlag1Index + i;
    tlagCmts[i] = i + 1;
  }

  events.AddLagTimes(parameters, tlagIndexes, tlagCmts);
  rates.MakeRates(events, nCmt);
  parameters.CompleteParameterHistory(events);

  Matrix<scalar, 1, Dynamic> zeros = Matrix<scalar, 1, Dynamic>::Zero(nCmt);
  Matrix<scalar, 1, Dynamic> init = zeros;

  // Construct the output matrix pred
  // The matrix needs to be a dynamically-sized matrix to be returned
  // by the function. The matrix is resized, but the resize function
  // of the eigen library is destructive. This means arbitrary values
  // are assigned to the elements in the matrix.
  // For the "COMPUTE PREDICTIONS" half of the function to run properly,
  // we need to set the values of each element to 0.
  Matrix<scalar, Dynamic, Dynamic>
    pred = Matrix<scalar, Dynamic, Dynamic>::Zero(nKeep, nCmt);

  /////////////////////////////////////////////////////////////////////////////
  // COMPUTE PREDICTIONS

  scalar dt, tprev = events.get_time(0);
  Matrix<scalar, 1, Dynamic> pred1;
  Event<scalar, scalar, scalar, scalar> event;
  ModelParameters<scalar, T_parameters, T_system> parameter;
  Rate<scalar, scalar> rate2;
  int iRate = 0;
  int ikeep = 0;

  for (int i = 0; i < events.get_size(); i++) {
    event = events.GetEvent(i);

    // Use index iRate instead of i to find rate at matching time, given there
    // is one rate per time, not per event.
    if (rates.get_time(iRate) != events.get_time(i)) iRate++;
    rate2 = rates.GetRate(iRate);

    for (int j = 0; j < nCmt; j++)
      rate2.rate[j] *= parameters.GetValue(i, F1Index + j);

    parameter = parameters.GetModelParameters(i);

    if ((event.get_evid() == 3) || (event.get_evid() == 4)) {  // reset events
      dt = 0;
      init = zeros;
    } else {
      dt = event.get_time() - tprev;
      pred1 = Pred1(dt, parameter, init, rate2.get_rate(), f);
      init = pred1;
    }

    if (((event.get_evid() == 1) || (event.get_evid() == 4))
      && (((event.get_ss() == 1) || (event.get_ss() == 2)) ||
      (event.get_ss() == 3))) {  // steady dose event
      pred1 = PredSS(parameter, parameters.GetValue(i, F1Index
        + event.get_cmt() - 1) * event.get_amt(), event.get_rate(),
        event.get_ii(), event.get_cmt(), f);

      if (event.get_ss() == 2) init += pred1;  // steady state without reset
      else
        init = pred1;  // steady state with reset (ss = 1)
    }

    if (((event.get_evid() == 1) || (event.get_evid() == 4)) &&
      (event.get_rate() == 0))  // bolus dose
      init(0, event.get_cmt() - 1) += parameters.GetValue(i, F1Index
      + event.get_cmt() - 1) * event.get_amt();

    if (event.get_keep()) {
      pred.row(ikeep) = init;
      ikeep++;
    }
  tprev = event.get_time();
  }
  return pred;
}

#endif
