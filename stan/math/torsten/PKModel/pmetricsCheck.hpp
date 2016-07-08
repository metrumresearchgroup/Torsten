// version 0.8

/**
 * Checks that the arguments the user inputs in the Torsten
 * functions are valid. 
 * 
 * @tparam T0 type of scalars for the model parameters.
 * @tparam T1 type of scalar for time of events. 
 * @tparam T2 type of scalar for amount at each event.
 * @tparam T3 type of scalar for rate at each event.
 * @tparam T4 type of scalar for inter-dose inteveral at each event.
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
 * @param[in] ss steady state approximation at each event (0: no, 1: yes)
 * @param[in] function The name of the function for which the check is being 
 *                     performed. 
 * @param[in] model object that contains basic structural information about 
 *                  a compartment model.
 * @return void
 * 
 */


#include <stan/model/model_header.hpp>
#include "PKModel_class.hpp"
#include "Event.hpp"
#include "Rate.hpp"
#include "ModelParameters.hpp"

using std::vector;
using Eigen::Dynamic;

template <typename T0, typename T1, typename T2, typename T3, typename T4> 
void pmetricsCheck(const vector< Matrix<T0, Dynamic, 1> >& pMatrix, 
                   const vector<T1>& time,
                   const vector<T2>& amt,
                   const vector<T3>& rate,
                   const vector<T4>& ii,
                   const vector<int>& evid,
                   const vector<int>& cmt,
                   const vector<int>& addl,
                   const vector<int>& ss,
                   const char* function,
                   PKModel& model) 
{
  
  // TEST ARGUMENTS FOR EVENT SCHEDULE
  
  if(!(time.size()>0)) invalid_argument(function,
     "length of time vector,", time.size(), "", "needs to be positive and greater than 0!");
  
  std::string message = ", but must be the same as the length of the time array: " 
  + boost::lexical_cast<string>(time.size()) + "!"; 
  const char* length_error = message.c_str();
  
  if (!(amt.size()==time.size())) invalid_argument(function,
      "the length of the amount (amt) array is", amt.size(), "",
      length_error);
  if (!(rate.size()==time.size())) invalid_argument(function,
      "the length of the rate array is", rate.size(), "",
      length_error);
  if (!(evid.size()==time.size())) invalid_argument(function,
      "the length of the event ID (evid) array is", evid.size(), "",
      length_error);
  if (!(cmt.size()==time.size())) invalid_argument(function,
      "the length of the compartment (cmt) array is", cmt.size(), "",
      length_error);  
  
  
  std::string message2 = ", but must be either 1 or the same as the length of the time array: " 
  + boost::lexical_cast<string>(time.size()) + "!"; 
  const char* length_error2 = message2.c_str();  
  
  if (!(ii.size()==time.size())||(ii.size()==1)) invalid_argument(function,
      "the length of the interdose interval (ii) array is", ii.size(), "",
      length_error2);
  if (!(addl.size()==time.size())||(addl.size()==1)) invalid_argument(function,
      "the length of the additional dosing (addl) array is", ii.size(), "",
      length_error2);
  if (!(ss.size()==time.size())||(ss.size()==1)) invalid_argument(function,
      "the length of the steady state approximation (ss) array is", ss.size(), "",
      length_error2);
  
  
  std::string message3 = ", but must be the same as the length of the additional dosing (addl) array: " 
  + boost::lexical_cast<string>(addl.size()) + "!"; 
  const char* length_error3 = message3.c_str();     
  if (!(ss.size()==time.size())||(ss.size()==1)) invalid_argument(function,
      "the length of steady state approximation (ss) array is", ss.size(), "",
      length_error3);
  
  
  
  
  // TEST ARGUMENTS FOR PARAMETERS
  
  if (!((pMatrix.size()==time.size())||(pMatrix.size()==1))) invalid_argument(function,
      "length of the parameter (2d) vector,", pMatrix.size(), "",
      length_error2);
  if(!(pMatrix[0].size()>0)) invalid_argument(function,
     "the number of parameters per event is", pMatrix[0].size(), "", "but must be greater than 0!");
  
  std::string message4 = ", but must equal the number of parameters in the model: " 
  + boost::lexical_cast<string>(model.GetNParameter()) + "!"; 
  const char* length_error4 = message4.c_str();    
  
  if (!(pMatrix[0].size() == model.GetNParameter())) invalid_argument(function,
      "The number of parameters per event (length of a vector in the first argument; reminder, the first argument is an array of vectors) is", pMatrix[0].size(), "", 
      length_error4);
  
}