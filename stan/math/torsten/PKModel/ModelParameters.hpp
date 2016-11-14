#ifndef STAN_MATH_TORSTEN_PKMODEL_MODELPARAMETERS_HPP
#define STAN_MATH_TORSTEN_PKMODEL_MODELPARAMETERS_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/ExtractVector.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>

// FIX ME: using statements should be declared 
// within the scope of functions
using Eigen::Matrix;
using Eigen::Dynamic;
using std::vector;
using stan::math::var;
using boost::math::tools::promote_args;

template<typename T_time, typename T_parameters, typename T_system> 
  class ModelParameterHistory; 

/**
 * The ModelParameters class defines objects that contain the parameters of a
 * compartment model at a given time. 
 */
template<typename T_time, typename T_parameters, typename T_system>
class ModelParameters {

private:
	T_time time;
	vector<T_parameters> RealParameters;
	Matrix<T_system, Dynamic, Dynamic> K;

public:
	ModelParameters() {
		vector<T_parameters> v(1, 0);
		Matrix<T_system, Dynamic, Dynamic> K_(0, 0);
		time = 0;
		RealParameters = v;
		K = K_;
	}

	ModelParameters(const T_time& p_time, 
	                const vector<T_parameters>& p_RealParameters,
	                const Matrix<T_system, Dynamic, Dynamic>& p_K) {
		time = p_time;
		RealParameters = p_RealParameters;
		K = p_K;
	}

	int CountParameters() const { 
	  return RealParameters.size();
	}
	
	Matrix<T_system, Dynamic, Dynamic> RateMatrix() const {
	  return K;
	}

	void Print() {
		std::cout << time << " ";
		for(int i = 0; i < RealParameters.size(); i++)
		  std::cout << RealParameters[i] << " ";
		if (K.rows() != 0) std::cout << K;
		std::cout << std::endl;
	}
	
	//declare friends
	friend class ModelParameterHistory<T_time, T_parameters, T_system>;
	
	// When befriending a template function, cannot do a partial instantiation.
	// Either all or none of the terms must be instantiate. Since we need to use 
	// template for Rate, we do not instantiate any of the template arguments. 	
	template<typename T_1, typename T_2, typename T_3, typename T_4> 
	friend Matrix<typename promote_args< T_1, T_2, T_3>::type, 1, Dynamic> 
	Pred1_one(const T_1&,
			  const ModelParameters<T_1, T_3, T_4>&,
			  const Matrix<typename promote_args<T_1, T_2, T_3>::type,
			    1, Dynamic>& init,
			  const vector<T_2>& rate);
	
	template<typename T_1, typename T_2, typename T_3, typename T_4> 
	friend Matrix<typename promote_args< T_1, T_2, T_3>::type, 1, Dynamic> 
	Pred1_two(const T_1&,
			  const ModelParameters<T_1, T_3, T_4>&,
			  const Matrix<typename promote_args<T_1, T_2, T_3>::type,
			    1, Dynamic>& init,
			  const vector<T_2>& rate);
			  
	template<typename T_1, typename T_2, typename T_3, typename T_4, typename F>
	friend
	Matrix<typename promote_args< T_1, T_2, T_3>::type, 1, Dynamic> 
	Pred1_general_solver(const T_1& dt,
		  				 const ModelParameters<T_1, T_3, T_4>& parameter, 
		 				 const Matrix<typename promote_args<T_1, T_2, T_3>::type, 1,
		 				   Dynamic>& init, 
		  				 const vector<T_2>& rate,
		  				 const F& f);
	
	template<typename T_1, typename T_2, typename T_3, typename T_4>
    friend Matrix<typename promote_args< T_1, T_3, T_4>::type, 1, Dynamic> 
	Pred1_linCpt(const T_1&,
			     const ModelParameters<T_1, T_2, T_3>&,
			     const Matrix<typename promote_args<T_1, T_3, T_4>::type, 1,
			       Dynamic>& init,
			     const vector<T_4>& rate);
	
	template<typename T_1, typename T_2, typename T_3, typename T_4, typename T_5,
	  typename T_6>
	friend Matrix<typename promote_args< T_1, T_2, T_3, 
	  typename promote_args<T_4, T_5>::type>::type, 1, Dynamic>
	PredSS_one(const ModelParameters<T_1, T_5, T_6>& parameter, 
               const T_2& amt, 
               const T_3& rate,
               const T_4& ii, 
               const int& cmt);		
	
	template<typename T_1, typename T_2, typename T_3, typename T_4, typename T_5,
	  typename T_6>
	friend Matrix<typename promote_args< T_1, T_2, T_3,
	  typename promote_args<T_4, T_5>::type>::type, 1, Dynamic>	   
	PredSS_two(const ModelParameters<T_1, T_5, T_6>& parameter, 
			   const T_2& amt, 
			   const T_3& rate,
			   const T_4& ii, 
		   	   const int& cmt);	

	template <typename T_0, typename T_1, typename T_2, typename T_3,
	  typename T_4, typename F, typename T_5>
	friend
	Matrix<typename promote_args<T_0, T_1, T_2, T_3,
	 typename promote_args<T_4, T_5>::type >::type, Dynamic, Dynamic>
	Pred(const vector< Matrix<T_0, Dynamic, 1> >& pMatrix,
     	 const vector<T_1>& time,
     	 const vector<T_2>& amt, 
     	 const vector<T_3>& rate,
    	 const vector<T_4>& ii,
  	   	 const vector<int>& evid,
    	 const vector<int>& cmt,
         const vector<int>& addl,
         const vector<int>& ss,
     	 PKModel model,
     	 const F& f,
     	 const std::vector<Matrix<T_5, Dynamic, Dynamic> >& system);
};


/**
 * The ModelParameterHistory class defines objects that contain a vector 
 * of ModelParameters, along with a series of functions that operate on
 * them.  
 */
template<typename T_time, typename T_parameters, typename T_system>
class ModelParameterHistory{

private:
	vector< ModelParameters<T_time, T_parameters, T_system> > MPV;
	
public:
	ModelParameterHistory(vector<T_time> p_time, 
						  vector< Matrix<T_parameters, Dynamic, 1> > p_RealParameters,
						  vector< Matrix<T_system, Dynamic, Dynamic> > p_K) {
		int nParameters = std::max(p_RealParameters.size(), p_K.size()); 
		vector<T_parameters> col;
		MPV.resize(nParameters);
		int j, k;
		for(int i = 0; i < nParameters; i++) {
		    (p_RealParameters.size() == 1) ? j = 0 : j = i;
			col = ExtractVector(p_RealParameters[j], 0, "col");
			(p_K.size() == 1) ? k = 0 : k = i;
			MPV[i] = ModelParameters<T_time, T_parameters, T_system>
		      (p_time[i], col, p_K[k]);
		}
	}

	ModelParameters<T_time, T_parameters, T_system> GetModelParameters(int i) {
		ModelParameters<T_time, T_parameters, T_system>
		  newPara(MPV[i].time, MPV[i].RealParameters, MPV[i].K);
		return newPara;
	}
	
	T_parameters GetValue(int iEvent, int iParameter) {
		// MPV.size gives us the number of events
		// MPV[0].RealParameters.size gives us the number of parameters for the first event. 
		// The code assumes this number is the same for all events
		assert((iEvent >= 0) && (iEvent < MPV.size()));
		assert((iParameter >= 0) && (iParameter < MPV[0].RealParameters.size()));
		return MPV[iEvent].RealParameters[iParameter];
	}	

	void InsertModelParameters(ModelParameters<T_time, T_parameters, T_system> p_M) {
	  MPV.push_back(p_M);
	}

	int get_size() {
	  return MPV.size();
	}

	struct by_time {
		bool operator()(ModelParameters<T_time, T_parameters, T_system> const &a, 
		                ModelParameters<T_time, T_parameters, T_system> const &b) {
			return a.time < b.time;
		}
	};

	void Sort() {
	  std::sort(MPV.begin(), MPV.end(), by_time());
	}

	bool Check() {
		// check that elements are in chronological order. 
		int i = MPV.size() - 1;
		bool ordered = true;
		
		while((i > 0) && (ordered)) {
			ordered = (MPV[i].time >= MPV[i-1].time);
			i--;
		}
		return ordered;
	}

	void Print(int j) {
		std::cout << MPV[j].time;
		for (int i = 0; i < MPV[j].RealParameters.size(); i++)
		  print(MPV[j].RealParameters[i], false);
		std::cout << std::endl;
	}	

	/**
 	 * COMPLETE MODEL PARAMETERS
 	 * 
 	 * Completes parameters so that it contains model parameters for each event in events. 
 	 * If parameters contains only one set of parameters (case where the parameters are 
 	 * constant), this set is replicated for each event in events.
 	 * Otherwise a new parameter vector is added for each new event (isnew = true). This
  	 * parameter vector is identical to the parameter vector at the subsequent event. If the
  	 * new event occurs at a time posterior to the time of the last event, than the 
  	 * new vector parameter equals the parameter vector of the last event. This amounts to
 	 * doing an LOCF (Last Observation Carried Forward). 
 	 *
 	 * Events and Parameters are sorted at the end of the procedure.
 	 *
 	 * @param[in] parameters at each event
 	 * @param[in] events elements (following NONMEM convention) at each event
 	 * @return - modified parameters and events.
 	 *
 	 */
	template<typename T_amt, typename T_rate, typename T_ii>
	void CompleteParameterHistory(EventHistory<T_time, T_amt, T_rate, T_ii>& events) {
		int nEvent, nParameters, len_Parameters;
	
		nEvent = events.get_size();
		nParameters = MPV[0].RealParameters.size(); // number of parameters per event
		len_Parameters = MPV.size(); // number of events for which
		                             // parameters are determined
		assert(nEvent > 0);
		assert(nParameters > 0);
		assert(len_Parameters > 0);
		if (!Check()) Sort();
		if (!events.Check()) events.Sort();

		MPV.resize(nEvent);
		int iEvent = 0;
		
		for (int i = 0; i < len_Parameters-1; i++) {
			while(events.Events[iEvent].isnew) iEvent++; // skip new events 
			assert(MPV[i].time == events.Events[iEvent].time); // compare time of "old" events 
															   // to time of parameters.  
			iEvent++;
		}
		
		if (len_Parameters == 1)  {
			for(int i = 0; i < nEvent; i++) {
				MPV[i].RealParameters = MPV[0].RealParameters;
				MPV[i].K = MPV[0].K;
				MPV[i].time = events.Events[i].time;
				events.Events[i].isnew = false;
			}
		}
		else {  // address the case where the parameters are event dependent.
			vector<T_time> times(nEvent, 0);
			for(int i = 0; i < nEvent; i++) times[i] = MPV[i].time;
			iEvent = 0;	
		
		int k, j = 0;
		ModelParameters<T_time, T_parameters, T_system> newParameter;
		
			for(int i = 0; i < nEvent; i++) {
				while(events.Events[iEvent].isnew) {
					/* Three cases:
					 * (a) The time of the new event is higher than the time of the last
					 *     parameter vector in parameters (k = len_parameters). 
					 *     Create a parameter vector at the the time of the new event, 
					 *     with the parameters of the last.
					 *     (Last Observation Carried Forward)
					 * (b) The time of the new event matches the time of a parameter vector
					 *     in parameters. In this case, this parameter vector is replicated.
					 * (c) (a) is not verified and no parameter vector occurs at the time 
					 *     of the new event. A new parameter vector is created at the time 
					 *     of the new event, and its parameters are equal to the parameters
					 *     of the subsequent parameter vector in parameters.
					 */
					 
					// find the index corresponding to the time of the new event in the 
					// times vector.  
					k = SearchReal(times, len_Parameters-1, events.Events[iEvent].time);
					
					if((k == len_Parameters)||(events.Events[iEvent].time == MPV[k-1].time)) 
									newParameter = GetModelParameters(k-1);
					else newParameter = GetModelParameters(k);
					
					newParameter.time = events.Events[iEvent].time;
					MPV[len_Parameters + j] = newParameter; // Include new element in 
															// MPV 
					events.Events[iEvent].isnew = false;
					if (iEvent < nEvent-1) iEvent++;
					j++;
				}
			
			if (iEvent < nEvent-1) iEvent++;
			}
		}
		if(!Check()) Sort();
	}

	// declare friends
	friend class ModelParameters<T_time, T_parameters, T_system>;
	template<typename T1, typename T2, typename T3, typename T4> friend class Events;
	
	template <typename T_0, typename T_1, typename T_2, typename T_3,
	  typename T_4, typename F, typename T_5>
	friend
	Matrix<typename promote_args<T_0, T_1, T_2, T_3,
	 typename promote_args<T_4, T_5>::type >::type, Dynamic, Dynamic>
	Pred(const vector< Matrix<T_0, Dynamic, 1> >& pMatrix,
     	 const vector<T_1>& time,
     	 const vector<T_2>& amt, 
     	 const vector<T_3>& rate,
    	 const vector<T_4>& ii,
  	   	 const vector<int>& evid,
    	 const vector<int>& cmt,
         const vector<int>& addl,
         const vector<int>& ss,
     	 PKModel model,
     	 const F& f,
     	 const std::vector<Matrix<T_5, Dynamic, Dynamic> >& system);
};
		
#endif
