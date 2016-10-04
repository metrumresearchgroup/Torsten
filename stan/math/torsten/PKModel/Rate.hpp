#ifndef STAN_MATH_TORSTEN_PKMODEL_RATE_HPP
#define STAN_MATH_TORSTEN_PKMODEL_RATE_HPP

#include <Eigen/Dense>
#include <stan/math/torsten/PKModel/functions.hpp>

// FIX ME: using statements should occur within
// the scope of functions
using std::vector;
using namespace Eigen;
using stan::math::var;
using boost::math::tools::promote_args;

//forward declaration
template<typename T_time, typename T_rate> class RateHistory; 

/**
 * The Rate class defines objects that contain the rate in each compartment
 * at each time of the event schedule (but not nescessarily at each event) 
 */
template<typename T_time, typename T_rate>
class Rate {
private:
	T_time time;
	vector<T_rate> rate; // rate for each compartment
public:
	Rate() {
		vector<T_rate> v(1,0);
		time = 0;
		rate = v;
	}

	Rate(T_time p_time, vector<T_rate> p_rate) {  //parameters: time and rate vectors
		time = p_time;
		rate = p_rate;
	}
	
	void Print() {
		int i;
		print(time, false);
		for(i=0;i<rate.size();i++){print(rate[i], false);}
		std::cout<<std::endl;
	}		
	
	friend class RateHistory<T_time, T_rate>;
	template <typename T_amt, typename T_ii> 
	friend void MakeRates(EventHistory<T_time, T_amt, T_rate, T_ii>&, 
	  RateHistory<T_time,T_rate>&);
	
	template <typename T0, typename T1, typename T2, typename T3,
	typename T4, typename F, typename T5> 
	friend
	Matrix<typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic, Dynamic> 
	Pred(const vector< Matrix<T0, Dynamic, 1> >& pMatrix,
     	 const vector<T1>& time,
     	 const vector<T2>& amt, 
     	 const vector<T3>& rate,
    	 const vector<T4>& ii,
  	     const vector<int>& evid,
    	 const vector<int>& cmt,
         const vector<int>& addl,
         const vector<int>& ss,
     	 PKModel model,
     	 const F& f,
     	 const Matrix<T5, Dynamic, Dynamic>& system);	
};

/**
 * The RateHistory class defines objects that contain a vector of rates,
 * along with a series of functions that operate on them.  
 */
template <typename T_time, typename T_rate>
class RateHistory {
private:
	vector<Rate<T_time, T_rate> > Rates;
	
public:	
	RateHistory() {
		Rate<T_time, T_rate> initRate;
		Rates.resize(1);
		Rates[0]=initRate;
	}

	template <typename T0, typename T1>
	RateHistory(vector<T0> p_time, vector<vector<T1> > p_rate) {
		int i, nRate = p_rate.size();
	
		Rates.resize(nRate);
		for(i=0;i<nRate;i++) {Rates[i] = Rate<T_time, T_rate>(p_time[i],p_rate[i]);}
	}
	
	RateHistory(int nEvent) {
		int i;
		Rate<T_time, T_rate> initRate;
		Rates.resize(nEvent);
		for(i=0;i<nEvent;i++) {Rates[i] = initRate; }
	}
	
	bool Check() {
		int i=Rates.size()-1;
		bool ordered=true;
		
		while((i>0)&&(ordered)) {
			ordered=(Rates[i].time >= Rates[i-1].time);
			i--;
		}
		
		return ordered;
	}
	
	Rate<T_time,T_rate> GetRate(int i) {
		Rate<T_time,T_rate> newRate(Rates[i].time, Rates[i].rate);
		return newRate;
	}
	
	void InsertRate(Rate<T_time,T_rate> p_Rate){Rates.push_back(p_Rate);}

	void RemoveRate(int i) {	
		assert(i >= 0);
		Rates.erase(Rates.begin()+i);
	}
	
	int Size(){return Rates.size();}
	
	void Print(int j) {
		int i, nCmt = Rates[j].rate.size();
		print(Rates[j].time, false);
		for(i=0;i<nCmt;i++) print(Rates[j].rate[i], false);
		std::cout << std::endl;
	}

	struct by_time {
		bool operator()(Rate<T_time,T_rate> const &a, Rate<T_time,T_rate> const &b) {
			return a.time < b.time;
		}
	};


	void Sort(){std::sort(Rates.begin(),Rates.end(),by_time());}

	template <typename T_amt, typename T_ii>
	void MakeRates(EventHistory<T_time, T_amt, T_rate, T_ii>& events, int nCmt) {	
 	  int i=0, j, iRate=0, k, l; 
	  T_time endTime, time_init = 0;
	  Event<T_time, T_amt, T_rate, T_ii> newEvent;
	  vector<T_time> EventTimes(events.get_size(), 0);
	  vector<T_rate> rate_init(nCmt, 0);
	  Rate<T_time, T_rate> newRate(time_init, rate_init);

	  if(!events.Check()) events.Sort();

	  for(j = 0; j < events.get_size(); j++) {
		  if(j == 0 || events.get_time(j) != events.get_time(j-1)) {
			  newRate.time = events.get_time(j);
			  InsertRate(newRate); 
	  	}
	  }

	  RemoveRate(0); //remove rate created by default constructor. 

	  if(!Check()) Sort();	

	  //Create time vector for rates 
	  vector<T_time> RateTimes(Size(), 0);
	  for(j = 0; j < Size(); j++) RateTimes[j] = Rates[j].time; 

	  //Create time vector for events 
	  for(j = 0; j < events.get_size(); j++) EventTimes[j] =
	    events.Events[j].time;
	  
	  while(i < events.get_size()) {	
		  if((events.Events[i].evid == 1 || events.Events[i].evid == 4)
			  && (events.Events[i].rate > 0 && events.Events[i].amt > 0)) {
			  
			  endTime = events.Events[i].time +
			    events.Events[i].amt/events.Events[i].rate;
			  newEvent = newEvent(endTime, 0, 0, 0, 2, events.Events[i].cmt,
			    0, 0, false, true);
			  events.InsertEvent(newEvent);
			  if(!events.Check()) events.Sort();
			  EventTimes.push_back(endTime);
			  std::sort (EventTimes.begin(), EventTimes.end());

			  //Only create a new Rate if endTime does not correspond to a time that is 
			  //already in RateHistory. - CHECK
			  if(!find_time(RateTimes, endTime)) {
				  newRate.time = endTime;
				  InsertRate(newRate);
				  if(!Check()) Sort();
				  RateTimes.push_back(endTime);
				  std::sort (RateTimes.begin(), RateTimes.end());
			  }			

			  //Find indexes at which time of event and endtime occur. 
			  l = SearchReal(RateTimes, events.get_size(), events.Events[i].time);
			  k = SearchReal(RateTimes, events.get_size(), endTime);

			  //Compute Rates for each element between the two times
			  //for(iRate = l - 1; iRate < k; iRate++)
			  for(iRate = l ; iRate < k; iRate++)
  				Rates[iRate].rate[events.Events[i].cmt-1] += events.Events[i].rate;
  		}
  		i++;
  		iRate = 0;
  	}
	
  	//Sort events and rates
  	if(!Check()) Sort();
  	if(!events.Check()) events.Sort();
  	
  }
	
	// declare friends 
	friend class Rate<T_time, T_rate>;
	template <typename T_amt, typename T_ii>
	friend void MakeRates(EventHistory<T_time, T_amt, T_rate, T_ii>&, 
													RateHistory<T_time, T_amt>&);
	
	template <typename T0, typename T1, typename T2, typename T3, typename T4,
	  typename F, typename T5> 
	friend
	Matrix<typename promote_args<T0, T1, T2, T3, T4>::type, Dynamic, Dynamic> 
	Pred(const vector< Matrix<T0, Dynamic, 1> >& pMatrix,
     	 const vector<T1>& time,
     	 const vector<T2>& amt, 
     	 const vector<T3>& rate,
    	 const vector<T4>& ii,
  	   	 const vector<int>& evid,
    	 const vector<int>& cmt,
         const vector<int>& addl,
         const vector<int>& ss,
     	 PKModel model,
     	 const F& f,
     	 const Matrix<T5, Dynamic, Dynamic>& system);	
};

#endif
