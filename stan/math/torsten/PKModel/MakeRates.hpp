// CHECK - this function doesn't seem to get used.

#ifndef STAN_MATH_TORSTEN_PKMODEL_MAKERATES_HPP
#define STAN_MATH_TORSTEN_PKMODEL_MAKERATES_HPP

#include <iostream>
#include <stan/math/torsten/PKModel/Event.hpp>
#include <stan/math/torsten/PKModel/Rate.hpp>
#include <stan/math/torsten/PKModel/SearchReal.hpp>

const int nCmtMax = 50; // FIX ME - find a way not to specify a limit on nCmt
bool find_time(dvector, double); //forward declare

/**
 * Make infusion rate history and augment event history with end-of-infusion times.
 * At the end of the procedure, both events and rates are sorted.  
 *
 * @param[in] events 
 * @param[in] rates
 * @return - modified events and rates. 
 */
void MakeRates(EventHistory& events, RateHistory& rates) {		
  using std::vector;
  
  if(!events.Check()) events.Sort(); 
  
  vector<double> rate_init(nCmtMax, 0);
  Rate newRate(0, rate_init);
  int nEvent = events.get_size();					    
  for(int j = 0; j < nEvent; j++)
    if((j == 0) || (events.get_time(j) != events.get_time(j - 1))) {
      newRate.time = events.get_time(j);
      rates.InsertRate(newRate); 
    }
	
  rates.RemoveRate(0); //remove rate created by default constructor. 
	
  if(!rates.Check()) rates.Sort();	
	
  //Create time vector for rates 
  vector<double> RateTimes(rates.get_size(), 0);
  for(int j = 0; j < rates.get_size(); j++) RateTimes[j] = rates.Rates[j].time; 

  //Create time vector for events.
  vector<double> EventTimes(nEvent, 0); 
  for(j = 0; j < events.get_size(); j++) EventTimes[j] = events.Events[j].time;

  int i = 0;
  double endTime;
  Event newEvent;
  while(i < nEvent) {	
    if(((events.Events[i].evid == 1) || (events.Events[i].evid == 4))
      &&((events.Events[i].rate>0)&&(events.Events[i].amt>0))) {
      endTime = events.Events[i].time + events.Events[i].amt / events.Events[i].rate;
	  newEvent = newEvent(endTime, 0, 0, 0, 2, events.Events[i].cmt, 0, 0, false, true);
      events.InsertEvent(newEvent);
      if (!events.Check()) events.Sort();
      EventTimes.push_back(endTime);
      std::sort(EventTimes.begin(), EventTimes.end());

      // Only create a new Rate if endTime does not correspond to a time that is 
      // already in RateHistory. 
      if (!find_time(RateTimes, endTime)) {
        newRate.time = endTime;
        rates.InsertRate(newRate);
        if (!rates.Check()) rates.Sort();
        RateTimes.push_back(endTime);
        std::sort(RateTimes.begin(), RateTimes.end());
      }			
			
      //Find indexes at which time of event and endtime occur. 
      int l = SearchReal(RateTimes, nEvent, events.Events[i].time),
        k = SearchReal(RateTimes, nEvent, endTime);

      //Compute Rates for each element between the two times
      for (int iRate = l - 1; iRate < k; iRate++)
        rates.Rates[iRate].rate[events.Events[i].cmt-1] += events.Events[i].rate;
    }
    i++;
  }
	
	//Sort events and rates
	if(!rates.Check()) rates.Sort();
	if(!events.Check()) events.Sort();
	
}

bool find_time(dvector v, double time)
{
	int k, size;
	bool found=false;
	size = v.size();
	k = SearchReal(v, size, time);
	if((k <= size) && (time == v[k-1])) found = true; 
	return found;
}
		
		
#endif