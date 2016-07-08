// version 0.8

#ifndef PKMODEL_MAKERATES_HPP
#define PKMODEL_MAKERATES_HPP

#include <iostream>
#include "Event.hpp"
#include "Rate.hpp"
#include "SearchReal.hpp"

using std::vector;
typedef vector<double> dvector;


/**
 * Make infusion rate history and augment event history with end-of-infusion times.
 * At the end of the procedure, both events and rates are sorted.  
 *
 * @param[in] events 
 * @param[in] rates
 * @return - modified events and rates. 
 */

const int nCmtMax = 50; // FIX ME - find a way not to specify a limit on nCmt

bool find_time(dvector, double); //forward declare
	
void MakeRates(EventHistory& events, RateHistory& rates)
{	
 	int i=0, j, iRate=0, k, l, nEvent = events.get_size(); 
	double endTime, time_init=0;
	Event newEvent;
	dvector EventTimes(nEvent, 0), rate_init(nCmtMax,0);
	Rate newRate(time_init,rate_init);
	
	if(!events.Check()) events.Sort(); // events need to be sorted, else the function will
									                   // not compute one rate per time, but most likely
									                   // one rate per event. 
							    
	for(j=0;j<nEvent;j++)
	{
		if((j==0)||(events.get_time(j) != events.get_time(j-1)))
		{
			newRate.time = events.get_time(j);
			rates.InsertRate(newRate); 
		}
	}
	
	rates.RemoveRate(0); //remove rate created by default constructor. 
	
	if(!rates.Check()) rates.Sort();	
	
	//Create time vector for rates 
	dvector RateTimes(rates.get_size(),0);
	for(j=0;j<rates.get_size();j++){RateTimes[j] = rates.Rates[j].time;} 

	//Create time vector for events. 
	for(j=0;j<events.get_size();j++){EventTimes[j]=events.Events[j].time;}
	
	while(i < nEvent)
	{	
		if(((events.Events[i].evid==1)||(events.Events[i].evid==4))
			&&((events.Events[i].rate>0)&&(events.Events[i].amt>0)))
		{
			endTime = events.Events[i].time + events.Events[i].amt/events.Events[i].rate;
			newEvent = newEvent(endTime, 0, 0, 0, 2, events.Events[i].cmt, 0, 0, false, true);
			events.InsertEvent(newEvent);
			if(!events.Check()) events.Sort();
			EventTimes.push_back(endTime);
			std::sort (EventTimes.begin(), EventTimes.end());
						
			//Only create a new Rate if endTime does not correspond to a time that is 
			//already in RateHistory. 
			if(!find_time(RateTimes, endTime))
			{
				newRate.time = endTime;
				rates.InsertRate(newRate);
				if(!rates.Check()) rates.Sort();
				RateTimes.push_back(endTime);
				std::sort (RateTimes.begin(), RateTimes.end());
			}			
			
			//Find indexes at which time of event and endtime occur. 
			l=SearchReal(RateTimes, nEvent, events.Events[i].time);
			k=SearchReal(RateTimes, nEvent, endTime);

			//Compute Rates for each element between the two times
			for(iRate=l-1; iRate < k; iRate++)
			{
				rates.Rates[iRate].rate[events.Events[i].cmt-1] += events.Events[i].rate;
			}
		}
		i++;
		iRate = 0;
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