// version 0.8

#ifndef PKMODEL_SEARCHREAL_HPP
#define PKMODEL_SEARCHREAL_HPP

#include <iostream>
using std::vector;
using stan::math::var;

//DEV: may need to overload functions to work on fvar.

int min(int a, int b); //forward declare 

/**
 * Returns the position of the largest value smaller or greater than srchNum
 * Assumes that v is sorted.
 * The numltm is used to set an upper limit for the index
 * the function can return. If numltm is greater than the size of v, 
 * there is no upper limit and the function searches the entire vector.
 *
 * @tparam: type of scalar in input vector 
 * @param[in]: v searched vector
 * @param[in]: numltm maximum index
 * @param[srchNum]: searched Number 
 * @return: index of largest value <= srchNum
 *
 */

template<typename T>
inline int SearchReal(std::vector<T> v, int numltm, T srchNum)
{

	int first=0, last, mid, real_limit;
	
	assert(numltm >= 0);
	real_limit = min(numltm, v.size());	//limit cannot exceed size of searched vector
	last = real_limit - 1;
	
	if(srchNum < v[first]) mid = first;
	else if(srchNum >= v[last]) mid = last+1;
	else 
	{
		while(first <= last)
		{
			mid = (first+last)/2; 
			if (srchNum < v[mid]) last = mid-1;
			else if (srchNum > v[mid]) first = mid+1;
			else first = last + 1;
		}
	
		while(srchNum >= v[mid]){mid += 1;}
	}
	
	return mid;
}


inline int min(int a, int b)
{
	if(a < b) return a;
	else return b;
}


#endif