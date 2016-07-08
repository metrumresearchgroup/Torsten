// version 0.8
// DEV - replace overload of function with template function

#ifndef PKMODEL_EXTRACTVECTOR_HPP
#define PKMODEL_EXTRACTVECTOR_HPP

#include <iostream>
#include <string>
#include <stan/model/model_header.hpp>
#include <Eigen/Dense>

using std::vector;
using std::string;
using namespace Eigen;
using stan::math::var;

/**
 * Extracts the nth column or row from a matrix and
 * returns it as an std::vector. Indexing starts at 0, meaning the 
 * first row or column is found for n=0. 
 *
 * @param[in]: matrix from which the row or column will be extracted
 * @param[in]: n the row or column number
 * @param[in]: str a string that specifies whether the function extracts
 *             a column or a row. 
 * @return: output dvector
 *
 */

vector<double> ExtractVector(Matrix<double, Dynamic, 1> matrix, int n, string str)
{
	int i, length;
	vector<double> output; 
	
	assert((str == "row")||(str == "col"));
	
	if(str == "col") //extract a column. 
	{ 
		length = matrix.rows();
		assert(n < matrix.cols());
		output.resize(length);
		for(i=0;i<length;i++){output[i] = matrix(i, n);}
	}
	
	else // extract a row
	{
		length = matrix.cols();
		assert(n < matrix.rows());
		output.resize(length);
		for(i=0;i<length;i++){output[i] = matrix(n, i);}
	}
	
	return output;
}

vector<var> ExtractVector(Matrix<var, Dynamic, 1> matrix, int n, string str)
{
	int i, length;
	vector<var> output; 
	
	assert((str == "row")||(str == "col"));
	
	if(str == "col") // extract a column. 
	{ 
		length = matrix.rows();
		assert(n < matrix.cols());
		output.resize(length);
		for(i=0; i<length; i++) {output[i] = matrix(i, n);}
	}
	
	else // extract a row
	{
		length = matrix.cols();
		assert(n < matrix.rows());
		output.resize(length);
		for(i=0; i<length; i++) {output[i] = matrix(n, i);}
	}
	
	return output;
}



#endif