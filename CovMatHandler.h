#ifndef __COVMATHANDLER__H
#define __COVMATHANDLER__H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <omp.h>

using namespace std;

/**
 * This class handles the covariance matrices
 * More specific: it loads, stores them and calculates errors out of them
 */
class CovMatHandler
{
public:

	//Constructor
	CovMatHandler();
	
	//Setter for the number of bins
	void set_num_points(size_t num_points);

	//reader of the covariance matrices
	void read_covmat(string filename, size_t index);
	
	//setter for @templatemat
	void set_templatemat(vector<vector<int>>& power);
	
	//calculates a parameter depending matrix based upon @templatemat and specific parameters
	void get_specific_mat(vector<double>& parameters, vector<vector<double>>& spec_mat);
	vector<vector<double>> get_specific_mat(vector<double>& parameters);
	
	//calculates the variance
	double get_var(const vector<vector<double>>& mat, size_t index);
	
	//setter for a diagonal matrix (aka use only the variances of the parameters)
	void set_diagmat(vector<double> vec, size_t index);
	
private:

	//adds two vectors
	vector<int> add_powers(const vector<int>& a, const vector<int>& b);
	
	//getter for an entry of a specific matrix
	inline double get_specific_element(size_t row, size_t col, const vector<vector<double>>& powers, const size_t npars);
	
	//getter for the highest occuring power in every parameter
	void set_highest_powers();
	
	//getter of the value of a specific parameter set to a power up to the highest occuring power
	vector<vector<double>> get_powers(vector<double>& parameters);

	//checks, if a matrix is symmetric
	void is_symm(size_t index);
	
	//transforms the covariance matrix so, that only a triangle matrix needs to be used
	void make_triangle(size_t index);

	//covariance matrix; structure: covmat[bin number][row][coloumn]                                                                                                                                 
	vector<vector<vector<double>>> covmat;
	
	//template matrix; structure: templatemat[row][coloumn][parameter number]
	vector<vector<vector<int>>> templatemat;
	
	//vector that contains the highest occuring power of every parameter
	vector<int> highest_powers;
	
	//flag that is true, if the matrix is a diagonal matrix
	int is_diag = 0;

};

#endif
