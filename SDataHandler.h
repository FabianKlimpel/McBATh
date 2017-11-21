#ifndef __SDATAHANDLER__H
#define __SDATAHANDLER__H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>

using namespace std;

/**
* class for reading the interpolation file
* the reading process is divided in 4 void functions
* the public parameters are for the external access by the BCModel constructor
*/
class SDataHandler
{
public:
	SDataHandler(string interpolation, int sameorder);
	
	size_t get_num_histos();

	size_t getnum_runs();
	vector<string> get_runs();
	size_t get_dimension();
	vector<string> getvar_names();
	vector<string> get_analysis();
	vector<string> get_observable();	
	vector<size_t> getchanging_histo();

	vector<double> getvar_minvalue();
	vector<double> getvar_maxvalue();

	vector<size_t> getfit_num_params();
	vector<vector<double>> getfit_params();
	vector<vector<double>> getfit_error();

	 
	
	

	
private:

	void parameter_reading();
	void var_min_max_reading();
	void fit_params_reading(int sameorder);
	void analysis_name();
	void set_runs();

	vector<size_t> fit_num_params;
	vector<vector<double>> fit_params, fit_error;
	size_t num_runs = 0, dimension = 0;
	vector<string> runs, var_names, analysis, observable;
	vector<size_t> changing_histo;
	vector<double> var_minvalue, var_maxvalue;

	size_t pos_begin, pos_end, length;
	string line, parameternames;
};

#endif
