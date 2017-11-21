#ifndef __REFHANDLER__H
#define __REFHANDLER__H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <limits>
#include "ConfigHandler.h"
#include <omp.h>
#include "SDataHandler.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TAxis.h>

/**
* This class handles the simulated reference data
*/
class RefHandler
{

public:

	RefHandler(ConfigHandler& ch, SDataHandler& sdh);

	vector<vector<double>> getref_values();
	vector<vector<double>> getref_values_errors();
	vector<size_t> getnum_bins();


private:
	
	/**
	* @_parameter_values: contains the sampled parameter values
	* @_bin_values: contains the simulated reference data
	* @_bin_values_err: contains the simulated reference data error
	* @_normalvectors: contains the normalvectors of the anchor points
	* @_num_bins: contains the number of bins in each observable
	*/
	vector<vector<double>> ref_values, ref_values_errors;
	vector<size_t> _num_bins;

	//reader for the parameter values and the simulated reference data (error)
	void readrefdata(ConfigHandler& ch, SDataHandler& sdh);

};

#endif

