#ifndef __BMPROF__H
#define __BMPROF__H

#include "Batman.h"

using namespace std;

class BMprof: public Batman
{
	
public:

	//Constructor
	BMprof(const string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh, int fit_order);
	
	//Destructor
	~BMprof();
	
	//LogLikelihood calculator
	double LogLikelihood(const vector<double>& pars);
	
	//Plotter
	void binlogprobdistribution();
	void chi2plot();
	void chi2distribution();
	
private:

	//Calculators of the fit value / variance at a certain point in the parameter space
	double fit_value(const vector<double>& parameters, const int& index);
	double fit_var(const vector<double>& anchors, const int& index);
	
	//number of fit parameters per polynomial function
	size_t _fit_num_params;

};

#endif
