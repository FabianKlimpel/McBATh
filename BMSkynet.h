#ifndef __BMSKYNET__H
#define __BMSKYNET__H

#include "Batman.h"
#include "CovMatHandler.h"

using namespace std;

class BMSkynet: public Batman
{
	
public:

	//Constructor
	BMSkynet(const string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh, ConfigHandler& ch);
	
	//Destructor
	~BMSkynet();
	
	//LogLikelihood calculator
	double LogLikelihood(const vector<double>& pars);
	
	//plotter
	void binlogprobdistribution();
	void chi2plot();
	void chi2distribution();
	
	//setter for the # of covariance matrices = # of chains
	void resize_covmats(size_t nchains);
	
	
private:

	//calculators for the interpolation value / variance
	double fit_value(const vector<double>& anchors, const int& index);
	double fit_var(const vector<vector<double>>& tmp_specific_mat, size_t index);
	
	//Container for the parameter depending matrices, that are used for faster variance computation. The size of the matrix is = max. # of fit parameters in use
	//structure: _specific_mat[chain number][row][coloumn]
	vector<vector<vector<double>>> _specific_mat;
	
	//storage for the function values / variances
	//structure: fval/fvar[chain number][bin number]
	vector<vector<double>> fval, fvar;
	
	//storage for the # of fit parameters in every bin
	vector<size_t> fit_num_params;
	
	//Handler of the covariance matrices
	CovMatHandler _cmh;
	size_t rejected_bins;
};

#endif
