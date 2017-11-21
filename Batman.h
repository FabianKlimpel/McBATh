// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BATMAN__H
#define __BAT__BATMAN__H

#include <BAT/BCModel.h>
#include <BAT/BCMath.h>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLine.h>
#include <omp.h>
#include <TMath.h>
#include "SDataHandler.h"
#include "RefHandler.h"
#include "WeightsHandler.h"
#include "Power.h"


// This is a Batman header file.
// Model source code is located in file Batman/Batman.cxx

// ---------------------------------------------------------
class Batman : public BCModel
{

public:

   	// Constructor
	Batman(const std::string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh, ConfigHandler& ch);
	Batman(const std::string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh);

   	// Destructor
    virtual ~Batman();

   	// Overloaded LogLikelihood to implement model
    virtual double LogLikelihood(const std::vector<double>& pars) = 0;

	//Tracker & Plotter
	void updateHistory();
	void updateHistory_prerun();
	virtual void binlogprobdistribution() = 0;
	virtual void chi2plot() = 0;
	void chainhistoryplot();
	void chainhistoryplot_prerun();
	virtual void chi2distribution() = 0;

//protected variables and a function for derivation
protected:

	//storage for the cumulative number of bins after a number of histograms
	std::vector<size_t> num_bins_cumulative;
	
	//minimum/maximum allowed parameter value
	std::vector<double> var_minvalue, var_maxvalue;
	
	//fit/refernce parameter values/errors
	std::vector<std::vector<double>> fit_params, fit_errors, ref_values, ref_values_errors;
	
	//parameter names
	std::vector<std::string> var_names;
	
	//weights per histogram
	std::vector<int> weights;
	
	//number of iterations in the pre-/main-run
	int iteration_number;
	
	//object for the creation of the powers of the parameters in the respective monomial
	Power _powlist;
	
	//list of powers
	std::vector<std::vector<int>> _power;
	
	//function for delivering the list of powers for @_power
	void windofchange(int max_power);

	//histories of the evolution of the parameters in every chain in the pre-/main-run 
	std::vector<std::vector<std::vector<double>>> MCMC_history, MCMC_history_prerun;
	
	//flag for plotting chain histories during the pre run
	int _prerunplots;
	
	//number of steps in the pre run between the history plots
	size_t _numprerunplotsteps;

	std::vector<std::vector<int>> weights_mat;

	//std::vector<std::vector<double>> history_storage; 

};
// ---------------------------------------------------------

#endif
