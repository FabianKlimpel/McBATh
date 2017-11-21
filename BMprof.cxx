#include "BMprof.h"

/**
 * Constructor
 * @model_name: name of the model
 * @sdh: container for fit parameter information
 * @rh: container of reference information
 * @wh: container of weight information
 * @fit_order: order of the polynomial
 */
BMprof::BMprof(const string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh, int fit_order) : Batman(model_name, sdh, rh, wh)
{
		
	//storing the number of fit parameters
	_fit_num_params = sdh.getfit_num_params()[0];
			
	//setting up the list of powers
	_powlist.setn(var_names.size());
	windofchange(fit_order);
	
}

// ---------------------------------------------------------
BMprof::~BMprof()
{
    // destructor
}

/**
 * This function is the heartpiece of the Bayesian Analysis. It calculates Loglikelihood of the parameters
 * @parameters: this is a vector with a length = #parameters in the model, containing a value for each parameter
 * @logprob: this parameter contains (in the end) the Loglikelihood of the model given @parameters
 * @tmp_parameters: [0,1]-mapped values of @parameters based upon the limits of the parameter values
 * @anchors: value of the @parameters - depending part of every monomial
 */ 
double BMprof::LogLikelihood(const vector<double>& parameters)
{
	double logprob = 0;
	
	//If the program is in its mainrun, the history of the evolution of the Markovchains will be saved.
	if(fMCMCPhase == BCEngineMCMC::kMainRun)
		updateHistory();
	
	//If the program is in its prerun, the history of the evolution of the Markovchains will be stored and printed, if wanted	.
	if(fMCMCPhase == BCEngineMCMC::kPreRun && _prerunplots)
		updateHistory_prerun();

	//setting up @tmp_parameters and calculate its value
	vector<double> tmp_parameters;
	tmp_parameters.resize(parameters.size());
	
	for(size_t i = 0; i < parameters.size(); i++)
		tmp_parameters[i] = (parameters[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);
	
	//calculate @anchors based on the information of @tmp_parameters and @_power
	std::vector<double> anchors;
	anchors.assign(_power.size(), 1);
	
	for(size_t i = 0; i < _power.size(); i++)
		for(size_t j = 0; j < _power[0].size(); j++)
			anchors[i] *= pow(tmp_parameters[j], _power[i][j]);
					
	#pragma omp parallel for reduction(+:logprob) schedule(guided)
	for(size_t i = 0; i < ref_values.size(); i++)
		for(size_t j = 0; j < ref_values[i].size(); j++)
		{
			//calculating the Loglikelihood by walking over every bin and getting the LogGauss of the interpolation and the reference values of the current bin.
			//the weights improve/reduce the importance of a certain bin
			logprob += weights[i] * BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], 
						sqrt(ref_values_errors[i][j] * ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j)));
		
		}
	
	//Return the log of the conditional probability p(data|pars).
    return logprob;
}

/**
 * This function calculates the function value of the interpolation function in a certain bin
 * @anchors: fit parameter independent part of every monimial
 * @index: # of the bin
 * @result: function value the will be returned
 */
double BMprof::fit_value(const vector<double>& anchors, const int& index){
		
	double result = 0;
	
	//walk over every monomial and multadd every part of the polynomial together
	for(size_t i = 0; i < _fit_num_params; i++)
		result += fit_params[index][i] * anchors[i];
	
	
	return result;
}

/**
 * This function calculates the function variance of the interpolation function in a certain bin
 * @anchors: fit parameter independent part of every monimial
 * @index: # of the bin
 * @result: function variance that will be returned
 */
double BMprof::fit_var(const vector<double>& anchors, const int& index){
	
	double result = 0;
	
	//walk over every monomial and multadd every part of the polynomial together
	for(size_t i = 0; i < _fit_num_params; i++)
		result += fit_errors[index][i] * anchors[i];
		
	//square the error and return the variance
	return result * result;
	
}


/**
 * This function plots the Loglikelihood of every single bin using the best fit parameters estimated by BAT
 * @mode: parameter for the mode of the parameters
 * @tmp_parameters: [0,1]-mapped parameter values based on the limits of allowed parameter values
 * @anchors: value of the @mode - depending part of every monomial
 * @maximum, @minimum: parameters for the maximum/minimum value of the binwise Loglikelihood, so that the plotranges can be set
 * @tc, @th1d: parameter for a storageframe and the histogram
**/ 
void BMprof::binlogprobdistribution(){

	//calculating the mode of the parameters using Minuit
	std::vector<double> mode = FindMode(GetBestFitParameters());
	
	//setting up @tmp_parameters and calculate its value	
	vector<double> tmp_parameters;
	tmp_parameters.resize(mode.size());
	
	for(size_t i = 0; i < mode.size(); i++)
		tmp_parameters[i] = (mode[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);

	//calculate @anchors based on the information of @tmp_parameters and @_power	
	std::vector<double> anchors;
	anchors.assign(_power.size(), 1);
	
	for(size_t i = 0; i < _power.size(); i++)
		for(size_t j = 0; j < _power[0].size(); j++)
			anchors[i] *= pow(tmp_parameters[j], _power[i][j]);
								
	//walk over all different bins, calculate the Loglikelihood of the bin and searching for its maximum/minimum, so that the plotrange can be set later
	double maximum = 0, minimum = 0;
	
	for(size_t i = 0; i < ref_values.size(); i++)
		for(size_t j = 0; j < ref_values[i].size(); j++)
		{
			if(BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j))) > maximum)
				maximum = BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j)));
			if(BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j))) < minimum)
				minimum = BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j)));
		}

	//setting up the graph and storage container
	TCanvas* tc;
	TH1D* th1d;
	
	tc = new TCanvas("binlogprobdistribution.png");
	tc->cd();
	th1d = new TH1D("binlogprogdistribution", "binlogprobdistribution", (int) (1.5 * sqrt(fit_params.size())), minimum, maximum);

	th1d->GetXaxis()->SetTitle("Log(p)");
	th1d->GetYaxis()->SetTitle("#");
	
	//walk over every bin and saving the calculated Loglikelihood values in @th1d
	for(size_t i = 0; i < ref_values.size(); i++)
			for(size_t j = 0; j < ref_values[i].size(); j++)
				th1d->Fill(BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j))));
		

	//draw and save
	th1d->Draw();
	
	gPad->Update();
	tc->Print("binlogprobdistribution.png");
	tc->Clear();
}


/**
 * This function plots the Chi^2 value of each bin using the best fit parameters estimated by BAT and compares them to the referencedate of the respective bin
 * @mode: parameter for the mode of the parameters
 * @tmp_parameters: [0,1]-mapped parameter values based on the limits of allowed parameter values
 * @anchors: value of the @mode - depending part of every monomial
 * @maximum: this is a calculation of the upper plotrange limit
 * @tc, @th1d: parameter for a storageframe and the histogram
**/ 
void BMprof::chi2plot(){

	//calculating the mode of the parameters using Minuit
	std::vector<double> mode = FindMode(GetBestFitParameters());		

	//setting up @tmp_parameters and calculate its value		
	vector<double> tmp_parameters;
	tmp_parameters.resize(mode.size());
	
	for(size_t i = 0; i < mode.size(); i++)
		tmp_parameters[i] = (mode[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);

	//calculate @anchors based on the information of @tmp_parameters and @_power	
	std::vector<double> anchors;
	anchors.assign(_power.size(), 1);
	
	for(size_t i = 0; i < _power.size(); i++)
		for(size_t j = 0; j < _power[0].size(); j++)
			anchors[i] *= pow(tmp_parameters[j], _power[i][j]);
				
	//walk over all different bins, calculate the Loglikelihood of the bin and searching for its maximum, so that the plotrange can be set later
	double maximum = 0;
	
	for(size_t i = 0; i < ref_values.size(); i++)
		for(size_t j = 0; j < ref_values[i].size(); j++)
			if((fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j)) > maximum)
				maximum = (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j));
	

	//setting up the graph and storage container
	TCanvas* tc;
	TH1D* th1d;
	
	tc = new TCanvas("chi2plot.png");
	tc->cd();
	th1d = new TH1D("chi2plot", "chi2plot", (int) (1.5 * sqrt(fit_params.size())), 0, maximum);

	th1d->GetXaxis()->SetTitle("binwise Chi^2");
	th1d->GetYaxis()->SetTitle("#");
	
	//walk over every bin and saving the calculated Chi^2 values in @th1d
	for(size_t i = 0; i < ref_values.size(); i++)
		for(size_t j = 0; j < ref_values[i].size(); j++)
			th1d->Fill((fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j)));

	//draw and save
	th1d->Draw();
	
	gPad->Update();
	tc->Print("chi2plot.png");
	tc->Clear();
}


/**
 * This function calculates the Chi^2 over all bins using the best fit parameters calculated from the BAT run and compares them with the referencedata
 * @tc, @tf, @tl: parameters for a storageframe, a plot of the pdf of Chi^2 and the actual resulting Chi^2 value, respectively.
 * @sum: parameter for storing the binwise calculated Chi^2 value
 * @mode: parameter for storing the mode for the parameters calculated by BAT
 * @tmp_parameters: [0,1]-mapped parameter values based on the limits of allowed parameter values
 * @anchors: value of the @mode - depending part of every monomial
 * @ndof: parameter for calculating the number of dof
**/ 
void BMprof::chi2distribution(){

	TCanvas* tc = new TCanvas("chi2distribution.png");
	tc->cd();

	double sum = 0;
	//calculating the mode of the parameters by using Minuit
	std::vector<double> mode = FindMode(GetBestFitParameters());

	//setting up @tmp_parameters and calculate its value	
	vector<double> tmp_parameters;
	tmp_parameters.resize(mode.size());
	
	for(size_t i = 0; i < mode.size(); i++)
		tmp_parameters[i] = (mode[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);

	//calculate @anchors based on the information of @tmp_parameters and @_power	
	std::vector<double> anchors;
	anchors.assign(_power.size(), 1);
	
	for(size_t i = 0; i < _power.size(); i++)
		for(size_t j = 0; j < _power[0].size(); j++)
			anchors[i] *= pow(tmp_parameters[j], _power[i][j]);
					
	//loop over all bins
	for(size_t i = 0; i < ref_values.size(); i++)
		for(size_t j = 0; j < ref_values[i].size(); j++)
		{
			//simple Chi^2 for every bin
			sum += (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + fit_var(anchors, num_bins_cumulative[i] + j));
		}

		

	//getting the ndof
	int ndof = fit_params.size() - GetNParameters();
	
	//setting up the graphs and plotting them
	TF1* tf1 = new TF1("chi2","ROOT::Math::chisquared_pdf(x,[0])", 0, 2 * sum);
	
	tf1->SetTitle("Chi^2-distribution;Chi^2;p");
	
	tf1->SetParameter(0, ndof);
	tf1->Draw();

	TLine* tl = new TLine(sum, 0, sum, 1);
	tl->Draw("same"); 

	gPad->Update();
	tc->Print("chi2distribution.png");
	tc->Clear();
	
}
