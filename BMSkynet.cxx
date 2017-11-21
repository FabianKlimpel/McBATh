#include "BMSkynet.h"

/**
 * Constructor
 * @model_name: name of the model
 * @sdh: Container of fit information
 * @rh: Container of reference information
 * @wh: Container of weights
 * @ch: Container of configurations
 */
BMSkynet::BMSkynet(const string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh, ConfigHandler& ch) : Batman(model_name, sdh, rh, wh, ch)
{

	//store the number of fit parameters per bin
	fit_num_params = sdh.getfit_num_params();
		
	//setting up the number of parameters for the monomial wise powers
	_powlist.setn(var_names.size());
		
	size_t max = 0, max_power = 0;
		
	//search for the bin with the highest number of fit parameters
	for(size_t i = 0; i < fit_num_params.size(); i++)
		if(fit_num_params[i] > max)
			max = fit_num_params[i];
	
	//get the list of powers for as long as the list is smaller than the highest number of fit parameters			
	while(_power.size() < max)
	{
		windofchange(max_power);
		max_power++;
	}
		
	//check, if the whole covariance matrices should be used or only the diagonal terms
	if(ch._covmat)
	{
		cout << "loading covariance matrices..." << flush;
		
		//setting the number of covariance matrices = number of bins
		_cmh.set_num_points(fit_num_params.size());
		
		//read every covariance matrix
		#pragma omp parallel for
		for(size_t i = 0; i < fit_num_params.size(); i++)
			_cmh.read_covmat((ch._covfile + to_string(i)).c_str(), i);
	
		//set up the template matrix for faster computation during the pre-/main-run
		_cmh.set_templatemat(_power);
		
		cout << "complete" << endl;
	}
	else
	{
		cout << "setting up errors..." << flush;
		
		//setting the number of covariance matrices = number of bins
		_cmh.set_num_points(fit_num_params.size());
		
		//setting only the diagonal terms into the covariance matrices for every bin based upon the interpolation file
		for(size_t i = 0; i < fit_errors.size(); i++)
			_cmh.set_diagmat(fit_errors[i], i);
		
		//set up the template matrix for faster computation during the pre-/main-run		
		_cmh.set_templatemat(_power);
		
		cout << "complete" << endl;
	}
	
	//the following resizes are made to reduce the time consumption in the case, that the memory would be allocated and deallocated all the time
	//set the number of the parameter set specific matrices = # of chains
	resize_covmats(ch._nchains);
	
	//set the number of function values / variances = # of chains & # of bins
	fval.resize(ch._nchains);
	fvar.resize(ch._nchains);
	for(size_t i = 0; i < fval.size(); i++)
	{
		fval[i].resize(fit_num_params.size());
		fvar[i].resize(fit_num_params.size());
	}		
	
	rejected_bins = 0;
	
}


// ---------------------------------------------------------
BMSkynet::~BMSkynet()
{
    // destructor
}

/**
 * This function is the heartpiece of the Bayesian Analysis. It calculates Loglikelihood of the parameters
 * @parameters: this is a vector with a length = #parameters in the modell, containing a value for each parameter
 * @logprob: this parameter contains (in the end) the Loglikelihood of the model given @parameters
 * @currentchain: Storage of the # of the chain that is worked on right now. Needed for access the right function values / variances etc.
 * @i, @j: used for for-loops
 * @tmp_parameters: [0,1]-mapping of the @parameters values based upon the limits of their values
 * @anchors: fit parameter independent part of every monomial
**/ 
double BMSkynet::LogLikelihood(const vector<double>& parameters)
{
	//~ if((iteration_number % 1000 == 0) && (GetCurrentChain() == 0))
	//~ {
		//~ for(size_t i = 0; i < parameters.size(); i++)
			//~ if(GetRValueParameters(i) > 1.1)
				//~ break;
			//~ else
				//~ if(i == parameters.size() - 1 && GetRValueParameters(i) > 1.1)
					//~ SetScaleFactorLowerLimit(0);
	//~ }
	vector<vector<double>> specific_mat;
	vector<double> fval, fvar;
	
	fval.resize(fit_num_params.size());
	fvar.resize(fit_num_params.size());

	double logprob = 0;
	const size_t currentchain = GetCurrentChain();
	size_t i, j;
			
	
	//if the program is in its mainrun, the history of the evolution of the Markovchains will be saved.
	if(fMCMCPhase == BCEngineMCMC::kMainRun)
		updateHistory();
		
	//if the program is in its prerun, the history of the evolution of the Markovchains will be saved, if wanted.
	if(fMCMCPhase == BCEngineMCMC::kPreRun && _prerunplots)
		updateHistory_prerun();

	vector<double> tmp_parameters, anchors;
	tmp_parameters.resize(parameters.size());
	anchors.assign(_power.size(), 1);
	
	/**
	 * Construction of the nested parallel region. In this region every time consuming variable will be calculated. So, the logprob is more or less a simple look-up of values.
	 * The order of the computations is based upon the highest efficiency due to synchronization and earlier computations. This is important due to the high amount
	 * of necessary tasks.
	 */
	//~ #pragma omp parallel default(shared) firstprivate(currentchain)
	//~ {

		//calculating the [0,1]-mapping of the parameter
		//~ #pragma omp for
		for(j = 0; j < _power[0].size(); j++)
			tmp_parameters[j] = (parameters[j] - var_minvalue[j]) / (var_maxvalue[j] - var_minvalue[j]);
			
		
		//calculating the parameter specific matrix
		//This matrix is needed for the variance calculation, so it is no barrier needed at that point.
		//The calculation can start, if @tmp_parameters are calculated
		//~ #pragma omp single nowait		
		_cmh.get_specific_mat(tmp_parameters, specific_mat);
	

		//calculating @anchors; a barrier is needed for @anchors and for @_specific_mat because they are both needed in the following for-loop
		//~ #pragma omp for	private(j) schedule(dynamic)
		for(i = 0; i < _power.size(); i++)	
			for(j = 0; j < _power[0].size(); j++)
				anchors[i] *= pow(tmp_parameters[j], _power[i][j]);

				
		//calculating the function value / variance for every bin
		//~ #pragma omp for private(j) schedule(dynamic)
		for(i = 0; i < ref_values.size(); i++)
			for(j = 0; j < ref_values[i].size(); j++)
			{
				if(weights_mat[i][j] == 0)
					continue;
					
				fvar[num_bins_cumulative[i] + j] = fit_var(specific_mat, num_bins_cumulative[i] + j);
				fval[num_bins_cumulative[i] + j] = fit_value(anchors, num_bins_cumulative[i] + j);
				
				if(fvar[num_bins_cumulative[i] + j] < 0 && GetCurrentIteration() > 0 && weights_mat[i][j] != 0 
					&& fabs(fvar[num_bins_cumulative[i] + j]) > ref_values_errors[i][j] * ref_values_errors[i][j])
				{
					#pragma omp critical (reset)
					{
						weights_mat[i][j] = 0;
						rejected_bins++;
						cout << "bin #" << num_bins_cumulative[i] + j << " set to 0 by chain #" << currentchain << " in iteration " << GetCurrentIteration() << ". Overall deleted " << rejected_bins << " bin(s)" << endl 
							<< "variance of function: " << fvar[num_bins_cumulative[i] + j] << ", variance of data: "
							<< ref_values_errors[i][j] * ref_values_errors[i][j] << endl << "parameter values: ";
						for(size_t index = 0; index < parameters.size(); index++)
							cout << parameters[index] << "\t";
						cout << endl << i << "\t" << j << endl << endl;
						
						//~ fMCMCprob.clear();
						//~ fMCMCLogLikelihood.clear();	
						//~ fMCMCRValueParameters.clear();
					}

				}
			}

	//~ }
	
	
	double asdf;
	//calculating the Loglikelihood by walking over every bin and getting the LogGaus of the interpolation and the reference values of the current bin.
	//the weights improve/reduce the importance of a certain bin
	for(i = 0; i < ref_values.size(); i++)
		for(j = 0; j < ref_values[i].size(); j++)
		{
			if(weights_mat[i][j] == 0)
				continue;
				
			asdf = weights_mat[i][j] * BCMath::LogGaus(fval[num_bins_cumulative[i] + j], ref_values[i][j], 
						sqrt(ref_values_errors[i][j] * ref_values_errors[i][j] + fvar[num_bins_cumulative[i] + j]));
			//~ asdf = - weights_mat[i][j] * (fval[num_bins_cumulative[i] + j] - ref_values[i][j]) * (fval[num_bins_cumulative[i] + j] - ref_values[i][j])
						//~ / (ref_values_errors[i][j] * ref_values_errors[i][j] + fvar[num_bins_cumulative[i] + j]);	//modified
						
						
			if(std::isnan(-asdf))
				cout << currentchain << "\t" << asdf << "\t" << weights_mat[i][j] << "\t" << fval[num_bins_cumulative[i] + j] << "\t" << ref_values[i][j]
				<< "\t" << ref_values_errors[i][j] * ref_values_errors[i][j] << "\t" << fvar[num_bins_cumulative[i] + j] << endl;
				
				
			//~ if(fMCMCStatistics[currentchain].efficiency[0] == 0 && fvar[currentchain][num_bins_cumulative[i] + j] == 0)
				//~ cout << currentchain << "\t" << fMCMCprob[currentchain] << "\t" << asdf << endl;
		
			//~ logprob += weights_mat[i][j] * BCMath::LogGaus(fval[currentchain][num_bins_cumulative[i] + j], ref_values[i][j], 
						//~ sqrt(ref_values_errors[i][j] * ref_values_errors[i][j] + fvar[currentchain][num_bins_cumulative[i] + j]));
			logprob += asdf;
		}

	
	//Return the log of the conditional probability p(data|pars).
    return logprob;
}

/**
 * This function calculates the function value of the interpolation function in a certain bin
 * @anchors: fit parameter independent part of every monomial
 * @index: bin number
 * @result: function value 
 */
double BMSkynet::fit_value(const vector<double>& anchors, const int& index){
	
	double result = 0;
	
	//walk over every monomial and multadd the fit parameter and the rest together
	for(size_t i = 0; i < fit_num_params[index]; i++)
		result += fit_params[index][i] * anchors[i];
	
	return result;
}

/**
 * This function calculates the variance of the interpolation function in a certain bin
 * @mat: matrix based upon the Jacobian matrix with a fixed parameter set
 * @index: bin number
 */
double BMSkynet::fit_var(const vector<vector<double>>& mat, size_t index)
{
	//pass the arguments to @_cmh, it can calculate it
	return _cmh.get_var(mat, index);
	
}

/**
 * This function plots the Loglikelihood of every single bin using the best fit parameters estimated by BAT
 * @mode: parameter for the mode of the parameters
 * @tmp_parameters: [0,1]-mapping of the @parameters values based upon the limits of their values
 * @tmp_specific_mat: parameter depending specific matrix for error computation
 * @anchors: fit parameter independent part of every monomial
 * @maximum, @minimum: parameters for the maximum/minimum value of the binwise Loglikelihood, so that the plotranges can be set
 * @tc, @th1d: parameter for a storageframe and the histogram
**/ 
void BMSkynet::binlogprobdistribution(){

	//calculating the mode of the parameters using Minuit
	std::vector<double> mode = FindMode(GetBestFitParameters());
	
	//calculating the [0,1]-mapping of the parameter
	std::vector<double> tmp_parameters;
	tmp_parameters.resize(mode.size());
	
	for(size_t i = 0; i < mode.size(); i++)
		tmp_parameters[i] = (mode[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);
	
	//calculating @tmp_specific_mat
	std::vector<std::vector<double>> tmp_specific_mat = _cmh.get_specific_mat(tmp_parameters);
	
	//calculating @anchors
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
			double tmp_var = _cmh.get_var(tmp_specific_mat, num_bins_cumulative[i] + j);
			if(BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + tmp_var)) > maximum)
				maximum = BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + tmp_var));
			if(BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + tmp_var)) < minimum)
				minimum = BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + tmp_var));
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
		{
			th1d->Fill(BCMath::LogGaus(fit_value(anchors, num_bins_cumulative[i] + j), ref_values[i][j], sqrt(ref_values_errors[i][j] 
					* ref_values_errors[i][j] + _cmh.get_var(tmp_specific_mat, num_bins_cumulative[i] + j))));
		}

	//draw and save
	th1d->Draw();
	
	gPad->Update();
	tc->Print("binlogprobdistribution.png");
	tc->Clear();
}

/**
 * This function plots the Chi^2 value of each bin using the best fit parameters estimated by BAT and compares them to the referencedate of the respective bin
 * @mode: parameter for the mode of the parameters
 * @tmp_parameters: [0,1]-mapping of the @parameters values based upon the limits of their values
 * @tmp_specific_mat: parameter depending specific matrix for error computation
 * @anchors: fit parameter independent part of every monomial
 * @maximum: this is a calculation of the upper plotrange limit
 * @tc, @th1d: parameter for a storageframe and the histogram
**/ 
void BMSkynet::chi2plot(){

	//calculating the mode of the parameters using Minuit
	std::vector<double> mode = FindMode(GetBestFitParameters());
	
	//calculating the [0,1]-mapping of the parameter
	std::vector<double> tmp_parameters;
	tmp_parameters.resize(mode.size());
	
	for(size_t i = 0; i < mode.size(); i++)
		tmp_parameters[i] = (mode[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);
	
	//calculating @tmp_specific_mat
	std::vector<std::vector<double>> tmp_specific_mat = _cmh.get_specific_mat(tmp_parameters);
	
	//calculating @anchors
	std::vector<double> anchors;
	anchors.assign(_power.size(), 1);
	
	for(size_t i = 0; i < _power.size(); i++)
		for(size_t j = 0; j < _power[0].size(); j++)
			anchors[i] *= pow(tmp_parameters[j], _power[i][j]);
		
		
	//walk over all different bins, calculate the Loglikelihood of the bin and searching for its maximum, so that the plotrange can be set later
	double maximum = 0;

	for(size_t i = 0; i < ref_values.size(); i++)
		for(size_t j = 0; j < ref_values[i].size(); j++)
		{
			double tmp_var = _cmh.get_var(tmp_specific_mat, num_bins_cumulative[i] + j);
			if((fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + tmp_var) > maximum)
				maximum = (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
					/ (ref_values_errors[i][j] * ref_values_errors[i][j] + tmp_var);
		}
	

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
		{
			th1d->Fill((fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) * (fit_value(anchors, num_bins_cumulative[i] + j) -  ref_values[i][j]) 
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + _cmh.get_var(tmp_specific_mat, num_bins_cumulative[i] + j)));
		}


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
 * @tmp_parameters: [0,1]-mapping of the @parameters values based upon the limits of their values
 * @tmp_specific_mat: parameter depending specific matrix for error computation
 * @anchors: fit parameter independent part of every monomial
 * @ndof: parameter for calculating the number of dof
**/ 
void BMSkynet::chi2distribution(){

	TCanvas* tc = new TCanvas("chi2distribution.png");
	tc->cd();

	double sum = 0;
	//calculating the mode of the parameters by using Minuit
	std::vector<double> mode = FindMode(GetBestFitParameters());

	//calculating the [0,1]-mapping of the parameter
	std::vector<double> tmp_parameters;
	tmp_parameters.resize(mode.size());
	
	for(size_t i = 0; i < mode.size(); i++)
		tmp_parameters[i] = (mode[i] - var_minvalue[i]) / (var_maxvalue[i] - var_minvalue[i]);
	
	//calculating @tmp_specific_mat
	std::vector<std::vector<double>> tmp_specific_mat = _cmh.get_specific_mat(tmp_parameters);
	
	//calculating @anchors
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
				/ (ref_values_errors[i][j] * ref_values_errors[i][j] + _cmh.get_var(tmp_specific_mat, num_bins_cumulative[i] + j));
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

/**
 * Setter for the # of @_specific_mat entries = # of chains
 */
void BMSkynet::resize_covmats(size_t nchains){
	
	//resizing the vector
	_specific_mat.resize(nchains);
	
}
