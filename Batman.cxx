// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "Batman.h"

/**
 * Constructor
 * @model_name: name of the model
 * @sdh: container for fit parameter information
 * @rh: container of reference information
 * @wh: container of weight information
 */
Batman::Batman(const std::string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh, ConfigHandler& ch) : BCModel(model_name){

	//storing the data
	var_names = sdh.getvar_names();
	var_minvalue = sdh.getvar_minvalue();
	var_maxvalue = sdh.getvar_maxvalue();
	fit_params = sdh.getfit_params();
	fit_errors = sdh.getfit_error();
	ref_values = rh.getref_values();
	ref_values_errors = rh.getref_values_errors();
	weights_mat = wh.get_weights_mat(ref_values, ch._weights, sdh);
	num_bins_cumulative = sdh.getchanging_histo();
	_numprerunplotsteps = ch._numprerunplotsteps;
	_prerunplots = ch._prerunplots;
	
	//adding the parameter names, minimum/maximum value to the object
	for(size_t i = 0; i < var_names.size(); i++)
	{
		//~ if(i == 2)
			//~ AddParameter(var_names[i], 1, 2);
		//~ else
			AddParameter(var_names[i], var_minvalue[i], var_maxvalue[i]);
    }
    
 

    //Setting the prior. In this case it should be flat for every parameter.
	SetPriorConstantAll();
	
	//~ GetParameter(2).Fix(1.887e-8);
	//~ GetParameter(5).Fix(0.331672);
	
	//Setting the starting value of the iterationnumber for the history of the Markovchains (s.f. @Batman::updateHistory() / @Batman::updateHistory_prerun())
	iteration_number = 0;
	
}

/**
 * Constructor
 * @model_name: name of the model
 * @sdh: container for fit parameter information
 * @rh: container of reference information
 * @wh: container of weight information
 */
Batman::Batman(const std::string& model_name, SDataHandler& sdh, RefHandler& rh, WeightsHandler& wh) : BCModel(model_name){
	
	//storing the data
	var_names = sdh.getvar_names();
	var_minvalue = sdh.getvar_minvalue();
	var_maxvalue = sdh.getvar_maxvalue();
	fit_params = sdh.getfit_params();
	fit_errors = sdh.getfit_error();
	ref_values = rh.getref_values();
	ref_values_errors = rh.getref_values_errors();
	weights = wh.getweights();
	num_bins_cumulative = sdh.getchanging_histo();
	
	//adding the parameter names, minimum/maximum value to the object
	for(size_t i = 0; i < var_names.size(); i++)
    	AddParameter(var_names[i], var_minvalue[i], var_maxvalue[i]);
    

    //Setting the prior. In this case it should be flat for every parameter.
	SetPriorConstantAll();
	
	//Setting the starting value of the iterationnumber for the history of the Markovchains (s.f. @Batman::updateHistory() / @Batman::updateHistory_prerun())
	iteration_number = 0;
	
}


// ---------------------------------------------------------
Batman::~Batman()
{
    // destructor
}

/**
 * This function sets a new point in the history of the evolution of the Markovchains in the main run. 
**/ 
void Batman::updateHistory(){
	//due to the outer functions are always done in parallel, a critical region is needed in order to not mess up the storage of the data
	#pragma omp critical (MCMC_history_update)
	{
	
		//Check, if a chains reachs this point in a new iteration. This simple check is possible due to an implicit barrier at the end of #pragma omp for, that occurs 
		//in @BCEngineMCMC::GetNewPointMetropolis().
		if(iteration_number != GetCurrentIteration())
		{
			//If it is a new iteration, the new diced parameterpoints will be added to the history.
			//The loop is used for the case, that in one or more iterations in a row no parameterset could be diced, that was inside the allowed hypercube.
			//If those values lied outside, there was (due to BAT) simply no new dicing, but just a skip of the iterationsstep.
			for(int i = 0; i < GetCurrentIteration() - iteration_number; i++)
				MCMC_history.push_back(Getx());
				
			//Setting the iterationnumber to the current iteration, leading to the point, that in this iteration no other thread would be able to enter this if-region.
			//Unless this value is updated, every other thread could enter and there could a write error to @MCMC_history occur or the same parameterset could be added multiple times.
			iteration_number = GetCurrentIteration();	
		}

	}
}

/**
 * This function sets a new point in the history of the evolution of the Markovchains in the pre run. 
**/ 
void Batman::updateHistory_prerun(){
	//due to the outer functions are always done in parallel, a critical region is needed in order to not mess up the storage of the data
	#pragma omp critical (MCMC_history_update_prerun)
	{
	
		//Check, if a chains reachs this point in a new iteration. This simple check is possible due to an implicit barrier at the end of #pragma omp for, that occurs 
		//in @BCEngineMCMC::GetNewPointMetropolis().
		if(iteration_number != GetCurrentIteration())
		{
			//If it is a new iteration, the new diced parameterpoints will be added to the history.
			//The loop is used for the case, that in one or more iterations in a row no parameterset could be diced, that was inside the allowed hypercube.
			//If those values lied outside, there was (due to BAT) simply no new dicing, but just a skip of the iterationsstep.
			for(int i = 0; i < GetCurrentIteration() - iteration_number; i++)
				MCMC_history_prerun.push_back(Getx());
				
			//Setting the iterationnumber to the current iteration, leading to the point, that in this iteration no other thread would be able to enter this if-region.
			//Unless this value is updated, every other thread could enter and there could a write error to @MCMC_history_prerun occur or the same parameterset could be added multiple times.
			iteration_number = GetCurrentIteration();
			
			//The number of iterations in the pre run is unknown. In order to track its behaviour, the history should be plotted after a certain amount of iterations.
			if(iteration_number % _numprerunplotsteps == 0)		
				chainhistoryplot_prerun();
				
		}		

	}
}

/**
 * This function creates a graph of the evolution of each Markovchain of each parameter over the iterationsteps in the mainrun
 * @tc, @tg: parameter for a storageframe and the actual graph
 * @xaxisname: string for the axislabel
**/ 
//~ void Batman::chainhistoryplot(){
//~ 
	//~ TCanvas* tc;
	//~ TGraph* tg;
	//~ std::string xaxisname;
//~ 
	//~ //loop over all Markovchains
	//~ for(size_t i = 0; i < GetNChains(); i++)
		//~ //loop over all parameters
		//~ for(size_t j = 0; j < var_names.size(); j++)
		//~ {
			//~ //setting up the graph and its container
			//~ tc = new TCanvas((var_names[j] + "_" + std::to_string(i) + "_history.png").c_str());
			//~ tc->cd();
			//~ tg = new TGraph(MCMC_history.size());
			//~ 
			//~ xaxisname = "Chainevolution;" + var_names[j] + ";Iterationstep";
			//~ tg->SetTitle(xaxisname.c_str());
			//~ 
			//~ //disable the plotpoints; just showing a line connecting the points
			//~ tg->SetMarkerSize(0);
			//~ 
			//~ //walk over all points in the respective history and add them to the graph
			//~ for(size_t k = 0; k < MCMC_history.size(); k++)
				//~ tg->SetPoint(k, MCMC_history[k][i][j], k);
				//~ 
			//~ //draw and save
			//~ tg->Draw();
			//~ gPad->Update();
			//~ tc->Print((var_names[j] + "_" + std::to_string(i) + "_history.png").c_str());
			//~ tc->Clear();
		//~ }
//~ }

void Batman::chainhistoryplot(){

	TCanvas* tc;
	TGraph* tg = new TGraph[GetNChains()];
	for(size_t i = 0; i < GetNChains(); i++)
		tg[i].Set(MCMC_history.size());
	
	std::string xaxisname;

	//loop over all parameters
	for(size_t j = 0; j < var_names.size(); j++)
	{
		tc = new TCanvas((var_names[j] + "_history.png").c_str());
		
		//loop over all Markovchains
		for(size_t i = 0; i < GetNChains(); i++)
		{
			//setting up the graph and its container	
			
			xaxisname = "Chainevolution;Iterationstep;" + var_names[j];
			tg[i].SetTitle(xaxisname.c_str());
			
			//disable the plotpoints; just showing a line connecting the points
			tg[i].SetMarkerSize(0);
			tg[i].SetLineColor(i + 1);
			
			//walk over all points in the respective history and add them to the graph
			for(size_t k = 0; k < MCMC_history.size(); k++)
				tg[i].SetPoint(k, k, MCMC_history[k][i][j]);
				
			//draw and save
			if(i == 0)
			{
				tg[0].Draw();
				//~ tg[0].GetYaxis()->SetRangeUser(var_minvalue[j], var_maxvalue[j]);
				tg[0].GetYaxis()->SetRangeUser(GetParameter(j).GetLowerLimit(), GetParameter(j).GetUpperLimit());
				tg[0].Draw();
			}
			else
				tg[i].Draw("same");

		}
		
		gPad->Update();
		tc->Print((var_names[j] + "_history.png").c_str());
		tc->Clear();
		
	}
	
	delete[] tg;
}

/**
 * This function creates a graph of the evolution of each Markovchain of each parameter over the iterationsteps in the mainrun
 * @tc, @tg: parameter for a storageframe and the actual graph
 * @xaxisname: string for the axislabel
**/ 
void Batman::chainhistoryplot_prerun(){

	TCanvas* tc;
	TGraph* tg = new TGraph[GetNChains()];
	for(size_t i = 0; i < GetNChains(); i++)
		tg[i].Set(MCMC_history_prerun.size());
	
	std::string xaxisname;

	//loop over all parameters
	for(size_t j = 0; j < var_names.size(); j++)
	{
		tc = new TCanvas((var_names[j] + "_history_prerun.png").c_str());
		
		//loop over all Markovchains
		for(size_t i = 0; i < GetNChains(); i++)
		{
			//setting up the graph and its container	
			
			xaxisname = "Chainevolution;Iterationstep;" + var_names[j];
			tg[i].SetTitle(xaxisname.c_str());
			
			//disable the plotpoints; just showing a line connecting the points
			tg[i].SetMarkerSize(0);
			tg[i].SetLineColor(i + 1);
			
			//walk over all points in the respective history and add them to the graph
			for(size_t k = 0; k < MCMC_history_prerun.size(); k++)
				tg[i].SetPoint(k, k, MCMC_history_prerun[k][i][j]);
				
			//draw and save
			if(i == 0)
			{
				tg[0].Draw();
				tg[0].GetYaxis()->SetRangeUser(GetParameter(j).GetLowerLimit(), GetParameter(j).GetUpperLimit());
				//~ tg[0].GetYaxis()->SetRangeUser(var_minvalue[j] - 0.1 * (var_maxvalue[i] - var_minvalue[i]), var_maxvalue[j] + 0.1 * (var_maxvalue[i] - var_minvalue[i]));
				tg[0].Draw();
			}
			else
				tg[i].Draw("same");

		}
		
		gPad->Update();
		tc->Print((var_names[j] + "_history_prerun.png").c_str());
		tc->Clear();
		
	}
	
	delete[] tg;
}

/**
 * This function calculates the new powers for every value in every monomial
 * @max: the maximum power
 */
void Batman::windofchange(int max_power){

	_power = _powlist.getpoweroforder(max_power);

}
