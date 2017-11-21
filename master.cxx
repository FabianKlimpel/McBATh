#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <omp.h>
#include <TCanvas.h>
#include "ConfigHandler.h"
#include "RefHandler.h"
#include "SDataHandler.h"
#include "Batman.h"
#include "WeightsHandler.h"
#include "BMprof.h"
#include "BMSkynet.h"
#include "PreTune.h"

using namespace std;


int main(int argc, char** argv) {


	//check if config file is given
	if(argc > 1)
	{
		//read the config file
		ConfigHandler ch(argv[1]);
		
		//read the fit data
		SDataHandler sdh(ch._ipolfile, ch._sameorder);
		
		//read the reference data
		RefHandler rh(ch, sdh);

		//read the weights and sort them
		WeightsHandler wh(ch._weights);
		wh.sort(sdh);

		
/**
 * the following parts create and use a BAT-Object based on the informations of the previous data
 */
		
		//setting up the omp enviroment
		//The chains are distributed in threads, the likelihood will be calculated in another parallel environment. Therefore its nested.
		//The dynamic setting needs to be disabled in order to force the creation of the threads.
		//The setting of the number of threads is machine dependent.	
		omp_set_dynamic(0);
		omp_set_nested(1);
		omp_set_num_threads(ch._nchains);

 		// set nicer style for drawing than the ROOT default
    	BCAux::SetStyle();
		
    	// open log file
    	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
			
		
			
    	// Create new Batman pointer. Depending on the configuration, the respective object will be created.
    	Batman* m;
		
			
		if(ch._sameorder)
			m = new BMprof("prof", sdh, rh, wh, ch._order);
		else
			m = new BMSkynet("Skynet", sdh, rh, wh, ch);
    	
    	m->SetPrecision(BCEngineMCMC::kMedium);
		m->SetNChains(ch._nchains);
		//~ m->SetMinimumEfficiency(0.05);
		//~ m->SetMaximumEfficiency(0.15);
	
    	//check, if a previously performed tune should be read
    	if(!ch._pretunefile.empty())
		{
			//Read the tune data & set it up as the initial vectors
			PreTune pt(sdh, ch);
			vector<vector<double>> inits;
			for(size_t i = 0; i < ch._nchains; i++)
				inits.push_back(pt.get_inits());
			m->SetInitialPositions(inits);
			
		}
			
		//~ m->SetProposeMultivariate(false);
		
		//the number of iterations should assure, that a convergence can be reached and therefore an optimal mainrun can be performed
		//~ m->SetNIterationsPreRunMax(1000 * m->GetNIterationsPreRunMax());
		m->SetNIterationsPreRunMax(1e6);
		//~ m->SetScaleFactorLowerLimit(1e-16);
		
    	BCLog::OutSummary("Test model created");

		//~ m->SetIntegrationMethod(BCIntegrate::kIntMonteCarlo);
		//~ m->Normalize();

		// run MCMC, marginalizing posterior
    	m->MarginalizeAll(BCIntegrate::kMargMetropolis);
		
    	// run mode finding; by default using Minuit
  		m->FindMode(m->GetBestFitParameters());
		
  		// draw all marginalized distributions into a PDF file
  		m->PrintAllMarginalized(m->GetSafeName() + "_plots.pdf");

  		// print summary plots
   		m->PrintParameterPlot(m->GetSafeName() + "_parameters.pdf");
   		m->PrintCorrelationPlot(m->GetSafeName() + "_correlation.pdf");
   		m->PrintCorrelationMatrix(m->GetSafeName() + "_correlationMatrix.pdf");
   		m->PrintKnowledgeUpdatePlots(m->GetSafeName() + "_update.pdf");

  		// print results of the analysis into a text file
  		m->PrintSummary();

  		// close log file
   		BCLog::OutSummary("Exiting");
  		BCLog::CloseLog();
		
		// print more summaries
		m->chainhistoryplot();
		m->binlogprobdistribution();
		m->chi2distribution();
		m->chi2plot();
		
		// clean up the Batman
		delete m;
	}
	else
		cout << "no filenames specified" << endl;

	return 0;
}
