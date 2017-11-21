#include "RefHandler.h"

using namespace std;

/**
 * Constructor that sets up the anchor points and simulated reference values
 * @ch: object that keeps information about the files that need to be read
 */
RefHandler::RefHandler(ConfigHandler& ch, SDataHandler& sdh){

	//reading the information
	readrefdata(ch, sdh);
	
}

/**
 * This function reads the referencedata provided in a .root format.
 * These files need to be in the same directory as the program is started.
 * @dr: storage of the informations, gathered from the interpolation file (e.g. the concerned Analysis)
 * @num_histos: this parameter stores the number of all observables, that are taken into account
 * @tf: this parameter is needed in order to simply open the needed datafile
 * @tgae: this parameter extracts the needed observable out of the datafile
**/ 
void RefHandler::readrefdata(ConfigHandler& ch, SDataHandler& sdh){

	size_t num_histos = sdh.get_num_histos();
	vector<double> tmp, tmp_err;
	_num_bins.resize(num_histos);
	
	ref_values.resize(num_histos);
	ref_values_errors.resize(num_histos);
	
	cout << "start reference data reading ..." << flush;
	
	//walk over all histograms named in @dr
	//~ #pragma omp parallel for private(tmp, tmp_err, ref_num_bins)
	for(size_t i = 0; i < num_histos; i++)
	{
		//reading the root file and converting it to a data type,
		//that allows reading the properties of the plots
		TFile tf((sdh.get_analysis()[sdh.getchanging_histo()[i]] + ".root").c_str());
		TGraphAsymmErrors* tgae = (TGraphAsymmErrors*) tf.FindObjectAny(sdh.get_observable()[sdh.getchanging_histo()[i]].c_str());

		//reading the number of bins and the content of each
		_num_bins[i] = tgae->GetN();
		//~ cout << (sdh.get_analysis()[sdh.getchanging_histo()[i]] + ".root").c_str() << "\t" << ref_num_bins <<  endl;
		//TODO:check, ob fehler stets symmetrisch sind
		//simple push_back on values and errors
		for(size_t j = 0; j < _num_bins[i]; j++)
		{
			tmp.push_back(tgae->GetY()[j]);
			tmp_err.push_back(tgae->GetErrorY(j));
			//~ cout << tmp.back() << "\t" << tmp_err.back() << endl;
		}
		
		ref_values[i] = tmp;
		tmp.clear();
		ref_values_errors[i] = tmp_err;
		tmp_err.clear();
		
		
	
	}
	cout << "complete" << endl;
	
}

/**
 * This function is a getter for the number of bins in each histogram
 */
vector<size_t> RefHandler::getnum_bins(){
	
	return _num_bins;

}

vector<vector<double>> RefHandler::getref_values(){
	
	return ref_values;
	
}

vector<vector<double>> RefHandler::getref_values_errors(){
	
	return ref_values_errors;
	
}
