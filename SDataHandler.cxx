#include "SDataHandler.h"

using namespace std;

/**
 * Constructor, reading through the interpolationfile
 * @interpolation: name of the interpolationfile
 * @ifile: reader for the file
**/ 
SDataHandler::SDataHandler(string interpolation, int sameorder){

	ifstream ifile(interpolation);

	if (ifile.is_open()) 
	{
		cout << "reading interpolation file ...";

		//linewise stepping through the file
		while (getline(ifile, line))
		{
			//checking and reading the necessary informations if they are in that line
			//the check is done by comparing the first some characters of each line with certain triggerwords
			//TODO:if abragen rausziehen
			length = line.length();
			
			parameter_reading();
			
			var_min_max_reading();
			
			fit_params_reading(sameorder);
			
			analysis_name();
			
			set_runs();
		}

		ifile.close();
		cout << "complete" << endl;
	}
	else
		cout << "unable to read interpolation file" << endl;
}

/**
 * This function extract the number and names of the parameters
**/ 
void SDataHandler::parameter_reading(){
	//just extracting the line with the parameter names
	//this will be used after the amount of parameters is known
	//this is done due to the formation of the file
	if (line.substr(0, 10) == "ParamNames")
	{
		parameternames = line;
	}

	//reading the dimension and the parameter names
	if (line.substr(0, 9) == "Dimension")
	{
		//splitting the string, so that only the number of parameters remain
		pos_begin = 10;
		pos_end = line.find(" ", pos_begin + 1);

		dimension = atoi(line.substr(pos_begin, pos_end - pos_begin).c_str());

		//splitting @parameternames, so that the parameternames can be extractet
		pos_begin = 11;
		while(var_names.size() < dimension)
		{
			
			pos_end = parameternames.find(" ", pos_begin + 1);
			var_names.push_back(parameternames.substr(pos_begin + 1, pos_end - pos_begin - 1));
			pos_begin = pos_end;
		}
	}

}

/**
 * This function extract the minimum and maximum value of the parameters that was sampled
**/ 
void SDataHandler::var_min_max_reading(){

	if (line.substr(0, 12) == "MinParamVals")
	{
		//splitting the string in order to get the minimum values
		pos_begin = 14;

		while(var_minvalue.size() < dimension)
		{
			pos_end = line.find(" ", pos_begin + 1);
			var_minvalue.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));
			pos_begin = pos_end;
		}
		
	}

	if (line.substr(0, 12) == "MaxParamVals")
	{
		//splitting the string in order to get the maximum values
		pos_begin = 14;

		while(var_maxvalue.size() < dimension)
		{
			pos_end = line.find(" ", pos_begin + 1);
			var_maxvalue.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));
			pos_begin = pos_end;
		}
	
  	}

}

/**
 * This function extract the interpolationparameters
**/ 
void SDataHandler::fit_params_reading(int sameorder){

	vector<double> tmp;
	
	if (line.substr(2, 3) == "val")
	{
		
		//reading the value of every parameter in every bin
		pos_begin = 10;

		while (pos_begin + 1 < length) 
		{
		
			pos_end = line.find(" ", pos_begin + 1);
			tmp.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));		
			pos_begin = pos_end;

		}

		if(!sameorder)
			while(tmp.back() == 0)
				tmp.pop_back();
		
		fit_params.push_back(tmp);
		tmp.clear();
		
		
		
		//calculating number of parameters in use
		fit_num_params.push_back(fit_params.back().size());
	}
	
	if (line.substr(2, 3) == "err")
	{
	
		//reading the binwise error of the interpolation
		pos_begin = 10;

		while (pos_begin + 1 < length) 
		{
			pos_end = line.find(" ", pos_begin + 1);
			tmp.push_back(atof(line.substr(pos_begin, pos_end - pos_begin).c_str()));		
			pos_begin = pos_end;
		}
		
		if(!sameorder)
			while(tmp.back() == 0)
				tmp.pop_back();
			
		fit_error.push_back(tmp);
	}
}

/**
 * This function extract names of the analysises and corresponding observables.
 * All of them will be stored in one list each and the same name will be stored several times, for each bin once.
**/ 
void SDataHandler::analysis_name(){

	if (line.substr(0, 1) == "/")
	{

		//reading the analysis name
		pos_begin = 1;
		
		pos_end = line.find("/", pos_begin + 1);
		analysis.push_back(line.substr(pos_begin, pos_end - pos_begin));

		//reading the observable name
		pos_begin = pos_end + 1;

		pos_end = line.find("#", pos_begin);
		//~ if(pos_end > line.find(".dat", pos_begin))
			//~ observable.push_back(line.substr(pos_begin, pos_end - pos_begin - 4));
		//~ else
		observable.push_back(line.substr(pos_begin, pos_end - pos_begin));
		
		//If there are no observables saved yet, the index of the first changing histogram will be set to 0.
		if(observable.size() == 1)
			changing_histo.push_back(0);

		//If a new observable was read, its index will be added to @changing_histo.
		if(observable.size() > 1)
			if(observable.back() != observable[observable.size() - 2] || analysis.back() != analysis[analysis.size() - 2])
				changing_histo.push_back(observable.size() - 1);
	}

}

/**
 * This function gets the number of different observables mentioned in the interpolationfile
 * @tmp_analysis, @tmp_observable: temporary storage of a name; indicator of the change while walking through @analysis and @observable.
**/ 
size_t SDataHandler::get_num_histos(){

	//if there are no analysises, it returns 0.
	if(analysis.size() == 0)
		return 0;

	string tmp_analysis = analysis[0];
	string tmp_observable = observable[0];
	size_t num_histos = 1;

	//walk over all entries in @analysis
	for(size_t i = 1; i < analysis.size(); i++)
	{
		//if there is a change in the name compared to the predecessor, the strings will be updated to the value of the last ones and the counter of observables will be increased.
		if(analysis[i] != tmp_analysis || observable[i] != tmp_observable)
		{
			tmp_analysis = analysis[i];
			tmp_observable = observable[i];
			num_histos++;
		}
	}

	return num_histos;

}

/**
 * This function extract the number of parametersets (runs) that were used to create the interpolation and the number, every run is identified with.
**/ 
void SDataHandler::set_runs(){
	
	//extract the total number of runs
	if(num_runs == 0 && line.substr(0, 9) == "NumInputs")
	{
		pos_begin = 10;
		pos_end = line.find(" ", pos_begin + 1);

		num_runs = atoi(line.substr(pos_begin, pos_end - pos_begin).c_str());

	}

	//extract the identifiernumber of each run
	if(line.substr(0, 4) == "Runs")
	{
		pos_begin = 5;
		for(size_t i = 0; i < num_runs; i++)
		{
			pos_end = line.find(" ", pos_begin + 1);
			
			runs.push_back(line.substr(pos_begin + 1, pos_end - pos_begin - 1));
			
			pos_begin = pos_end;
		}
		
	}


}

size_t SDataHandler::getnum_runs(){
	
	return num_runs;
	
}

vector<string> SDataHandler::get_runs(){
	
	return runs;
	
}

size_t SDataHandler::get_dimension(){
	
	return dimension;
	
}

vector<string> SDataHandler::getvar_names(){
	
	return var_names;
	
}

vector<string> SDataHandler::get_analysis(){
	
	return analysis;
	
}

vector<string> SDataHandler::get_observable(){
	
	return observable;
	
}

vector<size_t> SDataHandler::getchanging_histo(){
	
	return changing_histo;
	
}

vector<double> SDataHandler::getvar_minvalue(){
	
	return var_minvalue;
	
}

vector<double> SDataHandler::getvar_maxvalue(){
	
	return var_maxvalue;
	
}

vector<size_t> SDataHandler::getfit_num_params(){
	
	return fit_num_params;
	
}

vector<vector<double>> SDataHandler::getfit_params(){
	
	return fit_params;
	
}

vector<vector<double>> SDataHandler::getfit_error(){
	
	return fit_error;
	
}

