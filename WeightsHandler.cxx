#include "WeightsHandler.h"

using namespace std;

/**
 * Constructor, that reads a given weightfile
 * @weightsfile: name of the file, that contains the informations about the weights
 * @pos_begin, @pos_end: begin und end positions in a string for data extraction
 * @ifile: reader for the file
**/ 
WeightsHandler::WeightsHandler(string weightsfile){

	if(!weightsfile.empty())
	{
		size_t pos_begin, pos_end;
		string line;
		ifstream ifile(weightsfile);

		if (ifile.is_open()) 
		{
			cout << "reading weights file ...";
			getline(ifile, line);
			//linewise stepping through the file
			while (getline(ifile, line)) 
			{
				//"/" is an indicator for content in the line
				if (line.substr(0, 1) == "/")
				{
					//extracting the analysisname
					pos_begin = 1;
					pos_end = line.find("/", pos_begin + 1);
					_analysis.push_back(line.substr(pos_begin, pos_end - pos_begin));
					//extracting the observablename
					pos_begin = pos_end + 1;

					
					if(line.find("#", pos_begin + 1) < line.size())
					{
						pos_end = line.find("#", pos_begin + 1);
						_observable.push_back(line.substr(pos_begin, pos_end - pos_begin));
						pos_begin = line.find("	", pos_end + 1);	
					}
					else
					{
						pos_end = line.find("	", pos_begin + 1);
						_observable.push_back(line.substr(pos_begin, pos_end - pos_begin));					
						pos_begin = pos_end;
					}

					_weights.push_back(atoi(line.substr(pos_begin).c_str()));
					
				}
			}
				

			ifile.close();
			cout << "complete" << endl;
		}
		else
			cout << "unable to read weights file" << endl;
	}
}

/**
 * This function sorts the weights so, that it can be easier combined with the order in the ConfigHandler object.
 * Furthermore, the ConfigHandler object will be truncated if a weight is not mentioned in @_weights
 * @ch: this object contains the analyses and observables in the data folders
 * @tmp_weights: temporary storage of the weights
 * @to_delete: list that stores all observables that will be deleted
 * @i, @j: countingvariables
 * @notfound: flag for adding an index to @to_delete
**/ 
void WeightsHandler::sort(SDataHandler& sdh){

	cout << "sorting weights ...";
	
	
	if(!_weights.empty())
	{
		vector<int> tmp_weights;

		//walk over all histograms (aka @sdh.get_num_histos())
		for(size_t i = 0; i < sdh.get_num_histos(); i++)
		{
			//walk over all weights in this object
			for(size_t j = 0; j < _analysis.size(); j++)	
				//If the observable & analysis of ch and this object matches, add it to @tmp_weights. 
				//Therefore the weights will be in a list in the same order as the bins mentioned in the interpolationfile
				if(sdh.get_analysis()[sdh.getchanging_histo()[i]] == _analysis[j] && sdh.get_observable()[sdh.getchanging_histo()[i]] == _observable[j])
				{
					tmp_weights.push_back(_weights[j]);
					break;
				}				
		}
			
		_weights = tmp_weights;
		_analysis.clear();
		_observable.clear();
	}
	else
		for(size_t i = 0; i < sdh.get_num_histos(); i++)
			_weights.push_back(0);
	
	cout << "complete" << endl;			
}

/**
 * This function return @_weights
 */
vector<int> WeightsHandler::getweights(){
	
	return _weights;
	
}

vector<vector<int>> WeightsHandler::get_weights_mat(vector<vector<double>>& bins, string weightsfile, SDataHandler& sdh){
	
	vector<vector<int>> result;
	result.resize(bins.size());
	
	for(size_t i = 0; i < result.size(); i++)
		result[i].assign(bins[i].size(), 0);
		
			
	size_t pos_begin, pos_end, intervall_start, intervall_end;
	int weight;
	ifstream ifile;
	string line, analysis, observable;
	ifile.open(weightsfile);
	
	if (ifile.is_open()) 
		{
			cout << "reading valid bins ..." << flush;
			getline(ifile, line);

			//linewise stepping through the file
			while (getline(ifile, line)) 
			{
				
				//"/" is an indicator for content in the line
				if (line.substr(0, 1) == "/")
				{
					//extracting the analysisname
					pos_begin = 1;
					pos_end = line.find("/", pos_begin + 1);
			
					analysis = line.substr(pos_begin, pos_end - pos_begin);

					//extracting the observablename
					pos_begin = pos_end + 1;
					
					if(line.find("#", pos_begin + 1) < line.size())
					{
						pos_end = line.find("#", pos_begin + 1);
						observable = line.substr(pos_begin, pos_end - pos_begin);
						
						pos_begin = pos_end;
						pos_end = line.find(":", pos_begin + 1);
						intervall_start = atoi(line.substr(pos_begin + 1, pos_end - pos_begin - 1).c_str());
						
						pos_begin = pos_end;
						pos_end = line.find("	", pos_begin + 1);	
						intervall_end = atoi(line.substr(pos_begin + 1, pos_end - pos_begin - 1).c_str());
						
						weight = atoi(line.substr(pos_end + 1).c_str());	
						
						//walk over all histograms (aka @sdh.get_num_histos())
						for(size_t i = 0; i < sdh.get_num_histos(); i++)
							if(sdh.get_analysis()[sdh.getchanging_histo()[i]] == analysis && sdh.get_observable()[sdh.getchanging_histo()[i]] == observable)
							{
								for(size_t j = intervall_start; j <= intervall_end; j++)
									result[i][j] = weight;
								break;
							}
					}
					else
					{
						pos_end = line.find("	", pos_begin + 1);
						observable = line.substr(pos_begin, pos_end - pos_begin);
					
						weight = atoi(line.substr(pos_end + 1).c_str());	
						
						//walk over all histograms (aka @sdh.get_num_histos())
						for(size_t i = 0; i < sdh.get_num_histos(); i++)
							if(sdh.get_analysis()[sdh.getchanging_histo()[i]] == analysis && sdh.get_observable()[sdh.getchanging_histo()[i]] == observable)
							{
								for(size_t j = 0; j < result[i].size(); j++)
									result[i][j] = weight;
								break;
							}
						
					}
				}	
			}	
		}
	
	ifile.close();
	cout << "complete" << endl;
			
	return result;
}
