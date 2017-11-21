#include "ConfigHandler.h"

/**
 * Constructor that reads a config file and sets all member variables according to the setted ones
 * @configfile: name of the config file
 * @ifile: input stream for the config file
 * @line: string that gets the content of a line in the configfile
 */
ConfigHandler::ConfigHandler(char* configfile){
	
	cout << "reading config file ...";
	ifstream ifile;
	ifile.open(configfile);
	string line;
	
	if(ifile.is_open())
	{
		//linewise file reading
		while(getline(ifile, line))
		{
			//checking for signal words and calling the respective function
			if(line.substr(0, 8) == "ipolfile")
				read_ipolfile(line);
				
			if(line.substr(0, 7) == "weights")
				read_weights(line);
				
			if(line.substr(0, 6) == "covmat")
				read_covmat(line);
				
			if(line.substr(0, 9) == "sameorder")
				read_sameorder(line);
				
			if(line.substr(0, 5) == "order")
				read_order(line);
				
			if(line.substr(0, 7) == "covfile")
				read_covfile(line);
				
			if(line.substr(0, 11) == "pretunefile")
				read_pretunefile(line);
				
			if(line.substr(0, 11) == "prerunplots")
				read_prerunplots(line);
				
			if(line.substr(0, 7) == "prsteps")
				read_numprerunplotsteps(line);
				
			if(line.substr(0, 7) == "nchains")
				read_nchains(line);
		}		
	}

	cout << "complete" << endl;
	
}

/**
 * Setter of the parameters, if they are found in the config file
 */
void ConfigHandler::read_ipolfile(string line){

	_ipolfile = line.substr(9);
	
}

void ConfigHandler::read_weights(string line){
	
	_weights = line.substr(8);

}

void ConfigHandler::read_covmat(string line){
	
	_covmat = atoi(line.substr(7).c_str());
	
}

void ConfigHandler::read_sameorder(string line){
	
	_sameorder = atoi(line.substr(10).c_str());
	
}

void ConfigHandler::read_order(string line){
	
	_order = atoi(line.substr(6).c_str());
	
}

void ConfigHandler::read_covfile(string line){
	
	_covfile = line.substr(8);
	
}

void ConfigHandler::read_pretunefile(string line){
	
	_pretunefile = line.substr(12);
	
}

void ConfigHandler::read_prerunplots(string line){
	
	_prerunplots = atoi(line.substr(12).c_str());
	
}

void ConfigHandler::read_numprerunplotsteps(string line){
	
	_numprerunplotsteps = atoi(line.substr(8).c_str());
	
}

void ConfigHandler::read_nchains(string line){
	
	_nchains = atoi(line.substr(8).c_str());
	
}
