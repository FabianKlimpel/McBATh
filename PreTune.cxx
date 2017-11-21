#include "PreTune.h"

/**
 * Constructor
 * The Construction of the object also reads the parameter values
 * @sdh: Contains the names of the parameters to tune
 * @ch: Contains the path to the tune file
 * @ifile: object to read from file
 * @line: storage of a line from the file
 * @line2: part of line; needed to split @line for comparison
 * @pos_begin, @pos_end: helper to split @line right
 */
PreTune::PreTune(SDataHandler& sdh, ConfigHandler& ch){

	//Setting up the read out from file
	ifstream ifile(ch._pretunefile);
	string line, line2;
	size_t pos_begin, pos_end;
	
	//if the file is open ...
	if(ifile.is_open())
	{
		cout << "reading pretune file ..." << flush;
		
		//resizing @inits, needed to store the data without concern about the order of appearance
		inits.resize(sdh.getvar_names().size());
		
		//walk over every line of the file
		while (getline(ifile, line))
		{
			//split the name of the parameter (if it is in that line)
			pos_begin = line.find("\t", 0);
			line2 = line.substr(0, pos_begin);
			
			//throw away all ' ' in order to compare them with the parameter names stored in @sdh
			while(line2.back() == ' ')
				line2.pop_back();
				
			//the splitted string will be compared with every parameter name
			for(size_t i = 0; i < sdh.getvar_names().size(); i++)		
				if(line2 == sdh.getvar_names()[i])
				{	
					//if matched, the value will be stored
					pos_end = line.find("\t", pos_begin + 1);
					inits[i] = atof(line.substr(pos_begin + 1, pos_end - pos_begin - 1).c_str());
				}
	
		}
		cout << "complete" << endl;
			
		//if the tuned parameters are outside the allowed parameter value hypercube, they will be set upon the maximum/minimum value
		for(size_t i = 0; i < inits.size(); i++)
		{
			if(inits[i] > sdh.getvar_maxvalue()[i])
				inits[i] = sdh.getvar_maxvalue()[i];
			if(inits[i] < sdh.getvar_minvalue()[i])
				inits[i] = sdh.getvar_minvalue()[i];
		}
	
		//print the inital values
		cout << "Initial values: " << flush;
		for(size_t i = 0; i < inits.size(); i++)
			cout << inits[i] << "\t";
		cout << endl;
	}
	//else print, that the file could not be read
	else
		cout << "unable to read from " << ch._pretunefile << endl;
}

/**
 * Destructor
 */
PreTune::~PreTune()
{}
	
/**
 * Getter that returns the tuned parameter values
 */
vector<double> PreTune::get_inits(){
	
	return inits;
	
}
