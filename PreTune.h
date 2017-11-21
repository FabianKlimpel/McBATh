#ifndef __PRETUNE__H
#define __PRETUNE__H

#include "SDataHandler.h"
#include "ConfigHandler.h"
#include <string>
#include <fstream>
#include <vector>

using namespace std;

/**
 * This class reads the parameter values from a previously performed tune. The tune needs to be saved in a Professor layout.
 */
class PreTune
{
	
public:

	/**
	 * Constructor
	 * The Construction of the object also reads the parameter values
	 * @sdh: Contains the names of the parameters to tune
	 * @ch: Contains the path to the tune file
	 */
	PreTune(SDataHandler& sdh, ConfigHandler& ch);
	
	/**
	 * Destructor
	 */
	~PreTune();
	
	/**
	 * Getter of the tuned parameters
	 */
	vector<double> get_inits();
	
private:

	//This is a vector that stores the tuned parameters
	vector<double> inits;

};

#endif
