#ifndef __CONFIGHANDLER__H
#define __CONFIGHANDLER__H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#include <dirent.h>

using namespace std;

/**
* This class handles the config file provided at the start of the program.
* It handles the reading of the file and stores all necessary information as the source paths of the data.
*/
class ConfigHandler
{
public:
	ConfigHandler(char* configfile);

	/**
	 * @_ipolfile: name of the interpolation file
	 * @_weights: name of the weights file
	 * @_covfile: basic name of the covariance matrix files
	 * @_pretunefile: name of the Professor file that contains a previously performed tune
	 * @_covmat: flag for usage of covariance matrices
	 * @_sameorder: flag, if every bin is a polynomial of the same order
	 * @_order: order of the polynomial, in case of @_sameorder is true
	 * @_prerunplots: flag for plotting chain history plots during the pre run
	 * @_numprerunplotsteps: number of steps between 2 plots of the chain histories during the pre run
	 */
	string _ipolfile = "", _weights = "", _covfile = "", _pretunefile = "";
	int _covmat = 0, _sameorder = 0, _order = 0, _prerunplots = 0;
	size_t _numprerunplotsteps = 0, _nchains = 1;

	
private:

	//setter for the member variables
	void read_ipolfile(string line);
	void read_weights(string line);
	void read_covmat(string line);
	void read_sameorder(string line);
	void read_order(string line);
	void read_covfile(string line);
	void read_pretunefile(string line);
	void read_prerunplots(string line);
	void read_numprerunplotsteps(string line);
	void read_nchains(string line);
		
};

#endif
