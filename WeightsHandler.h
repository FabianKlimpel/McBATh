#ifndef __WEIGHTSHANDLER__H
#define __WEIGHTSHANDLER__H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "ConfigHandler.h"
#include "SDataHandler.h"

using namespace std;

/**
* This class handles a weight file
*/
class WeightsHandler
{
public:
	//constructor that takes the name of the file
	WeightsHandler(string weightsfile);

	//sorts the readed weights in a way, that they are in the same order as in the ConfigHandler
	void sort(SDataHandler& sdh);

	//getter for the weigths
	vector<int> getweights();

	vector<vector<int>> get_weights_mat(vector<vector<double>>& bins, string weightsfile, SDataHandler& sdh);
	
	
private:

	//datacontainer for the weights, the analyses and the observables
	vector<int> _weights;
	vector<string> _analysis, _observable; 

};

#endif
