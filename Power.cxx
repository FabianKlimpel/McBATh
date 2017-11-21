#include "Power.h"
#include <iostream>
#include <algorithm>


Power::Power(){}


void Power::setpoweroforder(int order)
{
	
	vector<vector<int>> tmp_power;
	const vector<int> zero(n, 0);
	tmp_power.push_back(zero);
	
	for (int i = 0; i <= order; ++i) 
	{
		Counter c(n, i);
		while (c.next(n-1)) 
			if (c.sum() == i) 
				tmp_power.push_back(c.data());     
	}
				
				
	_powerlist[order] = tmp_power;
	
	tmp_power.clear();	
}

vector<vector<int>> Power::getpoweroforder(int order)
{
	if(order > ((int) _powerlist.size()) - 1)
		_powerlist.resize(order + 1);
	
	if(_powerlist[order].empty())
		setpoweroforder(order);
	
	return _powerlist[order];
	
}

void Power::setn(size_t size){
	
	n = size;
	
}
