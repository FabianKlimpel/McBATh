#ifndef __POWER__H
#define __POWER__H


#include <vector>
#include "Counter.h"
#include <iostream>

using namespace std;

  class Power {
  public:

   	 ///
   	 Power();

	vector<vector<int>> getpoweroforder(int order);
	void setn(size_t size);
   	
  private:

	void setpoweroforder(int order);

	vector<vector<vector<int>>> _powerlist;
	size_t n;

  };
#endif
