#include "Counter.h"
#include <iostream>
#include <algorithm>

  Counter::~Counter(void) { }

  bool Counter::next(int index) {
    if (_data[index] == _maxval) {
      if (index==0) return false;

      _data[index]=0;
      return next(index - 1);
    }
    else {
      _data[index]++;
      return true;
    }
  }

  int Counter::sum() {
    int sum_v = 0;
    for (size_t n = 0; n < _data.size(); n++)
      sum_v += _data[n];
    return sum_v;
  }
