#ifndef __COUNTER__H
#define __COUNTER__H


#include <vector>

//QUELLE: prof2 -> counter.h / counter.cc / Ipol.cc->mkStructure()
using namespace std;

  class Counter {
  public:

    ///
    Counter(size_t dim, int maxval) {
      for (unsigned int i=0; i< dim;++i) _data.push_back(0);
      _maxval=maxval;
    };
    ~Counter(void);

    bool next(int index);

    int sum();

    vector<int> data() { return _data;}

    void print();

  private:
    int _maxval;
    vector<int> _data;
  };
#endif
