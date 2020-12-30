/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef OUTPUT_H
#define OUTPUT_H

#include "correlations.h"
#include "splitter.h"
#include "data.h"

class Output {
  string out_pref;

  class Sorter;

public:
  Output(const string& pref) : out_pref(pref) {}

  void write(const vector<Split>& break_points, const vector<pair<int,int> >& positions, GenoData& data);
  void write_metrics(const vector<double>& metrics);
}; 

class Output::Sorter {
  const vector<Split>& data;
  
  struct SortObj {
    Sorter &sorter;
    SortObj(Sorter &sorter) : sorter(sorter) {}
    bool operator () (const long &i, const long &j) {return sorter.data[i].offset < sorter.data[j].offset;}
  };

  
public:
  Sorter(const vector<Split>& data) : data(data) {}

  vector<int> run();
};


#endif /* OUTPUT_H */
