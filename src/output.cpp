/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include <iostream>
#include <algorithm> 
#include <cmath>

#include "output.h"

void Output::write(const vector<Split>& break_points, const vector<pair<int,int> >& positions, GenoData& data) {
  string out_name = out_pref + ".breaks";
  cout << "Writing break point output to file '" << out_name << "'" << endl; 

  vector<int> order = Sorter(break_points).run();
  pair<int,int> bounds = data.get_bounds();
  
  ofstream out(out_name.c_str());
  out << "RANK\tMETRIC\tMETRIC_MIN\tINDEX_FILT\tINDEX_ALL\tPOSITION\tPOS_LOWER\tPOS_UPPER" << endl;
  out << "0\t0\t0\t0\t0\t" << bounds.first << "\tNA\tNA" << endl;
  for (int i = 0; i < break_points.size(); i++) {
    Split curr = break_points[order[i]];

    out << (order[i]+1) << "\t" << curr.metric << "\t" << curr.metric_min << "\t" << curr.offset << "\t" << positions[curr.offset].second << "\t";    
    int lower = positions[curr.offset].first, upper = positions[curr.offset+1].first;
    if (curr.position < 0) curr.position = round((lower + upper)/2.0);
    out << curr.position << "\t" << lower << "\t" << upper << endl;
  }
  out << "0\t0\t0\t" << positions.size() << "\t" << data.get_nsnps() << "\t" << bounds.second << "\tNA\tNA" << endl;
}

void Output::write_metrics(const vector<double>& metrics) {
  string out_name = out_pref + ".metric";
  cout << "Writing base metric values to file '" << out_name << "'" << endl; 
  ofstream out(out_name.c_str());
  for (int i = 0; i < metrics.size(); i++) out << metrics[i] << endl;  
}

vector<int> Output::Sorter::run() {
  vector<int> index(data.size());
  for (int i = 0; i < index.size(); i++) index[i] = i;
  sort(index.begin(), index.end(), SortObj(*this));  
  return index;
}
