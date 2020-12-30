/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef SPLITTER_H
#define SPLITTER_H

#include <deque>

#include "correlations.h"
#include "data.h"

struct Split {
  int offset, position; double metric, metric_min;
  Split() : offset(-1), position(-1), metric(2), metric_min(2) {};
  void set(int o, double m) {offset = o; metric = m; metric_min = min(metric_min, m);}
};

class Splitter {
  int min_size; double min_prop;
  double metric_margin, metric_max;

  deque<CorrelationMatrix*> blocks;
  vector<Split> break_points;

  void insert_block(CorrelationMatrix* block);
  void clear();

public:
  Splitter(Settings& settings);
  ~Splitter() {clear();}

  int run(CorrelationMatrix* input);
  const vector<Split>& get_breaks() {return break_points;}
};

class Refiner {
  GenoData& data;
  vector<Split> breaks;
  vector<pair<int,int> > positions;

public:
  Refiner(GenoData& data) : data(data) {data.set_thresh(0);}

  void refine(Splitter& analysis, Correlations& corrs);

  const vector<Split>& get_breaks() {return breaks;}
  const vector<pair<int,int> >& get_positions() {return positions;}
};

#endif /* SPLITTER_H */
