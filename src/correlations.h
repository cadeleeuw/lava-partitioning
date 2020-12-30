/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef CORRELATIONS_H
#define CORRELATIONS_H

#include <utility>

#include "data.h"                  

struct MatrixRow {
  float *begin, *end;

  MatrixRow() : begin(0), end(0) {}
  MatrixRow(float* begin, float* end) : begin(begin), end(end) {}

  int length() {return end-begin;}
  float& value(int index) {return begin[index];}
};

// for now, assuming regular matrix such that length of each row is either same as previous row or exactly one more
// and hence also that max. length of rows equals length of last row
class CorrelationMatrix {
  vector<MatrixRow> rows;
  vector<double> means;
  int block_offset, max_depth;

  double block_mean(int row_index); // for split right after row_index
  void compute_metric(int from=0, int to=-1);  

  void init();
  CorrelationMatrix(vector<MatrixRow>& input, vector<double>& input_means, int offset, bool trim); 
public:
  CorrelationMatrix(vector<MatrixRow>& input, int offset); 

  int get_offset() {return block_offset;}
  int get_size() {return rows.size();}
  const vector<double>& get_metric() {return means;}
  
  CorrelationMatrix* split(int index); //split after SNP at index (NB: index range starts at 0, not at CM offset!)
};


class Correlations {
  int depth, storage_size;
  
  vector<Buffer<float>*> storage;
  vector<MatrixRow> rows;
  vector<pair<int,int> > positions; //base pair position, full data SNP index

  class DataIterator;

  double compute_correlation(float* v1, float* v2, int n);
  void clear_storage();
  
public:
  Correlations(int snp_depth) : depth(snp_depth), storage_size(10000) {}
  ~Correlations() {clear_storage();}

  void compute(GenoData& data_src);
  int compute_block(GenoData& data_src, int from, int to);  

  int get_size() {return rows.size();}
  int get_depth() {return depth;}
  const vector<pair<int,int> >& get_positions() {return positions;} 
  CorrelationMatrix* get_matrix() {return new CorrelationMatrix(rows, 0);}  
};

class Correlations::DataIterator {
  GenoData& data_src;
  int block_size;
  
  pair<Buffer<float>*,Buffer<float>*> data;
  vector<float*> snps;
  vector<pair<int,int> >& positions;
  int offset, curr_lead, curr_trail;



public:
  DataIterator(GenoData& gd, vector<pair<int,int> >& pos, int size) : data_src(gd), positions(pos), block_size(size+1), data(0,0) {positions.clear();}
  ~DataIterator() {delete data.first; delete data.second;}

  float* advance_lead();
  float* get_next();   
};


#endif /* CORRELATIONS_H */
