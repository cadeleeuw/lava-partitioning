/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "correlations.h"

CorrelationMatrix::CorrelationMatrix(vector<MatrixRow>& input, vector<double>& input_means, int offset, bool trim) : block_offset(offset) {
  rows.swap(input); means.swap(input_means);
  init();
  
  if (trim) {
    int trim_index = -1;
    for (int i = 0; i < rows.size(); i++) {
      if (rows[i].length() > i) {
        trim_index = i;
        rows[i].begin = rows[i].end - i;
      }
    }
    max_depth = rows.back().length();
    compute_metric(0, trim_index+1);
  }
}

CorrelationMatrix::CorrelationMatrix(vector<MatrixRow>& input, int offset) : rows(input), block_offset(offset) {init();}

void CorrelationMatrix::init() {
  if (block_offset < 0) error("invalid offset for CorrelationMatrix object");
  if (rows.empty()) error("input for CorrelationMatrix object is empty");
  
  max_depth = rows.back().length();
  if (means.empty()) compute_metric();
}

double CorrelationMatrix::block_mean(int row_index) {
  double sum = 0; int count = 0;
  int offset = row_index+1, stop = rows.size() - offset;

  for (int r = 0; r < stop; r++) {
    MatrixRow& curr = rows[offset+r]; int use = curr.length() - r;
    if (use > 0) {
      for (float *read = curr.begin, *end = curr.end - r; read < end; ++read) sum += *read;
      count += use;
    }
    else break;
  }
  return sum/count;
}


void CorrelationMatrix::compute_metric(int from, int to) {
  if (from < 0) from = 0;
  if (to < 0 || to >= rows.size()) to = rows.size() - 1;

  means.resize(rows.size()-1);
  for (int i = from; i < to; i++) means[i] = block_mean(i);
}

CorrelationMatrix* CorrelationMatrix::split(int index) {
  if (index < 0 || index >= (rows.size()-1)) error("invalid split index for correlation matrix");
  
  vector<MatrixRow> split_rows(rows.begin()+index+1, rows.end());
  vector<double> split_means(means.begin()+index+1, means.end());
  
  rows.resize(index+1); 
  compute_metric(index - max_depth + 1);
  max_depth = rows.back().length();

  return new CorrelationMatrix(split_rows, split_means, block_offset+index+1, true);
}


void Correlations::compute(GenoData& data_src) {
  clear_storage();
  DataIterator data(data_src, positions, depth);

  storage.push_back(new Buffer<float>(depth*(depth+1)/2.0));
  float *write = storage.back()->get_data(), *end = write + storage.back()->size();
  int N = data_src.get_nrow(); 
  while (float* lead = data.advance_lead()) {
    if (end-write < depth) {
      storage.push_back(new Buffer<float>(depth,storage_size));    
      write = storage.back()->get_data(), end = write + depth*storage_size;
    }

    MatrixRow row(write,0);
    while (float* trail = data.get_next()) {
      float r = compute_correlation(lead, trail, N);
      *(write++) = r*r;
    }

    row.end = write; 
    rows.push_back(row);      
  }
  if (rows.size() != positions.size()) error("number of SNP positions does not match size of correlation matrix");
}

int Correlations::compute_block(GenoData& data_src, int from, int to) {
  clear_storage();
  
  Buffer<float> data;
  pair<int,int> loaded = data_src.load_data(data, positions, from, to-from);
  storage.push_back(new Buffer<float>(depth*loaded.first));
  
  float *write = storage.back()->get_data(); 
  int N = data_src.get_nrow();   
  for (int lead = 0; lead < loaded.first; lead++) {
    MatrixRow row(write,0);
    for (int trail = max(lead-depth,0); trail < lead; trail++) {
      float r = compute_correlation(data.get_column(lead), data.get_column(trail), N);
      *(write++) = r*r;      
    }
    row.end = write; 
    rows.push_back(row);      
  }
  if (rows.size() != positions.size()) error("number of SNP positions does not match size of correlation matrix");
  return rows.size();
}

double Correlations::compute_correlation(float* v1, float* v2, int n) {
  double sum = 0; float* end = v1 + n;
  while (v1 < end) sum += *(v1++) * *(v2++);
  return sum / (n-1);
}

void Correlations::clear_storage() {
  for (int i = 0; i < storage.size(); i++) delete storage[i];
  storage.clear(); positions.clear(); rows.clear();
}


float* Correlations::DataIterator::advance_lead() {
  if (!data.first) {    
    data.first = new Buffer<float>();
    data.second = new Buffer<float>();    
    snps.assign(2*block_size, 0);
    
    pair<int,int> loaded = data_src.load_data(*data.first, positions, 0, block_size); offset = loaded.second;
    for (int i = 0; i < loaded.first; i++) snps[i] = data.first->get_column(i);
    
    loaded = data_src.load_data(*data.second, positions, offset, block_size); offset += loaded.second;
    for (int i = 0; i < loaded.first; i++) snps[block_size+i] = data.second->get_column(i);
    curr_lead = 0;
  } else {
    curr_lead++; 
    if (curr_lead >= 2*block_size) {
      for (int i = 0; i < block_size; i++) snps[i] = snps[block_size+i];
      swap(data.first, data.second);

      pair<int,int> loaded = data_src.load_data(*data.second, positions, offset, block_size); offset += loaded.second;
      for (int i = 0; i < block_size; i++) snps[block_size+i] = (i < loaded.first) ? data.second->get_column(i) : 0;
      
      curr_lead = block_size; 
    }
  } 

  curr_trail = max(-1, curr_lead - block_size);   
  return snps[curr_lead];
}

float* Correlations::DataIterator::get_next() {
  curr_trail++; 
  return curr_trail < curr_lead ? snps[curr_trail] : 0;
}
