/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef DATA_H
#define DATA_H

#include <vector>
#include <fstream>
#include <cstring>

#include "global.h"

namespace DataUtils {
  template<typename T>
  static bool check_zero_bits() {
    T zero = 0, bits; 
    memset(&bits, 0, sizeof(T));  

    return bits == zero;
  }

  template<typename T>
  static void set_zero(T* buffer, long no_elem) {static bool do_loop = !check_zero_bits<T>();
    if (do_loop) {
      T* end = buffer + no_elem;
      while (buffer < end) *(buffer++) = 0;
    } else memset(buffer, 0, no_elem*sizeof(T));      
  }

  template<typename T>    
  string to_string(T num) {return static_cast<ostringstream*>( &(ostringstream() << num) )->str();}
};


template<typename T>
class Buffer {
protected:
  T* content;
  long capacity, rows, cols;
 
  Buffer<T>& operator=(const Buffer<T>& other) {assign(other);} 

  void set_size(long r, long c) {
    if (r*c > 0) {
      delete content; content = new T[r*c];
      rows = r; cols = c;
      capacity = r*c;
    } else clear();
  }
  
public: 
  Buffer(long init_rows=0, long init_cols=1, bool zero=false) : content(0), capacity(0), rows(0), cols(0) {
    set_size(init_rows, init_cols); 
    if (zero && content) DataUtils::set_zero(content, rows*cols);
  }
  ~Buffer() {clear();}

  void clear() {delete[] content; content = 0; rows = cols = capacity = 0;}
  void resize(long r, long c) {
    if (r*c <= 0) clear();
    else if (r*c > capacity) set_size(r,c);
    else {rows = r; cols = c;}
  }

  int nrow() {return !empty() ? rows : 0;}
  int ncol() {return !empty() ? cols : 0;}
  long size() {return content != 0 ? rows*cols : 0;}
  bool empty() {return size() == 0;}

  T& get(const long& r, const long& c) {return content[r+c*rows];}  
  T& operator[](const long& index) {return content[index];}  
  T* get_column(const long& c) {return content+c*rows;}
  T* get_data() {return content;}
};

namespace DataUtils {
  template<typename T> 
  void print_buffer(Buffer<T>& buff, int precision=6) {
    int rows = buff.nrow(), cols = buff.ncol();
    int prec = cerr.precision();
    cerr.precision(precision);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < cols; c++) cerr << buff.get(r, c) << "\t";
      cerr << endl;
    }          
    cerr.precision(prec);    
  }
};


class GenoData {
  string prefix;
  float maf_thresh;

  ifstream bed_file; 
  unsigned long long block_count;
  Buffer<char> raw_buffer, geno_buffer;   
  unsigned char geno_index[256][4]; 

  int no_indiv, no_snps;
  vector<int> position; //set to zero to skip
  pair<int,int> pos_bounds;

  void read_fam();
  void read_bim();
  void prep_bed();
  
  bool process_snp(char* raw, float*& target);
  pair<int,int> load_block(float* target, vector<pair<int,int> >& pos_target, int offset, int total);
  
public:
  GenoData(const string& prefix, float maf_thresh);

  void set_thresh(float thresh) {maf_thresh = thresh;}
  pair<int,int> load_data(Buffer<float>& target, vector<pair<int,int> >& pos_target, int offset, int total);
  
  int get_nrow() {return no_indiv;}
  int get_nsnps() {return no_snps;}
  pair<int,int> get_bounds() {return pos_bounds;}
};


#endif /* DATA_H */

