/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include <cmath>

#include "data.h"

GenoData::GenoData(const string& prefix, float maf_thresh) : prefix(prefix), maf_thresh(maf_thresh) {
  read_fam();  
  read_bim();
  prep_bed();
}

void GenoData::read_fam() {
  string fname = prefix + ".fam", line;
  ifstream fam(fname.c_str(), ifstream::in);

  cout << "Reading " << fname << "... ";
  no_indiv = 0;
  while (getline(fam, line)) no_indiv++;
  cout << "found " << no_indiv << " individuals in data" << endl;
}

void GenoData::read_bim() {
  string fname = prefix + ".bim", line, value;
  ifstream bim(fname.c_str(), ifstream::in);
  istringstream extract, convert;

  cout << "Reading " << fname << "... ";
  int line_no = 0, valid = 0, pos;
  while (getline(bim, line)) {
    line_no++; extract.clear(); extract.str(line);  
    for (int i = 0; i < 4; i++) {if (!(extract >> value)) error(string("not enough values on line ") + DataUtils::to_string(line_no));}

    convert.clear(); convert.str(value); convert >> pos;
    if (convert.eof() && !convert.fail() && pos > 0) {position.push_back(pos); valid++;}
    else position.push_back(0);
  }
  no_snps = position.size();

  pos_bounds.first = 0; pos_bounds.second = 0;
  for (int i = 0; i < no_snps && (pos_bounds.first == 0); i++) {if (position[i] > 0) pos_bounds.first = position[i];}
  for (int i = no_snps-1; i >= 0 && (pos_bounds.second == 0); i--) {if (position[i] > 0) pos_bounds.second = position[i];}  
  
  cout << "found " << valid << " SNPs (out of " << no_snps << ")" << endl;
}

void GenoData::prep_bed() {
  string fname = prefix + ".bed";
  cout << "Preparing file " << fname << "..." << endl;

  block_count = (unsigned long long) ceil(no_indiv/4.0);
  bed_file.open(fname.c_str(), ios::in|ios::binary|ios::ate);
  unsigned long bed_size = bed_file.tellg(), exp_bed_size = block_count * no_snps + 3; ///for SNP-major format

  char buffer[3];
  bed_file.seekg(0, ios::beg); bed_file.read(buffer, 3);
  if (bed_file.fail() || ((unsigned short) buffer[0] != 108 || (unsigned short) buffer[1] != 27)) error("file is not a valid .bed file");
  if ((unsigned short) buffer[2] != 1) {
    if (buffer[2] == 0) error("file is in individual-major format");
    else error("file-format specifier is not valid");    
  }
  if (bed_size != exp_bed_size) error("size of .bed file is inconsistent with number of SNPs and individuals in .bim and .fam files");

           
  unsigned char value_index[] = {1,0,2,3}; //hom1, miss, het, hom2
  for (int i = 0; i < 256; i++) {
    for (int j = 0; j < 4; j++) geno_index[i][j] = value_index[(i >> (2*j)) & 3];
  }
  
  raw_buffer.resize(block_count, 1);
  geno_buffer.resize(no_indiv, 1); 
}

pair<int,int> GenoData::load_block(float* target, vector<pair<int,int> >& pos_target, int offset, int total) {
  if (offset < 0 || offset >= no_snps) return pair<int,int>(0,0);
  
  char *raw = raw_buffer.get_data();
  bed_file.seekg(block_count*offset + 3);

  int no_loaded = 0, no_read = 0;
  for (int curr = offset; curr < no_snps; curr++) {
    bed_file.read(raw, block_count); no_read++;
    if (position[curr] > 0 && process_snp(raw, target)) {pos_target.push_back(pair<int,int>(position[curr],curr)); no_loaded++;}  
    if (no_loaded >= total) break;    
  }
  return pair<int,int>(no_loaded, no_read);
}

bool GenoData::process_snp(char* raw, float*& target) {
  char* geno = geno_buffer.get_data(); unsigned char* sub_buffer = 0;
  int counts[4] = {0,0,0,0}; int sub_i = 0; 

  for (int i = 0; i < no_indiv; i++) {
    if (i % 4 == 0) {sub_i = 0; sub_buffer = geno_index[(unsigned char) raw[i/4]];}
    geno[i] = sub_buffer[sub_i++];
    counts[geno[i]]++;
  }
  
  float sum = 0, sq = 0, nonzero = 0;
  for (int i = 1; i < 4; i++) {sum += (i-1)*counts[i]; sq += (i-1)*(i-1)*counts[i]; nonzero += (counts[i] != 0);}
  float mean = sum / no_indiv, freq = sum / (2*(no_indiv - counts[0]));

  sq += counts[0]*mean*mean;
  float sd = (sq - mean*sum) / (no_indiv-1); sd = sd > 0 ? sqrt(sd) : 0;

  if (nonzero < 2 || min(freq, 1-freq) < maf_thresh || sd <= 0) return false;
  
  float values[4] = {0};
  for (int i = 1; i < 4; i++) values[i] = ((i-1) - mean) / sd;
  for (int i = 0; i < no_indiv; i++) *(target++) = values[geno[i]];

  return true;
}

pair<int,int> GenoData::load_data(Buffer<float>& target, vector<pair<int,int> >& pos_target, int offset, int total) {
  if (target.nrow() != no_indiv || target.ncol() != total) target.resize(no_indiv, total);
  return load_block(target.get_data(), pos_target, offset, total);
}

