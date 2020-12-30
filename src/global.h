/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <sstream> 
#include <sys/stat.h>

using namespace std;

static void error(const string& msg, int code=1) {
  cerr << endl << "ERROR: " << msg << endl;
  exit(code);  
}

class Settings {
  bool is_dir(const string& basename) {struct stat status; return stat(basename.c_str(), &status) == 0 && S_ISDIR(status.st_mode);}
  bool is_file(const string& filename) {
    struct stat status;
    int code = stat(filename.c_str(), &status);
    return code == 0 && S_ISREG(status.st_mode);
  }
  
  template<typename T>
  bool convert_num(const string& value, T& target) {
    istringstream istr(value); istr >> target;
    return !istr.fail() && istr.eof();
  } 


public:
  string input_pref, output_pref;
  double maf_thresh;
  int snp_window;
  
  int split_size; double split_prop;
  double metric_margin, metric_max;
  bool print_metric, refine;

  Settings(int argc, char* argv[]) : maf_thresh(0.01), snp_window(200), split_size(1000), split_prop(0.1), metric_margin(0.01), metric_max(0.25), print_metric(false), refine(true) {
    if (argc < 2) error("no arguments provided");
    if (is_dir(argv[1])) error("file prefix is a directory");

    string suffix[] = {".bed", ".bim", ".fam"};
    for (int i = 0; i < 3; i++) {
      string fname = string(argv[1]) + suffix[i];
      if (!is_file(fname)) error(string("file '") + fname + "' not found");                  
    }
    input_pref = argv[1];
    output_pref = "ldblock";
    
    for (int a = 2; a < argc; a++) {
      if (string(argv[a]) == "-frq") {
        if (argc <= a+1) error("no value specified for argument '-frq'");
        if (!convert_num(argv[++a], maf_thresh)) error("value for argument '-frq' is not a number");
        if (maf_thresh < 0 || maf_thresh > 0.40) error("value for argument '-frq' should be between 0 and 0.4");
      } else if (string(argv[a]) == "-win") {
        if (argc <= a+1) error("no value specified for argument '-win'");
        if (!convert_num(argv[++a], snp_window)) error("value for argument '-win' is not a (whole) number");
        if (snp_window < 1) error("value for argument '-win' should be at least 1");
      } else if (string(argv[a]) == "-min-size") {
        if (argc <= a+1) error("no value specified for argument '-min-size'");
        if (!convert_num(argv[++a], split_size)) error("value for argument '-min-size' is not a (whole) number");
        if (split_size < 50) error("value for argument '-min-size' should be at least 50");
      } else if (string(argv[a]) == "-split-prop") {
        if (argc <= a+1) error("no value specified for argument '-min-prop'");
        if (!convert_num(argv[++a], split_prop)) error("value for argument '-min-prop' is not a number");
        if (split_prop < 0 || split_prop > 0.4) error("value for argument '-min-prop' should be between 0 and 0.4");        
      } else if (string(argv[a]) == "-margin") {
        if (argc <= a+1) error("no value specified for argument '-margin'");
        if (!convert_num(argv[++a], metric_margin)) error("value for argument '-margin' is not a number");
        if (metric_margin < 0 || metric_margin > 0.10) error("value for argument '-margin' should be between 0 and 0.1");
      } else if (string(argv[a]) == "-max") {
        if (argc <= a+1) error("no value specified for argument '-max'");
        if (!convert_num(argv[++a], metric_max)) error("value for argument '-max' is not a number");
        if (metric_max <= 0) error("value for argument '-max' should be greater than 0");
      } else if (string(argv[a]) == "-out") {
        if (argc <= a+1) error("no value specified for argument '-out'");
        output_pref = argv[++a];
      } else if (string(argv[a]) == "-print-metric") {
        print_metric = true;
      } else if (string(argv[a]) == "-refine") {
        if (argc <= a+1) error("no value specified for argument '-refine'");
        string value = argv[++a];
        if (value == "1") refine = true;
        else if (value == "0") refine = false;
        else error("value for argument '-refine' should be either 0 or 1");
      } else error(string("unknown argument '") + argv[a] + "'");
    }
    if (maf_thresh == 0) refine = false;
  }
}; 

#endif /* GLOBAL_H */
