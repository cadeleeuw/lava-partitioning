/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "global.h"
#include "data.h"
#include "correlations.h"
#include "splitter.h"
#include "output.h"

int main(int argc, char* argv[]) {
  Settings settings(argc, argv);  
  Output out(settings.output_pref);

  GenoData data(settings.input_pref, settings.maf_thresh);
  cout << endl;

  cout << "Computing correlations..." << endl;
  cout << "\twindow = " << settings.snp_window << endl;
  cout << "\tMAF threshold = " << settings.maf_thresh << endl;

  Correlations corrs(settings.snp_window);
  corrs.compute(data);  
  cout << "\tretained " << corrs.get_size() << " SNPs after filtering" << endl; 
  cout << endl;
  CorrelationMatrix* cm = corrs.get_matrix(); //is deleted by Splitter

  if (settings.print_metric) {
    out.write_metrics(cm->get_metric());
    cout << endl;
  }

  Splitter analysis(settings);
  cout << "Computing break points..." << endl;
  cout << "\tminimum size = " << settings.split_size << endl;
  cout << "\tminimum proportion = " << settings.split_prop << endl;
  cout << "\tmetric margin = " << settings.metric_margin << endl;
  cout << "\tmetric maximum = " << min(settings.metric_max, 1.0) << endl;

  int breaks = analysis.run(cm);
  if (breaks <= 0) error("unable to find any break points with current settings");
  cout << "\tfound " << breaks << " break points" << endl;
  cout << endl;


  if (settings.refine) {
    Refiner refiner(data);
    cout << "Refining break points for unfiltered data..." << endl;
    refiner.refine(analysis, corrs);
    out.write(refiner.get_breaks(), refiner.get_positions(), data);
  } else out.write(analysis.get_breaks(), corrs.get_positions(), data);

  
  cout << endl;
  cout << "Analysis complete. Goodbye." << endl;
  cout << endl;

  return 0;
}
 
