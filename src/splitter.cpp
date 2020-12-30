/** Copyright (C) 2021 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include <cmath>

#include "splitter.h"

void Splitter::clear() {
  for (int i = 0; i < blocks.size(); i++) delete blocks[i];
  blocks.clear(); break_points.clear();
}

Splitter::Splitter(Settings& settings) {
  min_size = settings.split_size > 0 ? settings.split_size : 1;
  min_prop = settings.split_prop;
  metric_margin = settings.metric_margin;
  metric_max = settings.metric_max;
}

void Splitter::insert_block(CorrelationMatrix* block) {
  if (!blocks.empty()) {
    for (deque<CorrelationMatrix*>::iterator curr = blocks.begin(); curr != blocks.end(); ++curr) {
      if (block->get_size() >= (*curr)->get_size()) {blocks.insert(curr, block); return;}
    }
  }
  blocks.push_back(block);
}

int Splitter::run(CorrelationMatrix* input) {
  clear(); blocks.push_front(input);
  
  while (!blocks.empty()) {
    CorrelationMatrix *cm1 = blocks.front(), *cm2; blocks.pop_front();
    int margin = max(int(ceil(cm1->get_size() * min_prop)), min_size), curr_size = cm1->get_size();
    
    if (curr_size >= margin*2) {
      Split curr; const vector<double>& metric = cm1->get_metric(); 
      
      if (curr_size > margin*2) {       
        for (int i = margin; i < curr_size - margin; i++) {
          if (metric[i] < curr.metric) curr.set(i, metric[i]);
        }
       
        if (curr.offset < 0) error("unknown failure when splitting block");

        if (curr.metric < metric_max && metric_margin > 0) {
          int offset = min(curr.offset, int(curr_size - curr.offset)) + 1, mid = curr_size / 2, end = curr_size - offset;
          double thresh = min(curr.metric + metric_margin, metric_max);
          for (int i = offset; i < end; i++) {
            if (metric[i] < thresh) {
              if (i < mid) {
                curr.set(i, metric[i]);
                end = curr_size - i;
              } else {
                if (i < (end - 1) || metric[i] < curr.metric) curr.set(i, metric[i]);
                break;
              }
            }
          }
        }
      } else curr.set(margin, metric[margin]);
      
      if (curr.metric < metric_max) {
        cm2 = cm1->split(curr.offset);
        curr.offset += cm1->get_offset();
        break_points.push_back(curr);
  
        if (cm2->get_size() > cm1->get_size()) swap(cm1, cm2);
        insert_block(cm1); insert_block(cm2);
      } else delete cm1;
    } else {delete cm1; break;}
  }  
  
  return break_points.size();
}
 
  
void Refiner::refine(Splitter& analysis, Correlations& corrs) {
  breaks = analysis.get_breaks(); positions = corrs.get_positions();
  int depth = corrs.get_depth();

  for (int b = 0; b < breaks.size(); b++) {
    Split& curr = breaks[b];
    int lower_bound = positions[curr.offset].second, upper_bound = positions[curr.offset+1].second;
      
    if (upper_bound > (lower_bound+1)) {
      int low = max(lower_bound - depth,0), high = upper_bound + depth;

      corrs.compute_block(data, low, high); 
      const vector<pair<int,int> >& curr_pos = corrs.get_positions();

      int i_lower = -1, i_upper = -1;
      for (int i = 0; i < curr_pos.size(); i++) {
        if (curr_pos[i].second == lower_bound) i_lower = i;
        if (curr_pos[i].second == upper_bound) i_upper = i;
      }
      if (i_lower >= 0 && i_upper > i_lower) {
        CorrelationMatrix* cm = corrs.get_matrix();
        const vector<double>& metric = cm->get_metric();
        
        int i_min = i_lower; double min_value = metric[i_lower];
        for (int i = i_lower+1; i < i_upper; i++) {
          if (metric[i] < min_value) {i_min = i; min_value = metric[i];}
        }

        if (i_min >= 0) curr.position = (curr_pos[i_min].first + curr_pos[i_min+1].first) / 2.0;
        delete cm;
      }
    }
  }
}
  