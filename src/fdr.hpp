#ifndef __fdr_hpp__
#define __fdr_hpp__

#include <vector>

namespace fdr {
  // Simple code for hanlding false discovery rate calculations.
  using namespace std;

  // Computes the score threshold for a desired FDR. 
  double threshold(double PIT, double fdr_needed, vector<double> &decoy_scores, vector<double> &target_scores);

  // Computes the q-value using the scores.
  double PIT_qvalue(double PIT, vector<double> &decoy_scores, vector<double> &target_scores, double score);

  // Computes the pi0 for the given decoy/null and target scores.
  double estimate_PIT(vector<double> &decoy_scores, vector<double> &target_scores);
}

#endif  
