#include <iostream>
#include <algorithm>
#include <limits>

#include <math.h>

#include "util.hpp"
#include "fdr.hpp"

namespace fdr {
  // O(n log n) algorithm, n = |target scores| + |decoy scores|
  // For each target score estimate a FDR.
  double threshold(double PIT, double fdr_needed, vector<double> &decoy_scores, vector<double> &target_scores) {
    double N_decoy = decoy_scores.size();
    double N_target = target_scores.size();

    // O(n log n) algorithm, n = |target scores| + |decoy scores|
    // For each target score estimate a FDR.
    size_t a = 0;
    while(a < target_scores.size()) {
      // Number of accepted target scors.
      double target_accept = target_scores.size() - a; // -1 is that right?

      // Find the number of decoy scores above the current target score.
      vector<double>::iterator it = lower_bound(decoy_scores.begin(), decoy_scores.end(), target_scores[a]);
      double decoy_accept = distance(it, decoy_scores.end());

      double fdr = PIT * (N_target / N_decoy) * decoy_accept  / target_accept;

      // NOTE: To compute q-value. Binary search in both target_scores
      // and decoy_scores to find current q-value score. Compute
      // target_accept and decoy_accept using distance()... multipy by
      // PIT correction factor.

      // If FDR threshold is achieved, return target score threshold.
      if(fdr < fdr_needed) {
	cout << "decoy_scores.size()=" << decoy_scores.size() << endl;
	cout << "target_scores.size()=" << target_scores.size() << endl;
	cout << "cutoff=" << target_scores[a] << endl;
	cout << "PIT=" << PIT << endl;
	cout << fdr_needed << ">=" << fdr << endl;
	cout << decoy_accept << "/" << target_accept << endl;
	
	return target_scores[a];
      }
      a++; // Otherwise, advance to next target score.
    }
    // If FDR is not acheived, accept no scores at all.
    return numeric_limits<double>::max();
  }

  // Compute a PIT adjusted q-value for a given score.
  double PIT_qvalue(double PIT, vector<double> &decoy_scores, vector<double> &target_scores, double score) {
    double N_decoy = decoy_scores.size(),  N_target = target_scores.size();

    // Find the number of decoy scores above the current target score.
    vector<double>::iterator it_decoy = lower_bound(decoy_scores.begin(), decoy_scores.end(), score);
    double decoy_accept = distance(it_decoy, decoy_scores.end());

    // Find the number of target scores above the current target score.
    vector<double>::iterator it_target = lower_bound(target_scores.begin(), target_scores.end(), score);
    double target_accept = distance(it_target, target_scores.end());

    // Compute q-value that adjusts for difference in number of target
    // scores and decoy scores.
    double qvalue = PIT * (N_target / N_decoy) * decoy_accept / target_accept;

    if(qvalue > 1) return 1; else return qvalue;
  }

  double estimate_PIT(vector<double> &decoy_scores, vector<double> &target_scores) {
    // Total number of deocy and target scores.
    double N_decoy = decoy_scores.size(),  N_target = target_scores.size();

    // Here we slowly increase p-value at which target spectra are marked as obviously incorrect.
    double pvalue_delta = 0.0005; // increased by pvalue_delta
    double PIT_pvalue = pvalue_delta; // start with a small p-value
    double PIT = 1; 
    if(decoy_scores.size() < 100 || target_scores.size() < 100) return PIT;

    vector<double> PITs;
    while(PIT_pvalue < 1) {
      // Get the score cutoff for the current p-value.
      size_t PIT_decoy = (size_t)floor(double(decoy_scores.size()) * (1.0 - PIT_pvalue));
      // Skip cases where p-value can't be computed.
      if(PIT_decoy == decoy_scores.size()) break;
      if(PIT_decoy == 0) break; 

      // Get the threshold of the current p-value.
      double PIT_threshold = decoy_scores.at(PIT_decoy);

      // Count the number of target spectra below this threshold.
      vector<double>::iterator it_target = upper_bound(target_scores.begin(), target_scores.end(), PIT_threshold);
      double PIT_target = distance(target_scores.begin(), it_target);

      // The PIT is the ratio of these incorrect target spectra over
      // the decoy spectra.  Note the correction for the ratio
      // between target and decoy scores!!!
      double PIT_cur = (PIT_target / double(PIT_decoy)) * (N_decoy / N_target);

      // cout << PIT_pvalue << "\t" << PIT_cur << "\t" << PIT_target << "\t" << PIT_decoy << endl;
      if(PIT_cur < 1) { PITs.push_back(PIT_cur); }
      else break; // Reached too high a p-value, estimate has large variance.

      PIT_pvalue += pvalue_delta;
    }
    // Keep the median PIT vlaue as the PIT estimate.  This assumes
    // the PIT values are pretty much flat for all of the p-values.
    if(PITs.size() > 0) PIT = util::median_unsafe(PITs);

    return PIT;
  }
}
