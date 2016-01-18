#include <iostream> 
#include <math.h>
#include <fstream>
#include <set>

#include "lcms.hpp"
#include "util.hpp"
#include "xml.hpp"

// MSVC++ is still moving to C99. Fuck you Microsoft.
#ifdef _MSC_VER
#include "win_stdint.h"
#else
#include <stdint.h>
#endif

// TODO: Implement the FFT algorithm for calculating isotope
// distributions from:
// Rockwood, A. L., Van Orden, S. L., and Smith, R. D. (1995). Rapid
// Calculation of Isotope Distributions. Anal. Chem. 67:15, 2699-2704.

#include "fft.hpp"

namespace lcms {
  using namespace std;

  void irun::recalibrate(recalmode_t recalmode, double winSz, vector<ms2::fragmentDB *> &DB) {
    if(DB.size() == 0) return; // No DBs, quit.
    vector<xml::recal_t> recal_pts;
    // Get recalibration points.
    for(size_t i = 0; i < xics.size(); i++) {
      // Get MS/MS spectrum with the highest score.
      ms2peak *peak = xics[i].max_ms2();
      if(peak != NULL && peak->num_ids() > 0) { 
	// Get the corresponding database .
	ms2::fragmentDB *db = DB[dbLIGHT];
	if(peak->iso == ms2::isoHEAVY) db = DB[dbHEAVY]; 
	else if(peak->iso == ms2::isoMEDIUM) db = DB[dbMEDIUM];
	if(db == NULL) continue; // No database, quit.
	// Get mass shift.
	vector<ms2::fragment *> ids;  peak->max_ids(ids);
	vector<ms2::mods_t> mods;     peak->max_mods(mods);
	// Compute precursor mass using modification and ID.
	double seq_mass = db->precursor_mass(ids[0], mods[0]);
	// Compute the m/z shift. 
	double delta_mz = ((seq_mass + peak->charge * mass::proton) - peak->mass()) / peak->charge;
	// Save recalibration value.
	xml::recal_t r;
	r.time = peak->retentionTime;
	r.delta_mz = delta_mz;
	r.mz = peak->mz; 
	recal_pts.push_back(r);
      }
    }
    cout << "recal_pts.size()=" << recal_pts.size() << endl;
    if(recal_pts.size() > 0) {  // If recalibration points availible run recalibration.
      string newfilename = filename + "recal";
      if(recalmode == recalmodeTIME) newfilename += "_time";
      if(recalmode == recalmodeMZ) newfilename += "_mz";
      xml::recalibrate(recalmode, winSz, recal_pts, filename, newfilename);
    }
  }

  // Functors for sorting XIC pointers.
  struct mz_sort_incr { bool operator() (xic *a, xic *b) const {  return a->mz < b->mz; } };
  struct mz_sort_decr { bool operator() (xic *a, xic *b) const {  return a->mz < b->mz; } };
  struct quant_sort_decr { bool operator() (xic *a, xic *b) const {  return a->quant > b->quant; } };

  void irun::get_range(float &startTime, float &endTime, float retentionTime) {
    vector<float>::iterator it = 
      lower_bound(retentionTimes.begin(), retentionTimes.end(), retentionTime);
    int rt_idx = it - retentionTimes.begin();
    int start_idx = rt_idx - 2 - conf.xic_peaks_missed; 
    int end_idx = rt_idx + 2 + conf.xic_peaks_missed;
    startTime = endTime = 0;
    if(start_idx < 0) startTime = retentionTimes[0] - 1.0f;
    else startTime = 0.5 * (retentionTimes[start_idx] + retentionTimes[start_idx+1]); // 3, 2

    if(end_idx >= (int)retentionTimes.size()) endTime = retentionTimes[retentionTimes.size()-1] + 1.0f;
    else endTime = 0.5 * (retentionTimes[end_idx] + retentionTimes[end_idx-1]);
  }

  // Remove noise and contaminant peaks from centroided peak data.  If
  // blanks are availible, use these to filter peaks as well.
  void irun::filter_new(double log10I_thresh){
    if(retentionTimes.size() <= 3) { cerr << "retentionTimes.size() <= 3" << endl; return; }
    size_t keep_cnt = 0;
    vector<bool> keep(data.size(), false);
    vector<peak *> results;

    for(size_t i = 0; i < data.size(); i++){
      if(log10(data[i].intensity) < log10I_thresh) continue;
      
      float startTime, endTime;
      get_range(startTime, endTime, data[i].retentionTime);

      // Query a rectangular window around data point.
      results.clear();
      queryraw(results, // retention time winow
	       startTime,  
	       endTime,
	       data[i].mz - mass::ppm2da(data[i].mz, conf.dmz) / 2, // m/z window  // PPM THRESHOLD
	       data[i].mz + mass::ppm2da(data[i].mz, conf.dmz) / 2); 

      // If the window returns enough peaks, keep the point. 
      if((int)(results.size()) > 1) { keep[i] = true; keep_cnt++; }
    }

    // Keep peaks labeled with "keep" flag.
    filtered.reserve(keep_cnt);
    for(size_t i = 0; i < data.size(); i++) { if(keep[i]) { filtered.push_back(data[i]); } }

    // Keep all MS/MS peaks.
    filtered2 = data2;
  }


  // Functor is used to sort along the retention time coordinate.
  struct retentionTimecoord{ bool operator () (peak a, peak b){ return a.retentionTime < b.retentionTime;} };

  // Gets rid of duplicates peaks (two peaks with the same retention
  // time). Replaces a duplicate set with peak that has the maximum
  // intensity.  Note cleaned chromatogram is in sorted order.
  void irun::clean_chromatogram(xic &s) {
    if(s.chrom.size() <= 1) return;

    // Copy original peaks and sort along retention time coordinate.
    vector<peak> orig = s.chrom;
    sort(orig.begin(), orig.end(), retentionTimecoord());

    s.chrom.clear(); // Clear the old chromatogram.

    size_t i = 0;
    // Iterate through each sorted peak.
    while(i < orig.size()) {
      // Find a run that has a delta retention time that is small.
      size_t j = i;  
      while(j < orig.size() - 1 && fabs(orig[j+1].retentionTime - orig[j].retentionTime) < 0.00001) j++;

      // Find the maximum intensity peak in that run.
      double Imax = -10; peak pmax;
      for(size_t k = i; k <= j; k++) {
	if(orig[k].intensity > Imax) { Imax = orig[k].intensity; pmax = orig[k]; }
      }

      // Keep only the maximumn intensity peak.
      if(Imax > 0) s.chrom.push_back(pmax);
      i = j + 1; // Advance to next set of peaks. 
    }
  }

  // Smooth chromatogram by running average. 
  void irun::smooth_chromatogram(xic &s) {
    if(s.chrom.size() <= 1) return;
    vector<peak> &chrom = s.chrom;
    vector<peak> smoothed(chrom.size());
    // 5 point smoothing
    const int W = 2; // Use 2 points before, current points, and 2 points ahead.
    for(int i = 0; i < (int)chrom.size(); i++) {
      int N = 0; 
      float sumI = 0;
      for(int w = i - W; w <= i + W; w++) {
	if(w >= 0 && w < (int)chrom.size()) { sumI += chrom[w].intensity; }
	N++;	
      }
      smoothed[i] = chrom[i];
      smoothed[i].intensity = sumI / float(N);
    }
    chrom = smoothed;
  }

  // Computes the area under a XIC using trapezoidal integration.
  void irun::trapezoidal_integration(xic &s){ 
    // Perform trapezoidal integration 
    s.quant = 0;
    for(unsigned int i = 0; i < s.chrom.size() - 1; i++){
      double timeStep = s.chrom[i+1].retentionTime - s.chrom[i].retentionTime;
      s.quant += ( s.chrom[i+1].intensity + s.chrom[i].intensity ) * 0.5 * timeStep;
    }
  }

  void irun::out_edges_new( peak* p,  vector< peak* > &edges ) {
    float startTime, endTime;
    get_range(startTime, endTime, p->x());
    double dyPPM = conf.peak_dmz;
    region<peak> query( startTime, endTime, 
			p->y() - (p->y() * dyPPM * 1e-6) / 2, 
			p->y() + (p->y() * dyPPM * 1e-6) / 2 );
    root->Query(query, edges);
  }


  // DFS scheme for labeling connected components to form XICs.
  void irun::dfs_label(int label, peak *start, vector<peak *> &L) {
    L.clear();
    vector<peak*> edges;
    L.push_back(start);
    while(L.empty() == false) {
      // Get current node from stack.
      peak* cur = L.back(); L.pop_back();
      cur->label = label; // Label node.
      // Get edges by range query.
      edges.clear(); 
      //out_edges(cur, edges);
      out_edges_new(cur, edges);
      for(size_t j = 0; j < edges.size(); j++) {
	// Nodes that haven't been labeled, add them to the stack.
	if(edges[j]->label == 0) {
	  edges[j]->label = 1;
	  L.push_back(edges[j]);
	}
      }
    }
  }

  // Locate XICs in the LC-MS data. XICs are peaks that occur
  // at about the same m/z over a "long" period of time.
  void irun::findxics(){
    // Create graph on filtered peaks. Find connected components in this graph.
    size_t num_components = 0;
    for(size_t i = 0;  i < filtered.size(); i++) filtered[i].label = 0; // Initialize labels.
    vector<peak *> L; L.reserve(1024);
    for(size_t i = 0;  i < filtered.size(); i++) {
      if(filtered[i].label != 0) continue; // Already visited vertex, skip.
      // Apply DFS to label component.
      dfs_label(num_components + 2, &filtered[i], L);
      num_components++; // Advance to next component label.
    }

    // We know the number of connected components so we know the number of XICs.
    vector<xic> raw_xics(num_components);

    // Read off connected components to create XICs.
    for(size_t i = 0; i < filtered.size();  i++){
      // Add peak to XIC which is the same as a connected
      // component in the peak graph.
      peak p;
      p.mz = filtered[i].mz; p.intensity = filtered[i].intensity;
      p.retentionTime = filtered[i].retentionTime;
      raw_xics[filtered[i].label-2].chrom.push_back(p);
    }

    // Keep only those XICs that are long enough and not too long.
    for(unsigned int i = 0; i < raw_xics.size(); i++){
      if(raw_xics[i].chrom.size() == 0) continue;

      // Sort XICs peaks by retention time. 
      sort(raw_xics[i].chrom.begin(), raw_xics[i].chrom.end());

      // If there are multiple peaks at same retention time value, picks only the most intense peak.
      clean_chromatogram(raw_xics[i]);
      smooth_chromatogram(raw_xics[i]);

      // TODO: Median and average filter chromatograms.
      //       1. Sputtering (eliminate by smoothing).
      //       2. Tailing?
      //       3. merged peaks (how?)
      // TODO: Split chromatograms here (to get isomers).

      double start = raw_xics[i].chrom[0].retentionTime;
      double end = raw_xics[i].chrom[raw_xics[i].chrom.size() - 1].retentionTime;

      // Not too short and not too long.
      if((int)raw_xics[i].chrom.size() >= conf.min_xic_peaks) {
	if(end - start < conf.max_xiclen) xics.push_back(raw_xics[i]);
      }

    }

    // Free up memory used by raw XICs.
    raw_xics.clear(); vector<xic> (raw_xics).swap(raw_xics);

    // Compute some basic statistics on these XICs. 
    for(unsigned int i = 0; i < xics.size(); i++){
      xic &s = xics[i];

      // If there are multiple peaks at same retention time value, picks only the most intense peak.
      clean_chromatogram(s);

      // Record start and end time of the XIC. 
      s.start = s.chrom[0].retentionTime;
      s.end = s.chrom[s.chrom.size() - 1].retentionTime;

      // m/z assigned to a XIC is the intensity weighted m/z of the
      // peaks in the XIC
      double maxIrt = 0, maxI = 0, sumI = 0;
      for(unsigned int j = 0; j < s.chrom.size(); j++){
	sumI += s.chrom[j].intensity; // Sum intensities.
	// Retention time is position of maximum 
	// intensity peak.
	if(s.chrom[j].intensity > maxI) {
	  maxI = s.chrom[j].intensity;
	  maxIrt = s.chrom[j].retentionTime;
	}
      }
      s.retentionTime = maxIrt;  // Center at maximum intensity peak in chromatogram.
      trapezoidal_integration(s); // Compute quantiative data.

      // Assign XIC the intensity weighted m/z.
      s.mz = 0;
      for(unsigned int j = 0; j < s.chrom.size(); j++){ s.mz += (s.chrom[j].intensity / sumI) * s.chrom[j].mz; }

      // NOTE: critical to set replicate and condition index so XIC
      // can be cross referenced back in cross irun grouping.
      s.repIndex = repIndex; s.condIndex = condIndex;  s.irunIndex = irunIndex;
      // Fill default info for XIC.
      s.type = xicUNKNOWN;  s.charge = 0; s.isotope_grouped = false;
    }
    // Build an index on XICs. 
    xics_all = new _2dtree<xic>(xics);
  }


  // Assign MS/MS spectra to monoisotopic XICs.  Use monoisotopic XIC
  // to update precursor m/z and to assign charge.  NOTE: This code
  // assures that an MS/MS spectrum is assigned to one and only one
  // XIC. This XIC is the one that has an MS1 peak that is closest to
  // the precursor m/z and retention time of the MS/MS spectrum.
  struct ms2xic_d {
    ms2peak *peak;
    xic *assigned;
    double distsq;
    double x() const { return peak->retentionTime; } double y() const { return peak->mz; }
  };
  void irun::assign2xic() {
    // Build a 2d-tree on MS/MS peaks that have been augmented with an
    // assignment and a distance.
    vector<ms2xic_d> spectra(filtered2.size());
    for(size_t i = 0; i < filtered2.size(); i++) {
      spectra[i].peak = &filtered2[i];
      spectra[i].assigned = NULL;
      spectra[i].distsq = numeric_limits<double>::max();
    }
    _2dtree_inplace<ms2xic_d> tree2(spectra);

    vector<ms2xic_d*> res;
    for(size_t x = 0; x < xics.size(); x++) {
      vector<peak> &chrom = xics[x].chrom; 

      // Use each peak in the XIC to query for MS/MS spectra in range
      // query region used for XIC construction.
      for(size_t c = 0; c < chrom.size(); c++) {
	double rt = chrom[c].retentionTime, mz = chrom[c].mz;
	/*region<ms2xic_d> query(rt - conf.peak_dtime / 2, 
			       rt + conf.peak_dtime / 2,
			       mz - mass::ppm2da(mz, conf.peak_dmz) / 2, 
			       mz + mass::ppm2da(mz, conf.peak_dmz) / 2);*/
	float startTime, endTime;
	get_range(startTime, endTime, rt);
	region<ms2xic_d> query(startTime, 
			       endTime,
			       mz - mass::ppm2da(mz, conf.peak_dmz) / 2, 
			       mz + mass::ppm2da(mz, conf.peak_dmz) / 2);
	res.clear(); tree2.Query(query, res);

	// Examine each returned spectrum.
	for(size_t r = 0; r < res.size(); r++) {
	  double drt = res[r]->peak->retentionTime - rt;
	  double dmz = res[r]->peak->mz - mz;
	  double distsq_cur = drt * drt + dmz + dmz;
	  // Assign XIC to spectrum if the XIC's peak is closer to this spectrum.
	  if(distsq_cur < res[r]->distsq) {
	    res[r]->assigned = &xics[x]; res[r]->distsq = distsq_cur;
	  }
	}
      }        
    }
    // Now each MS/MS spectrum is assigned to the XIC with a peak that is closest in retention time and m/z.
    for(size_t s = 0; s < spectra.size(); s++) {
      if(spectra[s].assigned == NULL) continue; 
      ms2peak *updated_peak = spectra[s].peak;
      xic *x = spectra[s].assigned;
      // MS/MS peak is assigned to an isotope XIC.
      if(x->type == xicMONOISO) {
	// NOTE: The intensity weighted m/z value in mz_update is a more accurate precursor m/z.
	updated_peak->mz = x->mz;
	if(updated_peak->charge == 0) updated_peak->charge = x->charge;
	x->ms2peaks.push_back(*updated_peak);
      } 
      else if(x->type == xicUNKNOWN || x->type == xicISOTOPE) {
	// Trust the instrument. It uses a real isotope distribution model to handle these cases. 
	// TODO: Add an isotope distribution calculator to PVIEW.
	x->type = xicMONOISO;
	x->charge = updated_peak->charge;
	updated_peak->mz = x->mz;
	x->ms2peaks.push_back(*updated_peak);
      }
    }
    // Rebuild 2d-tree and ajusted MS/MS peaks.
    delete root2; root2 = new _2dtree<ms2peak>(filtered2); 
  }

  // Apply filtering to XICs. 
  void irun::filter_xics() {
    if(conf.quantmode != quantALIGN) return;
    vector<xic> xics2;
    // Filtering eliminates XICs that are close and are contained
    // within another XIC.  These are caused by dominant peaks in the
    // data.
    for(size_t s = 0; s < xics.size(); s++) {
      vector<xic*> results;

      double length = xics[s].end - xics[s].start;

      queryxics(results, xics[s].start - 2 * length, xics[s].end + 2 * length, 
		xics[s].mz - 2 * mass::ppm2da(xics[s].mz, conf.xic_width), 
		xics[s].mz + 2 * mass::ppm2da(xics[s].mz, conf.xic_width)); 

      // NOTE: This is twice what is used for grouping across
      // iruns!!!!  NOTE change to xic_width (removes one paramter).
      // TODO: Make this fatness window a paramter. Keep relative to a
      // fixed xic_width!!!

      bool filter = false;
      // Did you get a XIC that fully overlaps/contains the current
      // XIC? Then this XIC is a suprious XIC.
      for(unsigned int i = 0; i < results.size(); i++){
	xic *s2 = results[i];
	if( s2 == &xics[s]) continue; // Skip current XIC.
	// If yes, filter the current XIC.
	if( s2->start <= xics[s].start && xics[s].end <= s2->end  ) {
	  filter = true;
	}
      }
      if(filter == false) { xics2.push_back(xics[s]);  }
    }

    // Save and re-index these remaining XICs.
    xics = xics2;
    delete xics_all;
    xics_all = new _2dtree<xic>(xics);
  }

  // Move 2d-trees, adjust positions of XICs and peaks.
  void irun::translate(double translation) {
    cout << description << ": delta=" << translation << endl;
    // Translate raw data.
    for(size_t p = 0; p < data.size(); p++) data[p].retentionTime += translation;
    for(size_t p = 0; p < data2.size(); p++) data2[p].retentionTime += translation;

    // Now the filtered peaks.
    for(size_t p = 0; p < filtered.size(); p++) filtered[p].retentionTime += translation;
    for(size_t p = 0; p < filtered2.size(); p++) filtered2[p].retentionTime += translation;

    // Translate XICs, adjusting start and end retention times as
    // well.  
    for(size_t s = 0; s < xics.size(); s++) {
      xics[s].retentionTime += translation;   
      xics[s].start += translation;  
      xics[s].end += translation;
      for(size_t i = 0; i < xics[s].ms2peaks.size(); i++) xics[s].ms2peaks[i].retentionTime += translation;
      for(size_t i = 0; i < xics[s].chrom.size(); i++) xics[s].chrom[i].retentionTime += translation;
    }
    data_range(); // Update data range
  }

  inline double _nonlinear_shift(vector<shift_t> &shifts, double x) {
    if(shifts.size() == 0) return 0;
    shift_t n(x);
    vector<shift_t>::iterator it = lower_bound(shifts.begin(), shifts.end(), n);
    double dx = 0;
    // Use nearest neighbor interpolation to pick the shift for the
    // given x value.
    // TODO: Replace this with linear or cubic interpolation!!!
    if(it == shifts.begin()) dx = it->dx;
    else if(it == shifts.end()) dx = shifts[shifts.size()-1].dx;
    else {
      vector<shift_t>::iterator it2 = it - 1;
      if( it->x - x < x - it2->x ) dx = it->dx;
      else dx = it2->dx;
    }
    return dx;
  }

  // Applies non-linear transformation on the LC-MS/MS irun data.
  void irun::nonlinear(vector<shift_t> &shifts) {
    // Translate raw data.
    for(size_t p = 0; p < data.size(); p++) data[p].retentionTime += _nonlinear_shift(shifts, data[p].retentionTime);
    for(size_t p = 0; p < data2.size(); p++) data2[p].retentionTime += _nonlinear_shift(shifts, data2[p].retentionTime);

    // Now the filtered peaks.
    for(size_t p = 0; p < filtered.size(); p++) filtered[p].retentionTime += _nonlinear_shift(shifts, filtered[p].retentionTime);
    for(size_t p = 0; p < filtered2.size(); p++) filtered2[p].retentionTime += _nonlinear_shift(shifts, filtered2[p].retentionTime);

    // NOTE: Do not recompute quantification value for XIC. Nonlinear
    // shift may introduce noise.  TODO: Figure out if this is the
    // case. It might actually improve things!!!!

    // Translate XICs, adjusting start and end retention times as well.  
    for(size_t s = 0; s < xics.size(); s++) {
      xics[s].retentionTime += _nonlinear_shift(shifts,xics[s].retentionTime);
      xics[s].start += _nonlinear_shift(shifts, xics[s].start);
      xics[s].end += _nonlinear_shift(shifts, xics[s].end);

      for(size_t i = 0; i < xics[s].ms2peaks.size(); i++) 
	xics[s].ms2peaks[i].retentionTime += _nonlinear_shift(shifts, xics[s].ms2peaks[i].retentionTime);

      for(size_t i = 0; i < xics[s].chrom.size(); i++) {
	xics[s].chrom[i].retentionTime += _nonlinear_shift(shifts, xics[s].chrom[i].retentionTime);
      }
    }

    // Rebuild indexes since points transformed nonlinearly.
    if(rawroot != NULL) { delete rawroot; rawroot = new _2dtree_inplace<peak>(data); }
    if(rawroot2 != NULL) { delete rawroot2; rawroot2 = new _2dtree<ms2peak>(data2); }
    if(root != NULL) { delete root; root = new _2dtree_inplace<peak>(filtered); }
    if(root2 != NULL) { delete root2; root2 = new _2dtree<ms2peak>(filtered2); }
    if(xics_all != NULL) {delete xics_all; xics_all = new _2dtree<xic>(xics); }

    data_range(); // Update data range
  }


  // Use conf.isotolerence window to do a range query to get an
  // overlapped XIC. If multiple XICs are found, return only the most
  // intense XIC.
  xic* irun::get_overlapped(xic *x, double isotopic_spacing, double tol) {
    // Find XICS in the range [mz + dmz - tol, mz + dmz + tol ). 
    vector<xic *> res, res_check;
    double mz_query = x->mz + isotopic_spacing;
    // Use start and end of the XIC, use m/z window. 
    region<xic> query(x->start, x->end, 
		      mz_query - mass::ppm2da(mz_query, tol),  
		      mz_query + mass::ppm2da(mz_query, tol));
    xics_all->Query(query, res); 
    // Find most intense XIC in queried region. Make sure it's not the query XIC itself.
    for(size_t r = 0; r < res.size(); r++) if(res[r] != x) res_check.push_back(res[r]);
    if(res_check.size() == 0) return NULL;
    xic *most_intense = res_check[0];
    for(size_t i = 1; i < res_check.size(); i++) { 
      if(res_check[i]->quant > most_intense->quant) most_intense = res_check[i]; 
    }
    return most_intense;
  }

  xic *irun::find_isotope(xic *cur, double mass_shift, bool ignore_type) {
    if(cur->type != xicMONOISO) return NULL; // not monoisotopic.
    if(!(cur->charge > 0)) return NULL; // No charge  information.
    if(cur->isotope_grouped) return NULL;  // Grouped already

    // Compute actual isotope shift based on XIC charge.
    double shift = mass_shift / double(cur->charge);
    xic *cand = get_overlapped(cur, shift, conf.isotolerence);     
    if(cand == NULL) return NULL; // no candidate
    if(cand->isotope_grouped) return NULL; // candidate XIC is already grouped
    xic *recip = get_overlapped(cand, -shift, conf.isotolerence);
    if(recip == NULL) return NULL; // recipricol query returned nothing
    if(cur != recip) return NULL; // recipricol is not the current XIC

    // Ideal case, matching charge and candidate is moniso.
    if(ignore_type) return cand;
    else { if(cur->charge == cand->charge && cand->type == xicMONOISO) return cand; }
    return NULL;
  }


  void irun::add_triple(xic *light, xic *medium, xic *heavy, vector<int> &labels) {
    isotope_group t; 
    t.light = light; t.medium = medium; t.heavy = heavy;
    // Mark XICs as grouped and set the label in the fragmentation
    // spectrum according to triple membership.
    if(light != NULL) {
      for(size_t p = 0; p < light->ms2peaks.size(); p++) light->ms2peaks[p].iso = ms2::isoLIGHT;
      t.light->isotope_grouped = true;
    }
    if(medium != NULL) {
      for(size_t p = 0; p < medium->ms2peaks.size(); p++) medium->ms2peaks[p].iso = ms2::isoMEDIUM;
      t.medium->isotope_grouped = true;
    }
    if(heavy != NULL) {
      for(size_t p = 0; p < heavy->ms2peaks.size(); p++) heavy->ms2peaks[p].iso = ms2::isoHEAVY;
      t.heavy->isotope_grouped = true;
    }
    t.labels = labels;

    // Compute all possible log2 ratios.
    if(t.medium != NULL && t.heavy != NULL) 
      //t.log2ratioHM = util::Log2(t.heavy->quant) - util::Log2(t.medium->quant);
      t.log2ratioHM = compute_log2ratio(t.heavy, t.medium);
    if(t.heavy != NULL && t.light != NULL) 
      //t.log2ratioHL = util::Log2(t.heavy->quant) - util::Log2(t.light->quant);
      t.log2ratioHL = compute_log2ratio(t.heavy, t.light);
    if(t.medium != NULL && t.light != NULL)
      //t.log2ratioML = util::Log2(t.medium->quant) - util::Log2(t.light->quant);
      t.log2ratioML = compute_log2ratio(t.medium, t.light);;

    if(t.heavy != NULL)  t.log2xicH = util::Log2(t.heavy->quant);
    if(t.medium != NULL) t.log2xicM = util::Log2(t.medium->quant);
    if(t.light != NULL)  t.log2xicL = util::Log2(t.light->quant);

    // Make sure condition, run, and replicate indicies are stored.
    t.condIndex = condIndex; t.repIndex = repIndex; t.irunIndex = irunIndex;
    t.id = 0;
    isotope_data.push_back(t);
  }

  struct iratio_t {
    float sumI, heavyI, lightI;
    bool operator < (const iratio_t &b) const { return sumI < b.sumI; } 
  };
  double irun::compute_log2ratio(xic *heavy_xic, xic *light_xic) {
    vector<peak> &light = light_xic->chrom, &heavy = heavy_xic->chrom;
    sort(light.begin(), light.end());  sort(heavy.begin(), heavy.end());

    /*vector<double> log2HL_ratios;
    for(size_t l = 0; l < light.size(); l++) {
      for(size_t h = 0; h < heavy.size(); h++) {
	if(fabs(light[l].retentionTime - heavy[h].retentionTime) < 0.00001) {
	  log2HL_ratios.push_back(util::Log2(heavy[h].intensity) - util::Log2(light[l].intensity));
	  break;
	}
      }
    }
    if(log2HL_ratios.size() > 0) return util::median_unsafe(log2HL_ratios);
    else return util::Log2(heavy_xic->quant) - util::Log2(light_xic->quant);*/

    vector<iratio_t> iratios;
    for(size_t l = 0; l < light.size(); l++) {
      for(size_t h = 0; h < heavy.size(); h++) {
	if(fabs(light[l].retentionTime - heavy[h].retentionTime) < 0.00001) {
	  iratio_t iratio;
	  iratio.sumI = light[l].intensity + heavy[h].intensity;
	  iratio.heavyI = heavy[h].intensity;
	  iratio.lightI = light[l].intensity;
	  iratios.push_back(iratio);
	  break;
	}
      }
    }
    if(iratios.size() > 0)  {
      sort(iratios.begin(), iratios.end());
      size_t istart = 0.2 * iratios.size();
      double Z = 0;
      for(size_t i = istart; i < iratios.size(); i++) Z += iratios[i].sumI;
      double log2_mu = 0;
      for(size_t i = istart; i < iratios.size(); i++) {
	log2_mu += (iratios[i].sumI / Z) * (util::Log2(iratios[i].heavyI) - util::Log2(iratios[i].lightI));
      }
      return log2_mu;
    }
    else return util::Log2(heavy_xic->quant) - util::Log2(light_xic->quant);
  }

  // Adds an isotope group, computing and storing any relevant information.
  void irun::add_pair(xic *light, xic *heavy, vector<int> &labels, int nitrogens) {
    isotope_group p; p.light = light; p.heavy = heavy;
    if(labels_reversed) { // **Handle reversed labels**
      //p.log2ratioHL = util::Log2(p.light->quant) - util::Log2(p.heavy->quant);
      p.log2ratioHL = compute_log2ratio(p.light, p.heavy);
      p.log2ratioHL = util::Log2(p.light->quant) - util::Log2(p.heavy->quant);
      p.log2xicH = util::Log2(p.light->quant);
      p.log2xicL = util::Log2(p.heavy->quant);
    }
    else {
      //p.log2ratioHL = util::Log2(p.heavy->quant) - util::Log2(p.light->quant);
      p.log2ratioHL = compute_log2ratio(p.heavy, p.light);
      p.log2xicH = util::Log2(p.heavy->quant); 
      p.log2xicL = util::Log2(p.light->quant);
    }
    for(size_t p = 0; p < light->ms2peaks.size(); p++) light->ms2peaks[p].iso = ms2::isoLIGHT;
    for(size_t p = 0; p < heavy->ms2peaks.size(); p++) heavy->ms2peaks[p].iso = ms2::isoHEAVY;

    // Mark XICs as grouped.
    p.light->isotope_grouped = true;  p.heavy->isotope_grouped = true;
    p.labels = labels; // Save labels used to create this isotope group.
    p.nitrogens = nitrogens;
    p.condIndex = condIndex; p.repIndex = repIndex; p.irunIndex = irunIndex;
    p.id = 0;
    isotope_data.push_back(p);
  }
  
  void irun::add_group(vector<xic *> &xics) {
    isotope_group p; 
    p.xics = xics;
    for(size_t x = 0; x < p.xics.size(); x++) { if(p.xics[x] != NULL) p.xics[x]->isotope_grouped = true;}
    // Mark XICs as grouped.
    p.condIndex = condIndex; p.repIndex = repIndex; p.irunIndex = irunIndex;
    p.id = 0;
    isotope_data.push_back(p);
  }

  void irun::add_group_phase(vector<xic *> &xics, int phase) {
    isotope_group p; 
    p.xics = xics;
    for(size_t x = 0; x < p.xics.size(); x++) { if(p.xics[x] != NULL) p.xics[x]->isotope_grouped = true;}
    // Mark XICs as grouped.
    p.condIndex = condIndex; p.repIndex = repIndex; p.irunIndex = irunIndex;
    p.id = 0;
    p.grouping_phase = phase;
    isotope_data.push_back(p);
  }

  // Medium and light mass shifts. 
  struct shift_label {  
    double mass_shift, mass_shift2;
    vector<int> labels; // Vector integers that has indicies of members of lcms::iso_labels.
    bool operator < (const shift_label &b) const { return mass_shift < b.mass_shift; }
  };

  // TODO: Too much repeated code here.
  void irun::find_isotope_pairs_and_triples() {
    cout << "pairs and triples" << endl;
    for(size_t i = 0; i < xics.size(); i++) xics[i].isotope_grouped = false;      // Flag all XICs as ungrouped. 

    // TODO: Change this so XICs with fragmentation spectra are
    // processed first.
    vector<xic*> xics_mz;     // Get pointers to XICs that are monoisotopic.
    for(size_t i = 0; i < xics.size(); i++)  // NOTE: These *should* have charge info.
      if(xics[i].type == xicMONOISO) xics_mz.push_back(&xics[i]);

    sort(xics_mz.begin(), xics_mz.end(), mz_sort_incr());      // Sort monoisotopic XICs on m/z dimension.

    // Generate all multicombinations of active isotope labels.
    vector< vector<int> > mcs;  // NOTE: The zeroth label is always the no-label label.
    util::multi_comb gen(iso_labels.size(), conf.max_label_count, mcs);

    // This tries all of the combinations e.g. for SILAC, KR, KK, RR.
    vector<shift_label> shift_labels;
    for(size_t m = 0; m < mcs.size(); m++) {
      vector<int> &mc = mcs[m]; // Get multicombination of labels.
      shift_label s; 
      if(conf.quantmode == quantISOPAIR) {
	double mass_shift = 0;  // Compute isotope shift for the multicombination of labels.
	for(size_t i = 0; i < mc.size(); i++) mass_shift += iso_labels[mc[i]].shift_heavy;
	if(mass_shift == 0) continue; // Skip no-shift case.
	s.mass_shift = mass_shift; 
	s.mass_shift2 = 0;
      }
      else {
	double shift_medium = 0;  
	for(size_t i = 0; i < mc.size(); i++) shift_medium += iso_labels[mc[i]].shift_medium;
	double shift_heavy = 0;
	for(size_t i = 0; i < mc.size(); i++) shift_heavy += iso_labels[mc[i]].shift_heavy;
	if(shift_medium == 0 || shift_heavy == 0) continue; // Skip no-shift case.
	s.mass_shift = shift_medium; 
	s.mass_shift2 = shift_heavy;
      }      
      s.labels = mc;
      shift_labels.push_back(s);  // Save/sort shift + label set.
    }
    // Sort shifts from small to large.
    sort(shift_labels.begin(), shift_labels.end());

    // Phase 1, get those that match on charge w/ no ambiguity.
    if(conf.quantmode == quantISOPAIR) {
      vector<xic *> candxics(shift_labels.size());
      // First group XICs with fragmentation spectra followed by
      // grouping XICs with matching charge.
      // Use more agressive pairing if 18O ATP data.
      const int max_mode = 2;  
      for(int mode = 0; mode < max_mode; mode++) {
	// Mode 0 = increasing m/z must have fragmentation spectra.
	// Mode 1 = increasing m/z matched charge and type. 

	// NOTE: Mode 2 and mode 3 are not used.
	// Mode 2 = increasing m/z must have fragmentation spectra, no type match.
	// Mode 3 = decreasing m/z must have fragmentation spectra, no type match

	for(size_t i = 0; i < xics_mz.size(); i++) {
	  if(xics_mz[i]->isotope_grouped) continue;
	  if(mode == 0 && xics_mz[i]->ms2peaks.size() == 0) continue; // mode 0 must have spectra
	  if(mode == 2 && xics_mz[i]->ms2peaks.size() == 0) continue; // mode 2 must have spectra
	  if(mode == 3 && xics_mz[i]->ms2peaks.size() == 0) continue; // mode 3 must have spectra
	  // Clear candidates for each shift. 
	  for(size_t c = 0; c < candxics.size(); c++) candxics[c] = NULL;
	  bool ignore_type = false; if(mode == 2 || mode == 3) ignore_type = true;
	  double sign = 1; if(mode == 3) sign = -1;
	  // Get XICs for each possible shift.
	  for(size_t s = 0; s < shift_labels.size(); s++) 
	    candxics[s] = find_isotope(xics_mz[i], sign * shift_labels[s].mass_shift, ignore_type);
	  int candcount = 0, cand_idx = -1;
	  for(size_t c = 0; c < candxics.size(); c++) { if(candxics[c] != NULL) { cand_idx = c; candcount++; } }
	  // Make sure only one is found.

	  if(candcount == 1) add_pair(xics_mz[i], candxics[cand_idx], shift_labels[cand_idx].labels); 
	}
	if(mode == 2) reverse(xics_mz.begin(), xics_mz.end()); //going to mode 3, so reverse
      }
    }
    else if(conf.quantmode == quantISOTRIPLE) {
      vector<xic *> candxics(shift_labels.size()), candxics2(shift_labels.size());
      for(int mode = 0; mode < 2; mode++) { // Process XICs with spectra first.
	for(size_t i = 0; i < xics_mz.size(); i++) {
	  if(xics_mz[i]->isotope_grouped) continue; // skip XICs that have been already grouped
	  if(mode == 0 && xics_mz[i]->ms2peaks.size() == 0) continue; 

	  // Clear candidates for each shift. 
	  for(size_t s = 0; s < shift_labels.size(); s++) { candxics2[s] = candxics[s] = NULL; }

	  // Get XICs for each shift.
	  for(size_t s = 0; s < shift_labels.size(); s++) {
	    candxics[s] = find_isotope(xics_mz[i], shift_labels[s].mass_shift);
	    candxics2[s] = find_isotope(xics_mz[i], shift_labels[s].mass_shift2);
	  }
	  // Count the number of candidate XICs for each shift.
	  int candcount = 0;
	  for(size_t s = 0; s < shift_labels.size(); s++) 
	    if(candxics[s] != NULL && candxics2[s] != NULL) candcount++;

	  // Only one candidate XIC triple for each shift found.
	  if(candcount == 1) {
	    // Get that candidate medium and heavy XICs. 
	    int medium_idx = -1, heavy_idx = -1;
	    for(size_t s = 0; s < shift_labels.size(); s++) 
	      if(candxics[s] != NULL && candxics2[s] != NULL) { medium_idx = heavy_idx = s; break; }
	    // Require that same shift combination be used.
	    add_triple(xics_mz[i], candxics[medium_idx], candxics2[heavy_idx], shift_labels[medium_idx].labels); 
	  }
	}
	// MH case
	for(size_t i = 0; i < xics_mz.size(); i++) {
	  if(xics_mz[i]->isotope_grouped) continue; // skip XICs that have been already grouped
	  if(mode == 0 && xics_mz[i]->ms2peaks.size() == 0) continue; 

	  for(size_t c = 0; c < candxics.size(); c++) candxics[c] = NULL; 	// Clear candidates for each shift. 
	  // Get XICs for each shift.
	  for(size_t s = 0; s < shift_labels.size(); s++) 
	    candxics[s] = find_isotope(xics_mz[i], shift_labels[s].mass_shift2 - shift_labels[s].mass_shift );
	  // Count candidate XICs for each shift.
	  int candcount = 0; for(size_t c = 0; c < candxics.size(); c++) if(candxics[c] != NULL) candcount++;
	  // Only one candidate XIC for each shift found.
	  if(candcount == 1) {
	    // Get that candidate medium and heavy XICs. 
	    for(size_t c = 0; c < candxics.size(); c++) 
	      if(candxics[c] != NULL) { add_triple(NULL, xics_mz[i], candxics[c], shift_labels[c].labels); break;}
	  }
	}
	// LM case
	for(size_t i = 0; i < xics_mz.size(); i++) {
	  if(xics_mz[i]->isotope_grouped) continue; // skip XICs that have been already grouped
	  if(mode == 0 && xics_mz[i]->ms2peaks.size() == 0) continue; 

	  for(size_t c = 0; c < candxics.size(); c++) candxics[c] = NULL; 	// Clear candidates for each shift. 
	  // Get XICs for each shift.
	  for(size_t s = 0; s < shift_labels.size(); s++) 
	    candxics[s] = find_isotope(xics_mz[i], shift_labels[s].mass_shift );
	  // Count candidate XICs.
	  int candcount = 0; for(size_t c = 0; c < candxics.size(); c++) if(candxics[c] != NULL) candcount++;
	  // Only one candidate XIC for each shift tried.
	  if(candcount == 1) {
	    for(size_t c = 0; c < candxics.size(); c++) 
	      if(candxics[c] != NULL) { add_triple(xics_mz[i], candxics[c], NULL, shift_labels[c].labels); break;}
	  }
	}
	// LH case
	for(size_t i = 0; i < xics_mz.size(); i++) {	
	  if(xics_mz[i]->isotope_grouped) continue; // skip XICs that have been already grouped
	  for(size_t c = 0; c < candxics.size(); c++) candxics[c] = NULL; 	// Clear candidates for each shift. 
	  // Get XICs for each shift.
	  for(size_t s = 0; s < shift_labels.size(); s++) 
	    candxics[s] = find_isotope(xics_mz[i], shift_labels[s].mass_shift2);
	  // Count candidate XICs.
	  int candcount = 0; for(size_t c = 0; c < candxics.size(); c++) if(candxics[c] != NULL) candcount++;
	  // Only one candidate XIC for each shift tried.
	  if(candcount == 1) {
	    for(size_t c = 0; c < candxics.size(); c++) 
	      if(candxics[c] != NULL) { add_triple(xics_mz[i], NULL, candxics[c], shift_labels[c].labels); break;}
	  }
	}
      }
    }
    cout << description << ": Found " << isotope_data.size() << " isotope groups." << endl;
    isotope_idx = new _2dtree<isotope_group>(isotope_data); // Build 2d-tree index on groups. 
  }

  xic *irun::find_15Npair(int &candcount, int &cnt_15N, xic *cur, fragmentDB *ms2db, float direction) {
    candcount = 0;
    if(cur->type != xicMONOISO) return NULL; // not monoisotopic, skip
    if(!(cur->charge > 0)) return NULL;      // no charge information, skip
    if(cur->isotope_grouped) return NULL;     // grouped already, skip

    // Query database for given m/z and charge of the monoisotopic
    // XIC.  NOTE: This also accounts for any active PTMs.
    vector<fragment *> res; ms2db->query_accurate_mass(cur->mz, cur->charge, res);
    if(res.size() == 0) return NULL; // No fragments were found.

    // Sort the nitrogen counts of each returned fragment.
    vector<int> nitrogens(res.size());
    for(size_t r = 0; r < res.size(); r++) nitrogens[r] = res[r]->nitrogens;
    sort(nitrogens.begin(), nitrogens.end());

    // Get unique nitrogen shifts to check.
    vector<int> nitrogen_shifts;
    size_t n = 0;
    while(n < nitrogens.size()) { // Note: nitrogens.size() >= 1 since res.size() >= 1
      while(n < nitrogens.size() - 1 && nitrogens[n] == nitrogens[n+1]) n++;
      nitrogen_shifts.push_back(nitrogens[n]);
      n++;
    }

    // Check each nitrogen shift returned by database search.
    vector<xic *> candxics(nitrogen_shifts.size(), NULL);
    for(size_t s = 0; s < nitrogen_shifts.size(); s++) {
      candxics[s] = find_isotope(cur, direction * mass::delta15N * double(nitrogen_shifts[s]));
    }

    // Make sure only one such candidate was returned.
    candcount = 0;
    for(size_t c = 0; c < candxics.size(); c++) if(candxics[c] != NULL) candcount++;
    if(candcount == 1) {
      // Get that candidate XIC. 
      xic *cand = NULL;
      for(size_t c = 0; c < candxics.size(); c++) {
	if(candxics[c] != NULL) { cand = candxics[c]; cnt_15N = nitrogen_shifts[c]; break; }
      }
      return cand;
    }
    // NOTE: This set might be resolvable by MS/MS scoring.
    /*else if(candcount > 1 && cur->ms2peaks.size() > 0) {
      }*/
    return NULL;
  }

  void irun::find_15Npairs(vector<fragmentDB *> &DB) {
    // Set all XICs to ungrouped by default.
    for(size_t i = 0; i < xics.size(); i++) xics[i].isotope_grouped = false; 

    // Get pointers to XICs that are monoisotopic.
    vector<xic*> xics_mz;
    for(size_t i = 0; i < xics.size(); i++)  // NOTE: These *should* have charge info.
      if(xics[i].type == xicMONOISO) xics_mz.push_back(&xics[i]);

    // Sort monoisotopic XICs on their quantification values from 
    // most intense to least intense. 
    sort(xics_mz.begin(), xics_mz.end(), quant_sort_decr());   

    // TODO: Go through the set with fragmentation spectra first.
    // Resolve any ambigious cases using MS/MS scoring.
    
    // Then go through the rest by decreasing intensity.
    int ambig15N = 0;
    // Phase 1, look for obvious candidates where charge matches.
    for(int mode = 0; mode < 2; mode++) { // Go through pairs with fragmentation spectra first.
      for(size_t i = 0; i < xics_mz.size(); i++) {
	if(xics_mz[i]->isotope_grouped) continue;
	if(mode == 0 && xics_mz[i]->ms2peaks.size() == 0) continue; 
	// Query from light to heavy.
	int candcountL2H = 0;
	int cand_15N_cnt_L2H = 0;
	xic *candL2H = find_15Npair(candcountL2H, cand_15N_cnt_L2H, xics_mz[i], DB[dbLIGHT], 1.0);
	// Query also from heavy to light. 
	int candcountH2L = 0;
	int cand_15N_cnt_H2L = 0;
	xic *candH2L = find_15Npair(candcountH2L, cand_15N_cnt_H2L, xics_mz[i], DB[dbHEAVY], -1.0);

	// This set might also be resolvable by MS/MS fragmentation scoring.
	if((candcountH2L > 0 && candcountL2H > 0) || candcountL2H > 1 || candcountH2L > 1) ambig15N++;

	//if(candL2H == NULL && candH2L == NULL) continue; // Nothing found.
	//if(candL2H != NULL && candH2L != NULL) continue; // Ambigious in heavy and light direction.
	if((candcountL2H == 1 && candcountH2L == 0) || (candcountH2L == 1 && candcountL2H == 0) ) {
	  if(candL2H != NULL) {
	    if(candL2H->isotope_grouped) continue;
	    add_pair(xics_mz[i], candL2H, cand_15N_cnt_L2H);
	  }
	  else {
	    if(candH2L->isotope_grouped) continue;
	    add_pair(candH2L, xics_mz[i], cand_15N_cnt_H2L);
	  }
	}
      }
    }
    cout << description << ": ambigious 15N pairs: " << ambig15N << endl;
    // Label MS/MS peaks as heavy and light.
    for(size_t i = 0; i < isotope_data.size(); i++) {
      xic *light = isotope_data[i].light;
      for(size_t p = 0; p < light->ms2peaks.size(); p++) light->ms2peaks[p].iso = ms2::isoLIGHT;
      xic *heavy = isotope_data[i].heavy;
      for(size_t p = 0; p < heavy->ms2peaks.size(); p++) heavy->ms2peaks[p].iso = ms2::isoHEAVY;
    }
    cout << description << ": Found " << isotope_data.size() << " isotope groups." << endl;
    isotope_idx = new _2dtree<isotope_group>(isotope_data); // Build 2d-tree index on groups. 
  }

  int irun::find_nKs(double &max_score, ms2::fragmentDB *curDB, vector<ms2peak> &ms2peaks) {
    max_score = 0;
    // Pre-process peaks.
    for(size_t m = 0; m < ms2peaks.size(); m++) curDB->pre_process(&ms2peaks[m]);

    // Perform search. Get peptide with highest search score.
    ms2peak *best_match = NULL; 
    for(size_t m = 0; m < ms2peaks.size(); m++) {
      vector<isotope_t> iso_none;
      curDB->search(&ms2peaks[m], min_mz2, max_mz2, iso_none,0);
      double target_score = ms2peaks[m].max_target_score();
      if(target_score > max_score) {
	best_match = &ms2peaks[m]; max_score = target_score;
      }
    }
    if(best_match == NULL) return -1;
    best_match->score_cutoff(max_score);

    // Get number of K or #'s for best matches.
    set<int> K_cnts;
    vector<fragment *> ids; best_match->max_ids(ids);
    vector<mods_t> mods; best_match->max_mods(mods);
    for(size_t d = 0; d < ids.size(); d++) {
      int K_cnt = 0;
      string seq = curDB->fasta.get_seqN(ids[d], mods[d]);
      for(size_t s = 0; s < seq.length(); s++) { if(seq[s] == 'K' || seq[s] == '#') K_cnt++;}
      if(K_cnt > 0) K_cnts.insert(K_cnt);
    }
    // Make sure this number is unique since best match can map to multiple peptides.
    if(K_cnts.size() != 1) return -1; 
	
    return *(K_cnts.begin());
  }
  void irun::find_stoich_patterns(vector<fragmentDB *> &DB) {
    if(DB.size() == 0) return; if(DB[dbLIGHT] == NULL || DB[dbHEAVY] == NULL) return;

    // Set all XICs to ungrouped by default.
    for(size_t i = 0; i < xics.size(); i++) xics[i].isotope_grouped = false; 

    vector<xic*> xics_mz;     // Get pointers to XICs that are monoisotopic.
    for(size_t i = 0; i < xics.size(); i++)  // NOTE: These *should* have charge info.
      if(xics[i].type == xicMONOISO && xics[i].ms2peaks.size() > 0) xics_mz.push_back(&xics[i]);

    sort(xics_mz.begin(), xics_mz.end(), mz_sort_incr());      // Sort monoisotopic XICs on m/z dimension.

    // #K cnt = 1
    // phase 0, heavy to light match charge, single xic added
    // phase 1, heavy to light no match, single xic added
    // phase 2, light to heavy match charge, single xic added
    // phase 3, light to heavy no match, single xic added
    // #K cnt > 1
    // phase 4, heavy to light, 
    // phase 5, light to heavy
    for(int phase = 0; phase < 6; phase++) {
      for(size_t i = 0; i < xics_mz.size(); i++) {
	if(xics_mz[i]->isotope_grouped) continue; // Skip XICs that have already been grouped.

	// Make a copy of MS/MS spectrum.
	vector<ms2peak> ms2peaks = xics_mz[i]->ms2peaks;

	fragmentDB *curDB = NULL;
	if(phase == 0 || phase == 1 || phase == 4) curDB = DB[dbHEAVY]; 
	if(phase == 2 || phase == 3 || phase == 5) curDB = DB[dbLIGHT];

	// Pre-process peaks.
	for(size_t m = 0; m < ms2peaks.size(); m++) curDB->pre_process(&ms2peaks[m]);

	double max_score = 0;
	int nKs = find_nKs(max_score, curDB, ms2peaks);
	if(nKs <= 0) continue;
	int xics_added = 0;
	vector<xic *> xics; xics.push_back(xics_mz[i]);
	for(int Ks = 1; Ks <= nKs; Ks++) {
	  double mass_shift = 0;
	  if(conf.acetyl_stoich_plus5)  { using namespace mass; mass_shift = ((C12*2 + H*2 + O) - (C13*2 + H2*3 - H + O)) * double(Ks);  }
	  else { using namespace mass; mass_shift = ((C12*2 + H*2 + O) - (C12*2 + H2*3 - H + O)) * double(Ks);  }
	  if(phase == 2 || phase == 3 || phase == 5) mass_shift *= -1.0; // light to heavy shift
	  bool ignore_type = true;
	  if(phase == 0 || phase == 2) ignore_type = false; // match on XIC type/charge if in phase
	  xic *isotope_xic = find_isotope(xics_mz[i], mass_shift, ignore_type) ;
	  if(isotope_xic != NULL) xics_added++;
	  xics.push_back(isotope_xic);
	}
	// First handle peptides with single lysines.
	if(phase == 0 || phase == 1 || phase == 2 || phase == 3) {
	  if(nKs == 1 && xics_added == 1) {
	    // Make heavy to light.
	    if(phase == 2 || phase == 3) reverse(xics.begin(), xics.end());
	    for(size_t m = 0; m < xics[0]->ms2peaks.size(); m++) xics[0]->ms2peaks[m].iso = isoHEAVY;
	    for(size_t m = 0; m < xics[1]->ms2peaks.size(); m++) xics[1]->ms2peaks[m].iso = isoLIGHT;
	    add_group_phase(xics, phase);
	  }
	}
	// Now handle peptides with multiple lysines.
	if(nKs > 1 && xics_added >= 1) {
	  if(phase == 4) {
	    // Mark MS/MS spectra with light heavy label.
	    for(size_t m = 0; m < xics[0]->ms2peaks.size(); m++) xics[0]->ms2peaks[m].iso = isoHEAVY;
	    if(xics[xics.size()-1] != NULL) {
	      for(size_t m = 0; m < xics[xics.size()-1]->ms2peaks.size(); m++) 
		xics[xics.size()-1]->ms2peaks[m].iso = isoHEAVY;
	    }
	    add_group_phase(xics, phase);
	  }
	  if(phase == 5) { // light to heavy pairing.
	    for(size_t m = 0; m < xics[0]->ms2peaks.size(); m++) xics[0]->ms2peaks[m].iso = isoLIGHT;
	    if(xics[xics.size()-1] != NULL) {
	      for(size_t m = 0; m < xics[xics.size()-1]->ms2peaks.size(); m++) 
		xics[xics.size()-1]->ms2peaks[m].iso = isoHEAVY;
	    }
	    // Always keep XIC groups in order of heavy to light.
	    reverse(xics.begin(), xics.end());
	    add_group_phase(xics, phase);
	  }
	}
      }    
    }
    cout << description << ": Found " << isotope_data.size() << " isotope groups." << endl;

    // Now add XIC pairs with 100 or 0% stiochiometry, using MS/MS search to label.
    for(size_t i = 0; i < xics_mz.size(); i++) {
      if(xics_mz[i]->isotope_grouped) continue; // Skip XICs that have already been grouped.
      vector<ms2peak> ms2peaks = xics_mz[i]->ms2peaks;
      double score_light = 0;
      int nKs_light = find_nKs(score_light, DB[dbLIGHT], ms2peaks);
      double score_heavy = 0;
      int nKs_heavy = find_nKs(score_heavy, DB[dbHEAVY], ms2peaks);
      if(score_heavy > 0 && score_light > 0) {
	if(nKs_light > 0 && score_light >= score_heavy) {
	  for(size_t m = 0; m < xics_mz[i]->ms2peaks.size(); m++) xics_mz[i]->ms2peaks[m].iso = isoLIGHT;
	  vector<xic *> xics; 
	  xics.push_back(NULL); xics.push_back(xics_mz[i]);
	  add_group(xics);
	}
	else if(nKs_heavy > 0 && score_heavy > score_light) {
	  for(size_t m = 0; m < xics_mz[i]->ms2peaks.size(); m++) xics_mz[i]->ms2peaks[m].iso = isoHEAVY;
	  vector<xic *> xics; 
	  xics.push_back(xics_mz[i]); xics.push_back(NULL);
	  add_group(xics);
	}
      }
    }
    cout << description << ": Found " << isotope_data.size() << " isotope groups (after search)." << endl;
  }

  // Collect 13C shifts relative to a given XIC and charge.
  void irun::collect13C(vector<xic *> &carbon13, xic *x, int charge, double sign) {
    xic *cur_xic = x;
    while(true) {
      xic *xic_13C = get_overlapped(cur_xic, sign * mass::delta13C / double(charge), conf.isotolerence);
      if(xic_13C == NULL) break;
      carbon13.push_back(xic_13C);
      cur_xic = xic_13C;
    }
  }

  // Labels monoisotopic XICs + assigns them charge.
  void irun::monoiso_and_charge() {
    // Sort XICs on their m/z, increasing.
    vector<xic*> xics_mz(xics.size()); 
    for(size_t i = 0; i < xics.size(); i++) xics_mz[i] = &xics[i];
    sort(xics_mz.begin(), xics_mz.end(), mz_sort_incr());      

    // First find isotope XICs.
    for(size_t i = 0; i < xics_mz.size(); i++) {
      xic *x = xics_mz[i]; // Scan through XICs with increasing m/z.
      // Find XIC at lower m/z that has the highest intensity.
      xic *maxI = NULL; 
      for(int charge = conf.max_charge; charge >= 1; charge--) {
	// NOTE: This looks at a "small time range" 
	xic *xic_13C = get_overlapped(x, -mass::delta13C / double(charge), conf.isotolerence);
	if(xic_13C == NULL) continue;
	if(maxI == NULL) { maxI = xic_13C; continue; }
	if(xic_13C->quant > maxI->quant) maxI = xic_13C; 
      }
      if(maxI == NULL) continue;
      // Mark as isotope, if XIC at lower m/z is 40% of the current
      // XIC intensity or higher. Sort of implements an isotope
      // distribution model
      if(maxI->quant > conf.spurxictol * x->quant) { x->type = xicISOTOPE; }
    }
    // Iterate allowing less an less intense 13C XICs.
    float alpha = 0.6; 
    for(int iter = 0; iter < 3; iter++) {
      for(size_t i = 0; i < xics_mz.size(); i++) {
	xic *x = xics_mz[i]; // Scan through XICs with increasing m/z.
	if(x->type == xicISOTOPE) continue; // Skip isotopes.
	if(x->type == xicMONOISO) continue; // Skip those labeled as monoiso..
	for(int charge = conf.max_charge; charge >= 1; charge--) {
	  xic *xic_13C = get_overlapped(x, mass::delta13C / double(charge), conf.isotolerence);
	  if(xic_13C != NULL && xic_13C->quant > alpha * x->quant && xic_13C->type != xicMONOISO) {
	    x->type = xicMONOISO; 
	    x->charge = charge;
	    break;
	  }
	}
      }
      alpha *= 0.5;
    }
  }


  // Compute max m/z (MS1 and MS2) and retention time values.
  void irun::data_range() {
    min_mz = numeric_limits<float>::max();   max_mz   = -(numeric_limits<float>::max());
    min_time = numeric_limits<float>::max(); max_time = -(numeric_limits<float>::max());
    // Get data ranges.
    for(size_t i = 0; i < data.size(); i++) {
      if(data[i].mz > max_mz) max_mz = data[i].mz;
      if(data[i].mz < min_mz) min_mz = data[i].mz;
      if(data[i].retentionTime < min_time) min_time = data[i].retentionTime;
      if(data[i].retentionTime > max_time) max_time = data[i].retentionTime;
    }
    for(size_t i = 0; i < filtered.size(); i++) {
      if(filtered[i].mz > max_mz) max_mz = filtered[i].mz;
      if(filtered[i].mz < min_mz) min_mz = filtered[i].mz;
      if(filtered[i].retentionTime < min_time) min_time = filtered[i].retentionTime;
      if(filtered[i].retentionTime > max_time) max_time = filtered[i].retentionTime;
    }
    min_mz2 = numeric_limits<float>::max();  max_mz2 = -(numeric_limits<float>::max());
    for(size_t i = 0; i < data2.size(); i++) {
      float mz2_high = data2[i].max_mz2(); 
      float mz2_low = data2[i].min_mz2();
      if(mz2_high > max_mz2) max_mz2 = mz2_high; 
      if(mz2_low < min_mz2) min_mz2 = mz2_low;
    }
    for(size_t i = 0; i < filtered2.size(); i++) {
      float mz2_high = filtered2[i].max_mz2(); 
      float mz2_low = filtered2[i].min_mz2();
      if(mz2_high > max_mz2) max_mz2 = mz2_high; 
      if(mz2_low < min_mz2) min_mz2 = mz2_low;
    }
  }

  // Computes histogram of intensity values. Computes an intensity
  // threshold that removes conf.hist_logI_thresh fraction of the data
  // peaks.
  double irun::histogram_threshold() {
    if(data.size() == 0) return 0;
    // Get min and max intensity.
    float min_intensity = numeric_limits<float>::max(), max_intensity = -(numeric_limits<float>::max());
    for(size_t i = 0; i < data.size(); i++) {
      if(data[i].intensity < min_intensity) min_intensity = data[i].intensity;
      if(data[i].intensity > max_intensity) max_intensity = data[i].intensity;
    }

    // Compute a histogram of log10(intensity) values.
    const int bins = 51200;
    vector<long> loghist(bins+1,0);
    double min_logI = log10(min_intensity), max_logI = log10(max_intensity);
    double range = max_logI - min_logI;
    for(size_t i = 0; i < data.size(); i++) {
      // Normalize between 0 and 1.
      double alpha = (log10(data[i].intensity) - min_logI) / range; 
      loghist[int(bins * alpha)]++; // update histogram bin
    }    
    
    // Advance index until fmin fraction of peaks are accounted for.
    double Ndata = data.size(), cum_min = 0;
    int bmin = 0;
    double fmin = conf.hist_log10I_thresh;
    while(cum_min / Ndata < fmin && bmin < bins) cum_min += loghist[bmin++];
    bmin--;

    // Compute an intensity threshold for this fraction of peaks.
    return double(bmin) / (double)bins * range + min_logI;
  }

  // Compute histogram of MS1 intensity values, adjust display values
  // so display is more pleasing.
  void irun::histogram_adjustment() {
    logImin_disp = 3; logImax_disp = 6; // Set default values.
    if(data.size() == 0) return;
    // Get min and max intensity.
    float min_intensity = numeric_limits<float>::max(), max_intensity = -(numeric_limits<float>::max());
    for(size_t i = 0; i < data.size(); i++) {
      if(data[i].intensity < min_intensity) min_intensity = data[i].intensity;
      if(data[i].intensity > max_intensity) max_intensity = data[i].intensity;
    }
    // Compute a histogram of log10(intensity) values.
    const int bins = 5120;
    vector<long> loghist(bins+1,0);
    double min_logI = log10(min_intensity), max_logI = log10(max_intensity);
    double range = max_logI - min_logI;
    for(size_t i = 0; i < data.size(); i++) {
      // Normalize between 0 and 1.
      double alpha = (log10(data[i].intensity) - min_logI) / range; 
      loghist[int(bins * alpha)]++; // update histogram bin
    }
    // Given the histogram, compute the logImin value that elimates
    // fmin percent of the data and the logI max value that eliminates
    // fmax of the data.
    double Ndata = data.size(), cum_min = 0, cum_max = 0;
    int bmin = 0, bmax = bins - 1;
    const double fmin = 0.30, fmax = 0.001;
    while(cum_min / Ndata < fmin && bmin < bins) cum_min += loghist[bmin++];
    bmin--;
    while(cum_max / Ndata < fmax && bmax >= 0)  cum_max += loghist[bmax--];
    bmax++;
    logImin_disp = double(bmin) / (double)bins * range + min_logI;
    logImax_disp = double(bmax) / (double)bins * range + min_logI;
  }

  // Main LCMS Processing algorithm for finding xics.
  // 
  // This is the main driver of all of the algorithms. Start here when
  // trying to understand this code.  TODO: Add better error
  // handling. Throw an exception if these steps go wrong.
  void irun::process(vector<ms2::fragmentDB *> &DB) { 
    cout << "analysis_id=" << conf.quantmode << endl;
    if(processed_flag) return; else processed_flag = true;

    // Otherwise perform processing.
    if(xml::load_mzXML_MS1(filename, retentionTimes, data, conf.ms1_profile) == false) { 
      data.clear(); cerr << description << ": failed loading MS1 data" << endl; return; 
    }
    if(data.size() == 0) return; // No MS1 data, nothing to do.

    if(conf.load_ms2) {
      bool ok = xml::load_mzXML_MS2(filename, data2, conf.ms2_profile);       // Load MS2 spectra if requested.
      if(!ok) { data2.clear(); cerr << description << ":failed loading MS2 data" << endl; return; }
    }

    cout << description << ": Indexing " << data.size() << " peaks." << endl;
    indexraw();  // Index raw centroided data LC-MS data using a 2d-tree
    data_range(); // Compute range of the data
    histogram_adjustment(); // Histogram intensity values for GUI display.
    if(conf.run_algorithms == false) return; // Stop processing here if requested.

    cout << description << ": Filtering raw peak data." << endl;
    // Find peak intensity threshold by histogram percentage thresholding.
    double hist_thresh = histogram_threshold();
    cout << description << ": hist_thresh=" << hist_thresh << endl;
    filter_new(hist_thresh); // Run peak filtering algorithm.
    cout << description << ": Indexing " << filtered.size() << " filtered peaks." << endl; 
    indexfiltered();  // Index filtered peaks w/ a 2d-tree. 

    if(conf.keep_raw == false) {       // Free up memory associated with raw data if requested.
      delete rawroot; rawroot = NULL;  delete rawroot2; rawroot2 = NULL;  
      data.clear(); vector<peak> (data).swap(data); 
      data2.clear(); vector<ms2peak> (data2).swap(data2); 
    }

    // Get rid of raw MS/MS peaks if requested.
    if(conf.keep_ms2 == false) { for(size_t i = 0; i < data2.size(); i++) data2[i].clear_ms2peaks(); }

    cout << description << ": Finding XICs." << endl;
    findxics();   // Find XICs. 
    cout << description << ": Found " << xics.size() << " XICs." << endl;

    // Free up memory associated with filtered peak data if requested.
    if(conf.keep_filtered == false) {
      delete root; root = NULL;
      filtered.clear(); vector<peak> (filtered).swap(filtered); 
    }

    cout << description << ": Identifying mono-iso XICs + determining charge." << endl;
    monoiso_and_charge(); // Determine monoisotopc XICs and assign charge information.

    cout << description << ": Filtering " << xics.size() << " XICs." << endl;
    filter_xics(); // XIC filtering of "close" XICs.

    cout << description << ": Assigning MS/MS peaks to XICs." << endl;
    assign2xic();  // Assign MS/MS peaks to an XIC.

    int xic_stat_monoiso = 0, xic_stat_monoiso_frag = 0;
    for(size_t x = 0; x < xics.size(); x++) { 
      if(xics[x].type == xicMONOISO) xic_stat_monoiso++; 
      if(xics[x].type == xicMONOISO && xics[x].ms2peaks.size() > 0) xic_stat_monoiso_frag++;
    }

    cout << description << ": XICs found " << xics.size() << endl;
    cout << description << ": monoisotopic XICs found " << xic_stat_monoiso << endl;
    cout << description << ": monoisotopic XICs found with frag. spectra " << xic_stat_monoiso_frag << endl;

    // Get rid of filtered MS/MS peaks.
    if(conf.keep_ms2  == false) { for(size_t i = 0; i < filtered2.size(); i++) filtered2[i].clear_ms2peaks();  }
    
    if(conf.keep_filtered == false) {
      // Free up filtered MS/MS peaks if requested.
      delete root2; root2 = NULL;
      filtered2.clear(); vector<ms2peak> (filtered2).swap(filtered2);
    }

    if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE) {
      cout << description << ": Processing isotope data." << endl;
      if(conf.isotope_15N) find_15Npairs(DB); 
      else find_isotope_pairs_and_triples();
    }
    if(conf.quantmode == quantAcetylSTOICH) {  find_stoich_patterns(DB); }

    // Run MS/MS database search for the current irun.
    cout << description << ": Running MS/MS DB search." << endl;
    if( conf.quantmode == quantISOPAIR  && DB[dbLIGHT] != NULL && DB[dbHEAVY] != NULL) {
      // Collect all heavy and light spectra.
      vector<ms2::ms2peak *> ms2_light, ms2_heavy;
      vector<isotope_t> cur_iso;
      int id_cnt = 0;
      for(size_t i = 0; i < isotope_data.size(); i++) {
	ms2_light.clear(); ms2_heavy.clear();
	isotope_data[i].light->get_ms2(ms2_light); isotope_data[i].heavy->get_ms2(ms2_heavy);

	// Run pre-processing for light and heavy.
	DB[dbLIGHT]->pre_process(ms2_light); DB[dbHEAVY]->pre_process(ms2_heavy);

	if(ms2_light.size() + ms2_heavy.size() == 0) continue;
	id_cnt++;

	// Get active isotopes for this group.
	vector<int> &labels = isotope_data[i].labels;
	cur_iso.clear();
	for(size_t l = 0; l < labels.size(); l++) 
	  if(labels[l] != 0)  // NOTE: Don't use the 0th no isotope label in iso_labels.
	    cur_iso.push_back(iso_labels[labels[l]]);

	DB[dbLIGHT]->search(ms2_light, min_mz2, max_mz2, cur_iso, isotope_data[i].nitrogens); 
	DB[dbHEAVY]->search(ms2_heavy, min_mz2, max_mz2, cur_iso, isotope_data[i].nitrogens);
      }
      cout << description <<  ": isotope group id count = " << id_cnt << endl;      
      if(conf.isotope_by_search) {
	cout << "isotope_data.size()=" << isotope_data.size() << endl;
	delete isotope_idx;
	for(size_t i = 0; i < xics.size(); i++) {
	  xic &x = xics[i];
	  if(x.type == xicMONOISO && x.ms2peaks.size() > 0 && x.isotope_grouped == false) {
	    vector<ms2::ms2peak *> ms2data;
	    // Ungrouped XIC. Determine if it's from heavy or light.
	    x.get_ms2(ms2data);
	    DB[dbLIGHT]->pre_process(ms2data); 
	    vector<isotope_t> isoempty;
	    DB[dbLIGHT]->search(ms2data, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
	    double light_score = x.max_td_score();
	    x.clear_targdecoy();
	    DB[dbHEAVY]->search(ms2data, min_mz2, max_mz2, isoempty); // Run MS/MS database search.	  
	    double heavy_score = x.max_td_score();
	    x.clear_targdecoy();
	    if(light_score > 0 && heavy_score > 0) {
	      isotope_group p;
	      p.condIndex = condIndex; p.repIndex = repIndex; p.irunIndex = irunIndex;
	      p.by_search = true;
	      if(light_score >= heavy_score) {
		// This is a light XIC.
		DB[dbLIGHT]->search(ms2data, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
		for(size_t m = 0; m < ms2data.size(); m++) ms2data[m]->iso = isoLIGHT;
		p.light = &xics[i];
		p.log2xicL = util::Log2(xics[i].quant);
	      }
	      else {
		// This is a heavy XIC.
		DB[dbHEAVY]->search(ms2data, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
		for(size_t m = 0; m < ms2data.size(); m++) ms2data[m]->iso = isoHEAVY;
		p.heavy = &xics[i];
		p.log2xicH = util::Log2(xics[i].quant);
	      }
	      isotope_data.push_back(p);
	    }
	  }
	}
	cout << "isotope_data.size()=" << isotope_data.size() << endl;
	isotope_idx = new _2dtree<isotope_group>(isotope_data); // Build 2d-tree index on groups. 
      }

    }
    else if( conf.quantmode == quantISOTRIPLE && DB[dbLIGHT] != NULL && DB[dbMEDIUM] != NULL && DB[dbHEAVY] != NULL) {
      // Collect all heavy and light spectra.
      vector<ms2::ms2peak *> ms2_light, ms2_medium, ms2_heavy;
      vector<isotope_t> cur_iso;
      for(size_t i = 0; i < isotope_data.size(); i++) {
	ms2_light.clear(); ms2_heavy.clear(); ms2_medium.clear();
	if(isotope_data[i].light != NULL) isotope_data[i].light->get_ms2(ms2_light); 
	if(isotope_data[i].medium != NULL) isotope_data[i].medium->get_ms2(ms2_medium); 
	if(isotope_data[i].heavy != NULL) isotope_data[i].heavy->get_ms2(ms2_heavy);

	// Run pre-processing for light, heavy, and medium.
	DB[dbLIGHT]->pre_process(ms2_light);
	DB[dbMEDIUM]->pre_process(ms2_medium); 
	DB[dbHEAVY]->pre_process(ms2_heavy);

	// Get active isotopes for this group.
	vector<int> &labels = isotope_data[i].labels;
	cur_iso.clear();
	for(size_t l = 0; l < labels.size(); l++) if(labels[l] != 0) cur_iso.push_back(iso_labels[labels[l]]);

	// TODO: Look for cached data here. This will be a little tricky with 
	// isotope information (assume isotope filtering applied already?). 
	// Run database search for light and heavy.
	DB[dbLIGHT]->search(ms2_light, min_mz2, max_mz2, cur_iso); 
	DB[dbMEDIUM]->search(ms2_medium, min_mz2, max_mz2, cur_iso); 
	DB[dbHEAVY]->search(ms2_heavy, min_mz2, max_mz2, cur_iso);
      }
    }
    else if(conf.quantmode == quantAcetylSTOICH) {
      // Collect light and heavy MS/MS spectra from isotope groups.
      for(size_t i = 0; i < isotope_data.size(); i++) {
	vector<ms2::ms2peak *> ms2data;
	for(size_t x = 0; x < isotope_data[i].xics.size(); x++) {
	  if(isotope_data[i].xics[x] != NULL)  isotope_data[i].xics[x]->get_ms2(ms2data);
	}
	vector<ms2::ms2peak *> ms2heavy, ms2light;
	for(size_t m = 0; m < ms2data.size(); m++) {
	  if(ms2data[m]->iso == isoHEAVY) ms2heavy.push_back(ms2data[m]);
	  else ms2light.push_back(ms2data[m]);
	}
	DB[dbHEAVY]->pre_process(ms2heavy);  // Pre-process spectra.
	DB[dbLIGHT]->pre_process(ms2light);  // Pre-process spectra.
	vector<isotope_t> isoempty;
	DB[dbHEAVY]->search(ms2heavy, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
	DB[dbLIGHT]->search(ms2light, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
      }
      // Collect remaining MS/MS spectra and search with light database to improve
      // q-value estimation.
      vector<ms2peak *> ms2data;
      for(size_t i = 0; i < xics.size(); i++) {
	if(xics[i].type == xicMONOISO && xics[i].isotope_grouped == false) {
	  xics[i].get_ms2(ms2data);
	}
      }
      DB[dbLIGHT]->pre_process(ms2data);  // Pre-process spectra.
      vector<isotope_t> isoempty;
      DB[dbLIGHT]->search(ms2data, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
      isotope_idx = new _2dtree<isotope_group>(isotope_data); // Build 2d-tree index on groups. 
    }
    else if(DB[dbLIGHT] != NULL) {
      // Get all MS/MS spectra from all XICs.
      vector<ms2::ms2peak *> ms2data;

      // NOTE: Only run MS/MS search if XIC is monoisotopic. 
      for(size_t i = 0; i < xics.size(); i++) if(xics[i].type == xicMONOISO) xics[i].get_ms2(ms2data);

      DB[dbLIGHT]->pre_process(ms2data);  // Pre-process spectra.
      
      vector<isotope_t> isoempty;
      DB[dbLIGHT]->search(ms2data, min_mz2, max_mz2, isoempty); // Run MS/MS database search.
    }
    // Run MS/MS database search for the current irun.
    cout << description << ": *DONE* running MS/MS DB search." << endl;

    // Free up memory used by MS/MS spectra in XICs peaks, if requested.
    if(conf.keep_ms2 == false) { for(size_t i = 0; i < xics.size(); i++) xics[i].clear_ms2peaks(); }
  } 

  void irun::get_shifts(irun *ref, vector<shift_t> &shifts, double dx_max) {
    for(size_t i = 0; i < xics.size(); i++) {
      if(xics[i].type == xicUNKNOWN) continue; // Unknown XIC, skip.
      // 1. Find the nearest point in reference irun.
      xic *refxic = ref->xics_all->NearestX(xics[i].x(), xics[i].y(), dx_max, 
					    xics[i].y() * conf.xic_width * 1e-6); 
      
      if(refxic != NULL && refxic->type == xics[i].type) {
	// Perform recipricol query.
	xic *recip = xics_all->NearestX(refxic->x(), refxic->y(), dx_max, 
					refxic->y() * conf.xic_width * 1e-6); 
	if(recip == &xics[i])  // If we get the same XIC back store the shift.
	  shifts.push_back(shift_t(xics[i].x(), refxic->x() - xics[i].x()));	  
      }
    }
  }

  void irun::align(irun *ref) {
    vector<shift_t> shifts;
    if(conf.align_translation) {
      // Get recipricol nearest shifts to reference irun.
      shifts.clear(); get_shifts(ref, shifts, conf.align_dtime);
      if(shifts.size() > 0) {
	// Apply translation.
	vector<double> dx(shifts.size());
	for(size_t i = 0; i < shifts.size(); i++) dx[i] = shifts[i].dx;
	translate(util::median_unsafe(dx));
      }
    }
    if(conf.align_nonlinear) {
      cout << "nonlinear alignment" << endl;
      shifts.clear(); get_shifts(ref, shifts, conf.align_dtime);
      if(shifts.size() > 0) {
	sort(shifts.begin(), shifts.end());
	// Compute running median estimates.
	vector<double> delta(shifts.size());
	for(size_t i = 0; i < shifts.size(); i++) delta[i] = shifts[i].dx;
	vector<double> delta_med; util::runmed(delta, delta_med);
	// Generate output using smoothed running medians.
	for(size_t i = 0; i < delta_med.size(); i++) shifts[i].dx = delta_med[i];
	nonlinear(shifts);
      }
    }
  }

  // Adjust the log2quant value for each XIC. Used for normalization.
  void irun::adjust_log2_quant(double delta) { 
    for(size_t i = 0; i < xics.size(); i++) xics[i].quant = pow(2,util::Log2(xics[i].quant) + delta); 
  }

}

