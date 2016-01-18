#ifndef __lcms_hpp__
#define __lcms_hpp__

#include <vector>
#include <string>
#include <limits>

#include "msconfig.hpp"
#include "2dtree.hpp"
#include "ms2.hpp"

// LC-MS/MS, single instrument run processing.
// 1. Noise filtering algorithm.
// 2. XIC detection algorithm.
// 3. Charge det + isotoping algorithm.
// 4. Isotope XIC grouping/tripling/15N algorithm.
// 5. Alignment to a reference run.

namespace lcms{
  using namespace std;
  using namespace ms2;

  // Minimal MS peak data. 
  struct peak {
    peak() { retentionTime = mz = intensity = 0; label = 0; }
    float retentionTime, mz, intensity;  
    int label; // Label used for DFS based grouping of peaks into XICs.  0 = unvisited, 1 = visitied, >= 2 component label
    // Accessors for _2dtree<>. 
    double x() const { return retentionTime; } double y() const { return mz; }
    // Used to sort peaks in a chromatogram.
    bool operator < (const peak &p) const { return retentionTime < p.retentionTime; }
  };

  enum xic_type { xicUNKNOWN, xicISOTOPE, xicMONOISO };

  // Peaks that occur at nearly the same m/z value over time.
  struct xic {
    xic_type type; 
    vector<peak> chrom; // Actual chromatogram data. 
    double mz, retentionTime, start, end; // (m/z, retention time data)
    double quant; // quantification value
    int charge; // charge of the XIC (assigned or acquired for MS/MS spectra)
    int condIndex, repIndex, irunIndex;  // Index of condition, replicate, and run this XIC originated from.
    vector<ms2peak> ms2peaks; // MS/MS peaks, if any that belong to XIC.
    bool isotope_grouped; // flag indicates XIC has been grouped
    int irunID; bool grouped; // Used for grouping across across instrument runs.

    // Accessors for _2dtree. 
    double x() const { return retentionTime; }  double y() const { return mz; }

    // Returns length of XIC in retention time.
    double length() { return end - start; }

    // Returns number of MS/MS IDs in earch MS/MS peak.
    int num_ids() { int ids = 0; for(size_t i = 0; i < ms2peaks.size(); i++) ids += ms2peaks[i].num_ids(); return ids; }

    // Returns maximum score for MS/MS peaks. NOTE: Assumes score > 0.
    double max_score() { double s = 0; for(size_t i = 0; i < ms2peaks.size(); i++) s = max(s, ms2peaks[i].max_score()); return s; }
    double max_td_score() { double s = 0; for(size_t i = 0; i < ms2peaks.size(); i++) s = max(s, ms2peaks[i].max_td_score()); return s; }

    // Gets pointers to all MS/MS peaks in XIC.
    void get_ms2(vector<ms2peak *> &peaks) { for(size_t i = 0; i < ms2peaks.size(); i++) peaks.push_back(&ms2peaks[i]); }

    // Returns MS/MS peak with highest score. Returns NULL if there are none.
    ms2peak *max_ms2() {
      double score = max_score();
      for(size_t i = 0; i < ms2peaks.size(); i++) {
	// TODO: Is == 0 right?
	if(ms2peaks[i].num_ids() > 0 && score - ms2peaks[i].max_score() == 0.0) return &ms2peaks[i];
      }
      return NULL;
    }
    // Free MS/MS peak memory from XIC or MS/MS peaks.
    void clear_ms2peaks() { for(size_t i = 0; i < ms2peaks.size(); i++) ms2peaks[i].clear_ms2peaks(); }
    void clear_chrom() { chrom.clear(); vector<peak> (chrom).swap(chrom); }
    void clear_targdecoy() { for(size_t i = 0; i < ms2peaks.size(); i++) { ms2peaks[i].clear_target(); ms2peaks[i].clear_decoy(); } }
  };

  enum ratio_type { isoHL = 0, isoML = 1, isoHM = 2};

  // Isotope group for SILAC data or dueterium or O18 labels or
  // whatever.  Change this struct to isotope_group. It now can handle
  // 3 labels (in the future, possibly more).
  struct isotope_group {
    long id; // unique ID for this measurement 
    vector<xic *> xics;
    int condIndex, repIndex, irunIndex; // used top cross reference back to condition/replciate number
    bool by_search;

    xic *light, *heavy, *medium; 
    double log2xicH, log2xicM, log2xicL; // usually log2(area under heavy XIC), but these get normalized
    double log2ratioHM, log2ratioHL, log2ratioML;
    vector<int> labels; // Integers that index iso_labels in lcms::irun.
                        // They indicate which labels were used for this group.
    int nitrogens;  // number of nitrogens in shift if 15N labeling is used
    int grouping_phase; // algorithm phase/iteration group was made

    isotope_group() {
      id = 0;
      medium = light = heavy = NULL; 
      log2ratioHM = log2ratioHL = log2ratioML = log2xicM = log2xicH = log2xicL = 0;
      condIndex = repIndex = irunIndex = -1; nitrogens = 0;
      by_search = false;
      grouping_phase = -1;
    }
    double x() const { 
      if(xics.size() > 0) {
	for(size_t x = 0; x < xics.size(); x++) { if(xics[x] != NULL) return xics[x]->retentionTime; }
      }
      if(light  != NULL) return light->retentionTime;
      if(medium != NULL) return medium->retentionTime;
      if(heavy  != NULL) return heavy->retentionTime;
      cerr << "x(): group problem" << endl; return 0;
    } 
    double y() const { 
      //if(xics.size() > 0) return xics[0]->mz;
      if(xics.size() > 0) {
	for(size_t x = 0; x < xics.size(); x++) { if(xics[x] != NULL) return xics[x]->mz; }
      }
      if(light  != NULL) return light->mz;
      if(medium != NULL) return medium->mz; 
      if(heavy  != NULL) return heavy->mz;
      cerr << "y(): group problem" << endl; return 0;
    }
    // NOTE: Assumes charges are all the same.
    int charge() {
      for(size_t x = 0; x < xics.size(); x++) { if(xics[x]->charge > 0) return xics[x]->charge;  }
      if(light != NULL  && light->charge > 0)  return light->charge;
      if(medium != NULL && medium->charge > 0) return medium->charge;
      if(heavy != NULL  && heavy->charge > 0)  return heavy->charge;
      return 0;
    }
        
    // Total number of MS/MS IDs assigned to both groups.
    int num_ids() { 
      int id_cnt = 0;
      for(size_t x = 0; x < xics.size(); x++) { if(xics[x] != NULL) id_cnt += xics[x]->num_ids(); }
      if(light != NULL) id_cnt += light->num_ids();
      if(medium != NULL) id_cnt += medium->num_ids();
      if(heavy != NULL) id_cnt += heavy->num_ids();
      return id_cnt;
    }
    
    // Return MS/MS peak corresponding to highest scoring ID.  Returns
    // NULL, if none exists.
    ms2peak *max_ms2() {
      if(xics.size() > 0) {
	double best_score = 0; ms2peak *best_ms2peak = NULL;
	for(size_t x = 0; x < xics.size(); x++) {
	  if(xics[x] == NULL) continue;
	  double xic_max_score = xics[x]->max_score();
	  if(xic_max_score > best_score) {
	    best_score = xic_max_score;
	    best_ms2peak = xics[x]->max_ms2();
	  }
	}
	return best_ms2peak;
      }
      else {
	double light_score = 0, medium_score = 0, heavy_score = 0;
	if(light != NULL) light_score = light->max_score();
	if(medium != NULL) medium_score = medium->max_score();
	if(heavy != NULL) heavy_score = heavy->max_score();
	if(heavy_score > medium_score && heavy_score > light_score) return heavy->max_ms2();
	if(medium_score > heavy_score && medium_score > light_score) return medium->max_ms2();
	if(light_score > heavy_score && light_score > medium_score) return light->max_ms2();
	return NULL;
      }
    }
  };

  // Delta alignment between two groups of runs.  Uses nearest points in run A and in run B.
  struct shift_t { 
    double x, dx; 
    shift_t(double x_, double dx_) { x = x_; dx = dx_; }
    shift_t(double x_) { x = x_; dx = 0; }
    bool operator < (const shift_t &s) const { return  x < s.x; }
  };

  // Description + file information for an LC-MS/MS runs.
  struct irun_file { 
    string description, filename; 
    bool operator < (const irun_file &b) const { return description < b.description; } 
  };
  
  // Recalibration mode.
  enum recalmode_t { recalmodeTIME, recalmodeMZ } ;

  // LCMS algorithms: Filtered LC-MS run data. Handle processing and finding xics.  
  // TODO: Change run to run or something like that. Why did I call it scan? 
  //       This is what happens when you teach yourself. 
  struct irun {
    msconfig &conf; // "Global" configuration information
    int condIndex, repIndex, irunIndex; // condition, replicate, run numbers
    string description; // Short description of the run.
    string filename; // mzXML file name
    bool processed_flag; // Flag indicates processing has been completed.
    bool labels_reversed; // Isotope label is swaped.
    // Data set range.
    float min_mz, max_mz, min_time, max_time, min_mz2, max_mz2;
    float logImin_disp, logImax_disp; // log10(intensity) display range as determined by histogram adjustment

    // WARNING: There are a bunch of vector<> types here. The indicies
    // rely on the memory allocated by these vector<> types to never
    // change. Do not add or remove elements from these vector<> types
    // once they are filled.
    vector<peak> data; // Raw centroid data peaks 
    vector<float> retentionTimes; 
    _2dtree_inplace<peak> *rawroot; // Index for raw peak data. 

    vector<ms2peak> data2; // Raw MS/MS peaks.
    _2dtree<ms2peak> *rawroot2; // Index for raw MS/MS peak data. 

    // And those that did not survive grouping.
    vector<peak> filtered; // New peaks should not be added or removed to preserve pointers.
    _2dtree_inplace<peak> *root; // Index for filtered peak data.

    vector<ms2peak> filtered2; // MS/MS peaks after filtering (do not add or remove to preserve pointers).
    _2dtree<ms2peak> *root2; // Index for filtered MS/MS data.

    vector<xic> xics; // XICs = peaks that occur at the same m/z value over time
    _2dtree<xic> *xics_all; // (m/z) Index used for querying XICs.

    vector<isotope_t> iso_labels; // isotope labels used in this data set, 
    // these are indexed by labels in isotope_group and isotope_triple

    vector<isotope_group> isotope_data; // Grouped XICs for isotope data.
    _2dtree<isotope_group> *isotope_idx; // (m/z,rt) index for querying isotope groups. 

    // Build a 2d-tree, a spatial index of each of the peak values. 
    void indexraw() { rawroot = new _2dtree_inplace<peak>(data); rawroot2 = new _2dtree<ms2peak>(data2); }

    // Build a 2d-tree of filtered peak data and create a map from a
    // memory location to an index in the filtered peak data.
    void indexfiltered() { root = new _2dtree_inplace<peak>(filtered);  root2 = new _2dtree<ms2peak>(filtered2);  }

    // Query a rectangular region of raw centroided peak data.
    void queryraw(vector<peak *> &results, double rtmin, double rtmax, double mzmin, double mzmax) {
      if(rawroot == NULL) return;
      region<peak> query(rtmin, rtmax,mzmin, mzmax);
      rawroot->Query(query,results);
    }

    // Query a rectangular region of MS/MS centroided peak data.
    void queryraw2(vector<ms2peak *> &results, double rtmin, double rtmax, double mzmin, double mzmax) {
      if(rawroot2 == NULL) return;
      region<ms2peak> query(rtmin, rtmax,mzmin, mzmax);
      rawroot2->Query(query,results);
    }

    // Query a rectangular region of filtered peak data. 
    void queryfiltered(vector<peak *> &results, double rtmin, double rtmax, double mzmin, double mzmax) {
      if(root == NULL) return;
      region<peak> query(rtmin, rtmax,mzmin, mzmax);
      root->Query(query, results);
    }

    // Query a rectangular region of filtered peak data. 
    void queryfiltered2(vector<ms2peak *> &results, double rtmin, double rtmax, double mzmin, double mzmax) {
      if(root2 == NULL) return;
      region<ms2peak> query(rtmin, rtmax, mzmin, mzmax);
      root2->Query(query, results);
    }

    // Remove noise peaks from LC-MS data.
    double histogram_threshold();

    //void filter(double log10I_thresh);
    void get_range(float &startTime, float &endTime, float retentionTime);
    void filter_new(double log10I_thresh);

    // Performs DFS to label MS peaks into connected components/XICs.
    void dfs_label(int label, peak *start, vector<peak *> &L);

    double compute_log2ratio(xic *light_xic, xic *heavy_xic);
  
    // Returns edges adjacent to given point based on range query.
    //void out_edges( peak* p,  vector< peak* > &edges );
    void out_edges_new( peak* p,  vector< peak* > &edges );

    // Return edges from the peak graph that are in a given query region.
    // NOTE: This code will return duplicate edges.
    void querypeakgraph( vector< pair<peak*,peak*> > &edges, 
			 double rtmin, double rtmax, double mzmin, double mzmax) {
      if(root == NULL) return;
      // Get the points within query region.
      vector<peak *> results; 
      region<peak> query(rtmin, rtmax, mzmin, mzmax);

      root->Query(query, results);
      // Look at outgoing edges in each result peak and return all edges.
      vector <peak *> outE;
      for(size_t r = 0; r < results.size(); r++){
	outE.clear();
	out_edges_new(results[r], outE);
	// Note some repetative edges are returned, but this is OK.
	for(size_t j = 0; j < outE.size(); j++) edges.push_back(pair<peak*,peak*>(results[r], outE[j]));
      }
    }

    // Applies filtering algorithm to XICs.
    void filter_xics();

    // Generates a clean chromatogram from chromatogram peaks in XIC
    // (no duplicate peaks from same time point).
    void clean_chromatogram(xic &s);
    // Smooth chroamtogram (for better quantification). 
    void smooth_chromatogram(xic &s);

    // Compute connected components in peak graph to find XICs. 
    void findxics();

    // Return XICs within given query range. 
    void queryxics( vector<xic *> &mz_results, double rtmin, 
			  double rtmax, double mzmin, double mzmax ){
      if(xics_all == NULL) return;
      // Query XICs in window.
      region<xic> query(rtmin, rtmax,mzmin, mzmax);
      xics_all->Query(query, mz_results);
    }

    // Returns isotope groups in a given window.
    void queryIsotopeGroups(vector< isotope_group *> &groups, double rtmin, double rtmax, double mzmin, double mzmax) {
      if(isotope_idx == NULL) return;
      region<isotope_group> query(rtmin, rtmax,mzmin, mzmax);  
      isotope_idx->Query(query, groups);
    }

    // Perform trapezoidal integration to quantify a XIC.
    void trapezoidal_integration(xic &s);

    // Get overlapping XICs... ignores spurious overlaps. 
    xic *get_overlapped(xic *x, double isotopic_spacing, double tol);

    // Translate, in retention time, the entire LCMS run.
    void translate(double translation);

    // Apply nonlinear transformation to data.
    void nonlinear(vector<shift_t> &shifts);

    // Assigns MS/MS peak to an XIC, created using non-MS/MS peaks.
    void assign2xic();

    // Finds the maximum and minimum ranges for MS, MS/MS, retention time data.
    void data_range();

    // Computes the log10(intensity) range logImin_disp and logImax_disp by
    // histogram adjustment. These ranges which are used to set the coloring for displayed peaks.
    void histogram_adjustment();

    // Determines mono-iso XICs and their charge.
    void collect13C(vector<xic *> &carbon13, xic *x, int charge, double sign = 1);
    void monoiso_and_charge();

    // Apply isotope grouping algorithm.
    xic *find_isotope(xic *cur, double mass_shift, bool ignore_type = false);
    void add_triple(xic *light, xic *medium, xic *heavy, vector<int> &labels);
    void add_pair(xic *light, xic *heavy, vector<int> &labels, int nitrogens = 0);
    void add_pair(xic *light, xic *heavy, int nitrogens) { 
      vector<int> nolabels; add_pair(light, heavy, nolabels, nitrogens); 
    }
    void add_group(vector<xic *> &xics);
    void add_group_phase(vector<xic *> &xics, int phase);

    void find_isotope_pairs_and_triples(); 

    xic *find_15Npair(int &candcount, int &cnt_15N, xic *cur, fragmentDB *ms2db, float direction);
    void find_15Npairs(vector<ms2::fragmentDB *> &DB);

    int find_nKs(double &max_score, ms2::fragmentDB *curDB, vector<ms2peak> &ms2peaks);
    void find_stoich_patterns(vector<ms2::fragmentDB *> &DB);


    irun(int condIndex_, int repIndex_, int irunIndex_, irun_file &data, 
	 msconfig &conf_, bool labels_reversed_ = false) : conf(conf_) {
      description = data.description;
      filename = data.filename;
      condIndex = condIndex_; repIndex = repIndex_;  irunIndex = irunIndex_;
      labels_reversed = labels_reversed_;

      rawroot = NULL; root = NULL; root2 = NULL; rawroot2 = NULL; 
      xics_all = NULL;  isotope_idx = NULL;  
      processed_flag = false;

      if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE) {
	// Create list of possible isotope labels.
	iso_labels.push_back(isotope_t("none", 0)); // 0th label is the no-label label

	// Read through configuration information and load active isotope labels.
	for(size_t i = 0; i < conf.activeIsotopes.size(); i++) {
	  if(0 <= conf.activeIsotopes[i] && conf.activeIsotopes[i] < (int)conf.isotopes.size()) {
	    iso_labels.push_back(conf.isotopes[conf.activeIsotopes[i]]);
	  }
	}
      }      
    }

    // Main driver of all of the irun processing. Builds trees and finds features.
    void process(vector<ms2::fragmentDB *> &DB);

    // Recalibrate along time or m/z dimension using a given window size.
    void recalibrate(recalmode_t recalmode, double winSz, vector<ms2::fragmentDB *> &DB);

    // Free up memory used by the irun. 
    ~irun() { 
      if(rawroot != NULL) delete rawroot;  if(rawroot2 != NULL) delete rawroot2;
      if(root != NULL) delete root; if(root2 != NULL) delete root2;
      if(xics_all != NULL) delete xics_all;
      if(isotope_idx != NULL) delete isotope_idx;
    }

    // Returns pointers to all isotope groups in LC-MS/MS irun.
    void get_isotope_groups(vector<isotope_group *> &groups) {
      for(size_t i = 0; i < isotope_data.size(); i++) groups.push_back(&isotope_data[i]);
    }

    void get_monoiso_ms2(vector<ms2::ms2peak *> &peaks) { 
      for(size_t i = 0; i < xics.size(); i++) if(xics[i].type == xicMONOISO) xics[i].get_ms2(peaks); 
    }

    // Adjust the log2quant value for each XIC. Used for normalization.
    void adjust_log2_quant(double delta);

    // Get nearest recipricol shifts between current irun and reference.
    void get_shifts(irun *ref, vector<shift_t> &shifts, double dx_max);

    // Perform alignment to a given reference irun.
    void align(irun *ref);
  };
}

#endif  //  __lcms_hpp__
