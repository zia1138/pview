#ifndef __analysis_hpp__
#define __analysis_hpp__

#include "lcms.hpp"
#include "util.hpp"

namespace lcms {
  // Heirarchical file organization (see also irun_file in lcms.hpp).
  struct replicate_files {  // replicate sets contain many irun files
    string description;  
    // TODO: enzyme_t protease;  read off protease from directory?
    vector<irun_file> irunFiles;  // note each irun_file has a filename and and a description (see lcms.hpp) 
    bool label_reversed; // Flag designates isotope labels are reversed.
    replicate_files() { label_reversed = false; }
    bool operator < (const replicate_files &b) const { return description < b.description; }
    void sort_all() { sort(irunFiles.begin(), irunFiles.end()); }
  };
  struct condition_files {  // conditions contain many replicate sets
    string description;  vector<replicate_files> replicateFiles; 
    bool operator < (const condition_files &b) const { return description < b.description; }
    void sort_all() {
      sort(replicateFiles.begin(), replicateFiles.end());
      for(size_t i = 0; i < replicateFiles.size(); i++) replicateFiles[i].sort_all();
    }
  };
  struct analysis_files {  // an analysis covers many conditions
    string description; vector<condition_files> conditionFiles; 
    void sort_all() {
      sort(conditionFiles.begin(), conditionFiles.end());
      for(size_t i = 0; i < conditionFiles.size(); i++) conditionFiles[i].sort_all();
    }
  };
  
  // A single replicate consists of one or more iruns.  In isotope
  // mode each replicate is ratio normalized if requested.
  struct replicate {
    int condIndex, repIndex; 
    vector<irun *> iruns; // A replicate set has a several iruns.
    string description;
    bool processed_flag;
    msconfig &conf;
    double log2_ratio_median; // loading error for isotope data

    // In a replicate set, allocate slots for iruns.
    replicate(int condIndex_, int repIndex_, replicate_files &data, msconfig &conf_) : conf(conf_) {
      condIndex = condIndex_; repIndex = repIndex_; description = data.description;
      processed_flag = false;
      bool label_rev = data.label_reversed;
      for(int irunIndex = 0; irunIndex < (int)data.irunFiles.size(); irunIndex++) {
	lcms::irun *s = new lcms::irun(condIndex, repIndex, irunIndex, data.irunFiles[irunIndex], conf, label_rev);
	iruns.push_back(s);
      }
      log2_ratio_median = 0;
    }

    // Free up memory used by iruns.
    ~replicate(){ for(unsigned int s = 0; s < iruns.size(); s++) delete iruns[s]; }

    // Main driver of replicate set processing.
    void process(vector<ms2::fragmentDB *> &DB) {
      if(conf.run_algorithms == false) return; 
      if(processed_flag) return; else processed_flag = true;
      // Process each irun (if it hasn't been processed already) in the replicate set first.
      for(size_t s = 0; s < iruns.size(); s++) iruns[s]->process(DB);
      // If requested, adjust ratios to be centered around median value of zero.
      if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE) normalize_isotope_ratios(); 
    }

    // Normalizes H/L ratios so their log2 values are centered around zero.
    // TODO: Also normalize H/M and M/L?!?!?
    void normalize_isotope_ratios() {
      vector<isotope_group *> groups;
      for(size_t i = 0; i < iruns.size(); i++) iruns[i]->get_isotope_groups(groups);
      if(groups.size() == 0) return;

      // Isotope mode normalizes ratios to be centered around zero.
      vector<double> ratios; 
      for(size_t i = 0; i < groups.size(); i++) { 
	// Ignore groups obtained by search.
	if(groups[i]->by_search == false) ratios.push_back(groups[i]->log2ratioHL);
      }
      if(ratios.size() > 0) {
	if(conf.skip_normalize_ratios == false) {
	  log2_ratio_median = util::median(ratios);
	  cout << description << ":log2_ratio_median=" << log2_ratio_median << endl;
	  // Now, subtract median from all log-2 ratios, if requested.
	  for(size_t i = 0; i < groups.size(); i++) {
	    groups[i]->log2ratioHL -= log2_ratio_median;
	    groups[i]->log2xicH -= 0.5 * log2_ratio_median;  groups[i]->log2xicL += 0.5 * log2_ratio_median;
	  }
	}
      }
    }
    // Accessor functions that return and count iruns in the replicate set.
    void get_iruns(vector<irun *> &irun_set) { for(size_t i = 0; i < iruns.size(); i++) irun_set.push_back(iruns[i]); }
    irun *getIrun(int irunIndex) { return iruns.at(irunIndex); }
    int numIruns() { return iruns.size(); }
  };

  // A condition consists for one or more replicates sets. 
  struct condition {
    int condIndex; // Index in the vector of conditions in the analysis object in analysis.hpp.
    msconfig &conf; // Global configuration information.
    string description; // Description of replicate set.
    vector<replicate *> replicates; 
    bool processed_flag; // set to true when processing has been completed

    condition(int condIndex_, condition_files &data, msconfig &conf_) : conf(conf_) { 
      processed_flag = false;
      // Save index and description.
      condIndex = condIndex_; description = data.description; 
      // Allocate space for replicate.
      for(int repIndex = 0; repIndex < (int)data.replicateFiles.size(); repIndex++) {
	lcms::replicate *rep = new lcms::replicate(condIndex, repIndex, data.replicateFiles[repIndex], conf);
	replicates.push_back(rep);
      }
    }

    // Free up memory used by iruns.
    ~condition(){ for(size_t r = 0; r < replicates.size(); r++) delete replicates[r];  }

    void process(vector<ms2::fragmentDB *> &DB) {
      if(conf.run_algorithms == false) return;
      if(processed_flag) return; else processed_flag = true;
      // Process each replicate set (of iruns) if it hasn't been processed already.
      for(size_t r = 0; r < replicates.size(); r++) replicates[r]->process(DB);
    }

    // Accessor functions for replicate sets and iruns. 
    void get_iruns(vector<irun *> &irun_set) { for(size_t i = 0; i < replicates.size(); i++) replicates[i]->get_iruns(irun_set); }
    irun *getIrun(int repIndex, int irunIndex) { return replicates.at(repIndex)->getIrun(irunIndex);  }
    int numReps() { return replicates.size(); }
    replicate *getRep(int repIndex) { return replicates.at(repIndex); }
    int numIruns(int repIndex) { return replicates.at(repIndex)->numIruns(); }
  };


  // Isotope groups grouped based on their sequence information
  struct isop_seqgroup {
    string seq; int Nmods; // sequence + number of modifications
    vector<isotope_group *> groups; // actual groups with the same sequence
    // Returns ratios from sequence group.
    void get_ratios(ratio_type rt, vector<double> &res, int condIndex, int repIndex = -1) {
      for(size_t p = 0; p < groups.size(); p++) {
	isotope_group *group = groups[p];
	if(group->by_search) continue; 
	if( (repIndex == -1 || groups[p]->repIndex == repIndex) && groups[p]->condIndex == condIndex) {
	  if(rt == isoHL && group->heavy && group->light )  res.push_back(groups[p]->log2ratioHL); 
	  if(rt == isoHM && group->heavy && group->medium ) res.push_back(groups[p]->log2ratioHM); 
	  if(rt == isoML && group->medium && group->light ) res.push_back(groups[p]->log2ratioML); 
	}
      }
    }
  };

  // Protein group and associated isotope groups.
  struct isop_pgroup {
    vector<int> acc;      // Indexes in ms2::fragmentDB::protein_meta
			  // assigned to this group (protein names).
    vector<isotope_group *> groups; // Isotope groups in this protein group

    // Isotope groups organized by sequence. 
    vector<isop_seqgroup> seq_groups, seq_groups_shared;
    
    // Separates isotope groups based on their sequence.
    void build_seq_groups(fasta_data &fasta);
  };

  // XICs grouped across iruns. Used for nonlinear alignment label-free quantification.
  struct align_group {
    double mz, retentionTime, start, end;  // mean (m/z, retention time, start and end retenion times)
    vector<xic*> xics; // XICs grouped across replicates
    // Position using centroid of all XICs.
    void position() {
      mz = retentionTime = start = end = 0;
      double N = (double) xics.size();
      for(size_t i = 0; i < xics.size(); i++) {
	mz += xics[i]->mz / N;	retentionTime += xics[i]->retentionTime / N;
	start += xics[i]->start / N; end += xics[i]->end / N;
      }
    }
    // Accessors used by _2dtree<> 
    double x() const { return retentionTime; } double y() const { return mz; }
    // Return number of XICs in this group.
    int countxics() { return xics.size(); }
    // Return the total number of MS/MS IDs in all grouped XICs.
    int num_ids() { int ids = 0; for(size_t i = 0; i < xics.size(); i++) ids += xics[i]->num_ids(); return ids; }    
    int num_ms2() { int ms2cnt = 0; 
      for(size_t i = 0; i < xics.size(); i++) ms2cnt += xics[i]->ms2peaks.size();
      return ms2cnt; 
    }
    // Returns maximum score from MS/MS peaks. 0 if none.
    double max_score() {
      double s = 0; 
      for(size_t i = 0; i < xics.size(); i++) s = max(s, xics[i]->max_score());
      return s;
    }
    // Returns first MS/MS peak corresponding to the best score.
    xic *max_xic() {
      double score = max_score();
      for(size_t i = 0; i < xics.size(); i++) 
	if(xics[i]->num_ids() > 0 && score - xics[i]->max_score() == 0) return xics[i];       
      return NULL;
    }
    // Collects the XICs that have been grouped in this replicate set.
    void getxics(vector<xic *> &xics_output) { xics_output = xics; }
    // Returns the MS/MS peak with the maximum score.
    ms2peak *max_ms2() { xic *x = max_xic();  if(x == NULL) return NULL; else return x->max_ms2();  }
    // Count XICs with given condition index and replicate set index.
    int countIruns(int condIndex, int repIndex) {
      int count = 0;
      for(size_t i = 0; i < xics.size(); i++) 
	if(xics[i]->condIndex == condIndex && xics[i]->repIndex == repIndex) count++;
      return count;
    }
    // Returns most frequently occuring charge in the alignment group
    int charge() {
      int max_charge = 0;
      for(size_t i = 0; i < xics.size(); i++) if(xics[i]->charge > max_charge) max_charge = xics[i]->charge;
      vector<int> counts(max_charge+1, 0);
      for(size_t i = 0; i < xics.size(); i++) counts.at(xics[i]->charge)++;
      int freq_charge = 0, max_counts = 0;
      for(int cur_charge = 1; cur_charge < (int)counts.size(); cur_charge++) {
	if(counts[cur_charge] > max_counts) {
	  freq_charge = cur_charge;
	  max_counts = counts[cur_charge];
	}
      }
      return freq_charge;
    }
  };

  // Protein groups and corresponding XICs grouped across iruns after
  // nonlinear alignment and grouping xic_pgroup -> align_group -> xic
  struct align_pgroup {
    vector<int> acc;  // Indexes in ms2::fragmentDB::protein_met
		      // assigned to this protein group.

    // NOTE: each align group has the same sequence/charge.
    vector<align_group *> align_groups;     // Select "best" align group to report value for protein.
    vector<align_group *> shared_align_groups;
    // NOTE: Each of these align_groups maps to only one protein group. It is not shared.
  };

  // All of these XICs are grouped based on their MS/MS id.
  struct seq_group {
    string seq; 
    int charge; // charge state of all XICs
    double max_ms2score; // maximum MS/MS score
    vector<xic *> xics; 
    bool operator == (string &seqB) const { return seq == seqB; }
    // Used to sequence groups. Lower rank sequnce group is larger has higher max_ms2score.
    bool operator < (const seq_group &b) const { 
      if(xics.size() == b.xics.size()) return max_ms2score > b.max_ms2score; else return xics.size() > b.xics.size();
    }
    int numxics() { return xics.size(); }
    void getxics(vector<xic *> &xics_out) { for(size_t i = 0; i < xics.size(); i++) xics_out.push_back(xics[i]);  }
    // Get XIC with the maximum MS/MS score.
    xic *max_ms2_xic() {
      xic *max_xic = NULL;
      double max_score = 0;
      for(size_t i = 0; i < xics.size(); i++) {
	if(xics[i]->max_score() > max_score) {
	  max_score = xics[i]->max_score();
	  max_xic = xics[i];
	}
      }
      return max_xic;
    }

    // Add XIC to sequence group, avoding multiple XICs per irun.
    void addxic(xic *x) {
      bool found = false;
      // 
      for(size_t i = 0; i < xics.size(); i++) {
	if(xics[i]->condIndex == x->condIndex) {
	  // Keep XICs with the highest MS/MS score if there are
	  // multiple from the same condition.
	  if(xics[i]->max_score() < x->max_score()) {
	    xics[i] = x;
	  }
	  found = true;
	  break;
	}
      }
      if(!found) xics.push_back(x);
    }
    int numConds() {
      // Find all unique condition indexes and count.
      vector<int> idxs;
      for(size_t i = 0; i < xics.size(); i++) {
	int condIndex = xics[i]->condIndex;
	bool found = false;
	for(size_t k = 0; k < idxs.size(); k++) if(idxs[k] == condIndex) found = true;
	if(!found) idxs.push_back(condIndex);
      }
      return idxs.size();
    }
  };

  // For a protein group quantification pick best unique XIC.
  // xic_pgroup -> sequence and charge -> xic
  struct xic_pgroup {
    vector<int> acc; // NOTE: Computing coverage only make sense for a single accession.
    vector<xic *> xics;  // ID'd XICs in this protein group
    vector<xic *> xics_shared;
    // Protein group level quantification from the largest seq_group
    // with the highest median MS/MS score.
    vector<seq_group> seq_groups;  // XICs grouped based on their modified sequence.
    // TODO: quantify off seq_groups[0] after sorting
    vector<seq_group> seq_groups_shared;

    void build_seq_groups(fasta_data &fasta, int min_cond) {
      _build_seq_groups(fasta, seq_groups, xics, min_cond);
      _build_seq_groups(fasta, seq_groups_shared, xics_shared, min_cond);
    }
    void _build_seq_groups(fasta_data &fasta, vector<seq_group> &seq_groups, vector<xic *> &xics, int min_cond);
  };

  // Contains several conditions that make up the analysis. Aligns
  // conditions to one another and groups replicate observations.
  struct analysis {
    string description; // Description of analysis (directory name usually).
    msconfig &conf; // Global configuration information.
    bool processed_flag; // Set to true if processing has been completed.

    // Range of MS and MS/MS peaks.
    float min_mz, max_mz, min_time, max_time, min_mz2, max_mz2, logImin_disp, logImax_disp;

    vector<align_group> align_groups;  
    _2dtree <align_group> *root;

    vector<condition *> conditions; // Conditions
    vector<irun *> iruns; // Pointers to individual iruns in replicate/fraction sets.
                          // Allows multithreaded processing of each individual irun.

    // See the ms2.cpp: dbtype_t enum. There should be one for each or null if they
    // don't exist. 
    vector<ms2::fragmentDB *> DB; // MS/MS fragment databases

    // One global location for fasta data, as opposed to having one copy for each database.
    vector<string> fasta_descr, fasta_files;  // List of FASTA files containing amino acid sequence data.
    vector<string> pepxml; // List of PepXML files containing external MS/MS IDs.

    // Global sequence data and any external DB search information.
    ms2::fasta_data fasta;

    // Protein group data for isotope quantification.
    vector<isop_pgroup>   isop_pgroups_data; // Protein group data.
    vector<isop_pgroup *> isop_pgroups; // Protein groups sorted on log2ratio.

    // Protein groups for alignment based and XIC-based quantification.
    vector<align_pgroup> align_pgroups;  // Label-free nonlinear alignment protein group data.

    // Note: XIC-based quantification cross references on sequence.
    vector<xic_pgroup> xic_pgroups;  // XIC-based quantificiation protein groups.

    // Provide set of mzXML files, grouped with description information
    // and grouped by condition name.
    analysis(msconfig &conf_, analysis_files &data, vector<string> &fasta_descr_, 
	     vector<string> &fasta_files_, vector<string> &pepxml_);

    ~analysis() {
      for(size_t s = 0; s < conditions.size(); s++) delete conditions[s];
      for(size_t i = 0; i < DB.size(); i++) if(DB[i] != NULL) delete DB[i];
      if(root != NULL) delete root;
    }

    void query(vector<align_group *> &results, double rtmin, double rtmax, double mzmin, double mzmax) {
      if(root == NULL) return;
      region<align_group> query(rtmin, rtmax,mzmin, mzmax);
      root->Query(query,results);
    }

    // Loads FASTA files and builds protein MS/MS databases.
    void load_protein_db() {
      if(pepxml.size() == 0)  fasta.load_fasta(fasta_descr, fasta_files);
      DB[dbLIGHT] = new ms2::fragmentDB(conf, fasta, dbLIGHT); 
      DB[dbLIGHT]->build(); // Build fragment database.
      if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE || 
	 conf.quantmode == quantAcetylSTOICH) {
	DB[dbHEAVY] = new ms2::fragmentDB(conf, fasta, dbHEAVY);
	DB[dbHEAVY]->build(); 
      }
      if(conf.quantmode == quantISOTRIPLE) {
	DB[dbMEDIUM] = new ms2::fragmentDB(conf, fasta, dbMEDIUM);
	DB[dbMEDIUM]->build();
      }
      
    }

    void process();     // Apply full data set processing (e.g. protein grouping).
    int num_iruns() { return iruns.size(); }     // Get total number of iruns.

    // Return index of reference irun. This one has the largest number of XICs. 
    int selectAlignReference() {
      int ref = -1;
      size_t n_xics = 0;
      for(int s = 0; s < (int)iruns.size(); s++) {
	if(iruns[s]->xics.size() > n_xics) {
	  ref = s;
	  n_xics = iruns[s]->xics.size();
	}
      }
      return ref;
    }

    void align_irun(int irun, int ref) { iruns.at(irun)->align(iruns.at(ref)); }

    // Process an individual irun (from ungrouped set). Enables multi-threaded processing.
    void process_irun(int idx) {  iruns.at(idx)->process(DB); }

    // Recalibrates an instrument run using confident search IDs.
    void recalibrate_irun(recalmode_t recalmode, double winSz, int idx) {  
      iruns.at(idx)->recalibrate(recalmode, winSz,  DB); 
    }

    void saveMassErrorsMS1(string filename); // save mass error histogram for precursor mass MS1
    void saveMassErrorsMS2(string filename); // save mass error histogram for MS2 product ion mass b or y ions

    void save_ms2peak_header(ofstream &csv); // save MS/MS peak information
    void save_ms2peak(ofstream &csv, ms2peak *peak);

    void save_iso_header(ofstream &csv, char eol_chr = '\t');
    void save_iso(ofstream &csv, isotope_group *group, char eol_char = '\t');

    // Nonlinear alignment based quantification.
    void form_align_group(xic *x, align_group &group);
    bool is_valid_align_group(align_group &group);
    void align_mode_quant();  
    void saveAlignTable(string filename); // Nonlinear alignment based table.

    // TODO: Peptide level output needs to be improved. Need to show mapping to multiple protein groups.
    // TODO: Peptide level output needs to be improved. Need to show mapping to multiple protein groups.
    // ZK: As it stands it isn't quite right.
    void saveAlignTablePeptides(string filename); // per peptide output
    void saveAlignTableInternal(string filename); // internal ratio correlations for every replicate
    void saveAlignGroupHeader(ofstream &csv); // save header information for alignment grops
    void saveAlignGroup(ofstream &csv, align_group *group);

    // XIC-based Label-free quantification by cross referncing XICs.
    void xic_based_quant();
    void saveXICBasedTable(string filename);         // Save a per-protein table.
    // TODO: Peptide level output needs to be improved. Need to show mapping to multiple protein groups.
    // TODO: Peptide level output needs to be improved. Need to show mapping to multiple protein groups.
    // ZK: As it stands it is not quite right.
    void saveXICBasedTablePeptides(string filename); // Save a per-peptide table.
    void saveXICBasedTableInternal(string filename); // Save internal correlations between non-overlapping peptides.

    void saveXIC_allids(string filename);

    void saveSeqGroupHeader(ofstream &csv); 
    void saveSeqGroup(ofstream &csv, seq_group &group);

    // Output rows in label-free quantification mode.
    void saveLabelFreeHeader(ofstream &csv);
    void saveLabelFreeRow(ofstream &csv, vector<xic *> &xics);     

    // Label-free normalization.
    void label_free_normalize(); 

    // Use for biomarker discovery, MS1 only.
    void saveLabelFreeMS1(string filename);

    // Isotope group output (valid for protein based and peptide based quantification).
    void save_isotope_data(string filename);
    void save_isotope_groups(string filename);          // Per-protein isotope output.
    void save_isotope_groups_internal(string filename); // Generate internal correlations per condition and per replicate set. 
    void save_isotope_groups_internal_pep(string filename); // Generate interal correlations on a per-peptide basis comparing +2 to >=+3
    void save_isotope_groups_corr(string filename); // Generate per replicate correlations of the data set.
    void save_isotope_groups_table(string filename); // Generate tabular output of isotope ratio data (protein-based output)

    void save_stoich_data(string filename); 

    // TODO: Peptide level output needs to be improved. Need to show mapping to multiple protein groups.
    // TODO: Peptide level output needs to be improved. Need to show mapping to multiple protein groups.
    // ZK: As it stands it is not quite right.
    void save_isotope_groups_table_pep(string filename);  // Save per peptide quantification table.

    // Used for sequence group internal correlations.
    void save_internal_pep_seq_group(ofstream &, int, int , string &, vector<isop_seqgroup> &); 
    int isotope_table_header(ofstream &csv);
    void save_protein_group(ofstream &csv, vector<int> &acc, bool extra_data = false);
    void save_isotope_ratios(ofstream &csv, int Ncol, vector<isotope_group *> &groups);

    // Compute protein groups from ID'd isotope groups. Assign shared groups to largest group.
    void protein_grouping(vector<isotope_group *> &groups); // isotope label quantification
    void protein_grouping(vector<align_group *> &align_groups); // retention time alignment quantification
    void protein_grouping(vector<xic *> &xics); // xic sequence cross referecning quantification

    // Return pointers to all isotope groups in all data sets.
    void get_isotope_groups(vector<isotope_group *> &groups) { 
      for(size_t i = 0; i < iruns.size(); i++) iruns[i]->get_isotope_groups(groups); 
    }

    // FDR cutoff code. NOTE: split charge calls split length if
    // request if not does a simple FDR cutoff.
    void fdr_split_charge(vector<ms2peak *> &peaks); // <- call this first not length or cutoff!!!!!
    void fdr_cutoff(vector<ms2peak *> &peaks);

    // Clear IDs where maximum score points to two or more sequences.
    void clear_ambig_maxscore(vector<ms2peak *> &peaks);

    // Return total number of protein groups given mode.
    int numProteinGroups() { 
      if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE || conf.quantmode == quantAcetylSTOICH) return isop_pgroups.size();
      else if(conf.quantmode == quantALIGN) return align_pgroups.size();
      else if(conf.quantmode == quantXICBASED) return xic_pgroups.size();
      else return 0;
    }

    // Access different types of protein groups.
    isop_pgroup *isotopeGroupProteinGroup(int groupIdx) { return isop_pgroups.at(groupIdx); }
    align_pgroup *alignProteinGroup(int groupIdx) { return &align_pgroups.at(groupIdx); }
    xic_pgroup *xicProteinGroup(int groupIdx) { return &xic_pgroups.at(groupIdx); }

    // Return meta information for a specific protein index.
    string get_meta(int idx) { return fasta.protein_meta.at(idx); }
    string get_meta(ms2::fragment *f) { return get_meta(f->protein_idx); }

    // Return sequence information for an MS/MS peak.
    string get_seqS(ms2::ms2peak *p) {
      vector<fragment *> ids;  p->max_ids(ids);
      vector<mods_t> mods;     p->max_mods(mods);
      return fasta.get_seqS(ids[0], mods[0]);
    }

    // Returns sequence with site score.
    string get_seqS(ms2::fragment *f, mods_t &m) { return fasta.get_seqS(f, m); }
    // Returns sequence w/ modification. This is critical for grouping
    // and removing ambigious localization, etc.
    string get_seqN(ms2::fragment *f, mods_t &m) { return fasta.get_seqN(f, m); }

    // Accessor functions for conditions, replicates, and iruns.
    irun *getIrun(int condIndex, int repIndex, int irunIndex) { return conditions.at(condIndex)->getIrun(repIndex, irunIndex); }
    int numConds() { return conditions.size(); }
    int numReps(int condIndex) { return conditions.at(condIndex)->numReps(); }
    int numIruns(int condIndex, int repIndex) { return conditions.at(condIndex)->numIruns(repIndex); }
    int numIruns(int condIndex) {
      int iruns = 0; 
      for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) 
	iruns += conditions.at(condIndex)->numIruns(repIndex);
      return iruns;
    }
      
    condition *getCond(int condIndex) { return conditions.at(condIndex); }
    replicate *getRep(int condIndex, int repIndex) { return conditions.at(condIndex)->getRep(repIndex); }

    // TODO: Confirm this is OK, especially repIndex positioning.
    // Maps to column when irun index is unique for a replicate run.
    int toColumn(int condIndex, int repIndex, int irunIndex) {
      int column = 0;
      for(int c = 0; c < condIndex; c++) column += numIruns(c);
      for(int r = 0; r < repIndex; r++)  column += numIruns(condIndex, r);
      // Last, add irunIndex offset.
      return column + irunIndex;
    }

    // Column output for conditions + replicate sets.  Here all
    // irunIndex numbers belong to one replicate set (indexed by
    // repIndex).
    int toColumn(int condIndex, int repIndex) {
      int column = 0; for(int c = 0; c < condIndex; c++) column += numReps(c);  return column + repIndex;
    }
  };
}

#endif //  __analysis_hpp__
