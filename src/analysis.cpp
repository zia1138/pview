#include <algorithm>
#include <iostream> 
#include <iomanip>
#include <fstream>
#include <math.h>
#include <map>
#include <set>

#include "analysis.hpp"
#include "xml.hpp"
#include "util.hpp"

namespace lcms { 
  // Setup data set analysis.  
  analysis::analysis(msconfig &conf_, analysis_files &data, 
		     vector<string> &fasta_descr_, 
		     vector<string> &fasta_files_, vector<string> &pepxml_) : conf(conf_), fasta(conf_) { 
    description = data.description;
    fasta_descr = fasta_descr_; fasta_files = fasta_files_; pepxml = pepxml_;

    // Set data structure pointers to null.
    root = NULL; processed_flag = false;
    DB.resize(ms2::NumDBTypes, NULL);
    
    // Iterate through each condition allocating iruns and replicates.
    // Assign a condition/replicate number index, and use prefix of
    // mzXML file to assign name to irun.
    for(int condIndex = 0; condIndex < (int)data.conditionFiles.size(); condIndex++) {
      // Allocate space for replicate set.
      condition *r = new condition(condIndex, data.conditionFiles[condIndex], conf);
      // Save replicate set in analysis object.
      conditions.push_back(r);    
    }

    // Collect pointers to all iruns used by process_irun().
    for(size_t i = 0; i < conditions.size(); i++) conditions[i]->get_iruns(iruns);

    // Initialize data ranges for all of the contained LC-MS/MS runs. 
    min_mz = numeric_limits<float>::max();   max_mz   = -(numeric_limits<float>::max());
    min_time = numeric_limits<float>::max(); max_time = -(numeric_limits<float>::max());
    min_mz2 = numeric_limits<float>::max();  max_mz2 = -(numeric_limits<float>::max());
    logImin_disp = 3; logImax_disp = 6;
  }

  
  // Used to sort isotope protein groups by the number of measurements
  // in that group.
  struct groups_size_cmp { bool operator () (const isop_pgroup *a, const isop_pgroup *b){ return a->groups.size() > b->groups.size();} };

  // Applies protein grouping to isotope groups. Handles shared peptides
  // by assigning them to largest protein group. Sorts gorups by
  // median log2 ratio.
  void analysis::protein_grouping(vector<isotope_group *> &isogroups) {
    // Get all (highest scoring) MS/MS peaks present in the groups.
    // Note ms2peaks is in one-to-one correspondence with groups.
    // There is a NULL value if the group does not have a peak.
    vector<ms2::ms2peak *> ms2peaks(isogroups.size());
    for(size_t i = 0; i < isogroups.size(); i++) ms2peaks[i] = isogroups[i]->max_ms2();

    // Construct protein groups and assign MS/MS peaks to these groups.
    vector< vector<int> > pgroup2idx;
    ms2::build_protein_groups(fasta.Nproteins(), ms2peaks, pgroup2idx);

    // Computes mapping from a MS/MS spectrum to a list of groups it belongs to.
    vector< vector<int> > peak2groups(ms2peaks.size());
    ms2::assign_protein_groups(fasta.Nproteins(), ms2peaks, pgroup2idx, peak2groups);

    // Compute mapping from group number to a list of unique groups.
    // Leverages one-to-one correspondence between ms2peaks and groups.
    vector< vector<lcms::isotope_group *> > pgroup2group(pgroup2idx.size());
    for(size_t i = 0; i < ms2peaks.size(); i++) {
      if(ms2peaks[i] == NULL) continue; // skip groups that have no assigned MS/MS peaks
      vector<int> &groups = peak2groups[i];
      for(size_t g = 0; g < groups.size(); g++) {
	if(groups[g] < 0) continue; // To a non-conclusive group, skip.
	pgroup2group.at(groups[g]).push_back(isogroups[i]);
      }
    }
      
    // Organize protein group access and isotope group data.
    isop_pgroups_data.resize(pgroup2idx.size());
    for(size_t i = 0; i < pgroup2idx.size(); i++) {
      // Get accession number information and isotope groups for protein group.
      isop_pgroup &data = isop_pgroups_data[i];
      data.acc   = pgroup2idx[i]; 
      data.groups = pgroup2group[i];
    }

    // Segregate IDs in a protein group into smaller groups based on sequence.
    for(size_t i = 0; i < isop_pgroups_data.size(); i++) isop_pgroups_data[i].build_seq_groups(fasta);

    // Sort pointers to isotope group data by the number of groups in that group.
    isop_pgroups.resize(isop_pgroups_data.size()); // Actually sort memory pointers.
    for(size_t i = 0; i < isop_pgroups_data.size(); i++) isop_pgroups[i] = &isop_pgroups_data[i];
    sort(isop_pgroups.begin(), isop_pgroups.end(), groups_size_cmp());
  }

  struct sort_align_group {
    // Sort descending on size of aligned XIC groups.  If size is the
    // same, use maximum MS/MS score to sort.
    bool operator () (align_group *a, align_group *b) const { 
      int xicA = a->countxics(), xicB = b->countxics();
      if(xicA == xicB) return a->max_score() > b->max_score(); else return xicA > xicB;
    } 
  };

  // Applies protein grouping to aligned XIC groups (grouped across iruns in alignment mode).
  void analysis::protein_grouping(vector<align_group *> &align_groups) {
    // Get maximum scoring MS/MS ID in each XIC group.
    vector<ms2::ms2peak *> ms2peaks(align_groups.size());
    for(size_t i = 0; i < align_groups.size(); i++) ms2peaks[i] = align_groups[i]->max_ms2();

    // Designate protein groups.
    vector< vector<int> > pgroup2idx; // map from aprotein group to list of protein meta indexes in DB[dbLIGHT]
    ms2::build_protein_groups(fasta.Nproteins(), ms2peaks, pgroup2idx);

    // Use map from protein index to group to assign MS/MS peaks to groups.
    vector< vector<int> > peak2groups(ms2peaks.size());
    ms2::assign_protein_groups(fasta.Nproteins(), ms2peaks, pgroup2idx, peak2groups);

    // Compute mapping from group number to a list of unique
    // groups.  Leverages one-to-one correspondence between
    // ms2peaks and groups.
    align_pgroups.resize(pgroup2idx.size());
    for(size_t i = 0; i < pgroup2idx.size(); i++) align_pgroups[i].acc = pgroup2idx[i];

    // Assign align groups to their corresponding protein groups.
    for(size_t i = 0; i < ms2peaks.size(); i++) {
      if(ms2peaks[i] == NULL) continue;
      vector<int> &groups = peak2groups[i];
      // Keep only XIC groups that map to only one protein group.
      if(groups.size() == 1) 
	align_pgroups.at(groups[0]).align_groups.push_back(align_groups[i]);
      else if (groups.size() > 1) {
	for(size_t g = 0; g < groups.size(); g++) {
	  if(groups[g] < 0) continue; // To a non-conclusive group.
	  // Assign shared sequence
	  align_pgroups.at(groups[g]).shared_align_groups.push_back(align_groups[i]);
	}
      }
    }
    // Sort aligned XIC groups based on their size and MS/MS score.
    // NOTE: align_groups[0] is the best aligned XIC group and is used
    // for protein level quantification.
    for(size_t i = 0; i < align_pgroups.size(); i++) {
      sort(align_pgroups[i].align_groups.begin(), align_pgroups[i].align_groups.end(), sort_align_group());
      sort(align_pgroups[i].shared_align_groups.begin(), align_pgroups[i].shared_align_groups.end(), sort_align_group());
    }
  }

  // Apply protein grouping to XICs.
  void analysis::protein_grouping(vector<xic *> &xics) {
    vector<ms2::ms2peak *> ms2peaks(xics.size());
    for(size_t i = 0; i < xics.size(); i++) ms2peaks[i] = xics[i]->max_ms2();

    // Designate protein groups.
    vector< vector<int> > pgroup2idx; // map from a protein group to list of protein meta indexes in DB[dbLIGHT]
    ms2::build_protein_groups(fasta.Nproteins(), ms2peaks, pgroup2idx);
    
    // Set accessions of XIC protein groups.
    xic_pgroups.resize(pgroup2idx.size());
    for(size_t i = 0; i < pgroup2idx.size(); i++) xic_pgroups[i].acc = pgroup2idx[i];

    // Use map from protein index to group to assign MS/MS peaks to groups.
    vector< vector<int> > peak2groups(ms2peaks.size());
    ms2::assign_protein_groups(fasta.Nproteins(), ms2peaks, pgroup2idx, peak2groups);

    for(size_t i = 0; i < ms2peaks.size(); i++) {
      if(ms2peaks[i] == NULL) continue; // skip XICs that have no IDs.
      vector<int> &groups = peak2groups.at(i);
      if(groups.size() == 0) continue; // no group assignment, skip
      if(groups.size() == 1) {
	// XIC maps to a single group, assign.
	xic_pgroups.at(groups[0]).xics.push_back(xics.at(i));
      }
      else {
	// Assign shared XICs to XIC groups.
	for(size_t g = 0; g < groups.size(); g++) {
	  if(groups[g] < 0) continue; // non-conclusive group, skip
	  xic_pgroups.at(groups[g]).xics_shared.push_back(xics.at(i));
	}
      }
    }
  }

  // First split MS/MS peaks on charge.
  void analysis::fdr_split_charge(vector<ms2peak *> &peaks) {
    vector<ms2peak *> peaks2, peaks3;
    for(size_t i = 0; i < peaks.size(); i++) {
      // Skip peaks with NO IDs!!!!
      if(peaks[i]->num_ids() == 0) continue;
      if(peaks[i]->charge <= 2) peaks2.push_back(peaks[i]);
      else peaks3.push_back(peaks[i]);
    }
    cout << "charge 1 2" << endl;
    fdr_cutoff(peaks2);
    cout << "charge >= 3" << endl;
    fdr_cutoff(peaks3);
    cout << endl;
  }

  // Applies FDR cutoff to all of the data.
  void analysis::fdr_cutoff(vector<ms2peak *> &peaks) {
    if(peaks.size() == 0) return;
    double cutoff = ms2::fdr_threshold(peaks, conf.ms2_fdr);    
    for(size_t i = 0; i < peaks.size(); i++) peaks[i]->score_cutoff(cutoff);
  }

  // Assures maximum scoring sequence is not ambigious.  Some MS/MS
  // spectra may get the same maximum score for a given
  // sequence. Think luecine <-> isoluecine or a fragmentation
  // spectrum which has been filtered too much.
  void analysis::clear_ambig_maxscore(vector<ms2peak *> &peak) {
    int ambig_count = 0;
    for(size_t i = 0; i < peak.size(); i++) {
      ms2peak *p = peak[i]; 
      if(p == NULL) continue; if(p->num_ids() == 0) continue;

      // Get all maximum scoring IDs and modifications.
      vector<fragment *> ids; p->max_ids(ids); 
      vector<mods_t> mods; p->max_mods(mods);

      // Get 0th sequence with modification applied.  Note that using
      // get_seqN I use the precursor mass to filter out any possible
      // ambigious modification site localizations.
      string seq = get_seqN(ids[0], mods[0]);
      
      // Confirm there are no other ambigious maximum scoring sequences.
      bool ambig = false;
      for(size_t i = 1; i < ids.size(); i++) {
	// Get ith sequence to check.
	string seq_check = get_seqN(ids[i], mods[i]);
	if(seq_check != seq) { 
	  cout << seq_check << " ambig " << seq << endl;
	  ambig = true; break; 
	}
      }
      if(ambig) { 
	ambig_count++; p->clear_ms2id(); 
      }
    }
    cout << "ambig_count = " << ambig_count << endl;
  }

  // XIC-based quantification by cross referencing XICs based on sequnece and charge across data sets.
  void analysis::xic_based_quant() {
    // Collect pointers to XICs that have IDs.
    vector<xic *> id_xics;
    for(size_t s = 0; s < iruns.size(); s++) {
      // Get XICS with IDs from each irun.
      vector<xic> &xics = iruns[s]->xics;
      for(size_t i = 0; i < xics.size(); i++) { if(xics[i].num_ids() > 0) id_xics.push_back(&xics[i]);  }
    }
    // Apply protein grouping to XICs (using max MS/MS score per XIC).
    protein_grouping(id_xics);
    // Within each protein group assign XICs to sequence groups.
    for(size_t i = 0; i < xic_pgroups.size(); i++)  xic_pgroups[i].build_seq_groups(fasta, conf.min_conds);

    // Keep only protein groups that have assigned sequence groups
    // that occur in a minimum number of conditions.
    vector<xic_pgroup> valid_pgroups;
    for(size_t i = 0; i < xic_pgroups.size(); i++) {
      // TODO: Figure this out.  Should we repeat the protein grouping
      // with XICs in sufficient number of conditions?  Do we loose
      // protein groups here????
      if(xic_pgroups[i].seq_groups.size() > 0) valid_pgroups.push_back(xic_pgroups[i]);
    }
    xic_pgroups = valid_pgroups;
 }

  // Splits isotope groups up by sequence.
  void isop_pgroup::build_seq_groups(fasta_data &fasta) {
    for(size_t i = 0; i < groups.size(); i++) {
      ms2::ms2peak *peak = groups[i]->max_ms2(); // Get best MS/MS peak with best ID.
      // Is this a shared/razor peptide?
      vector<int> uniq_idxs; peak->max_idxs(uniq_idxs);
      bool ms2_shared = uniq_idxs.size() > acc.size();

      // Get all maximum scoring IDs and modifications.
      vector<fragment *> ids; peak->max_ids(ids);
      vector<mods_t> mods; peak->max_mods(mods);

      // Get sequence information.  Note we can use mods[0] because
      // clear_ambig should have made all of these the same.
      string seq = fasta.get_seqN(ids[0], mods[0]);

      bool found = false;
      if(ms2_shared == false) { 
	for(size_t s = 0; s < seq_groups.size(); s++) {
	  if(seq_groups[s].seq == seq) { seq_groups[s].groups.push_back(groups[i]);  found = true; }
	}
      }
      else {
	for(size_t s = 0; s < seq_groups_shared.size(); s++) {
	  if(seq_groups_shared[s].seq == seq) { seq_groups_shared[s].groups.push_back(groups[i]);  found = true; }
	}
      }
      if(found == false) {
	isop_seqgroup group;
	group.seq = seq; group.Nmods = mods[0].count();
	group.groups.push_back(groups[i]);
	if(ms2_shared) seq_groups_shared.push_back(group); else seq_groups.push_back(group);
      }
    }
  }

  // XIC-based protein group organization by sequence and charge.
  void xic_pgroup::_build_seq_groups(fasta_data &fasta, vector<seq_group> &output_groups, vector<xic *> &xics, int min_cond) {
    vector<seq_group> groups;
    for(size_t i = 0; i < xics.size(); i++) {
      ms2peak *p = xics[i]->max_ms2();       // Get maximum scoring MS/MS peak.

      // Get all maximum scoring IDs and modifications.
      vector<fragment *> ids; p->max_ids(ids);
      vector<mods_t> mods; p->max_mods(mods);

      // NOTE: We can use ids[0] and mods[0] because clear_ambig_maxscore() got rid of 
      // any ambiguities (cases where there are multiple sequences).
      string seq = fasta.get_seqN(ids[0], mods[0]);
      int charge = xics[i]->charge;

      // Add XIC to sequence group.  NOTE: We need to differentiate by
      // chrage here because the quantification values different
      // charge states are not comparable. In contrast, isotope ratios
      // are comparable even though charge differs for the ratios.
      bool found = false;
      for(size_t s = 0; s < groups.size(); s++) {
	if(groups[s].charge == charge && groups[s] == seq) { groups[s].addxic(xics[i]); found = true;  break; }
      }
      // Create a new group, if none is found.
      if(found == false) {
	seq_group group;
	group.seq = seq; group.charge = charge; // Index by sequence and charge
	group.addxic(xics[i]);
	groups.push_back(group);
      }
    }

    // Compute maximum MS/MS score of sequence group using all XICs in that group.
    for(size_t i = 0; i < groups.size(); i++) {
      vector<xic *> &xics = groups[i].xics;
      double max_score = 0;
      for(size_t x = 0; x < xics.size(); x++) max_score = max(max_score, xics[x]->max_score());
      groups[i].max_ms2score = max_score;
    }

    // Sort sequence groups on size and then median MS/MS score.
    for(size_t i = 0; i < groups.size(); i++) {
      if(groups[i].numConds() >= min_cond) output_groups.push_back(groups[i]);
    }
    sort(output_groups.begin(), output_groups.end());
  }

  // Starting from a seed XIC x group XICs across iruns using best
  // bidirectional nearest neighbor queries.
  void analysis::form_align_group(xic *x, align_group &group) {
    vector<xic*> &xics = group.xics;
    if(x->grouped) return;
    // Start with a seed XIC.
    x->grouped = true; // Mark the seed as grouped.
    xics.push_back(x);
    for(int irun = 0; irun < (int)iruns.size(); irun++) {
      if(irun == x->irunID) continue;  // For every other irun perform nearest neighbor queries.
      xic *p = iruns[irun]->xics_all->NearestX(x->retentionTime, x->mz, 
					       conf.group_dtime, x->mz * conf.xic_width * 1e-6);
      if(p != NULL && p->grouped == false) {
	// Perform recipricol query back to current seed irun.
	xic *recip = iruns[x->irunID]->xics_all->NearestX(p->retentionTime, p->mz, 
							  conf.group_dtime, p->mz * conf.xic_width * 1e-6);
	// If best-bidirectional hit, add to group.
	if(recip == x) { p->grouped = true; xics.push_back(p);}
      }
    }
    group.position(); // Compute centroid of group.
  }

  // Apply thresholds on the number of conditions/replicates/iruns
  // an alignment group must span.
  bool analysis::is_valid_align_group(align_group &group) {
    if(group.xics.size() == 0) return false; // No XICs.
    int good_conds = 0; // Number of conditions with enough replicate sets.
    for(int condIndex = 0; condIndex < numConds(); condIndex++ ) {
      int good_reps = 0; // Number of replicate sets with enough iruns.
      for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) {
	// Count the number of iruns in a replicate set.
	if(group.countIruns(condIndex, repIndex) >= conf.min_iruns) good_reps++;
      }
      if(good_reps >= conf.min_reps) good_conds++;
    }
    return good_conds >= conf.min_conds;
  }

  struct xic_s {
    xic *x; double s;
    xic_s() { x = NULL; s = 0; }
    xic_s(xic *x_, double s_) { x = x_; s = s_; }
    bool operator < (xic_s b) const { return s < b.s; }
  };

  // This algorithm differs from our PNAS paper. It uses a "seed" XIC approach.
  void analysis::align_mode_quant() {
    // Assign an XIC a irun number and set as ungrouped. 
    for(int s = 0; s < (int)iruns.size(); s++) {
      vector<lcms::xic> &xics = iruns[s]->xics;
      for(size_t i = 0; i < xics.size(); i++) { xics[i].irunID = s; xics[i].grouped = false; }
    }

    // Create two sets of alignment grouping seeds.
    vector<xic_s> id_seeds, mi_seeds;
    for(int s = 0; s < (int)iruns.size(); s++) {
      vector<lcms::xic> &xics = iruns[s]->xics;
      for(size_t i = 0; i < xics.size(); i++) {
	// 1. Seeds wih IDs.
	if(xics[i].num_ids() > 0) id_seeds.push_back(xic_s(&xics[i], xics[i].max_score()));
	// 2. Seeds that are monoisotopic
	else if(xics[i].type == xicMONOISO)  mi_seeds.push_back(xic_s(&xics[i], xics[i].length()));
      }
    }
    sort(id_seeds.begin(), id_seeds.end()); // Sort seed XICs from low MS/MS score to high MS/MS score

    // Seed grouping from highest MS/MS score to lowest (note FDR cutoff must be applied).
    for(size_t i = id_seeds.size(); i > 0; i--) {
      // Start with a seed XIC.
      align_group group; form_align_group(id_seeds[i-1].x, group);
      if(is_valid_align_group(group)) align_groups.push_back(group);
    }

    // Sort monoisotopic seeds on their length.
    sort(mi_seeds.begin(), mi_seeds.end());
    // Split them by the median length value.
    vector<xic_s> mi_seeds1, mi_seeds2;
    for(size_t i = 0; i < mi_seeds.size() / 2; i++) mi_seeds1.push_back(mi_seeds[i]);
    for(size_t i = mi_seeds.size() / 2; i < mi_seeds.size(); i++) mi_seeds2.push_back(mi_seeds[i]);
    
    // Advance through seeds closest to median length value to those that are furthest..
    size_t i1 = mi_seeds1.size(), i2 = 0;
    while(true) {
      if(i1 > 0) {  // This set is decreasing in length.
	align_group group1; form_align_group(mi_seeds1[i1-1].x, group1);
	if(is_valid_align_group(group1)) align_groups.push_back(group1);
	i1--;
      }
      if(i2 < mi_seeds2.size()) { // This set is increasing in length.
	align_group group2; form_align_group(mi_seeds2[i2].x, group2);
	if(is_valid_align_group(group2)) align_groups.push_back(group2);
	i2++;
      }
      if(i1 == 0 && i2 == mi_seeds2.size()) break; // Stop when both lists are finished.
    }
    root = new _2dtree<align_group>(align_groups);     // Build a 2d-tree on aligned groups.

    // Perform protein grouping on aligned XIC groups that have IDs.
    vector<align_group *> align_groups_p;  
    for(size_t i = 0; i < align_groups.size(); i++) {
      if(align_groups[i].num_ids() > 0)	align_groups_p.push_back(&align_groups[i]);
    }
    protein_grouping(align_groups_p);
  }

  struct quant_increasing { bool operator() (xic *a, xic *b) const {  return a->quant < b->quant;  } };

  void analysis::label_free_normalize() { 
    if(conf.skip_normalize_label_free) return; // Don't normalize if requested

    // Only get log2 quant values from aligned groups that have MS/MS IDs.
    vector< vector<double> > log2quant(num_iruns());

    if(conf.quantmode == quantALIGN) {
      // For each aligned group with an MS/MS ID collect log2
      // quantification values separated by column number in
      // label-free output.
      for(size_t a = 0; a < align_groups.size(); a++) {
	if(align_groups[a].num_ids() == 0) continue;
	vector<xic*> xics; align_groups[a].getxics(xics);
	for(size_t i = 0; i < xics.size(); i++) {
	  int idx = toColumn(xics[i]->condIndex, xics[i]->repIndex, xics[i]->irunIndex);
	  log2quant[idx].push_back(util::Log2(xics[i]->quant));
	}
      }
      
      // Determine which columns are still empty.
      vector<bool> empty_log2quant(num_iruns(), false);
      for(size_t idx = 0; idx < log2quant.size(); idx++) 
	empty_log2quant[idx] = log2quant[idx].size() == 0;

      // Now use align groups that don't have IDs for normalization.
      for(size_t a = 0; a < align_groups.size(); a++) {
	vector<xic*> xics; align_groups[a].getxics(xics);
	for(size_t i = 0; i < xics.size(); i++) {
	  int idx = toColumn(xics[i]->condIndex, xics[i]->repIndex, xics[i]->irunIndex);
	  if(empty_log2quant[idx]) 
	    log2quant[idx].push_back(util::Log2(xics[i]->quant));
	}
      }
    }
    else if(conf.quantmode == quantXICBASED) {
      // For XIC-based quantification, use XICs with MS/MS IDs for
      // normalization.
      for(int condIndex = 0; condIndex < numConds(); condIndex++ ) {
	for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) {
	  for(int irunIndex = 0; irunIndex < numIruns(condIndex, repIndex); irunIndex++) {
	    lcms::irun *irun = getIrun(condIndex, repIndex, irunIndex);
	    int idx = toColumn(condIndex, repIndex, irunIndex);
	    vector<xic> &xics = irun->xics;
	    for(size_t i = 0; i < xics.size(); i++) 
	      if(xics[i].num_ids() > 0) log2quant[idx].push_back(util::Log2(xics[i].quant));
	  }
	}
      }      
    }

    // Compute the medians of the collected log2 quant values.
    vector<double> medians;
    for(size_t idx = 0; idx < log2quant.size(); idx++) {
      if(log2quant[idx].size() > 0) medians.push_back(util::median(log2quant[idx]));
    }
    if(medians.size() == 0) return; // No medians, no normalization.

    double median_of_medians = util::median(medians);
    cout << "median_of_medians=" << median_of_medians << endl;
    // Perform normalization. 
    for(int condIndex = 0; condIndex < numConds(); condIndex++ ) {
      for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) {
	for(int irunIndex = 0; irunIndex < numIruns(condIndex, repIndex); irunIndex++) {
	  lcms::irun *irun = getIrun(condIndex, repIndex, irunIndex);
	  int idx = toColumn(condIndex, repIndex, irunIndex);
	  // Center irun around its median value then rescale to median of medians.
	  if(log2quant[idx].size() > 0) {
	    double irun_median = util::median(log2quant[idx]);
	    irun->adjust_log2_quant(median_of_medians - irun_median);
	  }
	}
      }
    }
  }

  void analysis::process() {
    // Get data set ranges.
    vector<double> logImins, logImaxs;
    for(size_t i = 0; i < iruns.size(); i++) {
      if(iruns[i]->max_mz > max_mz) max_mz = iruns[i]->max_mz; if(iruns[i]->min_mz < min_mz) min_mz = iruns[i]->min_mz;
      if(iruns[i]->max_time > max_time) max_time = iruns[i]->max_time;
      if(iruns[i]->min_time < min_time) min_time = iruns[i]->min_time;
      if(iruns[i]->max_mz2 > max_mz2) max_mz2 = iruns[i]->max_mz2;
      if(iruns[i]->min_mz2 < min_mz2) min_mz2 = iruns[i]->min_mz2;
      logImins.push_back(iruns[i]->logImin_disp);
      logImaxs.push_back(iruns[i]->logImax_disp);
    }
    // Compute median minimum intensity display values for entire data set.
    if(logImins.size() > 0 && logImaxs.size() > 0) {
      conf.logImin_disp = util::median(logImins); conf.logImax_disp = util::median(logImaxs);
    }
    if(conf.run_algorithms == false) return;
    if(processed_flag) return; else processed_flag = true;
    
    // Process all conditions if they haven't been processed already.
    for(size_t i = 0; i < conditions.size(); i++) conditions[i]->process(DB); 

    // TODO: Get rid of this shit. 
    if(pepxml.size() > 0) xml::load_pepxml(pepxml, iruns, fasta);

    if(DB[dbLIGHT] != NULL) {
      vector<ms2::ms2peak *> ms2fdr;
      // Go through each instrument run. 
      for(size_t s = 0; s < iruns.size(); s++) {
	vector<lcms::xic> &xics = iruns[s]->xics;
	// Collect pointers to all MS/MS spectra in XICs
	for(size_t i = 0; i < xics.size(); i++) if(xics[i].num_ids() > 0) xics[i].get_ms2(ms2fdr);
      }
      // Apply FDR cutoff to all MS/MS peaks.
      fdr_split_charge(ms2fdr); 
      clear_ambig_maxscore(ms2fdr); // Clear ambigious MS/MS IDs.
    }

    if(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE || conf.quantmode == quantAcetylSTOICH) {
      // Apply protein grouping + FDR cutoffs if there is an MS/MS DB loaded.
      vector<lcms::isotope_group *> groups; get_isotope_groups(groups);
      // Assign IDs to isotope groups (pairs, triples, quads?).
      for(long id = 0; id < (long)groups.size(); id++) groups[id]->id = id;
      // Group MS/MS IDs by protein or protein group.
      protein_grouping(groups);
    }
    else if(conf.quantmode == quantALIGN || conf.quantmode == quantXICBASED) {
      // Label free quantification modes.
      if(conf.quantmode == quantALIGN) align_mode_quant();  // Group monoiso-XICs after alignment
      else xic_based_quant(); // Use IDs to cross reference.
      label_free_normalize();
    }
  }

  // ms2.scanNum ms2.score ms2.charge ms2.protein ms2.seq
  void analysis::save_ms2peak_header(ofstream &csv) {
    csv << "ms2.scanNum\t"; 
    csv << "ms2.rt\t"; csv << "ms2.mz\t"; 
    csv << "ms2.score\t";   csv << "ms2.qvalue\t";
    csv << "ms2.charge\t";  
    csv << "ms2.protein\t"; csv << "ms2.protein.nfrags\t"; 
    csv << "ms2.seq\t";
    csv << "ms2.missedcleaves\t";
    csv << "ms2.ppm.error\t"; 
    csv << "mz.theory\t";
    csv << "ms2.Nmods\t";
    if(conf.quantmode == quantISOPAIR && conf.isotope_15N) csv << "ms2.nitrogens\t";
  }

  // ms2.scanNum ms2.score ms2.charge ms2.protein ms2.seq
  void analysis::save_ms2peak(ofstream &csv, ms2peak *peak) {
    if(peak == NULL || DB[dbLIGHT] == NULL) {
      // 15N needs one more tab.
      if(conf.quantmode == quantISOPAIR && conf.isotope_15N) csv << "\t\t\t\t\t\t\t\t\t\t";
      else csv << "\t\t\t\t\t\t\t\t\t"; // output nothing
      return; // return
    }
    csv << peak->scanNum << '\t'; // ms2.scanNum 
    csv << peak->retentionTime << '\t'; // ms2.rt 
    csv << fixed << setprecision(8) << peak->mz << '\t'; // ms2.mz
    csv << peak->max_score() << '\t';  // ms2.score
    csv << fixed << setprecision(9) << peak->min_qvalue() << '\t';  // ms2.qvalue
    csv << peak->charge << '\t'; // ms2.charge

    vector<ms2::fragment *> ids;  peak->max_ids(ids);
    vector<ms2::mods_t> mods;     peak->max_mods(mods);

    vector<int> acc(ids.size());
    for(size_t i = 0; i < ids.size(); i++) acc[i] = ids[i]->protein_idx;
    save_protein_group(csv, acc, true);

    ms2::fragmentDB *db = DB[dbLIGHT];
    if(peak->iso == ms2::isoHEAVY) db = DB[dbHEAVY]; 
    else if(peak->iso == ms2::isoMEDIUM) db = DB[dbMEDIUM];
    if(db == NULL) { csv << "\t\t\t"; return; }

    // NOTE: Each maximum scoring ID has the same sequence!!!!
    string seq = get_seqS(ids[0], mods[0]);
    csv << seq << "\t"; // ms2.seq
    csv << ids[0]->misses << "\t"; // ms2.missedcleaves

    double ppm_error = 0;
    double seq_mass = db->precursor_mass(ids[0], mods[0]);
    // NOTE: Use peak->mass() since we actually measure the protonated mass.
    ppm_error = mass::massError2ppm(peak->neutral_mass() - seq_mass, peak->mass());
    csv << ppm_error << "\t"; // ms2.ppm.error

    double mz_theory = (seq_mass + double(peak->charge) * mass::proton) / double(peak->charge);
    csv << fixed << setprecision(8) << mz_theory << "\t"; // mz.theory

    csv << mods[0].count() << "\t"; // ms2.Nmods
    if(conf.quantmode == quantISOPAIR && conf.isotope_15N) {
      csv << db->nitrogen_cnt(ids[0], mods[0]) << "\t";
    }
  }

  // Output column names for a label-free row of quantification data.
  void analysis::saveLabelFreeHeader(ofstream &csv) {
    // Outputs columns with condition, replicate, and irun information.
    for(int condIndex = 0; condIndex < numConds(); condIndex++ ) {
      // Get the condition and replicate information for column names.
      string condition = getCond(condIndex)->description;
      // Iterate through replicate sets.
      for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) {
	string replicate = getRep(condIndex, repIndex)->description;
	// Iterate through iruns.
	for(int irunIndex = 0; irunIndex < numIruns(condIndex, repIndex); irunIndex++) {
	  string irun = getIrun(condIndex, repIndex, irunIndex)->description;
	  csv << "cond." << condition << ".rep." << replicate << ".irun." << irun << '\t';
	}
      }
    }
    // Add columns for condition medians.
    for(int condIndex = 0; condIndex < numConds(); condIndex++) {
      string condition = getCond(condIndex)->description;
      csv << "cond." << condition << ".median";
      // \t or newline if last column
      if(condIndex == numConds() - 1) csv << endl; else csv << '\t';
    }
  }

  // Save a row of label-free quantification data given a set of XICs.
  void analysis::saveLabelFreeRow(ofstream &csv, vector<xic *> &xics) {
    // A row contains value for each irun and a median for each condition.
    vector<double> row(num_iruns() , -1); // if row value = -1, it implies NA
    // NOTE: We can use -1 here for NA because all log2quant values are > 0.

    for(size_t i = 0; i < xics.size(); i++) {
      int column = toColumn(xics[i]->condIndex, xics[i]->repIndex, xics[i]->irunIndex);
      row.at(column) = xics[i]->quant;
    }

    // Actually output the populated row. -1 imples NA.
    for(size_t j = 0; j < row.size(); j++) {
      if(row[j] < 0) csv << "NA"; else csv << fixed << setprecision(4) << row[j]; 
      csv << '\t';
    }

    // Segregate XICs by condition number.
    vector< vector<double> > by_cond(numConds());
    for(size_t j = 0; j < xics.size(); j++) 
      by_cond.at(xics[j]->condIndex).push_back(util::Log2(xics[j]->quant));

    // Output condition medians as the last set of columns.
    for(int condIndex = 0; condIndex < numConds(); condIndex++) {
      if(by_cond[condIndex].size() == 0) csv << "NA";
      else csv << fixed << setprecision(4) << pow(2, util::median(by_cond[condIndex]));
      if(condIndex == numConds() - 1) csv << endl;
      else csv << '\t';
    }
  }

  // Generate a file with mass errors for label-free data.
  void analysis::saveMassErrorsMS1(string filename) {
    if(DB[dbLIGHT] == NULL) return;

    ofstream csv(filename.c_str());
    csv << "ms2.mz\t"; // m/z of MS/MS peak
    csv << "ms2.da.error\t"; // error in daltons
    csv << "ms2.ppm.error"; // error in PPMs
    csv << endl;

    for(size_t s = 0; s < iruns.size(); s++) {
      vector<xic> &xics = iruns[s]->xics;
      for(size_t x = 0; x < xics.size(); x++) {
	// Irun through MS/MS peaks.
	vector<ms2peak *> peaks; xics[x].get_ms2(peaks);
	for(size_t i = 0; i < peaks.size(); i++) {
	  ms2peak *peak = peaks[i];
	  if(peak->num_ids() == 0) continue;


	  // For those that have MS/MS ids, compute the mass error.
	  vector<ms2::fragment *> ids;  peak->max_ids(ids); 
	  vector<ms2::mods_t> mods; peak->max_mods(mods);
	  
	  ms2::fragmentDB *db = DB[dbLIGHT]; 
	  if(peak->iso == ms2::isoHEAVY) db = DB[dbHEAVY]; 
	  else if(peak->iso == ms2::isoMEDIUM) db = DB[dbMEDIUM];

	  // Get the precursor mass error.
	  double seq_mass = db->precursor_mass(ids[0], mods[0]);
	  double da_error = peak->neutral_mass() - seq_mass;
	  
	  // NOTE: Use peak->mass() since we actually measure the protonated mass!!!!
	  double ppm_error = mass::massError2ppm(da_error, peak->mass());
	  csv << peak->mz << '\t'; csv << da_error << '\t'; csv << ppm_error; csv << endl; 
	}
      }
    }

  }

  // Generate a file with mass errors for label-free data.
  void analysis::saveMassErrorsMS2(string filename) {
    if(DB[dbLIGHT] == NULL) return;
    // Only one column with MS/MS errors.
    ofstream csv(filename.c_str()); 
    csv << "ms2.da.error\t"; // error in daltons
    csv << "ms2.ppm.error" << endl;

    // Iterate through each irun.
    for(size_t s = 0; s < iruns.size(); s++) {
      vector<xic> &xics = iruns[s]->xics;
      // Iterate through each XIC in the irun.
      for(size_t x = 0; x < xics.size(); x++) {
	vector<ms2peak *> peaks; xics[x].get_ms2(peaks);
	// Irun through MS/MS peaks in each XIC.
	for(size_t i = 0; i < peaks.size(); i++) {
	  // Compute PPM errors for charge +2 spectra.
	  ms2peak *peak = peaks[i];
	  vector<double> ppm_errors, da_errors;
	  ms2::fragmentDB *db = NULL; if(peak->iso == ms2::isoHEAVY) db = DB[dbHEAVY]; else db = DB[dbLIGHT];
	  db->ms2_mass_error(peak, da_errors, ppm_errors, iruns[s]->min_mz2, iruns[s]->max_mz2);
	  // Print errors in a single columns (in PPMs). 
	  for(size_t k = 0; k < ppm_errors.size(); k++)   csv << da_errors[k] << "\t" << ppm_errors[k] << endl;
	}
      }
    }
  }


  void analysis::saveAlignGroupHeader(ofstream &csv) {
    csv << "group.mz\t"; 
    csv << "group.retentionTime\t"; 
    csv << "group.size\t";
    csv << "xic.c\t"; csv << "xic.r\t"; csv << "xic.s\t";
    save_ms2peak_header(csv);
  }

  // Saves data in XICs that hvae been grouped across iruns.
  void analysis::saveAlignGroup(ofstream &csv, align_group *group) {
    csv << fixed << setprecision(6) << group->mz << '\t'; // group.mz
    csv << group->retentionTime << '\t'; // group.retentionTime
    csv << group->countxics() << '\t'; // group.size

    // Get XIC corresponding to best MS/MS ID.
    lcms::xic *max_xic = group->max_xic();
      
    // Output information that allows cross reference of MS/MS spectrum.
    csv << getCond(max_xic->condIndex)->description << '\t'; // xic.cond
    csv << getRep(max_xic->condIndex, max_xic->repIndex)->description << '\t'; // xic.rep
    csv << getIrun(max_xic->condIndex, max_xic->repIndex, max_xic->irunIndex)->description << '\t'; // xic.irun

    // Get best MS/MS peak with best ID.
    ms2::ms2peak *peak = max_xic->max_ms2();
    save_ms2peak(csv, peak);   // ms2.scanNum ms2.score ms2.charge ms2.protein ms2.seq

    // Get XICs and populate the row of the table with quantative data.  
    vector<xic *> xics; group->getxics(xics);
    saveLabelFreeRow(csv, xics);
  }

  // Creates a CSV file with ID and quantification data for alignment based quantification.
  void analysis::saveAlignTable(string filename) {
    if(iruns.size() == 0 || conf.quantmode != quantALIGN || DB[dbLIGHT] == NULL) return;

    ofstream csv(filename.c_str());
    // Output the header of the CSV file and count the number of iruns.
    csv << "group.id\t";
    csv << "protein.group\t";
    saveAlignGroupHeader(csv);  saveLabelFreeHeader(csv);
    
    // Done with writing column header information.
    for(size_t g = 0; g < align_pgroups.size(); g++) {
      csv << g + 1 << '\t'; // group.id
      save_protein_group(csv, align_pgroups[g].acc);  // protein.group

      // NOTE: align_groups is sorted based on size and then MS/MS
      // score.  0th entry is the best aligned XIC group!!!
      saveAlignGroup(csv, align_pgroups[g].align_groups[0]);
    }
  }

  void analysis::saveAlignTablePeptides(string filename) {
    if(iruns.size() == 0) return;
    if(conf.quantmode != quantALIGN) return;
    if(DB[dbLIGHT] == NULL) return;

    ofstream csv(filename.c_str());
    // Output the header of the CSV file and count the number of iruns.
    csv << "group.id\t";
    csv << "protein.group\t";
    saveAlignGroupHeader(csv);
    saveLabelFreeHeader(csv);
    
    // Done with writing column header information.
    for(size_t i = 0; i < align_pgroups.size(); i++) {
      for(size_t g = 0; g < align_pgroups[i].align_groups.size(); g++) {
	csv << i + 1 << '\t'; // group.id
	save_protein_group(csv, align_pgroups[i].acc);  // protein.group
	saveAlignGroup(csv, align_pgroups[i].align_groups[g]);
      }
      for(size_t g = 0; g < align_pgroups[i].shared_align_groups.size(); g++) {
	csv << i + 1 << '\t'; // group.id
	save_protein_group(csv, align_pgroups[i].acc);  // protein.group
	saveAlignGroup(csv, align_pgroups[i].shared_align_groups[g]);
      }
    }
  }

  void analysis::saveXIC_allids(string filename) {
    if(DB[dbLIGHT] == NULL) return;
    ofstream csv(filename.c_str());

    csv << "condition" << '\t';
    csv << "replicate" << '\t';
    csv << "irun" << '\t';
    save_ms2peak_header(csv);
    csv << "log2.quant" << endl;

    // Outputs columns with condition, replicate, and irun information.
    for(int condIndex = 0; condIndex < numConds(); condIndex++ ) {
      // Get the condition and replicate information for column names.
      string condition = getCond(condIndex)->description;
      // Iterate through replicate sets.
      for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) {
	string replicate = getRep(condIndex, repIndex)->description;
	// Iterate through iruns.
	for(int irunIndex = 0; irunIndex < numIruns(condIndex, repIndex); irunIndex++) {
	  irun *cur_irun = getIrun(condIndex, repIndex, irunIndex);
	  // Get XICs from current run.
	  vector<xic> &xics = cur_irun->xics;
	  for(size_t x = 0; x < xics.size(); x++) {
	    // Scan through XICs skip NO ID cases.
	    if(xics[x].num_ids() == 0) continue;
	    // Get highest scoring MS/MS ID and output.
	    ms2peak *peak = xics[x].max_ms2();
	    csv << condition << '\t';
	    csv << replicate << '\t';
	    csv << cur_irun->description << '\t';
	    save_ms2peak(csv, peak);   // ms2.scanNum ms2.score ms2.charge ms2.protein ms2.seq
	    csv << fixed << setprecision(8) << util::Log2(xics[x].quant) << endl;
	  }
	}
      }
    }

  }


  // Nonlinear alignment based internal correlations.
  void analysis::saveAlignTableInternal(string filename) {
    // Correlate conditions to one another. 
    ofstream csv(filename.c_str());

    // Output header for computing internal correlations.
    csv << "group.id\t"; // protein group ID
    csv << "description\t"; // name of internal correlation
    csv << "group1.log2.ratio\t"; // group one from that condition
    csv << "group1.N\t"; // size of group 1
    csv << "group2.log2.ratio\t"; // group two from that condition
    csv << "group2.N"; // size of group 2
    csv << endl;

    // Confirm tryptic fragment groups show the same ratio! 
    // Go throgh each of O(N^2) conditions.
    for(int condIndex1 = 0; condIndex1 < numConds(); condIndex1++) {
      string &condition1 = getCond(condIndex1)->description;
      for(int condIndex2 = condIndex1 + 1; condIndex2 < numConds(); condIndex2++) {
	string &condition2 = getCond(condIndex2)->description;

	// Iterate through each alignment protein group.
	for(size_t g = 0; g < align_pgroups.size(); g++) {
	  vector<align_group *> &frag = align_pgroups[g].align_groups;
	  vector<align_group *> &frag_shared = align_pgroups[g].shared_align_groups;

	  int frag_cnt = 0;
	  vector<double> group1, group2;

	  // Compute ratio for each tryptic fragment in the protein group.
	  for(size_t s = 0; s < frag.size(); s++) {
	    
	    vector<xic *> xics; frag[s]->getxics(xics);	    
	    vector<double> c1, c2;
	    // Collect all log2 quantification values from condition 1 and condition 2.
	    for(size_t i = 0; i < xics.size(); i++) {
	      if(xics[i]->condIndex == condIndex1) { c1.push_back(util::Log2(xics[i]->quant)); }
	      if(xics[i]->condIndex == condIndex2) { c2.push_back(util::Log2(xics[i]->quant)); }
	    }
	    // Store ratio between the conditions.
	    if(c1.size() > 0 && c2.size() > 0) {
	      if(frag_cnt % 2 == 0) group1.push_back(util::median(c1) - util::median(c2));
	      else group2.push_back(util::median(c1) - util::median(c2));
	      frag_cnt++;
	    }
	  }
	  // TODO: Get rid of this repeated code!!!

	  // Compute ratio for *shared* tryptic fragment.
	  for(size_t s = 0; s < frag_shared.size(); s++) {
	    vector<xic *> xics; frag_shared[s]->getxics(xics);
	    vector<double> c1, c2;
	    // Collect all log2 quantification values from condition 1 and condition 2.
	    for(size_t i = 0; i < xics.size(); i++) {
	      if(xics[i]->condIndex == condIndex1) { c1.push_back(util::Log2(xics[i]->quant)); }
	      if(xics[i]->condIndex == condIndex2) { c2.push_back(util::Log2(xics[i]->quant)); }
	    }
	    if(c1.size() > 0 && c2.size() > 0) {
	      if(frag_cnt % 2 == 0) group1.push_back(util::median(c1) - util::median(c2));
	      else group2.push_back(util::median(c1) - util::median(c2));
	      frag_cnt++;
	    }
	  }
	  if(group1.size() > 0 && group2.size() > 0) {
	    csv << g + 1 << "\t"; // group.id 
	    csv << condition1 << "_" << condition2 << "\t"; // description
	    csv << setprecision(8) << fixed << util::median(group1) << '\t';
	    csv << group1.size() << '\t';
	    csv << setprecision(8) << fixed << util::median(group2) << '\t';
	    csv << group2.size() << endl;
	  }
	}
      }
    }
  }


  void analysis::saveSeqGroupHeader(ofstream &csv) {
    csv << "seq\t" ; 
    csv << "charge\t"; 
    csv << "max.ms2score\t"; 
    csv << "group.size\t";
    csv << "file.name\t";
    save_ms2peak_header(csv);
    saveLabelFreeHeader(csv);
  }

  void analysis::saveSeqGroup(ofstream &csv, seq_group &group) {
    csv << group.seq << '\t'; // ms2.seq
    csv << group.charge << '\t'; // charge
    csv << group.max_ms2score << '\t'; // max.ms2score
    csv << group.numxics() << '\t'; // group.size
    // If there are multiple XICs across several fractions the one
    // with the highest search score is used for XIC based
    // quantification.
    xic *x = group.max_ms2_xic();
    if(x == NULL) csv << '\t';
    else csv << getIrun(x->condIndex, x->repIndex, x->irunIndex)->filename << '\t'; // file.name
    if(x == NULL) save_ms2peak(csv, NULL);  else save_ms2peak(csv, x->max_ms2());
    saveLabelFreeRow(csv, group.xics);
  }

  // Creates a table of XIC-based data.
  void analysis::saveXICBasedTable(string filename) {
    cout << "saveXICBasedTable():" << filename << endl;
    if(iruns.size() == 0 || conf.quantmode != quantXICBASED || DB[dbLIGHT] == NULL) return;
    ofstream csv(filename.c_str());
    csv << "group.id\t";
    csv << "protein.group\t";
    csv << "protein.group.masses\t"; 
    saveSeqGroupHeader(csv);

    for(size_t g = 0; g < xic_pgroups.size(); g++) {
      csv << g + 1 << '\t'; // group.id
      save_protein_group(csv, xic_pgroups[g].acc); // protein.group

      for(size_t i = 0; i < xic_pgroups[g].acc.size(); i++) {
	double mass = DB[dbLIGHT]->get_mass(xic_pgroups[g].acc[i]);
	csv << mass / 1e3;
	if(i != xic_pgroups[g].acc.size() - 1) csv << " ||||| ";
      }
      csv << '\t';
      saveSeqGroup(csv, xic_pgroups[g].seq_groups[0]); 
    }
  }


  void analysis::saveXICBasedTablePeptides(string filename) {
    cout << "saveXICBasedTablePeptides():" << filename << endl;
    if(iruns.size() == 0 || conf.quantmode != quantXICBASED || DB[dbLIGHT] == NULL) return;
    cout << "saveXICBasedTablePeptides():RUNNING" << endl;
    ofstream csv(filename.c_str());
    csv << "group.id\t";
    csv << "protein.group\t";
    saveSeqGroupHeader(csv);

    for(size_t g = 0; g < xic_pgroups.size(); g++) {
      xic_pgroup &pgroup = xic_pgroups[g];
      // Save sequence groups.
      for(size_t s = 0; s < pgroup.seq_groups.size(); s++) {
	csv << g + 1 << '\t'; // group.id
	save_protein_group(csv, pgroup.acc); // protein.group
	saveSeqGroup(csv, pgroup.seq_groups[s]);
      }
      // Save shared sequence groups.
      for(size_t s = 0; s < pgroup.seq_groups_shared.size(); s++) {
	csv << g + 1 << '\t'; // group.id
	save_protein_group(csv, pgroup.acc); // protein.group
	saveSeqGroup(csv, pgroup.seq_groups_shared[s]);
      }
    }
  }

  void analysis::saveXICBasedTableInternal(string filename) {
    ofstream csv(filename.c_str());
    // Output header for computing internal correlations.
    csv << "group.id\t"; // protein group ID
    csv << "description\t"; // name of internal correlation
    csv << "group1.log2.ratio\t"; // group one from that condition
    csv << "group1.N\t"; // size of group 1
    csv << "group2.log2.ratio\t"; // group two from that condition
    csv << "group2.N"; // size of group 2
    csv << endl;

    // Confirm tryptic fragment groups show the same ratio! 
    for(int condIndex1 = 0; condIndex1 < numConds(); condIndex1++) {
      string &condition1 = getCond(condIndex1)->description;
      for(int condIndex2 = condIndex1 + 1; condIndex2 < numConds(); condIndex2++) {
	string &condition2 = getCond(condIndex2)->description;

	// Iterate through each alignment protein group.
	for(size_t g = 0; g < xic_pgroups.size(); g++) {
	  vector<seq_group> &frag = xic_pgroups[g].seq_groups;
	  vector<seq_group> &frag_shared = xic_pgroups[g].seq_groups_shared;

	  int frag_cnt = 0;
	  vector<double> group1, group2;

	  // Compute ratio for each tryptic fragment in the protein group.
	  for(size_t s = 0; s < frag.size(); s++) {
	    vector<xic *> xics; frag[s].getxics(xics);
	    vector<double> c1, c2;
	    // Collect all log2 quantification values from condition 1 and condition 2.
	    for(size_t i = 0; i < xics.size(); i++) {
	      if(xics[i]->condIndex == condIndex1) { c1.push_back(util::Log2(xics[i]->quant)); }
	      if(xics[i]->condIndex == condIndex2) { c2.push_back(util::Log2(xics[i]->quant)); }
	    }
	    // Store ratio between the conditions.
	    if(c1.size() > 0 && c2.size() > 0) {
	      if(frag_cnt % 2 == 0) group1.push_back(util::median(c1) - util::median(c2));
	      else group2.push_back(util::median(c1) - util::median(c2));
	      frag_cnt++;
	    }
	  }

	  // Compute ratio for *shared* tryptic fragment.
	  for(size_t s = 0; s < frag_shared.size(); s++) {
	    vector<xic *> xics; frag_shared[s].getxics(xics);
	    vector<double> c1, c2;
	    // Collect all log2 quantification values from condition 1 and condition 2.
	    for(size_t i = 0; i < xics.size(); i++) {
	      if(xics[i]->condIndex == condIndex1) { c1.push_back(util::Log2(xics[i]->quant)); }
	      if(xics[i]->condIndex == condIndex2) { c2.push_back(util::Log2(xics[i]->quant)); }
	    }
	    if(c1.size() > 0 && c2.size() > 0) {
	      if(frag_cnt % 2 == 0) group1.push_back(util::median(c1) - util::median(c2));
	      else group2.push_back(util::median(c1) - util::median(c2));
	      frag_cnt++;
	    }
	  }
	  if(group1.size() > 0 && group2.size() > 0) {
	    csv << g + 1 << "\t"; // group.id 
	    csv << condition1 << "_" << condition2 << "\t"; // description
	    csv << setprecision(8) << fixed << util::median(group1) << '\t';
	    csv << group1.size() << '\t';
	    csv << setprecision(8) << fixed << util::median(group2) << '\t';
	    csv << group2.size() << endl;
	  }
	}
      }
    }
    
  }

  // Saves monoisotopic aligned XICs that have no IDs.
  void analysis::saveLabelFreeMS1(string filename) {
    if(iruns.size() == 0 || conf.quantmode != quantALIGN) return;

    ofstream csv(filename.c_str());
    csv << "group.mz\t";   csv << "group.retentionTime\t";  
    csv << "group.size\t"; csv << "group.charge\t";
    csv << "group.N.ms2\t";
    saveLabelFreeHeader(csv);

    for(size_t i = 0; i < align_groups.size(); i++) {
      if(align_groups[i].num_ids() > 0) continue; 
      csv << fixed << setprecision(6) << align_groups[i].mz << '\t'; // group.mz
      csv << align_groups[i].retentionTime << '\t'; // group.retentionTime
      csv << align_groups[i].countxics() << '\t';   // group.size
      csv << align_groups[i].charge() << '\t';      // group.charge
      csv << align_groups[i].num_ms2() << '\t';     // group.N.ms2
      saveLabelFreeRow(csv, align_groups[i].xics);
    }
  }

  // Use the accompoany R script to generate plots showing internal
  // group correlations on a per condition basis and a per replicate
  // basis. Correlates non-overlapping tryptic peptides. 
  void analysis::save_isotope_groups_internal(string filename) {
    if(DB[dbLIGHT] == NULL) return; // no database, nothing to output

    ofstream csv(filename.c_str());
    // Output header for computing internal correlations.
    csv << "group.id\t";                    // protein group ID
    csv << "description\t";                 // name of internal correlation
    csv << "group1.log2.ratio\tgroup1.N\t"; // group one from that condition
    csv << "group2.log2.ratio\tgroup2.N";   // group two from that condition
    csv << endl;

    // Geneates replicate set internal correlations.  Generates
    // per-replicate internal correlations. Iterate through each
    // condition and ratio type (for triple labels).
    for(int rtypeidx = 0; rtypeidx < 3; rtypeidx++) {
      ratio_type rtype = (ratio_type)rtypeidx;
      string rtypestr;
      switch(rtype) { case isoHL: rtypestr = "HL"; break; case isoML: rtypestr = "ML"; break; case isoHM: rtypestr = "HM"; break; }
      for(int condIndex = 0; condIndex < numConds(); condIndex++) {
	string condition = getCond(condIndex)->description;
	// Iterate through each replicate in that condition.
	for(int repIndex = 0; repIndex <= numReps(condIndex); repIndex++) {
	  string replicate;
	  // Note that if repIndex == numReps(condIndex) replicate will
	  // be the null string. This indicates we collect everything,
	  // regardless of the replicate index number.
	  if(repIndex < numReps(condIndex)) replicate = getRep(condIndex, repIndex)->description;
	
	  // Go through each protein group.
	  for(size_t g = 0; g < isop_pgroups.size(); g++) {
	    isop_pgroup &data = *isop_pgroups[g];
	    vector<double> group1, group2, g1_unique, g2_unique, g1_shared, g2_shared;
	    int g1g2idx = 0;

	    // Go through each sequence group and collect ratios for that sequence.
	    for(size_t s = 0; s < data.seq_groups.size(); s++) {
	      // Skip sequences that have modifications.
	      if(data.seq_groups[s].Nmods > 0) continue; 

	      // Handle index number and replicate index number.
	      vector<double> ratios; 
	      if(repIndex < numReps(condIndex))	data.seq_groups[s].get_ratios(rtype, ratios, condIndex, repIndex);
	      else data.seq_groups[s].get_ratios(rtype, ratios, condIndex); // Handle case were all reps are combined.

	      // Update group1 and group 2 ratios.
	      if(ratios.size() > 0) {
		if(g1g2idx % 2 == 0) { group1.push_back(util::median(ratios)); g1_unique.push_back(group1.back()); }
		else { group2.push_back(util::median(ratios)); g2_unique.push_back(group2.back()); }
		g1g2idx++;
	      }
	    }
	    // Repeat for shared peptides (probably introduces some noise).
	    for(size_t s = 0; s < data.seq_groups_shared.size(); s++) {
	      // Skip sequences with modifications.	      
	      if(data.seq_groups_shared[s].Nmods > 0) continue; 

	      vector<double> ratios; 
	      if(repIndex < numReps(condIndex))	data.seq_groups_shared[s].get_ratios(rtype, ratios, condIndex, repIndex);
	      else data.seq_groups_shared[s].get_ratios(rtype, ratios, condIndex);
	    
	      if(ratios.size() > 0) {
		if(g1g2idx % 2 == 0) { group1.push_back(util::median(ratios)); g1_shared.push_back(group1.back()); }
		else { group2.push_back(util::median(ratios)); g2_shared.push_back(group2.back()); }
		g1g2idx++;
	      }
	    }

	    string idstr = rtypestr;
	    if(replicate == "") idstr += "_" + condition; else idstr += "_" + condition + "_" + replicate;
	    // Only produce output if enough measurements. 
	    if(group1.size() > 0 && group2.size() > 0) {
	      csv << g + 1 << '\t'; // protein group.id for cross referencing
	      csv << idstr << '\t';
	      csv << util::median(group1) << '\t' << group1.size() << '\t';
	      csv << util::median(group2) << '\t' << group2.size() << endl;
	    }
	    if(g1_unique.size() > 0 && g2_unique.size() > 0) {
	      csv << g + 1 << '\t'; // protein group.id for cross referencing
	      csv << idstr << "_unique" << '\t';
	      csv << util::median(g1_unique) << '\t' << g1_unique.size() << '\t';
	      csv << util::median(g2_unique) << '\t' << g2_unique.size() << endl;
	    }
	    if(g1_shared.size() > 0 && g2_shared.size() > 0) {
	      csv << g + 1 << '\t'; // protein group.id for cross referencing
	      csv << idstr << "_shared" << '\t';
	      csv << util::median(g1_shared) << '\t' << g1_shared.size() << '\t';
	      csv << util::median(g2_shared) << '\t' << g2_shared.size() << endl;
	    }
	  }
	}
      }
    }
  }

  int analysis::isotope_table_header(ofstream &csv) {
    int Ncol = 0;
    // Outputs columns with condition, replicate, and irun information.
    for(int condIndex = 0; condIndex < numConds(); condIndex++ ) {
      // Get the condition and replicate information for column names.
      string condition = getCond(condIndex)->description;
      // Iterate through replicate sets.
      for(int repIndex = 0; repIndex < numReps(condIndex); repIndex++) {
	string replicate = getRep(condIndex, repIndex)->description;
	csv << "cond." << condition << ".rep." << replicate << '\t';
	Ncol++; // Keep track of number of replicate set columns.
      }
    }
    // Add columns for condition medians.
    for(int condIndex = 0; condIndex < numConds(); condIndex++) {
      string condition = getCond(condIndex)->description;
      csv << "cond." << condition << ".median";
      // \t or newline if last column
      if(condIndex == numConds() - 1) csv << endl; else csv << '\t';
    }
    return Ncol;
  }

  // Save isotope ratios by condition.
  void analysis::save_isotope_ratios(ofstream &csv, int Ncol, vector<isotope_group *> &groups) {
    vector< vector<double> > row(Ncol); 
    vector< vector<double> > by_cond(numConds());

    // Build row information.
    for(size_t i = 0; i < groups.size(); i++) {
      // By replicate set.
      if(groups[i]->by_search == false) { // Ignore cases obtained by search.
	row.at(toColumn(groups[i]->condIndex, groups[i]->repIndex)).push_back(groups[i]->log2ratioHL);
	// Also grouped by condition.
	by_cond.at(groups[i]->condIndex).push_back(groups[i]->log2ratioHL);
      }
    }

    // Output row for replicate sets.
    for(size_t i = 0; i < row.size(); i++) {
      if(row[i].size() > 0)  csv << fixed << setprecision(8) << util::median(row[i]) << '\t';
      else csv << "NA\t";
    }

    // Output condition medians as the last set of columns.
    for(int condIndex = 0; condIndex < numConds(); condIndex++) {
      if(by_cond[condIndex].size() == 0) csv << "NA";
      else csv << fixed << setprecision(8) << util::median(by_cond[condIndex]);
      if(condIndex == numConds() - 1) csv << endl;
      else csv << '\t';
    }
  }

  struct groupdata_t{
    int acc; string meta; int nfrags;
    bool operator < (const groupdata_t &B) const { return meta < B.meta; }
  };
  void analysis::save_protein_group(ofstream &csv, vector<int> &acc, bool extra_data) {
    // Get group meta information.
    vector<groupdata_t> groupdata;
    for(size_t i = 0; i < acc.size(); i++) {
      groupdata_t data;
      data.acc = acc[i];
      data.meta = DB[dbLIGHT]->get_meta(acc[i]);
      data.nfrags = DB[dbLIGHT]->get_nfrags(acc[i]);
      groupdata.push_back(data);
    }
    // Sort for uniform output.
    sort(groupdata.begin(), groupdata.end());
    // protein.group
    for(size_t i = 0; i < groupdata.size(); i++) {
      if(i != groupdata.size() - 1) csv << groupdata[i].meta << " ||||| "; 
      else csv << groupdata[i].meta;
    }
    csv << '\t'; 
    if(extra_data) {
      for(size_t i = 0; i < groupdata.size(); i++) {
	if(i != groupdata.size() - 1) csv << groupdata[i].nfrags << " ||||| "; 
	else csv << groupdata[i].nfrags;
      }
      csv << '\t';
    }
  }

  // Save sequence group information for peptide ratio/peptide ratio
  // correlation.
  void analysis::save_internal_pep_seq_group(ofstream &csv, int groupid, 
					     int condIndex, string &condition, 
					     vector<isop_seqgroup> &seq_groups) {
    // For each sequence group, get all groups from the current
    // condition.  All of these groups have the same sequence.
    for(size_t s = 0; s < seq_groups.size(); s++) {	
      vector<isotope_group *> &groups = seq_groups[s].groups;
      for(int rtypeidx = 0; rtypeidx < 3; rtypeidx++) {
	ratio_type rtype = (ratio_type)rtypeidx;
	string rtstr;
	switch(rtype) { case isoHL: rtstr = "HL"; break; case isoML: rtstr = "ML"; break; case isoHM: rtstr = "HM"; break;}
	vector<double> group1, group2;
	for(size_t p = 0; p < groups.size(); p++) {
	  isotope_group *group = groups[p];
	  // Skip groups not in current condition.
	  if(group->condIndex != condIndex) continue;
	  double log2ratio = 0; bool has_ratio = false;
	  if(rtype == isoHL && group->heavy && group->light )  { log2ratio = group->log2ratioHL; has_ratio = true; }
	  if(rtype == isoHM && group->heavy && group->medium ) { log2ratio = group->log2ratioHM; has_ratio = true; }
	  if(rtype == isoML && group->medium && group->light ) { log2ratio = group->log2ratioML; has_ratio = true; }
	  // Group1 = charge 1 ratios, group2 = charge >= 3 ratios.
	  if(has_ratio) {
	    if(group->charge() == 2) group1.push_back(log2ratio);
	    else if(group->charge() >= 3)  group2.push_back(log2ratio);
	  }
	}
	if(group1.size() > 0 && group2.size() > 0) {
	  // Output correlation for sequence..
	  csv << groupid << "\t" << rtstr << "_" << condition << "\t"; 
	  csv << seq_groups[s].seq << "\t";
	  csv << fixed << setprecision(8) << util::median(group1) << "\t"  << group1.size() << "\t";
	  csv << fixed << setprecision(8) << util::median(group2) << "\t"  << group2.size() << "\t";
	  csv << seq_groups[s].Nmods;
	  csv << endl;
	}
      }
    }
  }

  // Output internal correlations on a per-peptide basis.
  void analysis::save_isotope_groups_internal_pep(string filename) {
    if(DB[dbLIGHT] == NULL) return;
    ofstream csv(filename.c_str());
    csv << "group.id\tdescription\tseq\t"; 
    csv << "group1.log2.ratio\t"; // group one from that condition
    csv << "group1.N\t"; // size of group 1
    csv << "group2.log2.ratio\t"; // group two from that condition
    csv << "group2.N\t"; // size of group 2
    csv << "Nmods";
    csv << endl;

    // Handle each condition separately.
    for(int condIndex = 0; condIndex < numConds(); condIndex++) {
      string condition = getCond(condIndex)->description;
      condition += "_pep"; // indicate it is a peptide correlation
      // Go through each protein group.
      for(size_t g = 0; g < isop_pgroups.size(); g++) {
	// Get sequence groups for that protein group.
	save_internal_pep_seq_group(csv, g+1, condIndex, condition, isop_pgroups[g]->seq_groups);
	save_internal_pep_seq_group(csv, g+1, condIndex, condition, isop_pgroups[g]->seq_groups_shared);
      }
    }    
  }

  // Saves an isotope ratio table on the per-peptide level.
  void analysis::save_isotope_groups_table_pep(string filename) {
    if(DB[dbLIGHT] == NULL) return;
    ofstream csv(filename.c_str());

    csv << "group.id\t"; // identifier/index for protein group
    csv << "protein.group\t"; // protein group with each name sorted and  ||||| as a separator
    csv << "seq\t"; 
    csv << "seq.shared\t";

    // Save header with replicate and condition information.
    int Ncol = isotope_table_header(csv);

    // Go through each protein group.
    for(size_t g = 0; g < isop_pgroups.size(); g++) {
      // Get sequence group for that protein group.
      vector<isop_seqgroup> &seq_groups = isop_pgroups[g]->seq_groups;
      vector<isop_seqgroup> &seq_groups_shared = isop_pgroups[g]->seq_groups_shared;

      // For each sequence group, save group.
      for(size_t s = 0; s < seq_groups.size(); s++) {
	csv << g + 1 << '\t'; // group.id
	save_protein_group(csv, isop_pgroups[g]->acc); // protein.group
	csv << seq_groups[s].seq << "\t"; // seq
	csv << "FALSE\t"; // seq.shared
	save_isotope_ratios(csv, Ncol, seq_groups[s].groups);
      }

      for(size_t s = 0; s < seq_groups_shared.size(); s++) {
	csv << g + 1 << '\t'; // group.id
	save_protein_group(csv, isop_pgroups[g]->acc); // protein.group
	csv << seq_groups_shared[s].seq << "\t"; // seq
	csv << "TRUE\t"; // seq.shared
	save_isotope_ratios(csv, Ncol, seq_groups_shared[s].groups);
      }

    }
  }

  // Save isotope data in tabular format grouped by replicate set and by condition.
  void analysis::save_isotope_groups_table(string filename) {
    if(DB[dbLIGHT] == NULL) return;
    ofstream csv(filename.c_str());

    csv << "group.id\t"; // identifier/index for protein group
    csv << "protein.group\t"; // protein group with each name sorted and  ||||| as a separator

    int Ncol = isotope_table_header(csv);

    for(size_t g = 0; g < isop_pgroups.size(); g++) {
      csv << g + 1 << '\t'; // group.id
      save_protein_group(csv, isop_pgroups[g]->acc); // protein.group
      save_isotope_ratios(csv, Ncol, isop_pgroups[g]->groups);      
    }
  }

  // Saves information for correlation ratios across replicates runs.
  void analysis::save_isotope_groups_corr(string filename) {
    if(DB[dbLIGHT] == NULL) return;
    ofstream csv(filename.c_str());

    csv << "group.id\t"; // id/index for protein group.
    csv << "corr.id\t"; csv << "rep1.id\t"; csv << "rep2.id\t";
    csv << "log2.ratio1\t"; csv << "ratio1.N\t";
    csv << "log2.ratio2\t"; csv << "ratio2.N" << endl;

    for(int condIndex = 0; condIndex < numConds(); condIndex++) {
      if( numReps(condIndex) <= 1 ) continue; // Need multiple replicate sets to do correlation.
      string condition = getCond(condIndex)->description;

      // Go through O(N^2) replicates. 
      for(int repIndex1 = 0; repIndex1 < numReps(condIndex); repIndex1++) {
	string replicate1 = getRep(condIndex, repIndex1)->description;

	for(int repIndex2 = repIndex1 + 1; repIndex2 < numReps(condIndex); repIndex2++) {
	  string replicate2 = getRep(condIndex, repIndex2)->description;

	  // Iterate through protein groups.
	  for(size_t g = 0; g < isop_pgroups.size(); g++) {
	    vector<isotope_group *> &groups = isop_pgroups[g]->groups;
	    
	    vector<double> rep1, rep2;
	    for(size_t i = 0; i < groups.size(); i++) {
	      if(groups[i]->condIndex == condIndex) {
		if(groups[i]->repIndex == repIndex1 && !groups[i]->by_search) 
		  rep1.push_back(groups[i]->log2ratioHL);
		else if(groups[i]->repIndex == repIndex2 && !groups[i]->by_search) 
		  rep2.push_back(groups[i]->log2ratioHL);
	      }
	    }
	    // We can correlate these two replciates.
	    if(rep1.size() > 0 && rep2.size() > 0) {
	      double log2ratio1 = util::median(rep1);
	      double log2ratio2 = util::median(rep2);
	      csv << g + 1 << '\t'; // group.id
	      
	      // Save an ID for the condition and that replicate set group.
	      csv << condition + "_" + replicate1 + "_" + replicate2 << '\t'; // corr.id
	      csv << replicate1 << '\t'; // rep1.id
	      csv << replicate2 << '\t'; // rep2.id
	      
	      csv << fixed << setprecision(8) << log2ratio1 << '\t' << rep1.size() << '\t';
	      csv << fixed << setprecision(8) << log2ratio2 << '\t' << rep2.size() << endl; 
	    }
	    // Next protein group.
	  }
	  // Next replicate index #2.
	}
	// Next Replicate index #1.
      }
      // Next condition index.
    }
  }


  void analysis::save_iso_header(ofstream &csv, char eol_chr) {
    csv << "ratio.id\t"; // identifier/index for ratio
    csv << "group.c\t"; // condition name of the ratio origin
    csv << "group.r\t"; // replicate from which ratio originated
    csv << "group.i\t"; // instrument run the ratio originated
    if(conf.quantmode == quantISOPAIR) {
      csv << "L.mz\t"; // m/z of light XIC
      csv << "L.rt\t"; // retention time of most intense peak in light XIC
      csv << "H.mz\t"; // ditto heavy
      csv << "H.rt\t"; // ditto heavy
      csv << "log2.HL.ratio\t"; // actual ratio measured usually heavy over light 
      csv << "log2.xicH\t"; // log2(area) of heavy XIC
      csv << "log2.xicL\t"; // log2(area) of light XIC  (swapped if labels are reversed)
      if(conf.isotope_15N) csv << "group.nitrogens\t";
    }
    else if(conf.quantmode == quantISOTRIPLE) {
      csv << "log2.HM.ratio\t";
      csv << "log2.HL.ratio\t";
      csv << "log2.ML.ratio\t";
      csv << "log2.xicH\t";
      csv << "H.mz\t"; // ditto heavy
      csv << "H.rt\t"; // ditto heavy
      csv << "log2.xicM\t";
      csv << "M.mz\t"; // ditto medium
      csv << "M.rt\t"; // ditto medium
      csv << "log2.xicL\t";
      csv << "L.mz\t"; // m/z of light XIC
      csv << "L.rt\t"; // retention time of most intense peak in light XIC
    }
    save_ms2peak_header(csv); // MS/MS information
    csv << "ms2.origin"; // designates highest scoring MS/MS ID originated from heavy XIC
    csv << eol_chr;
  }

  void analysis::save_iso(ofstream &csv, isotope_group *group, char eol_chr) {
    csv << group->id + 1 << '\t'; // ratio.id
    csv << getCond(group->condIndex)->description << '\t'; // group.cond
    csv << getRep(group->condIndex, group->repIndex)->description << '\t'; // group.rep
    csv << getIrun(group->condIndex, group->repIndex, group->irunIndex)->description << '\t'; // group.irun

    if(conf.quantmode == quantISOPAIR) {
      if(group->by_search) {
	if(group->light != NULL) {
	  csv << group->light->mz << '\t'; // L.mz
	  csv << group->light->retentionTime << '\t'; // L.rt
	  csv << "NA" << '\t'; // H.mz
	  csv << "NA" << '\t'; // H.rt
	  csv << "LIGHT" << '\t'; // log2.HL.ratio
	  csv << "NA" << '\t'; // log2.xicH
	  csv << fixed << setprecision(8) << group->log2xicL << '\t'; // log2.xicL
	}
	if(group->heavy != NULL) {
	  csv << "NA" << '\t'; // L.mz
	  csv << "NA" << '\t'; // L.rt
	  csv << group->heavy->mz << '\t'; // H.mz
	  csv << group->heavy->retentionTime  << '\t'; // H.rt
	  csv << "HEAVY" << '\t'; // log2.HL.ratio
	  csv << fixed << setprecision(8) << group->log2xicH << '\t'; // log2.xicH
	  csv << "NA" << '\t'; // log2.xicL
	}
	if(conf.isotope_15N) csv << group->nitrogens << '\t';
      }
      else {
	csv << group->light->mz << '\t'; // L.mz
	csv << group->light->retentionTime << '\t'; // L.rt
	csv << group->heavy->mz << '\t'; // H.mz
	csv << group->heavy->retentionTime  << '\t'; // H.rt
	csv << fixed << setprecision(8) << group->log2ratioHL << '\t'; // log2.HL.ratio
	csv << fixed << setprecision(8) << group->log2xicH << '\t'; // log2.xicH
	csv << fixed << setprecision(8) << group->log2xicL << '\t'; // log2.xicL
	if(conf.isotope_15N) csv << group->nitrogens << '\t';
      }
    }
    else if(conf.quantmode == quantISOTRIPLE) {
      if(group->medium != NULL && group->heavy != NULL) 
	csv << fixed << setprecision(8) << group->log2ratioHM << '\t';
      else csv << "NA\t";
      if(group->heavy != NULL && group->light != NULL) 
	csv << fixed << setprecision(8) << group->log2ratioHL << '\t';
      else csv << "NA\t";
      if(group->medium != NULL && group->light != NULL)
	csv << fixed << setprecision(8) << group->log2ratioML << '\t';
      else csv << "NA\t";
      if(group->heavy != NULL) {
	csv << fixed << setprecision(8) << group->log2xicH  << '\t'; 
	csv << fixed << setprecision(8) << group->heavy->mz  << '\t'; 
	csv << fixed << setprecision(8) << group->heavy->retentionTime  << '\t'; 
      }
      else csv << "NA\tNA\tNA\t";
      if(group->medium != NULL) {
	csv << fixed << setprecision(8) << group->log2xicM  << '\t'; 
	csv << fixed << setprecision(8) << group->medium->mz  << '\t'; 
	csv << fixed << setprecision(8) << group->medium->retentionTime  << '\t'; 
      }
      else csv << "NA\tNA\tNA\t";
      if(group->light != NULL) {
	csv << fixed << setprecision(8) << group->log2xicL  << '\t'; 
	csv << fixed << setprecision(8) << group->light->mz  << '\t'; 
	csv << fixed << setprecision(8) << group->light->retentionTime  << '\t'; 
      }
      else csv << "NA\tNA\tNA\t";
    }
	
    // Get best MS/MS peak with best ID.
    ms2::ms2peak *peak = group->max_ms2(); save_ms2peak(csv, peak);   // ms2.scanNum ms2.score ms2.charge ms2.protein ms2.seq
	
    // ms2.origin
    if(peak->iso == ms2::isoLIGHT) csv << "LIGHT"; 
    else if(peak->iso == ms2::isoMEDIUM) csv << "MEDIUM";
    else if(peak->iso == ms2::isoHEAVY) csv << "HEAVY";
    else csv << "UNKNOWN";

    csv << eol_chr;
  }

  // Creates a CSV file with ID and quantification data. 
  void analysis::save_isotope_groups(string filename) {
    // || DB[dbLIGHT] == NULL || DB[dbHEAVY] == NULL
    if(!(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE)) return;

    ofstream csv(filename.c_str());
    // Output the header of the CSV file and count the number of iruns.
    csv << "group.id\t"; // identifier/index for protein group
    csv << "protein.group\t";
    csv << "protein.group.nfrags\t";
    csv << "group.N\t"; // number of measurements in protein group
    csv << "group.log2.HL.ratio\t"; // median of all of the ratios in protein group
    save_iso_header(csv);
    csv << "ms2.shared"; // shared peptide that maps to multiple protein groups, but was assigned to this one
    csv << endl;

    for(size_t groupId = 0; groupId < isop_pgroups.size(); groupId++) { 
      isop_pgroup *data = isop_pgroups[groupId]; // Output data for each protein group.

      // Get median log2 ratio for the protein group.
      vector<double> ratios;
      for(size_t p = 0; p < data->groups.size(); p++) {
	// NOT RIGHT WHEN BY SEARCH.
	if(data->groups[p]->by_search == false)	ratios.push_back(data->groups[p]->log2ratioHL);
      }
      double group_median = 0;
      if(ratios.size() > 0) group_median = util::median(ratios);

      // Now output all groups in protein group.
      for(size_t p = 0; p < data->groups.size(); p++) {
	isotope_group *group = data->groups[p];
	csv << groupId + 1 << '\t'; // group.id
	save_protein_group(csv, data->acc, true); // protein.group protein.group.nfrags
	csv << data->groups.size() << '\t'; // group.N
	if(ratios.size() > 0) csv << fixed << setprecision(8) << group_median << '\t'; // group.log2.HL.ratio
	else csv << fixed << setprecision(8) << "NA" << '\t'; // group.log2.HL.ratio
	save_iso(csv, group);
	// ms2.shared
	ms2::ms2peak *peak = group->max_ms2();
	vector<int> uniq_idxs; peak->max_idxs(uniq_idxs);
	if(uniq_idxs.size() > data->acc.size()) csv << "TRUE" << endl;  // NOTE: Last column so no tab!!!!!
	else csv << "FALSE" << endl;
      }
    }
  }

  void analysis::save_isotope_data(string filename) {
    if(!(conf.quantmode == quantISOPAIR || conf.quantmode == quantISOTRIPLE)) return;
    ofstream csv(filename.c_str());
    save_iso_header(csv, '\n');

    vector<lcms::isotope_group *> groups; 
    get_isotope_groups(groups);

    for(long idx = 0; idx < (long)groups.size(); idx++) {
      if(groups[idx]->num_ids() > 0) save_iso(csv, groups[idx], '\n');
    }
  }

  void analysis::save_stoich_data(string filename) {
    if(!(conf.quantmode == quantAcetylSTOICH)) return;
    ofstream csv(filename.c_str());

    csv << "group.c\t"; // condition name of the ratio origin
    csv << "group.r\t"; // replicate from which ratio originated
    csv << "group.i\t"; // instrument run the ratio originated
    save_ms2peak_header(csv); // MS/MS information
    csv << "xics.mz\t"; // m/z of xics
    csv << "xics.rt\t"; // retention time of XICs
    csv << "xics.log2.quant\t"; // log2(area) of sXIC
    csv << "grouping.phase" << endl;

    vector<lcms::isotope_group *> groups;  get_isotope_groups(groups);

    for(long idx = 0; idx < (long)groups.size(); idx++) {
      if(groups[idx]->num_ids() == 0) continue;

      csv << getCond(groups[idx]->condIndex)->description << '\t'; // group.cond
      csv << getRep(groups[idx]->condIndex, groups[idx]->repIndex)->description << '\t'; // group.rep
      csv << getIrun(groups[idx]->condIndex, groups[idx]->repIndex, groups[idx]->irunIndex)->description << '\t'; // group.irun
      
      ms2peak *peak = groups[idx]->max_ms2();
      save_ms2peak(csv, peak);

      vector<xic*> &xics = groups[idx]->xics;
      for(size_t x = 0; x < xics.size(); x++) { 
	if(xics[x] != NULL) csv << xics[x]->mz;
	else csv << "NA";
	if(x != xics.size() - 1) csv <<  " ||||| ";  else csv << "\t";
      }

      for(size_t x = 0; x < xics.size(); x++) { 
	if(xics[x] != NULL) csv << xics[x]->retentionTime;
	else csv << "NA";
	if(x != xics.size() - 1) csv <<  " ||||| ";  else csv << "\t";
      }

      for(size_t x = 0; x < xics.size(); x++) { 
	if(xics[x] != NULL) csv << fixed << setprecision(8) << util::Log2(xics[x]->quant);
	else csv << "NA";
	if(x != xics.size() - 1) csv <<  " ||||| ";  else csv << '\t';
      }
      csv << groups[idx]->grouping_phase << endl;
    }

  }
  

}

