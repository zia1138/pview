#ifndef __ms2_hpp__
#define __ms2_hpp__

#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <limits>

#include "msconfig.hpp"
#include "util.hpp"

// In this header file:
// 1. MS/MS ID algorithm
// 2. Protein grouping algorithm.
// 3. FDR cutoff selection algorithm.
// 4. Simple fasta file loader (no data validataion). 
namespace ms2 {
  using namespace std;

  // TODO: Add HCD activation. How the fuck does HCD work?
  enum activationMethod_t { CID_activation = 0, ETD_activation = 1, HCD_activation = 2}; 
  struct peak2 { float mz2, intensity2; };    // An MS/MS peak

  // Holds modified positions + modification at that position.
  struct modpos_t { 
    modpos_t() { pos = -1; prec_mass = frag_mass = 0; site_score = 0; mod_idx = fragmod_idx = -1; }
    modpos_t(int p, double mp, double mf) { 
      pos = p; prec_mass = mp; frag_mass = mf; site_score = 0; 
      mod_idx = fragmod_idx = -1;
    }
    int pos;  // Position of modified amino acid.
    double prec_mass, frag_mass; // Precursor and fragment mass of modified amino acid.
    double site_score; // a score for this site and corresponding q-value
    int mod_idx, fragmod_idx;  // modification ID
  };

  // Holds a set of modifications for a peptide.
  struct mods_t { 
    vector<modpos_t> sites; // Position + modification pairs.
    int count() { return sites.size(); }
    // Returns fragment mass modified amino acid at position p.
    double get_frag_mass(int p) {
      for(size_t k = 0; k < sites.size(); k++) 
	if(sites[k].pos == p) return sites[k].frag_mass; 
      return 0;
    }
    // Returns site scoreat position p.
    double site_score(int p) {
      for(size_t k = 0; k < sites.size(); k++) if(sites[k].pos == p) return sites[k].site_score; 
      return 0;
    }
    // Get precursor mass of modification at position p.
    double get_prec_mass(int p) {
      for(size_t k = 0; k < sites.size(); k++) if(sites[k].pos == p) return sites[k].prec_mass; 
      return 0;
    }
    // Returns true if position p is modified.
    bool is_modified(int p) {
      for(size_t k = 0; k < sites.size(); k++) if(sites[k].pos == p) return true;
      return false;
    }
    int get_mod_idx(int p) {
      for(size_t k = 0; k < sites.size(); k++) if(sites[k].pos == p) return sites[k].mod_idx; 
      return -1;
    }
    int get_fragmod_idx(int p) {
      for(size_t k = 0; k < sites.size(); k++) if(sites[k].pos == p) return sites[k].fragmod_idx; 
      return -1;
    }
  };


  // Peptide fragment with protein index and pointer to sequence data.
  struct fragment {
    int protein_idx; // protein description index in fasta_data::protein_meta.
    int seq_idx;     // index to sequence data in fasta_data::seq_data usually the same as protein_idx, 
                     // but different for pepxml external input
    int pos, len, misses;    // position, length, and number of missed cleavages in seq_data[seq_idx]
    double fragmass; // Each fragment is sorted on neutral mass.
    int nitrogens; // Number of nitrogens in this fragment.
    bool is_decoy;   // Flag indicates fragment originates from a the decoy protein sequence.

    fragment() { fragmass = -1;  pos = len = misses = 0; seq_idx = protein_idx = -1; is_decoy = false; }
    fragment(double m) { fragmass = m; }
    fragment(int i, int p, int l, double m, int miss, bool isd) {
      protein_idx = i; pos = p; len = l; fragmass = m; misses = miss; is_decoy = isd;
    }
    bool operator < (const fragment &b) const { return fragmass < b.fragmass; }
  };

  // Functor for sorting pointers to frag_data.
  struct cmp_frag_mass {
    bool operator () (const fragment *A, const fragment *B) const { return A->fragmass < B->fragmass; } 
  };

  // MS/MS ID and score assigned to an MS/MS peak.
  struct ms2id_t {
    ms2id_t() { score = 0; qvalue = 1; }
    fragment frag; // note that this is a copy of a fragment, not a pointer, allows data to be cached
    mods_t mods;    // any modified positions, 3 values position, precursor mass, fragment mass
    double score, qvalue; // assigned score and q-value of assigned score 
  };

  // Isotope status of an MS/MS spectrum.
  enum isotope_type { isoUNKNOWN, isoLIGHT, isoHEAVY, isoMEDIUM };

  // MS/MS peak that is in an LC-MS/MS instrument run.
  struct ms2peak { 
    ms2peak() { charge = 0; scanNum = -1; iso = isoUNKNOWN; activation = CID_activation; } 
    double retentionTime, mz; // precursor retention time, m/z, and intensity
    int charge;  // instrument determined charge for MS/MS id algorithm
    vector<peak2> ms2; // centroided MS/MS spectrum
    isotope_type iso; // for isotope pairing or tripling, isotope type
    long scanNum; // scan number in original mzXML file
    activationMethod_t activation; // CID or ETD (default is CID)
    vector<ms2id_t> results;       // MS/MS id information from DB search after FDR correction.
    vector<ms2id_t> target, decoy; // Top scoring search results to target and decoy databases.

    // Accessors for _2dtree<>. 
    double x() const { return retentionTime; } double y() const { return mz; }

    // Returns maximum m/z value in spectrum.
    float max_mz2() {
      float max_mz2 = -(numeric_limits<float>::max());
      for(size_t i = 0; i < ms2.size(); i++) if(ms2[i].mz2 > max_mz2) max_mz2 = ms2[i].mz2;
      return max_mz2;
    }
    // Returns minimum m/z value in spectrum.
    float min_mz2() {
      float min_mz2 = numeric_limits<float>::max();
      for(size_t i = 0; i < ms2.size(); i++) if(ms2[i].mz2 < min_mz2) min_mz2 = ms2[i].mz2;
      return min_mz2;
    }
    
    // Return number of MS/MS fragments assigned.
    int num_ids() { return max(results.size(), target.size() + decoy.size()); }

    // Applies a score cutoff on MS/MS DB search results.  
    void score_cutoff(double cutoff) {
      for(size_t i = 0; i < target.size(); i++) {
	if(target[i].score >= cutoff - 1e-7) results.push_back(target[i]);
      }      
      clear_target(); clear_decoy();
    }
    // Returns maximum score. Returns 0 if no maximum score.  Assumes scores are > 0.
    double max_score() { double s = 0; for(size_t i = 0; i < results.size(); i++) s = max(s,results[i].score); return s; }
    double max_target_score() { double s = 0; for(size_t i = 0; i < target.size(); i++) s = max(s, target[i].score); return s;  }
    double max_decoy_score() { double s = 0; for(size_t i = 0; i < decoy.size(); i++) s = max(s, decoy[i].score); return s;  }
    double max_td_score() { return max(max_target_score(), max_decoy_score());  }

    // Return minimum q-value.
    double min_qvalue() { 
      double s = 0, q = 1; 
      for(size_t i = 0; i < results.size(); i++) {
	if(results[i].score > s ) { s = results[i].score; q = results[i].qvalue; }
      }
      return q; 
    }

    // Return fragments that have the best score.
    void max_ids(vector<fragment *> &ids) {
      double score = max_score();
      for(size_t i = 0; i < results.size(); i++) if(fabs(score - results[i].score) < 1e-7) ids.push_back(&results[i].frag);
    }

    // Return their corresponding modifications.
    void max_mods(vector<mods_t> &mods) {
      double score = max_score();
      for(size_t i = 0; i < results.size(); i++) {
	if(fabs(score - results[i].score) < 1e-7) mods.push_back(results[i].mods);
      }
    }

    // For the highest scoring MS/MS ID returns the set of protein
    // indexes to which that ID maps. Used for forming protein groups.
    void max_idxs(vector<int> &pidxs) {
      vector<fragment *> ids; max_ids(ids); // Get maximum scoring IDs.
      for(size_t i = 0; i < ids.size(); i++) {
	int protein_idx = ids[i]->protein_idx; // Get protein index.
	bool found = false;
	// Has this protein index been added?
	for(size_t j = 0; j < pidxs.size(); j++) if(protein_idx == pidxs[j]) found = true;
	// Keep only unique protein indexes.
	if(found == false) pidxs.push_back(protein_idx);
      }
    }
    // Return the number of protein indexes in the highest scoring
    // MS/MS ID.  Also used for forming protein groups.
    int max_idx_cnt() { vector<int> idxs; max_idxs(idxs); return idxs.size(); }

    // Get mass of the fragment, including any protons (postive mode only). 
    double mass() { return double(charge) * mz; }
    // Returns neutral mass. POSTIVE MODE NEUTRAL MASS ONLY!!!!       
    double neutral_mass() { return double(charge) * (mz  - mass::proton);  }
    // Clears MS/MS db search results after FDR correction.
    void clear_ms2id() { results.clear(); vector<ms2id_t> (results).swap(results);   }
    // Free up memory used by MS/MS peaks.
    void clear_ms2peaks() {  ms2.clear(); vector<peak2> (ms2).swap(ms2); }
    // Clears MS/MS target/decoy search results.
    void clear_decoy() { decoy.clear(); vector<ms2id_t> (decoy).swap(decoy); }
    void clear_target() { target.clear(); vector<ms2id_t> (target).swap(target); }
  }; 

  // Loads fasta file. Prepends '#' N-term symbol and appends '*' as C-term symbol.
  bool load_fasta(string descr, string fasta, vector<string> &meta, vector<string> &seq);

  // Fasta sequence data with concatenated reversed decoy
  // database. All created fragment databases share this sequence
  // information.
  struct fasta_data {
    msconfig &conf;
    vector<string> protein_meta; // meta information for protein 
    vector<bool> is_decoy_protein; // flag indicating decoy or not, one-to-one with protein_meta
    vector<string> seq_data; // amino acid sequence data, in one-to-one with seq_data

    fasta_data(msconfig &conf_) : conf(conf_) { }
    // Total number of protein sequence present in the database.
    size_t Nproteins() { return protein_meta.size(); }
    bool load_fasta(vector<string> &descr, vector<string> &files);
    void rev_seq(string &seq); 
    // Given a fragment, pulls out sequence data w/o modification annotations.
    void get_seq_nomod(fragment *f, string &seq) {
      // Use '[' only if not the actual protein N-term.
      bool cleavage_N_term = f->pos == 0; if(cleavage_N_term == false) seq.push_back('[');
      // Collect sequence data. 
      string &source_seq = seq_data[f->seq_idx];
      for(int i = f->pos; i < f->pos + f->len; i++) seq.push_back(source_seq[i]);
      // Use ']' only if not the actual protein C-term.
      bool cleavage_C_term = (int)source_seq.length() == f->pos + f->len;
      if(cleavage_C_term == false) seq.push_back(']');
    }
    // Returns a sequence string the precursor mass shift of a
    // modification annotated.  This call is critical for dealing with
    // ambigious localization, grouping peptides, etc.
    string get_seqN(fragment *f, mods_t &m) {
      // NOTE: No SNP is applied here.
      string seq; get_seq_nomod(f, seq);
      string seq2;
      for(int i = 0; i < (int)seq.length(); i++) {
	seq2 += seq[i];
	if(m.is_modified(i)) {
	  seq2 += "[";
	  int mod_idx = m.get_mod_idx(i);
	  if(mod_idx >= 0) seq2 += conf.varMods.at(mod_idx).abbrev + "]";
	  else seq2 += util::toString(m.get_prec_mass(i)) + "]";
	}
      }
      return seq2;
    }  
    // Returns sequence string with precursor mass, site score and
    // fragment modification abbreviation.
    string get_seqS(fragment *f, mods_t &m) {
      string seq; get_seq_nomod(f, seq);
      string seq2;

      for(int i = 0; i < (int)seq.length(); i++) {
	seq2 += seq[i];
	if(m.is_modified(i)) {
	  seq2 += "[" + util::toString(m.get_frag_mass(i)) + ", ";
	  seq2 += util::toString(m.site_score(i));
	  int mod_idx = m.get_mod_idx(i), fragmod_idx = m.get_fragmod_idx(i);
	  if(mod_idx >= 0 && fragmod_idx >= 0) { 
	    seq2 += ", ";
	    seq2 += conf.varMods.at(mod_idx).fragmods.at(fragmod_idx).abbreviation;
	  }
	  seq2 += "]";
	}
      }
      return seq2;
    }      
  };

  const int NumDBTypes = 3;
  enum dbtype_t { dbLIGHT = 0, dbHEAVY = 1, dbMEDIUM = 2 };

  // Generates tryptic fragment database allowing a specified number of missed cleavages. 
  struct fragmentDB {
    msconfig &conf; // Global configuration information.
    vector<double> MW; // Molecular weight table.
    vector<int> nitrogens; // Amino acid to nitrogen count table.
    vector<varmod_t> MODs; // active, variable modification table, 0th mod is always the no-mod mod
    fasta_data &fasta; // Shared FASTA sequence data or imported MS/MS database results.
    vector<double> seq_mass; // mass of each accession in fasta file

    vector<int> fragment_count; // number of fragments 
    // Once created, do not modify this memory by adding/removing elements from vector<>.
    size_t frag_cnt;
    vector<fragment> frag_data;
    vector<fragment*> frags;  // Fragment list (sorted by mass). 

    vector<isotope_t> active_iso; // active isotopes
    dbtype_t dbtype; // database type, heavy, medium, light.

    // Constructor just stores reference to global configuration information.
    fragmentDB(msconfig &conf_, fasta_data &fasta_, dbtype_t type = dbLIGHT) : conf(conf_), fasta(fasta_) { 
      dbtype = type;
      // Read through configuration information and load active
      // isotope labels.  These are used to restrict fragment query
      // results for database search.
      for(size_t i = 0; i < conf.activeIsotopes.size(); i++) {
	if(0 <= conf.activeIsotopes[i] && conf.activeIsotopes[i] < (int)conf.isotopes.size()) {
	  active_iso.push_back(conf.isotopes[conf.activeIsotopes[i]]);
	}
      }
    }

    void build(); // Build fragment database. Call this before running anything else!!!
    void init_MWs(); // Build molecular weight table.

    // Compute the precursor mass of amino acid sequence.
    double precursor_mass(fragment *f, mods_t &m);
    int nitrogen_cnt(fragment *f, mods_t &m); 

    // Generate a, b, or c-ion and x,y,z-ion ladders.
    void abc_ion(double mwdelta, int charge, string &seq, mods_t &m, 
		 vector<double> &ladder, double mz1, double mz2, int cutoff = -1, 
		 bool exclude_proline = false);
    void xyz_ion(double mwdetal, int charge, string &seq, mods_t &m, 
		 vector<double> &ladder, double mz1, double mz2, int cutoff = -1,
		 bool exclude_proline = false);

    // Construct a ladder for a given peak. Different ladders
    // constructed based on charge and neutral mass.
    void build_ladder(activationMethod_t method, vector<double> &ladder, string &seq, 
		      mods_t &mods, int charge, double m, 
		      float min_mz2, float max_mz2, int cutoff = -1);

    // Compute MS/MS mass errors using b-ions and y-ioins.
    void ms2_mass_error(ms2peak *p, vector<double> &da_errors, 
			vector<double> &ppm_errors, float min_mz2, float max_mz2);

    // Returns protein meta information associated with tryptic fragment.
    const string &get_meta(fragment *f) { return fasta.protein_meta.at(f->protein_idx); }
    const string &get_meta(int idx) { return fasta.protein_meta.at(idx);  }

    // Return mass information for a protein sequence.
    double get_mass(int idx) { if(seq_mass.size() > 0) return seq_mass.at(idx); else return 0; }

    int get_nfrags(int idx) { if(fragment_count.size() > 0) return fragment_count.at(idx); else return -1; }

    // Digest sequence with given protein index. 
    void digest_all(); 
    void digest(size_t pidx, bool count = false);

    // Reverse a sequence, swapping cut site according to enzyme.
    void rev_seq(string &seq);

    // Applies an enzyme to digest sequence.
    bool enzyme(char aa1_, char aa2_) {  
      char aa1 = toupper(aa1_), aa2 = toupper(aa2_);
      if(aa2 == '*') return false; // Don't cut at C-term
      switch(conf.enzyme_id) {
      case ez_Trypsin: return trypsin(aa1, aa2);
      case ez_LysC: return lysc(aa1);
      case ez_GluC: return gluc(aa1); 
      case ez_ArgC: return argc(aa1, aa2);
      case ez_TrypsinP: return trypsinP(aa1);
      case ez_LysN: return lysn(aa2);
      case ez_ChymoTrypsin: return chymotrypsin(aa1,aa2);
      case ez_AspN: return aspn(aa2);
      default: return false;
      }
    }
    // TODO: Add all enzymes present in Mascot MS/MS search.
    // Glu-C cuts at the carboxyl side of glutamate
    bool gluc(char aa1) {  return aa1 == 'E'; }
    // Lys-C cuts at carboxyl side of lysines.
    bool lysc(char aa1) { return aa1 == 'K'; }
    // Trypsin cuts Lysine or arginine and not proline.
    bool trypsin(char aa1, char aa2) { return ((aa1 == 'K' || aa1 == 'R' ) && aa2 != 'P'); }
    // Version of trypsin that cuts at proline.
    bool trypsinP(char aa1) { return aa1 == 'K' || aa1 == 'R'; }
    bool argc(char aa1, char aa2) { return (aa1 == 'R' && aa2 != 'P'); }

    bool chymotrypsin(char aa1, char aa2) { return (aa1 == 'F' || aa1 == 'Y' || aa1 == 'W' || aa1 == 'L') && aa2 != 'P'; }
    bool aspn(char aa2) { return aa2 == 'D'; }
    bool lysn(char aa2) { return aa2 == 'K'; }

    // TODO: Add these enzymes:  AspN, ArgC, chymotrypsin.
    // TODO: Also don't forget to update rev_seq. 
    
    // Recursively calculate missed cleavages.
    int missed_cleave(size_t start, size_t p, size_t pidx, int misses, bool count);

    void query_accurate_mass(double mz, int charge, vector<fragment *> &res);
    // Range query MS/MS peak based on precursor mass.
    void precursor_query(vector<fragment *> &res, double measured_mass, double neutral_mass, double modmass = 0);

    // Range query MS/MS peak based on percuror mass, subtract
    // modifications to get putative modified peptides.
    void query_modmulticomb(ms2peak *p, vector< vector<int> > &mod_mcs, vector< vector<fragment *> > &res);
    void query_modmulticomb(double mass, double neutral_mass,
			    vector< vector<int> > &mod_mcs,
			    vector< vector<fragment *> > &res);
    void group_seq(vector<int> &mc, vector<fragment *> &seq_res, 
		   vector< vector<int> > &mod_mcs, vector< vector<fragment *> > &res);

    // Remove peaks that are likely to be caused by isotopes.
    void IsotopeFilter(vector<peak2> &ms2);
    // Filter out regions of spurious peaks. 
    void PeakFilter(vector<peak2> &ms2, double win = 100, size_t K = 8, size_t F = 100);

    // Apply MS/MS spectrum filters.
    void filter(ms2peak *p);
    void pre_process(ms2peak *p) { if(p == NULL) return; filter(p); sort_by_intensity(p);  }

    // NOTE: Call pre_process on peaks before calling search!!!!
    void pre_process(vector<ms2peak *> &ps) { for(size_t i = 0; i < ps.size(); i++) pre_process(ps[i]);  }

    // Sort peaks by decreasing intensitity.
    void sort_by_intensity(ms2peak *p);

    // Match peaks between filtered spectrum and theoretical peak ladder.
    void match(vector<peak2> &filtered, vector<int> &filteredF, vector<double> &ladder, vector<bool> &ladderF);

    // For given isotope labels, check to see if sequence contains the right content.
    bool check_nitro(string &seq, int nitro_cnt);
    bool check_isotope(string &seq_orig, vector<isotope_t> &iso);

    // Run database search. 
    void search(ms2peak *p, float min_mz2, float max_mz2, vector<isotope_t> &iso, int nitro_cnt = 0); 
    void search(vector<ms2peak *> &peaks, float min_mz2, float max_mz2, vector<isotope_t> &iso, int nitro_cnt = 0);

    // MS/MS scoring.
    double ms2score(vector<peak2> &filtered, vector<double> &ladder);

    // Given a modification combination (determined by precusor mass
    // queries) enumerate all possible aa specific modifications.
    bool enumerate_mods(string &seq, vector<int> &mod_multicomb, vector<mods_t> &mods);
  };

  // Creates protein groups, assigns them to MS/MS peaks given.
  // Non-conclusive protein groups have all shared peptides. These are
  // all shared with other protein groups.
  void build_protein_groups(size_t Nproteins, vector<ms2peak *> &peaks, 
			    vector< vector<int> > &group2idx, bool form_nonconclusive = false);

  // Assigns given MS/MS peaks to a protein groups.
  void assign_protein_groups(size_t Nproteins, vector<ms2peak *> &peaks, 
			     vector< vector<int> > &group2idx, vector< vector<int> > &peak2groups);

  // Find smallest score cutoff that achieves given "PIT adjusted FDR."
  double fdr_threshold(vector<ms2peak *> &peaks, double fdr_needed);
}

#endif // __ms2_hpp__
