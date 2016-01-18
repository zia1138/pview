#include <math.h>
#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <cstring>

#include <string.h>
#include <ctype.h>

#include "ms2.hpp"
#include "util.hpp"
#include "xml.hpp"

#include "fdr.hpp"

namespace ms2 {
  using namespace std;
  
  // Functor for peak to decreasing intensity sorting.
  struct cmp_decr_intensity2 { bool operator () (peak2 a, peak2 b) { return a.intensity2 > b.intensity2; } };
  // For increasing m/z. 
  struct cmp_incr_mz2 { bool operator () (peak2 a, peak2 b) { return a.mz2 < b.mz2; } };

  // Reverse a sequence, swapping cut site according to enzyme.
  void fasta_data::rev_seq(string &seq) {
    // Apply actual reversal of the sequence.
    std::reverse(seq.begin(), seq.end());
    seq[0] = '#'; seq[seq.size()-1] = '*'; 

    // Swap the cut site with the previous amino acid.
    // NOTE: DO NOT FORGET TO UPDATE THIS WHEN ADDING ENZYME.
    for(size_t k = 1; k < seq.size()-1; k++) {
      bool do_swap = false;
      char seqk = toupper(seq[k]); // , seqkm1 = toupper(seq[k-1]);
      switch(conf.enzyme_id) {
      case ez_Trypsin: case ez_TrypsinP:
	if(seqk == 'K' || seqk == 'R') { do_swap = true; }
	break;
      case ez_GluC: if(seqk == 'E') { do_swap = true; } break;
      case ez_LysC: case ez_LysN: if(seqk == 'K') { do_swap = true; } break;
      case ez_ArgC: if(seqk == 'R') { do_swap = true; } break;
      case ez_AspN: if(seqk == 'D') { do_swap = true; } break;
      case ez_ChymoTrypsin: 
	if(seqk == 'F' || seqk == 'Y' || seqk == 'W' || seqk == 'L') 
	  { do_swap = true; } break;
      default: break;
      }
      if(do_swap) { 
	char tmp = seq[k-1]; 
	seq[k-1] = seq[k]; 
	seq[k] = tmp; 
      }
    }
  }

  // Build fragment database by sorting fragments on molecular weight.
  void fragmentDB::build() {
    // Fill molecular weight table.
    MW.resize(128, 0); nitrogens.resize(128,0);  // note use 128 here since this is the number of ascii characters in char
    init_MWs(); // Pass heavy/light flag.

    // Compute the total mass of each protein sequence (TODO use only 1 SNP).
    seq_mass.resize(fasta.seq_data.size(), 0);
    for(size_t s = 0; s < fasta.seq_data.size(); s++) {
      double mass = 0;  
      string seq = fasta.seq_data[s];
      for(size_t i = 0; i < seq.length(); i++) mass += MW[seq[i]];
      seq_mass[s] = mass + mass::water;
    }
    fragment_count.resize(fasta.seq_data.size(), 0);

    digest_all();     // Digest all sequences.

    // Store pointers to fragments in frags.
    frags.resize(frag_data.size());
    for(size_t f = 0; f < frag_data.size(); f++) frags[f] = &frag_data[f];

    // Sort the DB on mass for fast range queries.
    sort(frags.begin(), frags.end(), cmp_frag_mass()); 

    // Increment fragment count for the protein.
    for(size_t f = 0; f < frag_data.size(); f++) {
      if(frag_data[f].misses == 0) fragment_count.at(frag_data[f].protein_idx)++;
    }
  }

  // Initializes a table of mapping from a character code to molecular
  // weights.  NOTE: The 'X' code is used to designate no modification
  // for variable modifications. TODO: Figure out if that is OK.
  void fragmentDB::init_MWs() {
    using namespace mass;
    MW['A'] = C*3 + H*5 + O + N;     // 71.037114  // alanine
    MW['C'] = C*3 + H*5 + O + N + S; // 103.009185 // cysteine
    MW['D'] = C*4 + H*5 + O*3 + N;   // 115.026943 // aspartic acid
    MW['E'] = C*5 + H*7 + O*3 + N;   // 129.042593 // glutamic acid
    MW['F'] = C*9 + H*9 + O + N;     // 147.068414 // phenylalanine
    MW['G'] = C*2 + H*3 + O + N;     // 57.021464  // glycine
    MW['H'] = C*6 + H*7 + O + N*3;   // 137.05891  // histidine
    MW['I'] = C*6 + H*11 + O + N;    // 113.08406  // isoleucine
    MW['K'] = C*6 + H*12 + O + N*2;  // 128.09496  // lysine
    MW['L'] = C*6 + H*11 + O + N;    // 113.08406  // leucine
    MW['M'] = C*5+ H*9 + O + N + S;  // 131.04049  // methionine
    MW['N'] = C*4 + H*6 + O*2 + N*2; // 114.04293  // asparagine
    MW['P'] = C*5 + H*7 + O + N;     // 97.052764  // proline
    MW['Q'] = C*5 + H*8 + O*2 + N*2; // 128.05858  // glutamine
    MW['R'] = C*6 + H*12 + O + N*4;  // 156.10111  // arginine
    MW['S'] = C*3 + H*5 + O*2 + N;   // 87.03203   // serine
    MW['T'] = C*4 + H*7 + O*2 + N;   // 101.04768  // threonine
    MW['V'] = C*5 + H*9 + O + N;     // 99.06841   // valine
    MW['W'] = C*11 + H*10 + O + N*2; // 186.07931  // tryptophan
    MW['Y'] = C*9 + H*9 + O*2 + N;   // 163.06333  // tyrosine

    // one nitrogen comes from the peptide bond for each amino acid
    nitrogens['A'] = 1; // alanine
    nitrogens['C'] = 1; // cysteine
    nitrogens['D'] = 1; // aspartic acid
    nitrogens['E'] = 1; // glutamic acid
    nitrogens['F'] = 1; // phenylalanine
    nitrogens['G'] = 1; // glycine
    nitrogens['H'] = 3; // histidine
    nitrogens['I'] = 1; // isoleucine
    nitrogens['K'] = 2; // lysine
    nitrogens['L'] = 1; // leucine
    nitrogens['M'] = 1; // methionine
    nitrogens['N'] = 2; // asparagine
    nitrogens['P'] = 1; // proline
    nitrogens['Q'] = 2; // glutamine
    nitrogens['R'] = 4; // arginine
    nitrogens['S'] = 1; // serine
    nitrogens['T'] = 1; // threonine
    nitrogens['V'] = 1; // valine
    nitrogens['W'] = 2; // tryptophan
    nitrogens['Y'] = 1; // tyrosine

    if( (dbtype == dbHEAVY && conf.isotope_15N) || conf.fixed_15N) {
      cout << "using 15N fixed modifications" << endl;
      // If DB type is N15, each amino acid is heavier by its peptide
      // bond plus any additional nitrogens.
      for(int i = 0; i < 128; i++) MW[i] += double(nitrogens[i]) * mass::delta15N;
    }

    MW['#'] = 0; // N-term
    MW['*'] = 0; // C-term are unmodified
    MW['['] = 0; // cleavage N-term
    MW[']'] = 0; // cleavage C-term

    // Apply fixed modifications.
    for(size_t i = 0; i < conf.fxMods.size(); i++) {
      // Shift mass of amino acids if fixed modification is active.
      fxmod_t &fx = conf.fxMods[i];
      if(fx.active) for(size_t a = 0; a < fx.AAs.size(); a++) MW[fx.AAs[a]] += fx.deltaMass;
    }

    // Apply weight of fixed modifications for heavy database.  Do not
    // apply fixed modifications if PepXML ID file has been imported.
    if(conf.quantmode == quantAcetylSTOICH) {
      using namespace mass;
      double acetate_shift = 0;
      if(conf.acetyl_stoich_plus5) {
	if(dbtype == dbLIGHT) acetate_shift = C12*2 + H*2 + O;
	if(dbtype == dbHEAVY) acetate_shift = C13*2 + H2*3 - H + O;
      }
      else {
	if(dbtype == dbLIGHT) acetate_shift = C12*2 + H*2 + O;
	if(dbtype == dbHEAVY) acetate_shift = C12*2 + H2*3 - H + O;
      }
      cout << "acetate_shift = " << acetate_shift << " ****************************************" << endl;
      MW['K'] += acetate_shift; MW['#'] += acetate_shift;
    }

    // Fixed modifications will show up in pepxml so you need to turn
    // off the fixed modification in the configuration.
    if(dbtype == dbHEAVY && conf.isotope_15N == false) { // Use heavy shifts.
      for(size_t i = 0; i < active_iso.size(); i++) {
	isotope_t &iso = active_iso[i];
	for(size_t a = 0; a < iso.AAs.size(); a++)  MW[iso.AAs[a]] += iso.shift_heavy;
      }
    }
    if(dbtype == dbMEDIUM) { // Use medium shifts.
      for(size_t i = 0; i < active_iso.size(); i++) {
	isotope_t &iso = active_iso[i];
	for(size_t a = 0; a < iso.AAs.size(); a++)  MW[iso.AAs[a]] += iso.shift_medium;
      }
    }
    // Save active variable modification set.
    MODs.push_back(varmod_t("", "X", 0)); // 0th mod is always the no modification modification

    // Add user specified modifications.
    for(size_t i = 0; i < conf.varMods.size(); i++) {
      varmod_t &M = conf.varMods[i]; 
      if(M.active) MODs.push_back(M);
    }
  }

  // Collect a, b, or c-ion ladders.
  void fragmentDB::abc_ion(double mwdelta, int charge, string &seq, mods_t &m, vector<double> &ladder, 
			   double mz1, double mz2, int cutoff, bool exclude_proline) {  // Compute b-ion ladder.
    if(cutoff < 0) cutoff = 0; // Cutoff eliminates abcions to the left of the position on the sequence.
    double mw = 0;
    for(int i = 0; i < (int)seq.length(); i++) {
      if(m.is_modified(i)) mw += m.get_frag_mass(i);
      else mw += MW[seq[i]];
      double m = (mw + mwdelta) + mass::proton; 
      if(charge == 2) m = ((mw + mwdelta) + 2.0 * mass::proton) / 2.0;
      // This checks to see if ions are within the m/z window [mz1,
      // mz2] and that thte position in the sequence is to the right
      // of the cutoff.  For a/b/c ions position to the right of the
      // cutoff are site determining.
      if(mz1 <= m && m <= mz2 && i >= cutoff) {
	if(exclude_proline) {
	  if(i + 1 < (int)seq.length() && seq[i+1] != 'P') ladder.push_back(m);
	}
	else ladder.push_back(m);
      }
    }    
  }

  // Collect x,y, or z-ion ladders.
  void fragmentDB::xyz_ion(double mwdelta, int charge, string &seq, mods_t &m, vector<double> &ladder, 
			   double mz1, double mz2, int cutoff, bool exclude_proline) {
    if(cutoff < 0) cutoff = seq.length(); // Eliminates x,y,or z-ions to the right of a given position on the sequence.
    double mw = 0;
    for(int i = 0; i < (int)seq.length(); i++) {
      int p = (int)seq.length() - i - 1;
      if(m.is_modified(p)) mw += m.get_frag_mass(p);
      else mw += MW[seq[p]];
      double m = (mw + mwdelta) + mass::proton + mass::water;
      if(charge == 2) m = ((mw + mwdelta) + 2.0 * mass::proton + mass::water) / 2.0;
      // Make sure in m/z range.  Keep positions to the left of the
      // cutoff site, For x/y/z ions these positions are site
      // determining.
      if(mz1 <= m && m <= mz2 && p <= cutoff) {
	if(exclude_proline) {
	  if(p > 0 && seq[p] != 'P') ladder.push_back(m); 
	}
	else ladder.push_back(m);      
      }
    }
  }

  /* Sequences use the standard N-term to C-term #CCTKTLKR* convention
     where # is the N-term and * is the C-term.

     Given a digested sequence. ] and [ designate the new C-term and N-term.

     #CCTK] [TLKR*   # N-term original and * C-term original
     ] new C-term cleavage and [ new N-term cleavage     
     
     Digest actually generates the fragments, accounting for missed cleavages.
  */
  void fragmentDB::digest_all() { 
    // First count the number of fragments generated for the given
    // enzyme and missed cleavage count.
    frag_cnt = 0;
    cout << "running digest count" << endl;
    for(size_t idx = 0; idx < fasta.seq_data.size(); idx++) digest(idx, true); 
    // Allocate memory for these fragments and generate.
    cout << "frag_cnt=" << frag_cnt << endl;
    frags.reserve(frag_cnt + 2);
    for(size_t idx = 0; idx < fasta.seq_data.size(); idx++) digest(idx); 
    cout << "*DONE* digest_all" << endl;
  }
  void fragmentDB::digest(size_t pidx, bool count) {
    size_t pos = 0; 
    while(pos < fasta.seq_data[pidx].length()) pos = missed_cleave(pos, pos, pidx, conf.missed_cleaves, count);
  }

  int fragmentDB::missed_cleave(size_t start, size_t p, size_t pidx, int misses, bool count) {
    string seq = fasta.seq_data[pidx];
    if(misses < 0 || p == seq.length()) return p;

    // Advance through the current sequence.
    while(p < seq.length()) {
      p++; // Advance through sequence.
      // Use specified enzyme to cut.
      if( p == seq.length() || enzyme(seq[p-1], seq[p]) ) {
	// Check length of the fragment.
	int aa_len = p - start;
	if(conf.min_aa <= aa_len && aa_len <= conf.max_aa) {
	  // Compute mass of the computed fragment. Handle cleavage N-term.
	  double m = 0;
	  bool badAA = false;
	  if(start != 0) m = MW['['];
	  for(size_t i = start; i < p; i++) {
	    if(isAAcode(seq[i]) == false) badAA = true;
	    m += MW[seq[i]];
	  }
	  if(p != seq.length() - 1) m += MW[']']; 	// Handle cleavage C-term.
	  m += mass::water; // Peptide bonds loose water. These are the ends of the peptide.
	  // If a bad amino acid code was found in this fragment, discard. 
	  if(badAA == false) {
	    if(count) frag_cnt++; 
	    else {
	      // Count the number of nitrogens in this fragment.
	      int nitro_cnt = 0; for(size_t i = start; i < p; i++) nitro_cnt += nitrogens[seq[i]];
	      fragment newfrag(pidx, start, aa_len, m, conf.missed_cleaves - misses, fasta.is_decoy_protein.at(pidx));
	      newfrag.seq_idx = pidx;
	      newfrag.nitrogens = nitro_cnt;
	      frag_data.push_back(newfrag);
	    }
	  }
	}
	// Found a digested fragment, stop.
	break; 
      }
    }

    // Address missed cleavages by continuing cut from position p, but
    // keep start position the same. Also decrease the number of
    // misses left.
    missed_cleave(start, p, pidx, misses - 1, count);

    return p; // Return position of first cleavage so algorithm can continue from that position.
  }

  // Query tryptic fragment database on precursor mass.
  //void fragmentDB::precursor_query(ms2peak *p, vector<fragment *> &res, double modmass) {
  void fragmentDB::precursor_query(vector<fragment *> &res, double measured_mass, double neutral_mass, double modmass) {
    // NOTE: We use p->mass() here to compute the neutral mass search
    // range because this is the actual experimental value
    // measured. It makes a slightly larger window than the neutral
    // mass (without proton mass).
    double tol = mass::ppm2massError(measured_mass, conf.mstol);

    // Subtract out the expected modification mass. Gives you the
    // unmodified mass.
    double mass = neutral_mass - modmass; 
    fragment querymass(mass - tol);
    // Binary search to find tryptic fragments.
    vector<fragment*>::iterator it = lower_bound(frags.begin(), frags.end(), &querymass, cmp_frag_mass());
    while(it != frags.end()) {
      // Collect fragments in query window.
      fragment *r  = *it;
      if(r->fragmass >= mass + tol) break;
      res.push_back(r);
      it++;
    }
  }

  // Peak with signal/nonsignal label.
  struct peak2S {
    peak2S() { signal = false; }
    peak2S(peak2 p_, bool signal_) { peak = p_; signal = signal_; }
    peak2 peak; bool signal;
  };

  // Sort signal/nonsignal peaks decreasing on intensity.
  struct cmp_intensity2S { bool operator() (peak2S *a, peak2S *b) { return a->peak.intensity2 > b->peak.intensity2; } };

  // Sort signal/nonsignal peaks increasing on m/z.
  struct cmp_mz2S { bool operator() (peak2S *a, peak2S *b) { return a->peak.mz2 < b->peak.mz2; } };

  // Isotope filter primarily optimized for ion-trap data. O(N log N).
  void fragmentDB::IsotopeFilter(vector<peak2> &ms2) {
    size_t N = ms2.size();

    // Keep track of peaks with a signal/non-signal label.
    vector<peak2S> specS; specS.reserve(N);
    for(size_t i = 0; i < N; i++) specS.push_back(peak2S(ms2[i], true));

    // Sort labeled peaks by intensity.
    vector<peak2S*> by_intensity(N);
    for(size_t i = 0; i < N; i++) by_intensity[i] = &specS[i];
    sort(by_intensity.begin(), by_intensity.end(), cmp_intensity2S());

    // Now sort labelled peaks by m/z. 
    vector<peak2S*> by_mz2 = by_intensity;
    sort(by_mz2.begin(), by_mz2.end(), cmp_mz2S());

    // Scan peaks from higest intensity to lowest.
    // NOTE: Everything is labelled as signal initially.
    for(size_t i = 0; i < N; i++) {
      // Skip already filtered peaks.
      if(by_intensity[i]->signal == false)  continue;

      double tol = mass::ppm2da(by_intensity[i]->peak.mz2, conf.ms2tol); // use PPMs by default
      if(conf.ms2tol_usedaltons) tol = conf.ms2tol; // use daltons if requested      

      for(int charge = 1; charge <= 2; charge++) {
	// Only get up to the second isotope. 
	for(int iso = 1; iso <= 2; iso++) {
	  double mz2iso = by_intensity[i]->peak.mz2 + double(iso) * mass::delta13C / double(charge);
	  double start_mz = mz2iso - tol;
	  double end_mz   = mz2iso + tol;
	  peak2S query; query.peak.mz2 = start_mz;
	  vector<peak2S*>::iterator it = lower_bound(by_mz2.begin(), by_mz2.end(), &query, cmp_mz2S() ); 
	  while(it != by_mz2.end()) {
	    peak2S *res = *it;
	    if(res->peak.mz2 > end_mz) break; // +2 + tol ]
	    // If the intensity is less than the current peak,  mark as non-signal.
	    if(res->peak.intensity2 < by_intensity[i]->peak.intensity2) { res->signal = false; }
	    it++; // Advance to next peak until out of range.
	  }
	}
      }
    }

    // Save signal peaks. 
    ms2.clear(); vector<peak2> (ms2).swap(ms2); // Free memory.
    for(size_t i = 0; i < N; i++) 
      if(by_intensity[i]->signal) ms2.push_back(by_intensity[i]->peak);  
  }


  // Keep peak if surrounded by less than K=8 peaks of higher
  // intensity in m/z window (default of 100 daltons) .  Do this until
  // threshold F = 100 peaks marked as keep or all peaks processed.
  // O( N log N) algorithm.  This algorithm cleans up ion trap data
  // pretty nicely.
  void fragmentDB::PeakFilter(vector<peak2> &ms2, double win, size_t K, size_t F) {
    size_t N = ms2.size();

    // Keep track of peaks with a signal/non-signal label.  Initially
    // all peaks are non-signal.
    vector<peak2S> specS; specS.reserve(N);
    for(size_t i = 0; i < N; i++) specS.push_back(peak2S(ms2[i], false));

    // Sort labeled peaks by intensity.
    vector<peak2S*> by_intensity(N);  for(size_t i = 0; i < N; i++) by_intensity[i] = &specS[i];
    sort(by_intensity.begin(), by_intensity.end(), cmp_intensity2S());
    
    // Now sort labelled peaks by m/z. 
    vector<peak2S*> by_mz2 = by_intensity;  sort(by_mz2.begin(), by_mz2.end(), cmp_mz2S());

    size_t keep = 0; // Keep nothing initially.
    // Scan peaks from higest intensity to lowest.
    for(size_t i = 0; i < N; i++) {
      // Query peaks in a window of size [-win / 2, +win / 2]
      peak2S query = *by_intensity[i];
      query.peak.mz2 -= win / 2;  // [ -100/2

      size_t high_count = 0;
      vector<peak2S*>::iterator it = lower_bound(by_mz2.begin(), by_mz2.end(), &query, cmp_mz2S() );
      while(it != by_mz2.end()) {
	// Use binary search to find peaks in a window.
	peak2S *res = *it;

	// Check to see if out of range.
	if(res->peak.mz2 > by_intensity[i]->peak.mz2 + win / 2) break; // +100/2]

	// NOTE: We'll get back the current peak, surely.  Count th
	// enumber of peaks in that window with higher intensity than
	// the current peak.
	if(res->peak.intensity2 > by_intensity[i]->peak.intensity2) { high_count++; }

	it++; // Advance to next peak until out of range.
      }
      if(high_count <= K) { 
	// If there are no more than K peaks of higher intensity
	// around the current peak, mark as signal.
	keep++; 
	by_intensity[i]->signal = true;
      }
      if(keep >= F) break; // F = 100 peaks marked as keep, then quit.
    }

    // Save signal peaks. 
    ms2.clear(); vector<peak2> (ms2).swap(ms2); // Free memory.
    for(size_t i = 0; i < N; i++) if(by_intensity[i]->signal) ms2.push_back(by_intensity[i]->peak);  
  }

  // Apply MS/MS spectrum filters.  
  void fragmentDB::filter(ms2peak *p) {
    vector<peak2> &filtered = p->ms2;
    IsotopeFilter(filtered); // Remove peaks likely caused by isotopes.
    PeakFilter(filtered);    // Remove dense regions of low intensity spurious peaks.
  }

  // Assumes ladder is sorted and filtered is decreasing on intensity!!!
  void fragmentDB::match(vector<peak2> &filtered, vector<int> &filteredF, 
			 vector<double> &ladder, vector<bool> &ladderF) {
    // All ladder peaks initially unmatched. 
    ladderF.resize(ladder.size(), false);

    // All peaks are initially unmatched. 
    filteredF.resize(filtered.size(), -1);

    // Iterates through peaks from most intense to least intense.
    for( size_t f = 0; f < filtered.size(); f++) {
      // Range query theoretical spectrum.
      peak2 &maxI = filtered[f];
      double mztol = mass::ppm2da(maxI.mz2, conf.ms2tol);
      if(conf.ms2tol_usedaltons)  mztol = conf.ms2tol;
      vector<double>::iterator it = lower_bound(ladder.begin(), ladder.end(),  maxI.mz2 - mztol);
      while(it != ladder.end()) {
	// Did the theory peak match a filtered peak?
	if(*it <= maxI.mz2 + mztol) {
	  // Yes. Then, don't match it again. This assures only most
	  // intense peak is assigned to a theoretical peak.
	  size_t L = distance(ladder.begin(), it);
	  if(ladderF[L] == false) {
	    ladderF[L] = true; 
	    // Mark current signal peak as matched and
	    filteredF[f] = L; // save which peak in the theoretical ladder this peak matched.
	  }
	}
	// Keep advancing through matched ladder peaks, marking them
	// as matched. 2 ladder peaks can be matched to a single
	// emprical peak.
	if(*it > maxI.mz2 + mztol) break;
	it++;
      }
    }
  }

  // Sorts decreasing by intensity. Note match() assumes that the
  // peaks have all been sorted by intesnity.
  void fragmentDB::sort_by_intensity(ms2peak *p) {
    vector<peak2> &filtered = p->ms2;
    // Filtered peaks so they are decreasing on intensity.
    sort(filtered.begin(), filtered.end(), cmp_decr_intensity2());
  }

  // Builds theoretical spectrum for database matching.  Ladders were
  // adapted from the OMSSA algorithm.

  void fragmentDB::build_ladder(activationMethod_t method, vector<double> &ladder, string &seq, mods_t &mods, int charge, 
				double m, float min_mz2, float max_mz2, int cutoff) {
    double abc_shift = 0, xyz_shift = 0;
    bool exclude_proline = false;
    // TODO: For ETD we need to do something more complex. Fuck. 
    // - Add the tricks in this paper for improved +2 and >+3 
    //    performance on ETD and ETCaD (http://pubs.acs.org/doi/abs/10.1021/pr100648r)
    if(method == ETD_activation) {
      abc_shift = mass::N + 3*mass::H;
      xyz_shift = -(mass::N + 3*mass::H) + mass::proton; // zradical ion
      exclude_proline = true;
    }
    if(charge < 3) { 
      // +1 +2 build (a,b,c)-ion (x,y,z)-ion ladders
      abc_ion(abc_shift, 1, seq, mods, ladder, min_mz2, max_mz2, cutoff, exclude_proline);  
      xyz_ion(xyz_shift, 1, seq, mods, ladder, min_mz2, max_mz2, cutoff, exclude_proline);
      if(method == ETD_activation) {
	// throw in some y-ions	(this definitely helps)
	xyz_ion(0, 1, seq, mods, ladder, min_mz2, max_mz2, cutoff, exclude_proline);
      }
    } 
    else {
      // Add b2+ ions and y2+ ions using OMSSA strategy.
      abc_ion(abc_shift, 1, seq, mods, ladder, min_mz2, max_mz2, cutoff, exclude_proline); 
      xyz_ion(xyz_shift, 1, seq, mods, ladder, min_mz2, max_mz2, cutoff, exclude_proline);
      abc_ion(abc_shift, 2, seq, mods, ladder, min_mz2, m/2, cutoff, exclude_proline); 
      xyz_ion(xyz_shift, 2, seq, mods, ladder, min_mz2, m/2, cutoff, exclude_proline);
    }
    // Sort ladder from smallest m/z to largest m/z.
    sort(ladder.begin(), ladder.end());
  }

  // O(N^2) algorithm groups redudant sequences returned by the range query.
  void fragmentDB::group_seq(vector<int> &mc, vector<fragment *> &seq_res, 
			      vector< vector<int> > &mod_mcs, vector< vector<fragment *> > &res) {
    string seq, seq2;
    vector<bool> grouped(seq_res.size(), false);
    for(size_t i = 0; i < seq_res.size(); i++) {
      if(grouped[i]) continue; // Skip processed fragment
      grouped[i] = true;
      // Get sequence.
      seq.clear(); fasta.get_seq_nomod(seq_res[i], seq); 

      // Add current sequence to set of results.
      vector<fragment *> res_m; res_m.push_back(seq_res[i]);

      // Find duplicate sequences in the results.
      for(size_t j = 0; j < seq_res.size(); j++) {
	if(grouped[j]) continue;
	// Are they the same sequence(?)
	seq2.clear(); fasta.get_seq_nomod(seq_res[j], seq2); 
	if(seq == seq2) {   // Same sequence, group together.
	  res_m.push_back(seq_res[j]);
	  grouped[j] = true;
	} 
      }
      // Now save multicombination w/ modification types + fragment/sequence set. 
      mod_mcs.push_back(mc);
      res.push_back(res_m);
    }
  }

  // Returns fragments querried just on accurate mass.
  // NOTE: Computes neutral (positive mode) mass!!!!
  void fragmentDB::query_accurate_mass(double mz, int charge, vector<fragment *> &res) {
    if(charge <= 0) return;
    double mass = double(charge) * mz;
    double neutral_mass = double(charge) * (mz  - mass::proton);
    vector< vector<int> > query_mcs;
    vector< vector<fragment *> > query_res;
    query_modmulticomb(mass, neutral_mass, query_mcs, query_res);
    for(size_t q = 0; q < query_res.size(); q++) {
      for(size_t r = 0; r < query_res[q].size(); r++) {
	res.push_back(query_res[q][r]);
      }
    }
  }

  // DB range query that subtracts out modification mass and gets
  // putative unmodified fragments.  Returns multicombination of
  // modifications for each multicombination it returns all fragments
  // with same sequence content.  One vector<int> entry in mod_mcs has
  // a multicombination of modifications from MODs.  and the
  // corresponding vector<fragment *> entry in res has all of the
  // redundant fragments this mulitcombination applies to.
  void fragmentDB::query_modmulticomb(double mass, double neutral_mass,
				      vector< vector<int> > &mod_mcs,
				      vector< vector<fragment *> > &res) {
    // Generate a multicombination of modification types. From a menu
    // of MODs.size() active modifications (which includes the 0th
    // no-mod mod) choose up to max_mod_aa mods w/ repetition.
    vector< vector<int> > multi_combs; 
    util::multi_comb gen(MODs.size(), conf.max_mod_aa, multi_combs);

    vector<int> modCounts(MODs.size(), 0); 
    vector< fragment *> query_res, seq_res;
    string seq;

    for(size_t m = 0; m < multi_combs.size(); m++) {
      vector<int> &mc = multi_combs[m];   // Get a generated multicombination.

      // Count the number of each type of modification present in the multicombination.
      for(size_t i = 0; i < modCounts.size(); i++) modCounts[i] = 0;
      for(size_t i = 0; i < mc.size(); i++) modCounts[mc[i]]++;

      // Skip the multicombination if we reached a hard limit on a specific modification type.
      bool skip_multicomb = false;
      for(size_t i = 0; i < modCounts.size(); i++) {
	if(MODs[i].max_count != 0 && modCounts[i] > MODs[i].max_count) { skip_multicomb = true; break;}
      }
      if(skip_multicomb) continue;

      // Compute the total modificaton mass (modification to the precursor mass).
      double delta_mass = 0; for(size_t k = 0; k < mc.size(); k++) delta_mass += MODs[mc[k]].d_prec_mass;

      // Query fragments with this shifted mass.
      query_res.clear(); precursor_query(query_res, mass, neutral_mass, delta_mass);
      
      // Examine the sequence content of each returned result.
      seq_res.clear();
      for(size_t q = 0; q < query_res.size(); q++) {
	// Get sequence information for the fragment.
	seq.clear(); fasta.get_seq_nomod(query_res[q], seq);
	bool keep_res = true;
	for(size_t k = 0; k < mc.size(); k++) {
	  if(mc[k] == 0) continue; // skip the no-mod mod

	  // For each modification, search the sequence for a 
	  // matching amino acid. 
	  bool found = false;
	  for(size_t s = 0; s < seq.length(); s++) {
	    if(MODs[mc[k]].is_modifiedAA(seq[s])) {
	      // If it is found, replace with a space and advance to next modification.
	      seq[s] = ' '; found = true; break;
	    }
	  }
	  // If a modifiable amino acid is not found, don't keep this
	  // sequence result.
	  if(!found) { keep_res = false; break; }
	}
	if(keep_res) seq_res.push_back(query_res[q]);
      }
      // If sequences found for this modification multicombination,
      // keep the results.
      if(seq_res.size() > 0 ) {
	mod_mcs.push_back(mc);
	res.push_back(seq_res);
      }
    }    
  }

  // Performs precursor mass query. Groups redudant sequences.
  // mod_mcs has a indicies in ms2::MODs and for each index set, res
  // has a list of fragments.
  void fragmentDB::query_modmulticomb(ms2peak *p, vector< vector<int> > &mod_mcs, vector< vector<fragment *> > &res) {
    vector< vector<int> > query_mcs;
    vector< vector<fragment *> > query_res;
    query_modmulticomb(p->mass(), p->neutral_mass(), query_mcs, query_res);

    // Group redudant sequences to avoid additional enumeration at fragment level.
    for(size_t q = 0; q < query_mcs.size(); q++) 
      group_seq(query_mcs[q], query_res[q], mod_mcs, res);
  }
  
  // Implements a simple MS/MS scoring scheme.
  double fragmentDB::ms2score(vector<peak2> &filtered, vector<double> &ladder) {
    if(filtered.size() == 0) return 0;

    // Match peaks between filtered spectrum and theoretical spectrum.
    vector<bool> theory_flag; vector<int> is_matched; 
    match(filtered, is_matched, ladder, theory_flag);

    // Count the number of theoretical peaks matched.
    int theory_matched = 0;
    for(size_t i = 0; i < theory_flag.size(); i++) if(theory_flag[i]) theory_matched++;

    // Compute the logI intensity of matched peaks.
    double logImatched = 0;    
    for(size_t i = 0; i < is_matched.size(); i++) 
      if(is_matched[i] >= 0) logImatched += util::Log2(filtered[i].intensity2);

    return logImatched * (double(theory_matched) / double(ladder.size()));
  }


  // Enumerates modifications.  NOTE: This code is tricky. Are there
  // any ways to simplify and limit enumeration?
  // Returns false if "too many" modifications enumerated.
  bool fragmentDB::enumerate_mods(string &seq, vector<int> &mod_multicomb, vector<mods_t> &mods) {       
    // Compute the cardinality of each modification type.
    vector<int> mod_card(MODs.size(), 0); 
    int total_mods = 0; // Also count the total number of modifications.
    for(size_t k = 0; k < mod_multicomb.size(); k++) {
      // Skip the no-mod modification (0th modification)
      if(mod_multicomb[k] != 0) { mod_card[mod_multicomb[k]]++; total_mods++; } 
    }
    // If no modifications, return (no modification).
    if(total_mods == 0) { mods.push_back(mods_t()); return true; }

    // Get pointers to active modifications. 
    vector<varmod_t*> active; 
    vector<bool> backref;
    for(size_t i = 0; i < MODs.size(); i++) {
      if(mod_card[i] > 0) {
	backref.push_back(false);
	active.push_back(&MODs[i]);
	for(int card = 1; card < mod_card[i]; card++) {
	  backref.push_back(true); // When backref is set to true, the cross product
	  // starts with the value of the previous slot. Mixes combinations with crosses.
	  // This allows proper handling repeats of the same modification type.
	  active.push_back(&MODs[i]);
	  // active will have multiple references to a modification we
	  // need to do this because we we allow only one modification
	  // per amino acid.
	}
      }
    }
    // Get the positions each active modification can modify.
    vector< vector<int> > positions(active.size()); 
    for(size_t pos = 0; pos < seq.length(); pos++) {
      for(size_t a = 0; a < active.size(); a++) {
	// Can this active mod modify position pos? If yes, save for that active modification.
	// Note that the ordering is critical here for the backref to work correctly.
	if(active[a]->is_modifiedAA(seq[pos])) positions[a].push_back(pos);
      }
    }
    // Compute (constrained) cross product of all modified positions.
    // We use this constrained approach because 
    vector<int> poscnts(positions.size());
    for(size_t a = 0; a < positions.size(); a++) poscnts[a] = positions[a].size();
    // Back ref will prevent the same modification from being applied
    // twice to the same amino acid (e.g. phospho 1 to position 1 and
    // phospho 2 to position 2 then phospho 1 to position 2 and
    // phospho 2 to position 1 (which are both equivalent).
    vector< vector<int> > posidxs; util::cross c(poscnts, backref, posidxs);

    // Cross each fragment modification type from each active
    // modification.  This correclty handles neutral losses because an
    // active modification can have several different fragment mass
    // shifts.
    vector<int> fragmodcnts(active.size());
    for(size_t a = 0; a < active.size(); a++) fragmodcnts[a] = active[a]->fragmods.size(); 
    vector< vector<int> > fragmodidxs; util::cross(fragmodcnts, fragmodidxs);

    vector<int> modpos;
    vector<modpos_t> mps(active.size());
    for(size_t p = 0; p < posidxs.size(); p++) {
      vector<int> &posidx = posidxs[p]; // Get position indexes.

      // Collect the actual modification positions in the sequence.
      bool doublemod = false;
      modpos.clear();
      for(size_t a = 0; a < posidx.size(); a++) {
	int pos = positions[a][posidx[a]];
	// Look for positions that are modified twice by two different active modifications.
	bool foundpos = false;	
	for(size_t m = 0; m < modpos.size(); m++) if(modpos[m] == pos) { foundpos = true; break; }
	if(!foundpos) modpos.push_back(pos);
	else { doublemod = true; break; } 
      }
      // This check assures we have a perfect matching.
      if(doublemod) continue;  // Skip this position set since a single position is double modified.

      // modpos now has a set of unique positions that are modified, save those positions.
      for(size_t a = 0; a < modpos.size(); a++) mps[a].pos = modpos[a];

      // Iterate through cross of fragment level modifications.
      for(size_t f = 0; f < fragmodidxs.size(); f++) {
	vector<int> &fragmodidx = fragmodidxs[f];

	// Use fragment modification to check amino acid content.
	bool badamino = false;
	for(size_t a = 0; a < fragmodidx.size(); a++) {
	  fragmod_t &fragmod = active[a]->fragmods[fragmodidx[a]];
	  if(fragmod.is_modifiedAA(seq[modpos[a]]) == false) { badamino = true; break; }
	  else {
	    // Save precursor mass and identifier.
	    mps[a].mod_idx = active[a]->idx;
	    mps[a].prec_mass = MW[seq[modpos[a]]] + active[a]->d_prec_mass;
	    // Save fragment mass and identifier.
	    mps[a].fragmod_idx = fragmod.idx;
	    mps[a].frag_mass = MW[seq[modpos[a]]] + fragmod.dfrag_mass;
	  }
	}
	if(!badamino) {
	  // Everything checks out, save modification set.
	  mods_t m; m.sites = mps;
	  mods.push_back(m);
	  // TODO: This is a total kludge. Figure out what to do when
	  // enumerate goes crazy!!!  This is not the best strategy
	  // for dealing with lots of mods and mod types.
	  if(mods.size() > conf.max_mod_enumerate) {
	    cout << "max_mod_enumerate" << endl;
	    mods.clear();
	    return false;
	  }
	}
      }
    }
    return true;
  }

  bool fragmentDB::check_nitro(string &seq, int nitro_cnt) {
    int seq_nitro = 0;
    for(size_t i = 0; i < seq.length(); i++) seq_nitro += nitrogens[seq[i]];
    return nitro_cnt == seq_nitro;
  }

  // Checks sequence for correct number of isotope labeled amino acids. 
  bool fragmentDB::check_isotope(string &seq_orig, vector<isotope_t> &iso) {
    string seq = seq_orig;
    // Account for all amino acids as determined by isotope pairing. 
    for(size_t s = 0; s < iso.size(); s++) {
      // For each active isotope find the corresponding amino acid.
      bool found = false;
      if(iso[s].AAs.size() == 0) found = true; // Handle case where no AAs for this shift.
      else {
	// Replace it with a space so subsequent isotopes won't match.
	for(size_t i = 0; i < seq.length(); i++) 
	  if(iso[s].find_aa(seq[i])) { found = true; seq[i] = ' ';  break; }
      }
      // Did not find a matching amino acid in this sequence, return false.
      if(!found) return false; 
    }
    for(size_t s = 0; s < iso.size(); s++) {
      if(iso[s].AAs.size() == 0) continue; // Handle case where there are no AAs for this shift.
      // Look for additional instances of the isotope, there should be none.
      for(size_t i = 0; i < seq.length(); i++) if(iso[s].find_aa(seq[i])) return false;
    }
    // Otherwise, everything is OK.
    return true;
  }


  // NOTE: Assumes MS/MS peaks in ms2peak *p are sorted by increasing intensity!!!!
  // CALL pre_process() first!!!!!
  void fragmentDB::search(ms2peak *p, float min_mz2, float max_mz2, vector<isotope_t> &iso, int nitro_cnt) {
    if(p->charge <= 0) return; // Skip invalid charge.

    // Do range query for multi-combinations of modification types.
    vector< vector<int> > mod_mcs;   // multicombinations of precursor modifications
    vector< vector<fragment *> > grouped_res;  // for each modification, a set of the same fragments by sequence (not by protein)
    query_modmulticomb(p, mod_mcs, grouped_res);

    vector<double> ladder;  vector<mods_t> mods; 
    string seq;
    vector<fragment *> res_target, res_decoy;

    // Iterate through each modification multicombination + sequences that match.
    for(size_t t = 0; t < mod_mcs.size(); t++) {
      // For the multicombination get grouped fragments (with same sequence).
      vector<fragment *> &res = grouped_res[t]; 

      // Get unmodified sequence data for the 0th fragment.
      seq.clear(); fasta.get_seq_nomod(res[0], seq); 

      // Check nitrogen content.
      if(nitro_cnt > 0 && check_nitro(seq, nitro_cnt) == false) continue;

      // Check sequence content for given isotopes (for SILAC enough Ks and Rs based on spacing between XICs).
      if(iso.size() > 0 && check_isotope(seq, iso) == false) continue;

      // Enumerate combinations and cross products of possible
      // modifications using selected multi-combination.
      mods.clear(); 
      if(enumerate_mods(seq, mod_mcs[t], mods) == false) { p->clear_ms2id(); break; }

      // Split results into decoy set and target set. See scoring below.
      res_target.clear(); res_decoy.clear();
      for(size_t i = 0; i < res.size(); i++) {
	if(res[i]->is_decoy) res_decoy.push_back(res[i]);
	else res_target.push_back(res[i]);
      }

      for(size_t m = 0; m < mods.size(); m++) {
	// Calculate mass ladders of returned tryptic fragment.
	ladder.clear(); 
	// TODO: Scan mod_mcs and see if we have to add a precursor peak.
	build_ladder(p->activation, ladder, seq, mods[m], p->charge, p->neutral_mass(), min_mz2, max_mz2);
	double score = ms2score(p->ms2, ladder); // Compute MS/MS score.
	// Compare current score to current target and decoy scores.
	double dtarget_score = score - p->max_target_score(), ddecoy_score = score - p->max_decoy_score();

	// Scores are strictly greater, clear out previous set.
	if(dtarget_score > 1e-7) p->clear_target(); if(ddecoy_score > 1e-7) p->clear_decoy();

	// We need only keep the maximum scoring decoy and target hits
	// for correct FDR calculation. This ends up saving a lot of
	// memory.  Need to check for -1e-7 since no such thing as
	// quality for floating point numbers.
	if(dtarget_score >= -1e-7) { // Allows through equal to and greater
	  for(size_t i = 0; i < res_target.size(); i++) {
	    ms2id_t id; // Save MS/MS ID.
	    id.score = score; id.frag = *(res_target[i]); id.mods = mods[m];
	    p->target.push_back(id);
	  }
	}
	if(ddecoy_score >= -1e-7) { // Allows through equal to or greater with epsilson of +-1e-7
	  for(size_t i = 0; i < res_decoy.size(); i++) {
	    ms2id_t id; // Save MS/MS ID.
	    id.score = score; id.frag = *(res_decoy[i]); id.mods = mods[m];
	    p->decoy.push_back(id);
	  }
	}
	// Advance to next modification set.
      }
    }
    
    // TODO: Change this to what fraction of the search score is
    // determined by site determining ions.
    // Currently setting the site score to 0.
    for(size_t t = 0; t < p->target.size(); t++) {
      vector<modpos_t> &sites = p->target[t].mods.sites;
      for(size_t i = 0; i < sites.size(); i++) sites[i].site_score = 0;
    }
    for(size_t d = 0; d < p->decoy.size(); d++) {
      vector<modpos_t> &sites = p->decoy[d].mods.sites;
      for(size_t i = 0; i < sites.size(); i++) sites[i].site_score = 0;
    }    
    // Now score each modification position.    
    /*for(size_t t = 0; t < p->target.size(); t++) {
      vector<modpos_t> &sites = p->target[t].mods.sites;
      seq.clear(); fasta.get_seq_nomod(&p->target[t].frag, seq);
      for(size_t i = 0; i < sites.size(); i++) {
	// Build ladder with only site determining product ions. 
	ladder.clear();
	build_ladder(p->activation, ladder, seq, 
		     p->target[t].mods, p->charge, p->neutral_mass(), 
		     min_mz2, max_mz2, sites[i].pos);
	// If there are two sites have about the same score, then it's
	// ambigious. 
	// Score those ions.
	sites[i].site_score = ms2score(p->ms2, ladder);
      }
      }*/
    // Repeat for decoy sites. 
    /*for(size_t d = 0; d < p->decoy.size(); d++) {
      vector<modpos_t> &sites = p->decoy[d].mods.sites;
      seq.clear(); fasta.get_seq_nomod(&p->decoy[d].frag, seq);
      for(size_t i = 0; i < sites.size(); i++) {
	ladder.clear();
	build_ladder(p->activation, ladder, seq, 
		     p->decoy[d].mods, p->charge, p->neutral_mass(), 
		     min_mz2, max_mz2, sites[i].pos);
	sites[i].site_score = ms2score(p->ms2, ladder);
      }
      }*/
    
  }

  // Compute the precursor mass of the fragment + the given modifications.
  double fragmentDB::precursor_mass(fragment *f, mods_t &m){
    double mass = 0;
    string seq; fasta.get_seq_nomod(f, seq); 
    for(size_t i = 0; i < seq.size(); i++) {
      if(m.is_modified(i)) mass += m.get_prec_mass(i);
      else mass += MW[seq[i]];
    }
    return mass + mass::water;
  }

  int fragmentDB::nitrogen_cnt(fragment *f, mods_t & /*m*/ ) {
    int nitro_cnt = 0;
    // TODO: Account for nitrogens in modifications.
    string seq; fasta.get_seq_nomod(f, seq); 
    for(size_t i = 0; i < seq.size(); i++) {
      nitro_cnt += nitrogens[seq[i]];
    }
    return nitro_cnt;
  }

  // Only do this for charge +2 peptides!!!  Compute mass errors in MS/MS spectra.
  void fragmentDB::ms2_mass_error(ms2peak *p, vector<double> &da_errors, 
				  vector<double> &ppm_errors, float min_mz2, float max_mz2) {
    if(p->charge != 2 || p->num_ids() == 0) return; // Must be charge 2 and must have IDs.

    // Get highest scoring MS/MS ID. 
    vector<ms2::fragment *> ids;  p->max_ids(ids); 
    vector<ms2::mods_t> mods; p->max_mods(mods);

    // Compute modified sequence information.
    string seq; fasta.get_seq_nomod(ids[0], seq); 

    // Calculate mass ladders of returned tryptic fragment.
    vector<double> ladder; 
    build_ladder(p->activation, ladder, seq, mods[0], p->charge, p->neutral_mass(), min_mz2, max_mz2);

    vector<peak2> &filtered = p->ms2; // Get peaks that survived filtering in pre_process.

    // Match the peaks.
    vector<bool> theory_flag; vector<int> is_matched; 
    match(filtered, is_matched, ladder, theory_flag);    

    // Compute mass errors for each of the matched peaks.
    for(size_t i = 0; i < is_matched.size(); i++) {
      if(is_matched[i] >= 0) {
	int L = is_matched[i];
	double da_error = filtered[i].mz2 - ladder[L];
	double ppm_error = mass::massError2ppm(da_error, filtered[i].mz2);
	da_errors.push_back(da_error);
	ppm_errors.push_back(ppm_error);
      }
    }
  }

  // Perform database search on all MS/MS spectra.
  void fragmentDB::search(vector<ms2peak *> &peaks, float min_mz2, float max_mz2, vector<isotope_t> &iso, int nitro_cnt) {
    for(size_t i = 0; i < peaks.size(); i++) search(peaks[i], min_mz2, max_mz2, iso, nitro_cnt);
  }

  // Used to keep track of minimum index count map to protein index
  struct cntidx_t {
    int min_idx_cnt, idx; 
    bool operator < ( cntidx_t b) const { return min_idx_cnt < b.min_idx_cnt; }
  };

  // Creates a mapping from group number, which ranges from 0 to
  // group2idx.size() - 1, to an index in protein_meta ranging from 0
  // to protein_meta.size() - 1 (or 0 to Nproteins). Note that MS/MS
  // peaks have member functions that count and collect the number of
  // protein indexes corresponding to highest scoring MS/MS ids.
  // NOTE: All peaks must support is the function ->max_idxs which
  // returns all of the unique protein indexes for that spectrum.
  void build_protein_groups(size_t Nproteins, vector<ms2peak *> &peaks, 
			    vector< vector<int> > &group2idx, bool form_nonconclusive) {
    // Map from protein index to a list of spectra.
    vector< vector<ms2peak *> > idx2peaks(Nproteins);

    // Using maximum scoring IDs to assign peak to a protein index.
    for(size_t i = 0; i < peaks.size(); i++) {
      if(peaks[i] == NULL) continue; // Skip no MS/MS peak case.
      vector<int> idxs; peaks[i]->max_idxs(idxs); // Get maximum scoring protein indexes.
      for(size_t j = 0; j < idxs.size(); j++)  // Create map from protein index to MS/MS spectra.
	idx2peaks[idxs[j]].push_back(peaks[i]);
    }

    // Create a map from a protein index to a minimum ID count.
    vector<cntidx_t> min_idx_cnt2idx; // Minimum ID count -> protein index map.
    vector<int> idx2min_idx_cnt(Nproteins);  // protein index -> minimum ID count

    for(size_t idx = 0; idx < idx2peaks.size(); idx++) {
      vector<ms2peak *> &peaks = idx2peaks[idx]; // Get peaks assigned to this protein index.
      if(peaks.size() == 0) continue; // No MS/MS peaks with an ID.

      // For a given protein, determine the most differentiating MS/MS
      // spectrum. The most differentiating spectrum maps to the
      // fewest number of protein indexes.
      cntidx_t c;
      c.min_idx_cnt = peaks[0]->max_idx_cnt();
      for(size_t j = 1; j < peaks.size(); j++) {
	// Remember the minimum number of proteins these spectra allow
	// us to distinguish.
	if((int)peaks[j]->max_idx_cnt() < c.min_idx_cnt) c.min_idx_cnt = peaks[j]->max_idx_cnt();
      }

      // Save minimum protein index count for that protein.
      idx2min_idx_cnt[idx] = c.min_idx_cnt;
      
      // Save a mapping form the count to the protein index.
      c.idx = idx;
      min_idx_cnt2idx.push_back(c);
    }

    // Start with proteins with the smallest most differentiating spectra.
    sort(min_idx_cnt2idx.begin(), min_idx_cnt2idx.end());
    util::ugraph_t G(Nproteins); // One node for each protein.

    for(size_t i = 0; i < min_idx_cnt2idx.size(); i++) {
      // Smaller min_idx_cnt protein indexes are visited first.
      int min_idx_cnt = min_idx_cnt2idx[i].min_idx_cnt;
      int idx = min_idx_cnt2idx[i].idx;

      // Get MS/MS peaks for that protein index.
      vector<ms2peak *> &peaks = idx2peaks[idx];      
      for(size_t j = 0; j < peaks.size(); j++) {
	// Process if the IDs in the MS/MS spectrum matches the most differentiating amount.
	if(peaks[j]->max_idx_cnt() == min_idx_cnt) {
	  vector<int> idxs; peaks[j]->max_idxs(idxs);	
	  for(size_t k = 0; k < idxs.size(); k++) {
	    // Only connect them in the graph to less differentiating proteins.
	    if(idx2min_idx_cnt[idxs[k]] >= min_idx_cnt)
	      G.add_edge(idx, idxs[k]);
	  }
	}
      }
      // Advance to next minimum index count and protein index.
    }

    // Compute connected components in protein graph.
    vector<int> labels; int num_components = G.connected_components(labels);

    // Group protein indexes by component.
    vector< vector<int> > component2idx(num_components);
    for(size_t idx = 0; idx < labels.size(); idx++) component2idx[labels[idx]].push_back(idx);
    
    // Get indexes assigned to component.
    for(size_t c = 0; c < component2idx.size(); c++) {
      vector<int> &idx = component2idx[c];

      // Count the number of spectra from this protein group.
      size_t peak_cnt = 0;
      for(size_t i = 0; i < idx.size(); i++) peak_cnt += idx2peaks[idx[i]].size();
      if(peak_cnt == 0) continue;       // No group if no supporting spectra are found.

      bool form_group = false; 
      // Examine all protein indexes in connected component.
      for(size_t i = 0; i < idx.size(); i++) {
	// Get MS/MS peaks assigned to a protein index.
	vector<ms2peak *> &peaks = idx2peaks[idx[i]];

	// Determine if the protein group is nonconclusive.
	size_t shared_cnt = 0;
	for(size_t p = 0; p < peaks.size(); p++) {
	  vector<int> idxs; peaks[p]->max_idxs(idxs);

	  // Get the unique indexes (avoid corner case where repeats within a protein).
	  set<int> unique_idxs; 
	  for(size_t k = 0; k < idxs.size(); k++) unique_idxs.insert(idxs[k]); 

	  // Remove all the indexes represented by the current group.
	  for(size_t x = 0; x < idx.size(); x++) unique_idxs.erase(idx[x]);

	  // If there are other indexes left, then this peak is shared.
	  if(unique_idxs.size() > 0) shared_cnt++;
	}

	// If any one of the spectra are not shared, then we can go
	// ahead and form the group.
	if(shared_cnt < peaks.size()) { form_group = true; break; }
	// Otherwise, if all peaks are shared, this is a
	// non-conclusive group. 
      }
      if(form_group || form_nonconclusive)
	group2idx.push_back(idx); // Form group if there is a non-shared spectrum.
    }
  }

  // Assigns give MS/MS peaks to a protein group using a mappign from
  // group number to index in protein_meta.  The end result is a
  // mapping from an MS/MS peak to the protein groups to which the
  // MS/MS peak belongs.

  // Note peak2groups[i] might contain a -1 for an assignment to a non-conclusive group.
  // group2idx comes from build_protein_groups. 
  void assign_protein_groups(size_t Nproteins,
			     vector<ms2peak *> &peaks, 
			     vector< vector<int> > &group2idx, vector< vector<int> > &peak2groups) {

    // Construct a reverse mapping from index in protein_meta to a protein group number
    vector<int> idx2group(Nproteins, -1); 
    // NOTE: -1 is critical here it assures that -1 group ends up in peak2groups.
    // for non-conclusive groups.
    for(size_t g = 0; g < group2idx.size(); g++) 
      for(size_t k = 0; k < group2idx[g].size(); k++) 
	idx2group[group2idx[g][k]] = g;

    // Assign MS/MS peak a set of groups.
    peak2groups.resize(peaks.size());
    for(size_t i = 0; i < peaks.size(); i++) {
      if(peaks[i] == NULL) continue; // Skip NULL peaks.
      if(peaks[i]->num_ids() == 0) continue; // Skip peaks with no MS/MS ID information.
      vector<int> idxs; peaks[i]->max_idxs(idxs);

      // Go through maximum scoring IDs and get groups this MS/MS spectrum is assigned to.
      vector<int> groups;
      for(size_t j = 0; j < idxs.size(); j++) {
	// Cross reference the group number by protein index.
	int group_num = idx2group[idxs[j]];

	// Do not add repeats to peak.
	bool found = false;
	for(size_t k = 0; k < groups.size(); k++) { if(groups[k] == group_num) { found = true; break; }	}
	// Add group, if not added already.
	if(found == false) groups.push_back(group_num); 
      }
      if(groups.size() == 0) continue;       // No groups assigned.

      // Rare case, but happens. Additional shared spectrum, on top of most discriminating spectra.
      if(groups.size() == 1 && groups[0] == -1) continue; 

      // Otherwise, save assigned groups.
      peak2groups[i] = groups;
      // NOTE: If more than one group then this is a razor peptide.
      // This additional group might be a -1 group which is a
      // nonconclusive group.
    }
  }

  // TODO: Make more robust to errors in FASTA file.  NOTE: This load
  // function appends an N-terminal symbol "#" and a C-terminal symbol
  // "*"!!!!
  bool load_fasta_file(string descr, string filename, vector<string> &metadata, vector<string> &seqdata) {
    string seq, line, line_seq, meta;
    ifstream fasta(filename.c_str());
    if(!fasta.is_open()) return false; // Unable to open file.
    descr = "[" + descr + "] "; // Add brackets around description.
    seq = '#'; // pre-pend N-term symbol
    while(!fasta.eof()) {
      line.clear();  getline(fasta, line, '\n'); // Load one line at a time.
      if(line.length() == 0) continue;
      if(line[line.length() - 1] == '\r') line = line.substr(0, line.length() - 1);
      
      if(line[0] == '>') {
	if(meta != "") {
	  // Append C-term symbol and save meta data and sequence data.
	  if(seq[seq.length() - 1] != '*') seq += '*';
	  metadata.push_back(descr + meta);
	  seqdata.push_back(seq);
	  seq = '#'; // Start with N-term symbol.
	}
	meta = line;
	meta = meta.substr(1);
	if(meta.size() > 0 && meta[0] == ' ') meta = meta.substr(1);
      }
      else { 
	line_seq.clear();
	// Remove any spaces in a sequence line.
	for(size_t i = 0; i < line.length(); i++) if(!isspace(line[i])) line_seq += line[i];
	// Then, append to sequence.
	seq += line_seq; line.clear();
      }
    }
    if(meta != "") { // Get very last protein sequence.
      // Append C-term symbol.
      if(seq[seq.length() - 1] != '*') seq += '*';
      metadata.push_back(descr + meta); 
      seqdata.push_back(seq);
    }

    return true; // TODO: return false if error message.
  }

  // Loads fasta files and constructs concatenated reverse decoy database.
  bool fasta_data::load_fasta(vector<string> &descr, vector<string> &files) {
    for(size_t i = 0; i < files.size(); i++) {
      cout << "Loading fasta file: " << files[i] << endl;
      // TODO: Handle errors here better.
      load_fasta_file(descr[i], files[i], protein_meta, seq_data); 
    }

    // All proteins are non-decoy proteins.
    is_decoy_protein.reserve(protein_meta.size() * 2);
    is_decoy_protein.resize(protein_meta.size());
    for(size_t i = 0; i < is_decoy_protein.size(); i++) is_decoy_protein[i] = false;

    // Concatenate (blank) protein meta information.
    size_t protein_meta_orig_size = protein_meta.size();

    for(size_t i = 0; i < protein_meta_orig_size; i++) {
      protein_meta.push_back("");       //protein_meta.push_back("REV_" + protein_meta[i]);
      is_decoy_protein.push_back(true); // Mark flag indicating decoy protein.
      // Reverse sequence.
      string seq = seq_data[i];
      rev_seq(seq); // Reverse sequence (accounting for cut site).
      seq_data.push_back(seq); // and append to sequence and SNP data.
    }
    cout << "Done loading FASTA file." << endl;
    return true;
  }

  // Finds the smallest score cutoff that achives given FDR.
  double fdr_threshold(vector<ms2peak *> &peaks, double fdr_needed) {
    vector<double> decoy_scores, target_scores;     // Collect decoy scores and target scores.

    // Separate scores based on is_decoy flag in ID.
    for(size_t i = 0; i < peaks.size(); i++) {
      // Skip NULL peaks or peaks with no IDs.
      if(peaks[i] == NULL) continue;  if(peaks[i]->num_ids() == 0) continue; 
      // Get maximum decoy and target scores.
      double max_target = peaks[i]->max_target_score(), max_decoy = peaks[i]->max_decoy_score();
      // If target is higher, keep targe score.
      if(max_target >= max_decoy) target_scores.push_back(max_target); else decoy_scores.push_back(max_decoy); 
    }
    if(decoy_scores.size() == 0) return 0; // No need to threshold, accept everything.

    // Sort scores in ascending order.
    sort(decoy_scores.begin(), decoy_scores.end());  sort(target_scores.begin(), target_scores.end());

    // Estimate \pi_0 or PIT.
    double PIT = fdr::estimate_PIT(decoy_scores, target_scores);

    // Map each score to a PIT adjusted q-value.
    for(size_t i = 0; i < peaks.size(); i++) {
      // Skip peaks w/ no IDs.
      if(peaks[i] == NULL) continue; if(peaks[i]->num_ids() == 0) continue;
      // Assign q-value to highest target scores.
      vector<ms2id_t> &target = peaks[i]->target;
      for(size_t k = 0; k < target.size(); k++) 
	target[k].qvalue = fdr::PIT_qvalue(PIT, decoy_scores, target_scores, target[k].score);
    }
    // Return score threshold.
    cout << "estimating search score" << endl;
    return fdr::threshold(PIT, fdr_needed, decoy_scores, target_scores);
  }
}
