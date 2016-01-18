#ifndef __msconfig_hpp__
#define __msconfig_hpp__

#include <vector>
#include <string>
#include <limits>
#include <ctype.h>

// Some useful masses.
namespace mass {
  const double neutron = 1.0086649156;
  const double proton =  1.00727646677;
  const double electron = 5.4857990943e-4;

  const double H = 1.0078250320710; // proton + electron 
  const double O = 15.9949146195616;
  const double N = 14.00307400486;
                   
  const double S = 31.97207069;
  const double P = 30.97376151;

  const double water = 2 * H + O;
  const double ammonia = N + H * 3;

  const double C12 = 12.000000;
  const double C = C12;
  const double C13 = 13.003354837810;
                     
  const double N14 = N;
  const double N15 = 15.00010889827;

  const double O16 = O;
  const double O17 = 16.9991317012;

  const double O18 = 17.99916107;
  
  const double H1 = H;
  const double H2 = 2.01410177784;
                    
  // Isotope label mass shifts. 
  const double delta15N = N15 - N14; 
  const double delta13C = C13 - C12; 
  const double delta18O = O18 - O16;
  const double delta17O = O17 - O16;
  const double delta2H  = H2 - H1;

  // Converts PPMs to Daltons.
  inline double ppm2da(double mz, double ppm) { return mz * ppm * 1e-6; }

  // Computes mass error in PPMs.
  // ppm = (experimental mass - calculated mass) / (experimental mass) x 10^6.
  inline double massError2ppm(double deltaMass, double mass) { return deltaMass / mass * 1e6; }

  // Computes mass error expected for a give experimental mass. 
  // (experimental mass - calculated mass) = ppm * (experimental mass) x 10^6.
  inline double ppm2massError(double mass, double ppm) { return ppm * mass * 1e-6; }
}

// Fixed modification.
struct fxmod_t {
  fxmod_t() { deltaMass = 0; active = false; system = false; }
  fxmod_t(std::string d, char aa, double dm, bool s = false) {  
    deltaMass = dm;  AAs.push_back(aa);  description = d;  
    active = false; system = s;
  }
  fxmod_t(std::string d, char aa0, char aa1, double dm, bool s = false) { 
    description = d; deltaMass = dm; AAs.resize(2); AAs[0] = aa0; AAs[1]= aa1;
    active = false; system = s;
  }
  fxmod_t(std::string d, char aa0, char aa1, char aa2, double dm, bool s = false) { 
    description = d; deltaMass = dm; AAs.resize(3); AAs[0] = aa0; AAs[1] = aa1; AAs[2]= aa2;
    active = false; system = s;
  }
  fxmod_t(std::string d, char aa0, char aa1, char aa2, char aa3, double dm, bool s = false) { 
    description = d; deltaMass = dm; AAs.resize(4); AAs[0] = aa0; AAs[1] = aa1; AAs[2] = aa2; AAs[3] = aa3;
    active = false; system = s;
  }
  double deltaMass; // Shift in precursor mass.
  std::string AAs; // Amino acid codes to which this fixed modification applies.
  std::string description; // Descriptive name.
  bool active, system; // Flags designating system active/inactive or system provided.
};

// Fragment modification. Stores changes in fragment masses.
struct fragmod_t {
  fragmod_t() { idx = -1; dfrag_mass = 0; }
  fragmod_t(std::string abbr, double dm, std::string _AAs) {idx = -1; abbreviation = abbr; dfrag_mass = dm; AAs = _AAs;}
  fragmod_t(std::string abbr, char a, double dm) { idx = -1; abbreviation = abbr; AAs += a; dfrag_mass = dm; }
  fragmod_t(std::string abbr, char a0, char a1, double dm) { 
    idx = -1;
    abbreviation = abbr;
    AAs += a0; AAs += a1; 
    dfrag_mass = dm; 
  }
  fragmod_t(std::string abbr, char a0, char a1, char a2, double dm) { 
    idx = -1;
    abbreviation = abbr; 
    AAs += a0; AAs += a1; AAs += a2;
    dfrag_mass = dm; 
  }
  int idx; // identifier position varmod_t::fragmods
  std::string abbreviation; // abbreviation S-phos, K3me, etc. Kac
  std::string AAs; // Amino acid codes with this fragment modification
  double dfrag_mass; // shift in fragment mass S, S98, S80, T, T98, T80, Y so repeat character codes
  bool is_modifiedAA(char aa) {
    for(size_t i = 0; i < AAs.size(); i++) if(aa == AAs[i]) return true;
    return false;
  }
};

// Designates a variable modification type. Stores the changes in
// precursor mass.
struct varmod_t {
  varmod_t() { idx = -1; active = false; d_prec_mass = 0; system = false; max_count = 0; }
  varmod_t(std::string d, std::string a, double dm1, bool sys = false, int mc = 0) { 
    description = d;   // modification description
    abbrev = a;   // also abbreviation
    d_prec_mass = dm1; // change in precursor mass
    max_count = mc;    // maximum cardinality for this modification.
    active = false;    // flag designates if active for use or not
    system = sys;      // system modifictation, don't delete these!
  }
  int idx; // numeric identifier (position in Varmods)
  std::string description; // long description of mod
  std::string abbrev; // mod abbreviation
  double d_prec_mass; // shift in precursor mass
  std::vector<fragmod_t> fragmods; // amino acid code + shift in fragment mass
  int max_count; // threshold on maximum number of this modification type (0 = unlimited)  
  bool active, system; // flag designates active/inactive... system designates provided for user (do not delete)
  // Returns true if the given amino acid code is modified by this variable mod.
  bool is_modifiedAA(char aa) {
    // Code is found by searching all fragment modifications.
    for(size_t i = 0; i < fragmods.size(); i++) if(fragmods[i].is_modifiedAA(aa)) return true;
    return false;
  }
  void updateFragModIdx() { for(int idx = 0; idx < (int)fragmods.size(); idx++) fragmods[idx].idx = idx; }
};

// An isotope label that induces a given dalton shift.  Treat these as
// fixed modifications for heavy == true in ms2.cpp
struct isotope_t {
  // multiple constructors for convenience
  isotope_t() { idx = -1; system = false; shift_medium = shift_heavy = 0; }
  isotope_t(std::string d, double s, bool sys = false) { 
    idx = -1; description = d; shift_medium = 0; shift_heavy = s; system = sys;
  }
  isotope_t(std::string d, char aa, double s, bool sys = false) { 
    idx = -1; description = d; AAs.push_back(aa); shift_medium = 0; shift_heavy = s; system = sys; 
  }
  isotope_t(std::string d, char aa1, char aa2, double s, bool sys = false) { 
    idx = -1; description = d; AAs.push_back(aa1); AAs.push_back(aa2); 
    shift_medium = 0; shift_heavy = s; system = sys;
  }
  isotope_t(std::string d, char aa, double sH, double sM, bool sys = false) {
    idx = -1; description = d; AAs.push_back(aa);
    shift_heavy = sH, shift_medium = sM; system = sys;
  }
  int idx;
  std::string AAs; // amino acid codes that induces this shift
  std::string description; // descriptive name for modification
  double shift_heavy, shift_medium;
  bool system;
  bool find_aa(char aa) {   // Find the amino acid code in set.
    for(size_t i = 0; i < AAs.size(); i++) if(AAs[i] == aa) return true;
    return false;
  }
  
};


// Returns true if the given character is a valid amino acid code. 
inline bool isAAcode(char c) {
  switch(toupper(c)) {
  case 'A': case 'C': case 'D': case 'E': case 'F': case 'G': case 'H': case 'I':
  case 'K': case 'L': case 'M': case 'N': case 'P': case 'Q': case 'R': case 'S':
  case 'T': case 'V': case 'W': case 'Y': case '#': case '*': case '[': case ']':
    return true;
  default: return false;
  }
}


typedef int quantmode_t;
const int quantXICBASED = 0;
const int quantISOPAIR = 1; 
const int quantALIGN = 2;
const int quantISOTRIPLE = 3;
const int quantAcetylSTOICH = 4;

// Supported enzymes. 
enum enzyme_t { ez_Trypsin = 0, 
		ez_LysC = 1, 
		ez_GluC = 2,
		ez_ArgC = 3, 
		ez_TrypsinP = 4,
		ez_LysN = 5,
		ez_ChymoTrypsin = 6,
		ez_AspN = 7,
		ez_NumEnzymes = 8};

// Provides a single structure where all paramaters are stored.  NOTE:
// When adding a paramter, make sure to update GUI components and
// pview.xml configuration file output (in xml.cpp).
struct msconfig{
  // Processing/memory usage related flags.
  quantmode_t quantmode;
  bool run_algorithms; // true = run a peak grouping algorithms
  bool keep_raw; // false, free up memory used by raw data in GUI
  bool keep_filtered; // false, free up filtered peak memory
  bool keep_ms2; // false, free memory used by MS/MS peaks
  bool load_ms2; // Load MS/MS spectra or not.
  bool ms1_profile, ms2_profile; // profile mode flags 

  int threads;   // Number of threads used to load data.

  // XIC finding parameters (in XIC tab in GUI)
  // double dtime; // Filtering range query size, retention time.
  double dmz; // Filtering range query size, m/z, in PPMs.
  //int peakthresh; // Used for filtering (number of peaks that need to be present). 
  double hist_log10I_thresh; // Fraction of peaks discarded based on log10I histogram.

  // double peak_dtime; // Query parameters for constructing peak graphs.  
  double peak_dmz; // mz part of of range for peak graph, in PPMs.
  int xic_peaks_missed;
  int min_xic_peaks; // minimum number of peaks in XIC
  double max_xiclen; // MAXimum XIC length in seconds.

  // Label-free specific paramters 
  double align_dtime; // Used to find nearest XIC for alignment. 

  int min_iruns; // number of iruns in which an XIC must appear
  int min_reps; // number of replicate set in which XIC groups must appear
  int min_conds; // numer of conditions in which XIC groups must appear

  double xic_width; // "width" of XIC query box in PPMs.
  double group_dtime; // stretch in time dimension for grouping across iruns

  bool align_translation; // retention time adjustment, simple translation to reference
  bool align_nonlinear; // nonlinear alignment using running medians algorithm
  bool skip_normalize_label_free;  // skip label-free median of median normalization

  // Isotope specific parameters. 
  double isotolerence; // Error tolerence for isotope spacing.
  int max_label_count; // Allow multiple labels in a tryptic fragment (missed cleavages)
  bool skip_normalize_ratios; // Skip normalizing isotope ratios so their log2 median is zero. 
  std::vector<int> activeIsotopes; // among the availible isotopes, these are the active set
  bool isotope_15N; // use integrative approach to 15N quantification
  bool isotope_by_search; // Assign ungrouped isotopes by database search.

  // MS/MS protein database search parameters
  double ms2tol; // fragmentation spectrum tolerence in PPMs
  bool ms2tol_usedaltons; // use daltons instead
  double mstol; // precursor mass tolerence in PPMs
  unsigned int max_mod_enumerate; // How many modification positions to enumerate before giving up.
  int max_aa; // Maximum number of AAs in fragment
  int min_aa; // Minimum number of AAs in fragment
  int missed_cleaves; // Number of missed cleavages.
  double ms2_fdr; // false discovery rate threshold for MS/MS ids.
  bool fixed_15N; // true, then apply fixed 15N mods

  int enzyme_id; // ID of the enzyme to use, see ms2.hpp/ms2.cpp

  // MS/MS modification related flags (new tab in config dialog)
  int max_mod_aa; // maximum number of modified amino acids
  
  std::vector<fxmod_t> fxMods; // menu of possible fixed modifications.

  std::vector<varmod_t> varMods; // menu of possible variable modifications
  std::vector<isotope_t> isotopes; // menu of possibe isotopes (including system isotopes)

  // Some paramters (not exposed to user).
  int max_charge; float spurxictol; // see lcms.cpp - used for monoisotopic xic determination
  double logImin_disp, logImax_disp;
  bool filter_param_changed;
  // Stoichiometry params
  bool acetyl_stoich_plus5;

  // Set default configuration parameters.
  msconfig() {
    // Algorithm activation and memory usage.
    quantmode = quantXICBASED;
    run_algorithms = true; 
    keep_raw = true;
    keep_filtered = true;
    keep_ms2 = load_ms2 = true;
    ms1_profile = ms2_profile = false;
    threads = 1;

    // XIC comptuation parameters.
    //dtime = 32; 
    dmz = 14; // in PPMs
    //peakthresh = 5; 
    hist_log10I_thresh = 0.05; 
    
    //peak_dtime = 15; 
    xic_peaks_missed = 1;
    min_xic_peaks = 5;
    peak_dmz = 14;  // in PPMs
    max_xiclen = 400; 
    isotolerence = 3.5;  // in PPMs

    // Grouping/alignment rectanges
    min_conds = min_reps = min_iruns = 1;
    align_dtime = 400;
    xic_width = 36; // in PPMs
    group_dtime = 100;
    align_nonlinear = align_translation = true;
    skip_normalize_label_free = false;

    // Isotope pairing specific parameters.
    max_label_count = 1;
    skip_normalize_ratios = false;
    isotope_15N = false;
    isotope_by_search = false;

    // MS/MS protein DB search paramters
    ms2tol = 0.5; // in Daltons
    ms2tol_usedaltons = true;
    mstol = 20;  // PPM
    missed_cleaves = 0; // default is usually 2
    max_mod_aa = 0;
    ms2_fdr = 0.01; // 1% default FDR
    fixed_15N = false;

    // Use trypsin by default. TODO: Replace with an enum ? 
    enzyme_id = ez_Trypsin;
	
    resetVarMods(); resetMods();  resetIsotopes();
 
    // Fixed paramters (not exposed to user).
    max_charge = 8; spurxictol = 0.4;
    min_aa = 6;  // Hard coded: Do not allow user to adjust this. 
    max_aa = 40; // Possibly have different set for ETD mode!
    max_mod_enumerate = 100000; 
    logImin_disp = 3; logImax_disp = 6;  // logI display range, set within lcms.cpp and analysis.cpp

    acetyl_stoich_plus5 = false;
  }

  static const int SysVarModPhospho = 0;

  void resetVarMods() {
    using namespace mass;
    varMods.clear();
    double HPO3 = H + P + O * 3, H3PO4 = H * 3 + P + O * 4;
    varmod_t phospho("STY phosphorylation", "HPO3", HPO3, true);
    phospho.fragmods.push_back(fragmod_t("HPO3", 'S', 'T', 'Y', HPO3));
    fragmod_t ST_loss("HPO3-H3PO4", 'S', 'T', HPO3 - H3PO4);
    phospho.fragmods.push_back(ST_loss);
    varMods.push_back(phospho);
    
    varmod_t metox("methionine oxidation", "oxid", O, true);
    metox.fragmods.push_back(fragmod_t("oxid", 'M', O)); 
    varMods.push_back(metox); 

    varmod_t acylNterm("acylated protein N-term", "ac", C * 2 + 2 * H + O, true);
    acylNterm.fragmods.push_back(fragmod_t("ac", '#', C * 2 + 2 * H + O));
    varMods.push_back(acylNterm);

    varmod_t phospho_noloss("STY phosphorylation (no loss)", "HPO3", HPO3, true);
    phospho_noloss.fragmods.push_back(fragmod_t("HPO3", 'S', 'T', 'Y', HPO3));
    varMods.push_back(phospho_noloss);

    varmod_t pyroglu_gln("Pyro-glu from Gln (Q)", "Pyro" , mass::O, true);
    pyroglu_gln.fragmods.push_back(fragmod_t("Pyro-glu", 'Q', -H * 3 + -N));
    varMods.push_back(pyroglu_gln);

    varmod_t deamidation("deamidation of Asn (N)", "Deamid", -H + -N + O, true);
    deamidation.fragmods.push_back(fragmod_t("Deamid", 'N', -H + -N + O));
    varMods.push_back(deamidation);

    updateVarModIdx();
  }

  // Updates indicies to variable modifications, used for various lookup operations. 
  // NOTE: This *MUST* be called if a new variable modification is added.
  void updateVarModIdx() {
    for(int idx = 0; idx < (int)varMods.size(); idx++) {
      varMods[idx].idx = idx;
      varMods[idx].updateFragModIdx();
    }
  }
  void resetMods() { 
    fxMods.clear();
    // Fixed modification interface.
    using namespace mass;
    double H3C2NO = 3 * H + C * 2 + N + O;
    fxMods.push_back(fxmod_t("carboxyamidomethylation of cysteine", 'C', H3C2NO, true));
    // For phosphopepide enrichment the carboxylate groups (COO-) are methylated sometimes.
    fxMods.push_back(fxmod_t("carboxylate groups to methyl esters", ']', '*', 'E', 'D', 14.01565, true));
    fxMods.push_back(fxmod_t("MMTS, methylthio", 'C', 45.987721, true));
  }

  void resetIsotopes() {
    activeIsotopes.clear(); isotopes.clear();
    isotopes.push_back(isotope_t("Lys-13C6-15N2", 'K', 6 * mass::delta13C + 2 * mass::delta15N, true));
    isotopes.push_back(isotope_t("Arg-13C6-15N4", 'R', 6 * mass::delta13C + 4 * mass::delta15N, true));
    isotopes.push_back(isotope_t("Arg-13C6/Lys-13C6", 'K','R', 6 * mass::delta13C, true));
    isotopes.push_back(isotope_t("Arg-13C6", 'R', 6 * mass::delta13C, true));
    isotopes.push_back(isotope_t("Lys-13C6", 'K', 6 * mass::delta13C, true));
    isotopes.push_back(isotope_t("Arg-15N4", 'R', 4 * mass::delta15N, true));
    isotopes.push_back(isotope_t("D5, N-term", '[', '#', 5 * mass::delta2H, true));
    isotopes.push_back(isotope_t("Lys-2H4", 'K', 4 * mass::delta2H, true));
    // Triple labels. 
    isotopes.push_back(isotope_t("Lys-13C6-15N2 and Lys-2H4", 'K', 
				 6 * mass::delta13C + 2 * mass::delta15N, 4 * mass::delta2H, true));
    isotopes.push_back(isotope_t("Arg-13C6-15N4 and Arg-13C6",'R', 
				 6 * mass::delta13C + 4 * mass::delta15N, 6 * mass::delta13C,true));
    updateIsotopeIdx();
  }
  void updateIsotopeIdx() {
    for(int idx = 0; idx < (int)isotopes.size(); idx++) isotopes[idx].idx = idx;
  }
};


#endif // __msconfig_hpp__

