#ifndef __xml_hpp__
#define __xml_hpp__

#include <string>
#include <vector>
#include <map>

#include "lcms.hpp"

namespace xml {
  using namespace std;
  using namespace lcms;

  // Saves/loads pview.xml files.
  bool save_config(string filename, msconfig &conf);
  bool load_config(string filename, msconfig &conf);

  // TODO: Get rid of this PepXML crap at some point.
  bool load_pepxml(vector<string> &pepxml, vector<irun*> &iruns, fasta_data &fasta);

  // Parser for MS1 and MS2 data.
  bool load_mzXML_MS1(string filename, vector<float> &retentionTimes, vector<peak> &data, bool is_profile = false);
  bool load_mzXML_MS2(string filename, vector<ms2peak> &data2, bool is_profile = false);

  // Recalibration points with retention time, m/z, and m/z shift.
  struct recal_t { 
    double time, delta_mz, mz; 
    recal_t() { time = delta_mz = mz = 0; }
    recal_t(double time, double delta_mz, double mz) : time(time), delta_mz(delta_mz), mz(mz) {}
  };
  void recalibrate(recalmode_t recalmode, double winSz, vector<recal_t> &recal, string src, string dst);
}

#endif  
