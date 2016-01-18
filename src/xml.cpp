// Code for parsing XML formats relevant to MS data. Uses the blazingly fast expat XML parser.

// UNDEFINE IF COMPILING EXPAT STATICALLY!!!!
// #define XML_STATIC
#include <expat.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <limits>

#include <ctype.h>
#include <string.h>

// MSVC++ is still moving to C99. Fuck you Microsoft.
#ifdef _MSC_VER
#include "win_stdint.h"
#else
#include <stdint.h>
#endif

#include "xml.hpp"
#include "b64.h"
#include "centroider.hpp"

namespace xml {
  using namespace lcms;
  using namespace std;

  // Static callback handlers, declaration.
  static void startElementCallback(void *data, const XML_Char *el, const XML_Char **attr);
  static void endElementCallback(void *data, const XML_Char *el);
  static void charactersCallback(void *data, const XML_Char *s, int len);

  // TODO: Get this working so that an exceptional handler can handle errors.
  struct XMLParseException {
    int line; string error;
    XMLParseException(int line_, const char *errorStr) { line = line_; error = errorStr;}
  };
  struct XMLHandlerException {
    string error; 
    XMLHandlerException(string errorStr) { error = errorStr;}
  };
  // Generic base-class hanlder for XML files.  Wraps expat functionality. 
  class XmlHandler {
    bool okflag;
    XML_Parser parser;
  public:
    XmlHandler() {
      parser = XML_ParserCreate(NULL);
      okflag = true;
      XML_SetUserData(parser, this);
      XML_SetElementHandler(parser, startElementCallback, endElementCallback);
      XML_SetCharacterDataHandler(parser, charactersCallback);
    }
    virtual ~XmlHandler() { XML_ParserFree(parser); }
    virtual void startElement(const XML_Char *el, const XML_Char **attr) = 0;
    virtual void endElement(const XML_Char *el) = 0;
    virtual void characters(const XML_Char *s, int len) = 0;
    bool ok() { return okflag; }
    bool isElement(const char *n1, const XML_Char *n2) { return (strcmp(n1, n2) == 0); }
    int getNumVal(const XML_Char **attr) { // Counts the number of attributes.
      int n_attr = 0;
      int i = 0;
      while(attr[i]) { n_attr++; i += 2; }
      return n_attr;
    }
    void getAllPairs(vector< pair<string, string> > &nameval, const XML_Char **attr) {
      int i = 0;
      while(attr[i]) {
	nameval.push_back( make_pair( string(attr[i]), string(attr[i+1]) ) );
	i += 2;
      }
    }
    string getValue(const char* name, const XML_Char **attr) {
      for (int i = 0; attr[i]; i += 2) { if (isElement(name, attr[i])) return string(attr[i + 1]); }
      return string();
    }
    // Functions for accessing and type conversion of attributes with
    // default vals in case the attribute is not found.
    double getValue(const char *name, const XML_Char **attr, double val) {
      bool ok;
      double valStr = util::fromString<double>(getValue(name, attr), &ok);
      if(ok) return valStr; else return val;
    }
    int getValue(const char *name, const XML_Char **attr, int val) {
      bool ok; int valStr = util::fromString<int>(getValue(name, attr), &ok);
      if(ok) return valStr; else return val;
    }
    long getValue(const char *name, const XML_Char **attr, long val) {
      bool ok; long valStr = util::fromString<long>(getValue(name, attr), &ok);
      if(ok) return valStr; else return val;
    }
    bool getValue(const char *name, const XML_Char **attr, bool val) {
      string valStr = getValue(name, attr);
      bool ok; int valInt = util::fromString<int>(valStr, &ok);
      if(ok) return valInt != 0; else return val;
    }
    void parse(string filename) {
      char Buff[16384];
      FILE* fp = fopen(filename.c_str(), "r");
      if(fp == NULL) return; // TODO Handle these errors with dialog boxes or exceptions.
      int done = 0, len = 0;
      while(done == 0) {
	len = fread(Buff, 1, 16384, fp);
	if (ferror(fp)) {
	  fprintf(stderr, "Read error\n"); 
	  okflag = false;
	  break; // TODO: Handle this error with an exception.
	}
	done = feof(fp);
	XML_Status status = XML_Parse(parser, Buff, len, done);
	if (done == 0 && status == XML_STATUS_ERROR ) {
	  // TODO handle with dialog boxes.
	  cerr << "Parse error at line " << XML_GetCurrentLineNumber(parser) << ": " << XML_ErrorString(XML_GetErrorCode(parser));
	  okflag = false;
	  //TODO: throw XMLParseException(XML_GetCurrentLineNumber(parser), XML_ErrorString(XML_GetErrorCode(parser)));
	  break; // TODO: Handle this error with an exception.
	}
      }
      fclose(fp);
    }
  };
  static void startElementCallback(void *data, const XML_Char *el, const XML_Char **attr) { ((XmlHandler*) data)->startElement(el, attr);}
  static void endElementCallback(void *data, const XML_Char *el) { ((XmlHandler*) data)->endElement(el);}
  static void charactersCallback(void *data, const XML_Char *s, int len) { ((XmlHandler*) data)->characters(s, len); }



  // State space of XML SAX parser: mzIGNORE = waiting for a scan, SCAN =
  // processing scan, and PEAKS = looking for peak data, SCAN2 = looking
  // for MS2 data, PEAKS2 = looking for MS2 peaks, PRECURSOR = looking
  // for precursor data, PRECURSOR2 = found precursor data
  enum mzxml_sax_state { mzIGNORE, SCAN, PEAKS, SCAN2, PEAKS2, PRECURSOR, PRECURSOR2 }; 

  inline uint32_t SwapBytes(uint32_t n) { 
    return ((n&0xff)<<24) | ((n&0xff00)<<8) | ((n&0xff0000)>>8) | ((n&0xff000000)>>24); 
  }

  inline uint64_t SwapBytes(uint64_t n) {
    return ((n&0x00000000000000ffll)<<56) | 
           ((n&0x000000000000ff00ll)<<40) | 
           ((n&0x0000000000ff0000ll)<<24) | 
           ((n&0x00000000ff000000ll)<<8)  |
           ((n&0x000000ff00000000ll)>>8)  | 
           ((n&0x0000ff0000000000ll)>>24) |
           ((n&0x00ff000000000000ll)>>40) | 
           ((n&0xff00000000000000ll)>>56);
  }
  
  template<typename uint, typename real> 
  void decodemzint(vector<float> &mz, vector<float> &intensity, string &buffer, vector<unsigned char> &decoded) {
    if(buffer.size() == 0) return;
    decoded.resize(B64_NAMESPACE::b64_decode(buffer.c_str(), buffer.length(), NULL, 0));
    if(decoded.size() == 0) return;
    unsigned int decodedLength = 
      B64_NAMESPACE::b64_decode(buffer.c_str(), buffer.length(), &decoded[0], decoded.size());

    // Convert decoded to data to host byte order from network byte order.
    uint *fixed = (uint*)&decoded[0];
    for(unsigned int i = 0; i < decodedLength / sizeof(uint); i++) fixed[i] = SwapBytes(fixed[i]);
    real *mz_int = (real*)fixed; 
    mz.clear(); intensity.clear();
    for(unsigned int i = 0; i < decodedLength / sizeof(real); i+=2){
      mz.push_back(mz_int[i]);
      intensity.push_back(mz_int[i+1]);
    }
  }


  // Rough picture of the simple XML parser's state diagram.
  // For MS1 data:
  // mzIGNORE -> SCAN -> PEAKS -> **mzIGNORE**
  // For MS2 data:
  // mzIGNORE SCAN2 -> PRECURSOR -> PRECURSOR2 -> PEAKS2 -> **mzIGNORE**
  class MzXMLHandlerMS1 : public XmlHandler {  
    mzxml_sax_state state;
    vector<float> &retentionTimes;
    vector<peak> &data; bool is_profile;
    vector<unsigned char> decoded; 
    vector<float> mz, intensity;
    vector<centroider::peak> profile, centroid;
    string buffer;
    double retentionTime;
    int precision; // 32 or 64 
    long scanNum;
  public:
    MzXMLHandlerMS1(vector<float> &retentionTimes, vector<peak>  &data, bool is_profile) : 
      retentionTimes(retentionTimes), data(data), is_profile(is_profile) { state = mzIGNORE; }
    void startElement ( const XML_Char *el, const XML_Char **attr ) {
      // Start processing a scan. 
      if(isElement("scan", el)) {
	string msLevel = getValue("msLevel", attr);

	// First, get the retention time.
	string retentionStr = getValue("retentionTime", attr);
	for(int i = 0; i < (int)retentionStr.length(); i++) 
	  if(isalpha(retentionStr[i])) retentionStr[i] = ' ';

	retentionTime = util::fromString<double>(retentionStr);

	// Save scan number.
	scanNum = getValue("num", attr, -1L); // TODO: Handle this error somehow!!!

	// Set the state so that we now look for peaks.
	if(msLevel == "1" || msLevel == "0") state = SCAN; else state = mzIGNORE;
      }
      else if(state == SCAN && isElement("peaks", el) ) {
	state = PEAKS;
	precision = getValue("precision", attr, (int)32);
      }
    }
    // The peak and precursor data is gathered here.
    void characters(const XML_Char *s, int len){ if(state == PEAKS) buffer.append(s, len); }
    void endElement(const XML_Char *el) { 
      if(state == PEAKS && isElement("peaks", el)) {
	if(precision == 64) decodemzint<uint64_t, double>(mz, intensity, buffer, decoded);
	else decodemzint<uint32_t, float>(mz, intensity, buffer, decoded);
	if(is_profile) {
	  profile.clear(); centroid.clear();
	  for(size_t i = 0; i < mz.size(); i++) profile.push_back(centroider::peak(mz[i],intensity[i]));
	  centroider::gausfit(profile, centroid);
	  for(size_t i = 0; i < centroid.size(); i++) {
	    peak p; p.retentionTime = retentionTime; 
	    p.mz = centroid[i].mz; p.intensity = centroid[i].intensity;
	    if(p.intensity > 0 && p.retentionTime > 0 && p.mz > 0) data.push_back(p); // Store peak.
	    else {
	      cout << p.mz << " " << p.intensity << endl;
	    }
	  }
	}
	else {
	  for(size_t i = 0; i < mz.size(); i++) {
	    peak p; p.retentionTime = retentionTime; 
	    p.mz = mz[i]; p.intensity = intensity[i];
	    if(p.intensity > 0 && p.retentionTime > 0 && p.mz > 0) data.push_back(p); // Store peak.
	  }
	}
	retentionTimes.push_back(retentionTime);
	state = mzIGNORE; // Reset the state to look for the next MS scan. 
      }
      buffer.clear();
    }
  };

  class MzXMLHandlerMS2 : public XmlHandler {  
    mzxml_sax_state state;
    vector<unsigned char> decoded;
    vector<float> mz, intensity;
    vector<ms2peak> &data2; bool is_profile;
    vector<centroider::peak> profile, centroid;
    string buffer;
    double retentionTime, precursorMz;
    long scanNum;
    int charge, precision; 
    ms2::activationMethod_t activation;
  public:
    MzXMLHandlerMS2(vector<ms2peak>  &data2, bool is_profile) : data2(data2), is_profile(is_profile) { state = mzIGNORE; }
    void startElement ( const XML_Char *el, const XML_Char **attr ) {
      // Start processing a scan. 
      if(isElement("scan", el)) {
	string msLevel = getValue("msLevel", attr);

	// First, get the retention time.
	string retentionStr = getValue("retentionTime", attr);
	for(int i = 0; i < (int)retentionStr.length(); i++) 
	  if(isalpha(retentionStr[i]) ) retentionStr[i] = ' ';
	retentionTime = util::fromString<double>(retentionStr);

	// Save scan number.
	scanNum = getValue("num", attr, -1L); // TODO: Handle this error somehow!!!

	// Set the state so that we now look for peaks.
	if(msLevel == "2") state = SCAN2; else state = mzIGNORE;
      }
      else if(state == SCAN2 && isElement("precursorMz", el)) {
	// Read activation method.
	activation = ms2::CID_activation; // Assume CID by default.

	string methodStr = getValue("activationMethod", attr);
	if(methodStr == "CID") activation = ms2::CID_activation;
	else if (methodStr == "ETD") activation = ms2::ETD_activation;
	else if (methodStr == "HCD") activation = ms2::HCD_activation;
	  
	charge = getValue("precursorCharge", attr, (int)0);
	state = PRECURSOR;
      }
      else if(state == PRECURSOR2 && isElement("peaks", el)) {
	state = PEAKS2;  // We have reached MS2 peak data.
	precision = getValue("precision", attr, (int)32);
      }
    }
    // The peak and precursor data is gathered here.
    void characters(const XML_Char *s, int len) { if(state == PRECURSOR || state == PEAKS2) buffer.append(s, len); }
    void endElement(const XML_Char *el) { 
      if(state == PRECURSOR && isElement("precursorMz", el)) {
	bool okprec; precursorMz = util::fromString<double>(buffer, &okprec);
	// TODO: Check ok flag!!!
	state = PRECURSOR2;
      }
      else if(state == PEAKS2) {
	if(precursorMz > 0) {
	  if(precision == 64) decodemzint<uint64_t, double>(mz, intensity, buffer, decoded);
	  else decodemzint<uint32_t, float>(mz, intensity, buffer, decoded);
	  ms2peak p; 
	  p.scanNum = scanNum; p.retentionTime = retentionTime; p.mz = precursorMz;
	  p.activation = activation;
	  p.charge = charge;
	  if(is_profile) {
	    profile.clear(); centroid.clear();
	    for(size_t i = 0; i < mz.size(); i++) profile.push_back(centroider::peak(mz[i],intensity[i]));
	    centroider::gausfit(profile, centroid);
	    for(size_t i = 0; i < centroid.size(); i++) {
	      peak2 p2; p2.mz2 = centroid[i].mz; p2.intensity2 = centroid[i].intensity;
	      if(p2.intensity2 > 0) p.ms2.push_back(p2); // Store peak.
	    }
	  }
	  else {
	    for(size_t i = 0; i < mz.size(); i++) {
	      peak2 p2;
	      p2.mz2 = mz[i]; p2.intensity2 = intensity[i];
	      if(p2.intensity2 > 0)  p.ms2.push_back(p2);
	    }
	  }
	  data2.push_back(p); 
	}
	state = mzIGNORE; // Reset the state to look for the next MS scan. 
      }
      buffer.clear();
    }
  };  

  // TODO: Get copy handler to preserve meta data in mzXML file.
  struct CopyHandler : public XmlHandler {  
    ofstream copy;
    string buffer;
    vector< pair<string, string> > nameVal;
    vector<string> estack;
    CopyHandler(const string &filename) : copy(filename.c_str()) {  
      copy << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
    }
    void startElement ( const XML_Char *el, const XML_Char **attr ) {
      string element = el;

      nameVal.clear(); getAllPairs(nameVal, attr);

      copy << string(2 * estack.size(), ' ');
      copy << "<" << element;
      if(nameVal.size() > 0) {
	copy << " " << nameVal[0].first << "=\"" << nameVal[0].second << "\"";
	if(nameVal.size() > 1) copy << endl;
	for(size_t i = 1; i < nameVal.size(); i++) {
	  copy << string(2 * estack.size() + (element.length() + 2), ' ');
	  copy << nameVal[i].first << "=\"" << nameVal[i].second << "\"";
	  if(i != nameVal.size() - 1) copy << endl;
	}
      }
      copy << ">" << endl;
      estack.push_back(element);
    }
    ~CopyHandler() {}
    void characters(const XML_Char *s, int len) { buffer.append(s, len);  }
    void endElement(const XML_Char *el) { 
      string element = el;
      if(buffer.size() > 0) copy << buffer;
      copy << "</" << element << ">" << endl;
      estack.pop_back();
      buffer.clear();
    }
  };

  // Used to sort recalibration points on time.
  struct rtime_cmp { bool operator () (const recal_t &a, const recal_t &b) const { return a.time < b.time; } };
  // Use to sort recalibration points on m/z.
  struct rmz_cmp { bool operator () (const recal_t &a, const recal_t &b) const { return a.mz < b.mz; } };

  // XML handler for recalibration. 
  // TODO: Need to merge with the copy handler so meta data is preserved.
  struct RecalibrateHandler : public XmlHandler {  
    recalmode_t recalmode;
    mzxml_sax_state state;
    string buffer;
    vector<unsigned char> decoded; 
    vector<float> mz, intensity, mz_int;
    double retentionTime, winSz;
    int precision; // 32 or 64 
    ofstream mzXML;
    double calibration;
    vector<recal_t> &recal;
    vector<double> deltamz; 
    RecalibrateHandler(recalmode_t recalmode, double wSz, string dst, vector<recal_t> &recal) : 
      recalmode(recalmode), mzXML(dst.c_str()), recal(recal) { 
      winSz = wSz;
      // Sort on retention time or on m/z based on recalibration mode.
      switch(recalmode) {
      case recalmodeTIME: sort(recal.begin(), recal.end(), rtime_cmp()); break;
      case recalmodeMZ:   sort(recal.begin(), recal.end(), rmz_cmp());   break;
      }
      state = mzIGNORE;  
      mzXML << "<mzXML>" << endl; 
      calibration = 0;
    }
    ~RecalibrateHandler() { mzXML << "</mzXML>" << endl; }
    void startElement ( const XML_Char *el, const XML_Char **attr ) {
      // Start processing a scan. 
      if(isElement("scan", el)) {
	mzXML << "\t<scan num=\"" <<  getValue("num", attr) << "\" "; 
	mzXML << "retentionTime=\"" <<  getValue("retentionTime", attr) << "\" "; 
	mzXML << "msLevel=\"" << getValue("msLevel", attr) << "\" ";
	mzXML << ">" << endl;

	// First, get the retention time.
	string retentionStr = getValue("retentionTime", attr);
	for(int i = 0; i < (int)retentionStr.length(); i++) 
	  if(isalpha(retentionStr[i])) retentionStr[i] = ' ';
	retentionTime = util::fromString<double>(retentionStr);

	string msLevel = getValue( "msLevel", attr);
	// Set the state so that we now look for peaks.
	if(msLevel == "1" || msLevel == "0") state = SCAN; 
	else if(msLevel == "2") state = SCAN2; 
	else state = mzIGNORE;
      }
      else if(state == SCAN && isElement("peaks", el) ) {
	state = PEAKS;
	precision = getValue("precision", attr, (int)32);
      }
      else if(state == SCAN2 && isElement("precursorMz", el)) {
	mzXML << "\t<precursorMz activationMethod=\"" << getValue("activationMethod", attr) << "\" ";
	mzXML << "precursorCharge=\"" << getValue("precursorCharge", attr) << "\" >";
	state = PRECURSOR;
      }
      else if(state == PRECURSOR2 && isElement("peaks", el)) {
	state = PEAKS2;  // We have reached MS2 peak data.
	precision = getValue("precision", attr, (int)32);
      }
    }
    // The peak and precursor data is gathered here.
    void characters(const XML_Char *s, int len) { 
      if(state == PRECURSOR || state == PEAKS2 || state == PEAKS) buffer.append(s, len); 
    }
    double recalibrateMz(double mz) {
      if(recalmode != recalmodeMZ) return 0; // No shift if in time recalibration mode.
      recal_t query; query.mz = mz - winSz;
      double mz_end = mz + winSz;
      vector<recal_t>::iterator it = lower_bound(recal.begin(), recal.end(), query, rmz_cmp());      
      deltamz.clear();
      while(it != recal.end()) {
	deltamz.push_back(it->delta_mz);
	if(it->mz > mz_end) break;
	it++;
      }
      if(deltamz.size() == 0) return 0; // No recalibration points, found.
      // Use the trimmed mean to update calibration.
      sort(deltamz.begin(), deltamz.end());
      size_t trim = 0;
      if(deltamz.size() > 10) trim = deltamz.size() * 0.1;
      double sum = 0, N = 0;
      for(size_t i = trim; i < deltamz.size() - trim; i++) {
	sum += deltamz[i];
	N += 1;
      }
      return sum / N;
    }
    // Time recalibration.
    void recalibrateTime() {
      // Only update recalibration if recalibration mode is correct.
      if(recalmode != recalmodeTIME) return;
      // Use current retention time to find recalibration points that
      // are in [retentionTime-winSz, retentionTime+winSz].
      recal_t query; query.time = retentionTime - winSz;
      vector<recal_t>::iterator it = lower_bound(recal.begin(), recal.end(), query, rtime_cmp());
      deltamz.clear();
      while(it != recal.end()) {
	deltamz.push_back(it->delta_mz);
	if(it->time > retentionTime + winSz) break;
	it++;
      }
      // If there are more than 10 calibration points,
      if(deltamz.size() >= 10) {
	// Use the trimmed mean to update calibration.
	sort(deltamz.begin(), deltamz.end());
	size_t trim = deltamz.size() * 0.1;
	double sum = 0, N = 0;
	for(size_t i = trim; i < deltamz.size() - trim; i++) {
	  sum += deltamz[i]; N += 1;
	}
	// Update the calibration.
	calibration = sum / N;
      }
    }
    void endElement(const XML_Char *el) { 
      if(isElement("scan", el)) { mzXML << "\t</scan>" << endl; }
      if(state == PRECURSOR && isElement("precursorMz", el)) {
	bool okprec; double precursorMz = util::fromString<float>(buffer, &okprec);
	if(!okprec) mzXML << buffer << "</precursorMz>" << endl;
	else {
	  recalibrateTime();
	  precursorMz += recalibrateMz(precursorMz);
	  mzXML << fixed << setprecision(15) << precursorMz << "</precursorMz>" << endl;
	}
	state = PRECURSOR2;
      }
      if(state == PEAKS2) { // Leave MS/MS spectra alone.
	//mzXML << "\t<peaks precision=\"" << precision << "\">" << buffer << "</peaks>" << endl;
	mzXML << "\t<peaks>" << buffer << "</peaks>" << endl;
	state = mzIGNORE;
      }
      if(state == PEAKS && isElement("peaks", el)) {
	if(precision == 64) decodemzint<uint64_t, double>(mz, intensity, buffer, decoded);
	else decodemzint<uint32_t, float>(mz, intensity, buffer, decoded);
	recalibrateTime();
	for(size_t i = 0; i < mz.size(); i++) mz[i] += calibration;
	mz_int.resize(mz.size() * 2);
	for(size_t i = 0, j = 0; i < mz.size(); i++) {
	  mz_int[j] = mz[i] + recalibrateMz(mz[i]);
	  mz_int[j+1] = intensity[i];
	  j += 2;
	}
	size_t destLen = B64_NAMESPACE::b64_encode(NULL, sizeof(float) * mz_int.size(), NULL, 0);
	buffer.resize(destLen);
	uint32_t *to_htonl = (uint32_t *)&mz_int[0];
	for(unsigned int i = 0; i < mz_int.size() * sizeof(float) / sizeof(uint32_t); i++) 
	  to_htonl[i] = SwapBytes(to_htonl[i]); 

	B64_NAMESPACE::b64_encode(&mz_int[0], mz_int.size() * sizeof(float), &buffer[0], destLen);
	mzXML << "\t<peaks>" << buffer << "</peaks>" << endl;

	state = mzIGNORE; // Reset the state to look for the next MS scan. 
      }
      buffer.clear();
    }
  };


  bool load_mzXML_MS1(string filename, vector<float> &retentionTimes, vector<peak> &data, bool is_profile) {
    MzXMLHandlerMS1 *handler = new MzXMLHandlerMS1(retentionTimes, data, is_profile);
    cout << "Loading " << filename << endl;
    handler->parse(filename);
    cout << filename << ": DONE parsing XML data." << endl;
    cout << filename << ": total # of LC-MS peaks processed = " << data.size() << endl;    
    bool ok = handler->ok();
    delete handler;
    return ok;
  }

  bool load_mzXML_MS2(string filename, vector<ms2peak> &data2, bool is_profile) {
    MzXMLHandlerMS2 *handler = new MzXMLHandlerMS2(data2, is_profile);
    cout << "Loading " << filename << endl;
    handler->parse(filename);
    cout << filename << ": DONE parsing XML data." << endl;
    cout << filename << ": total # of MS/MS peaks loaded = " << data2.size() << endl;
    bool ok = handler->ok();
    delete handler;
    return ok;
  }

  void recalibrate(recalmode_t recalmode, double winSz, vector<recal_t> &recal, string src, string dst) {
    RecalibrateHandler *handler = new RecalibrateHandler(recalmode, winSz, dst, recal);
    //CopyHandler *handler = new CopyHandler(dst);
    handler->parse(src);
    delete handler;
  }

  class ConfigXMLHandler : public XmlHandler {  
    msconfig &conf;
  public:
    ConfigXMLHandler(msconfig &conf_) : conf(conf_) {    

    }
    void startElement(const XML_Char *el, const XML_Char **attr) {
      if(isElement("data", el)) {
	conf.quantmode = (quantmode_t)getValue("quantmode", attr, 0);
	conf.threads = getValue("threads", attr, conf.threads);
	conf.keep_raw = getValue("keep_raw", attr, conf.keep_raw);
	conf.keep_filtered = getValue("keep_filtered", attr, conf.keep_filtered);
	conf.load_ms2 = getValue("load_ms2", attr, conf.load_ms2);
	conf.keep_ms2 = getValue("keep_ms2", attr, conf.keep_ms2);
	conf.run_algorithms = getValue("run_algorithms", attr, conf.run_algorithms);
      }  
      else if(isElement("xics", el)) {
	conf.dmz = getValue("dmz", attr, conf.dmz);
	//conf.dtime = getValue("dtime", attr, conf.dtime);
	//conf.peakthresh = getValue("peakthresh", attr, conf.peakthresh);
	conf.hist_log10I_thresh = getValue("hist_log10I_thresh", attr, conf.hist_log10I_thresh);
	//conf.peak_dtime = getValue("peak_dtime", attr, conf.peak_dtime);
	conf.peak_dmz = getValue("peak_dmz", attr, conf.peak_dmz);
	conf.max_xiclen = getValue("max_xiclen", attr, conf.max_xiclen);
	conf.isotolerence = getValue("isotolerence", attr, conf.isotolerence);
	conf.xic_peaks_missed = getValue("xic_peaks_missed", attr, conf.xic_peaks_missed);
	conf.min_xic_peaks = getValue("min_xic_peaks", attr, conf.min_xic_peaks);
      }
      else if(isElement("free", el)) {
	conf.align_dtime = getValue("align_dtime", attr, conf.align_dtime);
	conf.min_conds = getValue("min_conds", attr, conf.min_conds);
	conf.min_reps = getValue("min_reps", attr, conf.min_reps);
	conf.min_iruns = getValue("min_iruns", attr, conf.min_iruns);
	conf.xic_width = getValue("xic_width", attr, conf.xic_width);
	conf.group_dtime = getValue("group_dtime", attr, conf.group_dtime);
	conf.align_translation = getValue("align_translation", attr, conf.align_translation);
	conf.align_nonlinear = getValue("align_nonlinear", attr, conf.align_nonlinear);
	conf.skip_normalize_label_free = getValue("skip_normalize_label_free", attr, conf.skip_normalize_label_free);
      }
      else if(isElement("isotope", el)) {
	conf.max_label_count = getValue("max_label_count", attr, conf.max_label_count);
	conf.skip_normalize_ratios = getValue("skip_normalize_ratios", attr, conf.skip_normalize_ratios);
	conf.isotope_15N = getValue("isotope_15N", attr, conf.isotope_15N);
	conf.isotope_by_search = getValue("isotope_by_search", attr, conf.isotope_by_search);
      }
      else if(isElement("ms2", el)) {
	conf.ms2tol = getValue("ms2tol", attr, conf.ms2tol);
	conf.ms2tol_usedaltons = getValue("ms2tol_usedaltons", attr, conf.ms2tol_usedaltons);
	conf.mstol = getValue("mstol", attr, conf.mstol);
	conf.missed_cleaves = getValue("missed_cleaves", attr, conf.missed_cleaves);
	conf.ms2_fdr = getValue("ms2_fdr", attr, conf.ms2_fdr);
	conf.enzyme_id = getValue("enzyme_id", attr, conf.enzyme_id);
	conf.max_mod_aa = getValue("max_mod_aa", attr, conf.max_mod_aa);
	conf.fixed_15N = getValue("fixed_15N", attr, conf.fixed_15N);
      }
      else if(isElement("sysfixedmod", el)) {
	int idx = getValue("activeindex", attr, -1); // -1 is OK here!
	if(0 <= idx && idx < (int)conf.fxMods.size()) conf.fxMods[idx].active = true;
      }
      else if(isElement("userfixedmod", el)) { // Found a fixed modificaiton.
	fxmod_t mod;
	mod.description = getValue("description", attr);
	mod.deltaMass = getValue("delta_mass", attr, (double)0.0);
	mod.active = getValue("active", attr, false); // OK
	mod.system = false;
	mod.AAs = getValue("aa", attr);
	conf.fxMods.push_back(mod);
      }
      else if(isElement("userisotope", el)) { // Read user defined isotopes first. 
	isotope_t iso;
	iso.description = getValue("description", attr); 
	iso.shift_heavy = getValue("shift_heavy", attr, (double)0.0);
	iso.shift_medium = getValue("shift_medium", attr, (double)0.0);
	iso.AAs = getValue("aa", attr);
	conf.isotopes.push_back(iso);
	conf.updateIsotopeIdx();
      }
      else if(isElement("active_isotope", el)) {
	int idx = getValue("index", attr, -1);
	if(0 <= idx && idx < (int)conf.isotopes.size()) conf.activeIsotopes.push_back(idx);
      }
      else if(isElement("sysvarmod", el)) {
	int idx = getValue("activeindex", attr, -1); // OK
	if(0 <= idx && idx < (int)conf.varMods.size()) conf.varMods[idx].active = true;
      }
      else if(isElement("uservarmod", el)) {
	varmod_t vr;
	vr.description = getValue("description", attr);
	vr.abbrev = getValue("abbrev", attr);
	vr.d_prec_mass = getValue("d_prec_mass", attr, (double)0.0);
	vr.max_count = getValue("max_count", attr, 0); // OK
	vr.active = getValue("active", attr, false); // OK
	string fragmodStr = getValue("fragmods", attr);

	vector<string> splitStr;
	size_t f = 0; 
	while(f < fragmodStr.size()) {
	  string cur;
	  while( f < fragmodStr.size() && fragmodStr[f] != '|') { cur += fragmodStr[f]; f++; }
	  splitStr.push_back(cur);
	  f++;
	}
	if(splitStr.size() % 3 == 0) { // should be a multiple of 3
	  for(size_t i = 0; i < splitStr.size(); i += 3) {
	    cout << splitStr[i] << endl;
	    cout << splitStr[i+1] << endl;
	    cout << splitStr[i+2] << endl;
	    cout << "*NEXT*" << endl;
	    double dm = util::fromString<double>(splitStr[i]);
	    string &AAs = splitStr[i+1];
	    string &abbrev = splitStr[i+2];
	    vr.fragmods.push_back(fragmod_t(abbrev, dm, AAs));	    
	  }
	}
        vr.system = false; // not a system provided variable modification
	conf.varMods.push_back(vr);
	conf.updateVarModIdx();
      }
    }
    void endElement(const XML_Char *) { }
    void characters(const XML_Char *, int ) { }
  };


  bool load_config(string filename, msconfig &conf) {
    // Re-populates fixed modifications, isotopes, and variable modifications. 
    conf.resetMods(); 
    conf.resetIsotopes();
    conf.resetVarMods();

    ConfigXMLHandler *handler = new ConfigXMLHandler(conf);
    cout << "Loading " << filename << endl;
    handler->parse(filename);
    bool ok = handler->ok();
    delete handler;
    return ok;
  }

  // Default handling of the output of a binary paramter. 
  void B(ofstream &config, string key, bool value) { if(value) config << key << "=\"1\" ";  else config << key << "=\"0\" "; }

  // Generate an XML configuration file. 
  bool save_config(string filename, msconfig &C) {
    ofstream config(filename.c_str());
    if(config.bad()) return false;

    // Data tab.
    config << "<pview>" << endl;
    config << "<data "; 
    config << "quantmode=\"" << (int)C.quantmode << "\" ";
    config << "threads=\"" << C.threads << "\" ";
    B(config, "keep_raw", C.keep_raw);
    B(config, "keep_filtered", C.keep_filtered);
    config << "load_ms2=\""; if(C.load_ms2) config << "1"; else config << "0"; config << "\" ";
    config << "keep_ms2=\""; if(C.keep_ms2) config << "1"; else config << "0"; config << "\" ";
    config << "run_algorithms=\""; if(C.run_algorithms) config << "1"; else config << "0"; config << "\" ";
    config << "/>" << endl << endl;;

    // XICs tab
    config << "<xics ";
    config << "dmz=\"" << C.dmz << "\" ";
    //config << "dtime=\"" << C.dtime << "\" ";
    //config << "peakthresh=\"" << C.peakthresh << "\" ";
    config << "hist_log10I_thresh=\"" << C.hist_log10I_thresh << "\" ";
    //config << "peak_dtime=\"" << C.peak_dtime << "\" ";
    config << "peak_dmz=\"" << C.peak_dmz << "\" ";
    config << "max_xiclen=\"" << C.max_xiclen << "\" ";
    config << "isotolerence=\"" << C.isotolerence << "\" ";
    config << "xic_peaks_missed=\"" << C.xic_peaks_missed << "\" ";
    config << "min_xic_peaks=\"" << C.min_xic_peaks << "\" ";
    config << "/>" << endl << endl;

    // (label)-free tab
    config << "<free ";
    config << "align_dtime=\"" << C.align_dtime << "\" ";
    config << "min_conds=\"" << C.min_conds << "\" ";
    config << "min_reps=\"" << C.min_reps << "\" ";
    config << "min_iruns=\"" << C.min_iruns << "\" ";
    config << "xic_width=\"" << C.xic_width << "\" ";
    config << "group_dtime=\"" << C.group_dtime << "\" ";
    B(config, "align_translation", C.align_translation);
    B(config, "align_nonlinear", C.align_nonlinear);
    B(config, "skip_normalize_label_free", C.skip_normalize_label_free);
    config << "/>" << endl << endl;


    // Isotope-(labelled) tab.
    config << "<isotope ";
    config << "max_label_count=\"" << C.max_label_count << "\" ";
    B(config, "skip_normalize_ratios", C.skip_normalize_ratios);
    B(config, "isotope_15N", C.isotope_15N);
    B(config, "isotope_by_search", C.isotope_by_search);
    B(config, "acetyl_stoich_plus5",C.acetyl_stoich_plus5); 
    config << "/>" << endl << endl;

    // MS2 tab.
    config << "<ms2 ";
    config << "ms2tol=\"" << C.ms2tol << "\" ";
    B(config, "ms2tol_usedaltons", C.ms2tol_usedaltons);
    config << "mstol=\"" << C.mstol << "\" " ; 
    config << "missed_cleaves=\"" << C.missed_cleaves << "\" "; 
    config << "ms2_fdr=\"" << C.ms2_fdr <<  "\" "; 
    config << "enzyme_id=\"" << C.enzyme_id << "\" ";
    config << "max_mod_aa=\"" << C.max_mod_aa << "\" ";
    B(config, "fixed_15N", C.fixed_15N);
    config << "/>" << endl << endl;

    // Output fixed mods.
    for(size_t i = 0; i < C.fxMods.size(); i++) {
      fxmod_t &fx = C.fxMods[i];
      if(fx.system) { // If a system modification, binary active or not.
	if(C.fxMods[i].active) config << "<sysfixedmod activeindex=\"" << i << "\" />" << endl << endl;;
      }
      else {
	// Otherise, save all of the data from a user-defined fixed modification.
	config << "<userfixedmod ";
	config << "description=\"" << fx.description << "\" ";
	config << "delta_mass=\"" << fixed << setprecision(10) << fx.deltaMass << "\" ";
	B(config, "active", fx.active);
	config << "aa=\"" << fx.AAs << "\" ";
	config << "/>" << endl << endl;
      }
    }

    // Output variable modifications.
    for(size_t i = 0; i < C.varMods.size(); i++) {
      varmod_t &vr = C.varMods[i];
      if(vr.system) {
	if(C.varMods[i].active) config << "<sysvarmod activeindex=\"" << i << "\" />" << endl << endl;;
      }
      else {
	config << "<uservarmod ";
	config << "d_prec_mass=\"" << fixed << setprecision(10) << vr.d_prec_mass << "\" ";
	config << "max_count=\"" << vr.max_count << "\" ";
	config << "description=\"" << vr.description << "\" ";
	B(config, "active", vr.active);
	ostringstream fragmodStr;
	for(size_t i = 0; i < vr.fragmods.size(); i++) {
	  fragmodStr << fixed << setprecision(10) << vr.fragmods[i].dfrag_mass << "|";
	  fragmodStr << vr.fragmods[i].AAs << "|";
	  fragmodStr << vr.fragmods[i].abbreviation; 
	  if( i != vr.fragmods.size() - 1) fragmodStr << "|";
	}
	config << "fragmods=\"" << fragmodStr.str() << "\" ";
	config << "/>" << endl << endl;
      }
    }

    // Output user defined isotopes
    for(size_t i = 0; i < C.isotopes.size(); i++) {
      isotope_t &iso = C.isotopes[i];
      if(iso.system == false) {
	config << "<userisotope ";
	config << "description=\"" << iso.description << "\" ";
	config << "shift_heavy=\"" << fixed << setprecision(10) << iso.shift_heavy << "\" ";
	config << "shift_medium=\"" << fixed << setprecision(10) << iso.shift_medium << "\" ";
	config << "aa=\"" << iso.AAs << "\" ";
	config << "/>" << endl << endl;
      }
    }

    // Output active heavy isotopes.
    for(size_t i = 0; i < C.activeIsotopes.size(); i++) 
      config << "<active_isotope index=\"" << C.activeIsotopes[i] << "\" />" << endl << endl;

    config << "</pview>" << endl;

    config.close();
    return true;
  }


  // Instrument runs through PepXML file and creates a mapping from
  // protein identifiers and protein descriptions.
  class PepXMLHandler1 : public XmlHandler {  
    map<string, string> &protein2descr;
  public:
    PepXMLHandler1(map<string, string> &p2d) : protein2descr(p2d) { }
    void startElement(const XML_Char *el, const XML_Char **attr) {
      if(isElement("search_hit", el) || isElement("alternative_protein", el)) {
	// Save mapping from a protein ID in tag protein to a
	// description in protein_descr.
	string protein = getValue("protein", attr);
	string protein_descr = getValue("protein_descr", attr);
	if(protein2descr.find(protein) == protein2descr.end()) protein2descr[protein] = protein_descr;
      }
    }
    void endElement(const XML_Char *) { }
    void characters(const XML_Char *, int ) { }
  };

  struct ms2peak_incr { bool operator() (ms2peak *a, ms2peak *b) const {  return a->scanNum < b->scanNum;  } };

  enum pepXMLstate { pepIGNORE, pepQUERY, pepHIT, pepHITDONE };
  class PepXMLHandler2 : public XmlHandler {  
    string spectrum, peptide;
    vector<string> proteins;
    pepXMLstate state;
    vector<ms2::modpos_t> sites;
    long scanNum;
    double score; bool foundScore;
    vector<string> &irun_descr, &seq_data;
    map<string, int> &protein2idx;
    vector< vector<ms2peak *> > &peaks;
  public:
    PepXMLHandler2(vector<string> &descr, vector< vector<ms2peak *> > &peaks, 
		   map<string, int> &p2idx, vector<string> &seqd) : 
      irun_descr(descr), seq_data(seqd), protein2idx(p2idx), peaks(peaks) {
      score = 1e-9;  foundScore = false;
      state = pepIGNORE; 
    }
    void startElement(const XML_Char *el, const XML_Char **attr) {
      if(state == pepIGNORE && isElement("spectrum_query", el)) {
	spectrum = getValue("spectrum", attr);
	scanNum = getValue("start_scan", attr, -1L);
	state = pepQUERY;
      }
      if(state == pepQUERY && isElement("search_hit", el)) {
	int hit_rank = getValue("hit_rank", attr, 0);
	if(hit_rank == 1) {
	  peptide = getValue("peptide", attr); // Get sequence
	  proteins.push_back(getValue("protein", attr)); // Collect protein IDs
	  state = pepHIT;
	  string peptide_prev_aa = getValue("peptide_prev_aa", attr);
	  string peptide_next_aa = getValue("peptide_next_aa", attr);

	  // Use N-term and C-term conventions in PVIEW
	  char Nterm = '[', Cterm = ']';
	  if(peptide_prev_aa.length() == 1 && peptide_prev_aa[0] == '-') Nterm = '#';
	  if(peptide_next_aa.length() == 1 && peptide_next_aa[0] == '-') Cterm = '*';
	  peptide = Nterm + peptide + Cterm;
	}
      }
      if(state == pepHIT && isElement("alternative_protein", el)) {
	string alt_protein = getValue("protein", attr);
	proteins.push_back(alt_protein);
      }
      if(state == pepHIT && isElement("mod_aminoacid_mass", el)) {
	ms2::modpos_t m;
	m.pos = getValue("position", attr, 0); // Note implicit subtract 1, but N-term convention used! See "search_hit" above.
	m.frag_mass = getValue("mass", attr, 0.0);
	m.prec_mass = m.frag_mass;  // This is always true, but PepXML doesn't let you get the precursor amino acid mass.
	m.mod_idx = m.fragmod_idx = -1;
	sites.push_back(m);
      }
      if(state == pepHIT && isElement("search_score", el)) {
	if(foundScore == false) {
	  score = max(0.0, getValue("value", attr, 0.0)); // Scores must be positive.
	  foundScore = true;
	}
      }
    }
    void endElement(const XML_Char *el) { 
      if(state == pepHIT && isElement("search_hit", el)) { state = pepHITDONE; }
      else if(isElement("spectrum_query", el)) {
	if(state == pepHITDONE) {
	  bool foundScan = false;  int scanIdx = 0;
	  for(int s = 0; s < (int)irun_descr.size(); s++) { 
	    if(spectrum.substr(0, spectrum.find_first_of('.')) == irun_descr[s]) {
	      foundScan = true; scanIdx = s;  break; 
	    }
	  }
	  if(foundScan) {
	    // Look for peak with current scan number.
	    ms2peak *respeak = NULL;
	    ms2peak query; query.scanNum = scanNum;
	    vector<ms2peak *>::iterator it = 
	      lower_bound(peaks[scanIdx].begin(), peaks[scanIdx].end(), &query, ms2peak_incr()); 
	    if(it != peaks[scanIdx].end()) {
	      ms2peak *q = *it;
	      if(q->scanNum == scanNum) respeak = q; // Found peak with scan number.
	    }
	    if(respeak != NULL) {
	      // Use scan number of find MS/MS peak.
	      int seq_idx = (int)seq_data.size(); 
	      seq_data.push_back(peptide); // "each entry in seq corresponds to an input sequence.
	      // Fill information for this MS/MS database hit.
	      for(size_t i = 0; i < proteins.size(); i++) {	    
		int pidx = protein2idx[proteins[i]]; // Map to index in vector<string> protein_meta;
		ms2::fragment f;
		f.protein_idx = pidx; f.seq_idx = seq_idx;
		f.pos = 0; f.len = peptide.length();  f.misses = 0; 
		f.fragmass = 0; 
		f.is_decoy = false; // TODO: Fill this data?
		
		ms2id_t id; // Save MS/MS ID for this peak. 
		id.frag = f;
		id.mods.sites = sites;
		id.score = score; id.qvalue = 1;
		respeak->results.push_back(id);
	      }
	    }
	  }
	}
	spectrum = peptide = ""; proteins.clear(); scanNum = -1; sites.clear(); score = 1e-9; foundScore = false;
	state = pepIGNORE;
      }
    }
    void characters(const XML_Char *, int ) { }
  };


  bool load_pepxml(vector<string> &pepxml, vector<irun*> &iruns, fasta_data &fasta) {
    // Create mapping from protein ID to protein description.
    map<string, string> protein2descr;
    for(size_t i = 0; i < pepxml.size(); i++) {
      cout << "Scanning proteins in " << pepxml[i] << endl;
      PepXMLHandler1 *handler1 = new PepXMLHandler1(protein2descr);
      handler1->parse(pepxml[i]);
      delete handler1;
    }

    // Populate protein meta and create map from protein ID to index
    // in protein meta.
    map<string, int> protein2idx;
    map<string, string>::iterator it = protein2descr.begin();
    int metaIdx = 0;
    while(it != protein2descr.end()) {
      fasta.protein_meta.push_back(it->first + " " + it->second);
      protein2idx[it->first] = metaIdx;
      metaIdx++;  it++;
    }

    vector<string> descriptions(iruns.size());
    vector< vector<ms2peak *> > peaks(iruns.size());
    for(size_t s = 0; s < iruns.size(); s++) {
      iruns[s]->get_monoiso_ms2(peaks[s]);
      sort(peaks[s].begin(), peaks[s].end(), ms2peak_incr()); // Sort pointers on scan number.
      descriptions[s] = iruns[s]->description;
    }

    for(size_t i = 0; i < pepxml.size(); i++) {
      PepXMLHandler2 *handler = new PepXMLHandler2(descriptions, peaks, protein2idx, fasta.seq_data);
      cout << "Loading " << pepxml[i] << endl;
      handler->parse(pepxml[i]);
      delete handler;
    }
    cout << "fasta.seq_data=" << fasta.seq_data.size() << endl;

    return true;
  }

}

