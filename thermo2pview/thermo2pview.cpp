 #include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

#include "thermo2pview.hpp"
#include "centroider.hpp"
#include "b64.h"

#include <QtEndian>
#include <QMessageBox>
#include <QFileDialog>
#include <QProgressDialog>

#ifdef _MSC_VER
#include "win_stdint.h"
#import "C:\\Program Files\\Thermo\\MSFileReader\\XRawFile2_x64.dll"
#endif

ThermoWin::ThermoWin(QMainWindow *parent) : QMainWindow(parent) {  
  setupUi(this);  
}

struct rawpeak { double mass; double intensity; };

struct PrecursorInfo {
  double dIsolationMass; double dMonoIsoMass; long nChargeState; long nScanNumber;
};

void ThermoWin::on_testButton_clicked(bool) {
  QString dir = QFileDialog::getExistingDirectory(this, "Open Directory");
  if(dir == "") return;

  vector<string> rawfiles, mzXMLfiles;
  QDirIterator it_dir(dir);
  while(it_dir.hasNext()) {
    it_dir.next();
    // Now only get files.
    QFileInfo info_file = it_dir.fileInfo(); 
    if(info_file.isDir()) continue;
    // Save all RAW files that were found.
    if(QString::compare(info_file.suffix(), "RAW", Qt::CaseInsensitive) == 0) {
      rawfiles.push_back(info_file.filePath().toStdString());
      QString mzXMLname = dir + "/" + info_file.baseName() + ".mzXML";
      mzXMLfiles.push_back(mzXMLname.toStdString());
    }
  }
  if(rawfiles.size() == 0) {
    QMessageBox::information(this, "Not found", "Did not find any Thermo .RAW files in " + dir + ".");
    return;
  }
  

#ifdef _MSC_VER 
  CoInitialize(NULL);

  using namespace MSFileReaderLib;
  IXRawfile4Ptr raw4(NULL);
  raw4.CreateInstance("MSFileReader.XRawfile.1", NULL, CLSCTX_INPROC_HANDLER | CLSCTX_INPROC_SERVER);

  for(size_t i = 0; i < rawfiles.size(); i++) {
    long nRet = raw4->Open(rawfiles[i].c_str());
    if(nRet != 0) QMessageBox::critical(this, "Error", "Unable to open .RAW file.");

    ofstream mzXML(mzXMLfiles[i].c_str());
    mzXML << "<mzXML>" << endl;

    raw4->SetCurrentController(0, 1);

    long first = 0; raw4->GetFirstSpectrumNumber(&first);
    long last = 0;  raw4->GetLastSpectrumNumber(&last);

    QString message = "Converting " + QString(rawfiles[i].c_str()) + ".";
    QProgressDialog progress(message, 0, first, last, this);
    progress.setWindowTitle("Conversion Progress");
    progress.setWindowModality(Qt::WindowModal);
    progress.setMinimumDuration(0);
    QCoreApplication::processEvents();
    for(long scanNum = first;  scanNum <= last; scanNum++) {
      progress.setValue(scanNum);
      long MSOrder; raw4->GetMSOrderForScanNum(scanNum, &MSOrder); 

      if((MSOrder == 0 || MSOrder == 1 || MSOrder == 2) == false) continue;

      double retentionTime; raw4->RTFromScanNum(scanNum, &retentionTime);
      retentionTime *= 60.0; // Convert to seconds.

      mzXML << "\t<scan num=\"" <<  scanNum << "\" "; 
      mzXML << "retentionTime=\"P" << fixed << setprecision(15) << retentionTime << "S\" "; 
      mzXML << "msLevel=\"" << MSOrder << "\" ";
      mzXML << ">" << endl;

      if(MSOrder == 2) {
	long ActivationType;raw4->GetActivationTypeForScanNum(scanNum, MSOrder, &ActivationType);
	string ActivationTypeStr;
	switch(ActivationType) {
	case 0: ActivationTypeStr = "CID"; break;
	case 1: ActivationTypeStr = "MPD"; break;
	case 2: ActivationTypeStr = "ECD"; break;
	case 3: ActivationTypeStr = "PQD"; break;
	case 4: ActivationTypeStr = "ETD"; break;
	case 5: ActivationTypeStr = "HCD"; break;
	case 6: ActivationTypeStr = "Any activation type"; break;
	case 7: ActivationTypeStr = "SA"; break;
	case 8: ActivationTypeStr = "PTR"; break;
	case 9: ActivationTypeStr = "NETD"; break;
	case 10: ActivationTypeStr = "NPTR"; break;
	default: break;
	}

	// TODO: Use the filterLine for the scan number to determine if supplemental activation was used. 
	// You need to search for "sa" in the filter line.

	// This function has trouble getting the monoisotopic m/z.
	double precursorMass; raw4->GetPrecursorMassForScanNum(scanNum, MSOrder, &precursorMass);

	// This call used by ReAdW produces a more accurate monoisotopic m/z. I'm not sure
	// why. Need to ask the Thermo people. 
	VARIANT varPrecursor; VariantInit(&varPrecursor);
	raw4->GetTrailerExtraValueForScanNum(scanNum, "Monoisotopic M/Z:" , &varPrecursor);
	if( varPrecursor.vt == VT_R4 ) precursorMass = varPrecursor.fltVal;
	else if( varPrecursor.vt == VT_R8 ) precursorMass = varPrecursor.dblVal;

	int precursorCharge = 0;
	VARIANT varCharge; VariantInit(&varCharge);
	raw4->GetTrailerExtraValueForScanNum(scanNum, "Charge State:" , &varCharge);

	switch (varCharge.vt) {
	case VT_I1: precursorCharge = varCharge.cVal; break;
	case VT_UI1: precursorCharge = varCharge.bVal; break;
	case VT_I2: precursorCharge = varCharge.iVal; break;
	case VT_UI2: precursorCharge = varCharge.uiVal; break;
	case VT_I4: precursorCharge = varCharge.lVal; break;
	case VT_UI4: precursorCharge = varCharge.ulVal; break;
	case VT_INT: precursorCharge = varCharge.intVal; break;
	case VT_UINT: precursorCharge = varCharge.uintVal; break;
	default: break;
	}	
	

	mzXML << "\t<precursorMz activationMethod=\"" << ActivationTypeStr << "\" ";
	mzXML << "precursorCharge=\"" << precursorCharge<< "\" >";
	mzXML << fixed << setprecision(15) << precursorMass;
	mzXML << "</precursorMz>" << endl;	
      }


      long isProfile; raw4->IsProfileScanForScanNum(scanNum, &isProfile);
      
      VARIANT massList; VariantInit(&massList);
      VARIANT peakFlags; VariantInit(&peakFlags);
      long arraySize;
      double peakWidth = 0;
      raw4->GetMassListFromScanNum(&scanNum, "", 0, 0, 0,  
				   (long)0, 
				   &peakWidth, 
				   &massList, &peakFlags, &arraySize);

      if(arraySize > 0) {
	using namespace centroider;
	// Get a pointer to the SafeArray
	SAFEARRAY FAR* psa = massList.parray;
	rawpeak* rawpeaks = NULL; SafeArrayAccessData( psa, (void**)(&rawpeaks) );
	vector<peak> profile(arraySize);
	for( long j = 0; j < arraySize; j++ )  {
	    profile[j].mz = rawpeaks[j].mass;
	    profile[j].intensity = rawpeaks[j].intensity;
	}
	// Release the data handle
	SafeArrayUnaccessData( psa );

	vector<peak> centroid;
	if(isProfile) centroider::gausfit(profile, centroid); else centroid = profile;
	
	if(centroid.size() == 0) mzXML << "\t<peaks precision=\"32\"></peaks>" << endl;
	else {
	  size_t destLen = B64_NAMESPACE::b64_encode(NULL, sizeof(peak) * centroid.size(), NULL, 0);
	  string buffer;  buffer.resize(destLen);
      
	  uint32_t *to_htonl = (uint32_t *)&centroid[0];
	  for(unsigned int i = 0; i < centroid.size() * sizeof(peak) / sizeof(uint32_t); i++) to_htonl[i] = qToBigEndian(to_htonl[i]); 

	  B64_NAMESPACE::b64_encode(&centroid[0], centroid.size() * sizeof(peak), &buffer[0], destLen);
	  mzXML << "\t<peaks precision=\"32\">" << buffer << "</peaks>" << endl;
	}
      }
      else mzXML << "\t<peaks precision=\"32\"></peaks>" << endl;
      
      if(massList.vt != VT_EMPTY) {
	SAFEARRAY FAR* psa = massList.parray;
	massList.parray = NULL;
	SafeArrayDestroy(psa);
      }
      if(peakFlags.vt != VT_EMPTY) {
	SAFEARRAY FAR* psa = peakFlags.parray;
	peakFlags.parray = NULL;
	SafeArrayDestroy(psa);
      }
      mzXML << "\t</scan>" << endl;
    }
    mzXML << "</mzXML>" << endl;
    raw4->Close();
  }

  CoUninitialize();
#endif

}

// Application start point.
int main(int argc, char *argv[]) {
  // Launch main application.
  QApplication app(argc, argv);
  ThermoWin thermo;
  thermo.show();
  return app.exec();
}
