#ifndef __CENTROIDER_HPP__
#define __CENTROIDER_HPP__

#include <vector>
#include <math.h>

namespace centroider {
  using namespace std;

  struct peak { 
    peak() { }
    peak(float mz_, float intensity_) { mz = mz_; intensity = intensity_; }
    float mz, intensity; 
    bool operator < (const peak p) const { return mz < p.mz; }
  };
  

  // 3x3 Cholesky matrix factorization.
  void Cholesky3x3(double& l11,  
		   double& l21, double& l22, 
		   double& l31, double& l32, double& l33,
		   double  a11, 
		   double  a21, double  a22,
		   double  a31, double  a32, double  a33) { 
    l11 = sqrt(a11);  
    l21 = a21 / l11; l22 = sqrt(a22 - l21 * l21);
    l31 = a31 / l11; l32 = (a32 - l31 * l21) / l22; l33 = sqrt(a33 - l31 * l31 - l32 * l32);
  }

  // 3x3 solver for upper triangler system.
  void BackSubstitution3x3(double &x1, double &x2, double &x3,
			   double u11, double u12, double u13,
			   double u22, double u23,
			   double u33,
			   double b1,  double b2,  double b3) 
  {
    x3 = b3 / u33;
    x2 = (b2 - x3 * u23) / u22;
    x1 = (b1 - x2 * u12 - x3 * u13) / u11;
  }

  // 3x3 solver for lower trianglar system
  void ForwardSubstitution3x3(double &x1, double &x2, double &x3,
			      double l11,
			      double l21, double l22,
			      double l31, double l32, double l33,
			      double b1,  double b2,  double b3) {
    x1 = b1 / l11;
    x2 = (b2 - l21 * x1) / l22;
    x3 = (b3 - l31 * x1 - l32 * x2) / l33;
  }

  // Efficiently computes xi^k where k is an integer.
  double powk(double xi, int k) {
    double r = xi;
    for(int i = 1; i <= k - 1; i++) r *= xi;
    return r;
  }

  // Fits a polynomial of degree 2 given a vector of (x,y) pairs.

  // TODO: Figure out how to get residual.
  void polyfit3(double &a, double &b, double &c, vector<double> &x, vector<double> &y) {
    size_t N = x.size();

    // Space for the lower triangle of the Vandermonde 3x3 matrix.
    double v11;
    double v21, v22;
    double v31, v32, v33;

    // Fill the Vandermonde matrix.
    v11 = (double)x.size();
    v21 = 0; 
    for(size_t i = 0; i < N; i++) v21 += x[i];
    v31 = 0;
    for(size_t i = 0; i < N; i++) v31 += powk(x[i], 2);
    v22 = v31;
    v32 = 0;
    for(size_t i = 0; i < N; i++) v32 += powk(x[i], 3);
    v33 = 0;
    for(size_t i = 0; i < N; i++) v33 += powk(x[i], 4);

    // Form b for normal equations.
    double b1, b2, b3;
    b1 = 0; for(size_t i = 0; i < N; i++) b1 += y[i];
    b2 = 0; for(size_t i = 0; i < N; i++) b2 += x[i] * y[i];
    b3 = 0; for(size_t i = 0; i < N; i++) b3 += powk(x[i], 2) * y[i];


    // 1. Compute the Cholesky factorization of the Vandermonde matrix.
    double l11;
    double l21, l22;
    double l31, l32, l33;

    Cholesky3x3(l11,
		l21, l22,
		l31, l32, l33,
		v11,
		v21, v22, 
		v31, v32, v33);

    // 2. Solve the lower trianglar system for w.
    double w1, w2, w3;
    ForwardSubstitution3x3(w1, w2, w3,
			   l11, 
			   l21, l22,
			   l31, l32, l33,
			   b1,   b2,  b3);

    double u11 = l11, u12 = l21, u13 = l31;
    double            u22 = l22, u23 = l32;
    double                       u33 = l33;

    // Solver the upper trianglar system for x.
    BackSubstitution3x3(a,     b,   c,
			u11, u12, u13,
			u22, u23,
			u33,  
			w1,   w2,  w3);  
  }

  bool IsMax(float C, float L1, float R1, float L2, float R2) {
    if (C > L1 && C > R1) return true;
    if (C > L2 && C == L1 && C > R1) return true;
    if (C > L1 && C == R1 && C > R2) return true;
    return false;
  }

  bool IsMaxStrict(float C, float L1, float R1, float L2, float R2) {
    return C > L1 && C > R1 && L1 > L2 && R1 > R2;
  }

  int CalcMinPeakIndex(vector<peak> &profile, int ind) {
    while (ind > 0 && 
	   profile[ind - 1].intensity > 0 &&  // stop at zero intensity
	   profile[ind - 1].intensity < profile[ind].intensity) { // decreasing intensity
      ind--;
    }
    return ind;
  }

  int CalcMaxPeakIndex(vector<peak> &profile, int ind) {
    while (ind < (int)profile.size() - 1 && 
	   profile[ind + 1].intensity > 0 && // stop at zero intensity
	   profile[ind + 1].intensity < profile[ind].intensity) { // decreasing intensity
      ind++;
    }
    return ind;
  }
  template<typename T> inline bool isnan_portable(T x) { return x != x; }
  template<typename T> inline bool isinf_portable(T x) { return !isnan_portable(x) && isnan_portable(x - x); }
  bool isfin_portable(float v) { return !isnan_portable(v) && !isinf_portable(v);  }

  // Centroider adapted from Maxquant by Cox and Mann. 
  void gausfit(vector<peak> &profile, vector<peak> &centroid) {
    if(profile.size() < 2) return;  
    vector<double> mz; mz.reserve(30);
    vector<double> logI; logI.reserve(30);

    //uint32_t failed_fit = 0;

    for (int i = 2; i < (int)profile.size() - 2; i++) {
      float L2 = profile[i - 2].intensity;
      float L1 = profile[i - 1].intensity;
      float C  = profile[i].intensity;
      float R1 = profile[i + 1].intensity;
      float R2 = profile[i + 2].intensity;
      // Skip if no peak is found.
      if(IsMax(C, L1, R1, L2, R2) == false) continue;

      peak p;

      // [minInd, maxInd] the interval is **inclusive**
      int minInd = CalcMinPeakIndex(profile, i);
      int maxInd = CalcMaxPeakIndex(profile, i);

      int len = maxInd - minInd + 1;
      if(len == 1) p = profile[minInd]; // assigns both m/z and intensity
      else if (len == 2) {
	p.intensity = max(profile[minInd].intensity, profile[maxInd].intensity);
	double i1 = profile[minInd].intensity, i2 = profile[maxInd].intensity;
	double mz1 = profile[minInd].mz, mz2 = profile[maxInd].mz;
	p.mz = (i1 * mz1 + i2 * mz2) / (i1 + i2);
      }
      else {
	// Get maximum intensity value.
	float maxI = 0;
	for(int k = minInd; k <= maxInd; k++) 
	  if(profile[k].intensity > maxI) maxI = profile[k].intensity;

	// Center m/z points in the middle of the range.  Prevents
	// numeric problems in the Vendermonde matrix calculation.
	double mzCent = profile[i].mz;

	// Collect data points, compute log of intensity.
	mz.clear(); logI.clear();
	double maxShift = 0;
	for(int k = max(minInd, i - 2); k <= min(maxInd, i + 2); k++) {
	  double dmz = profile[k].mz - mzCent;
	  if(fabs(dmz) > maxShift) maxShift = fabs(dmz); // limit the m/z adjustment.
	  mz.push_back(dmz);
	  logI.push_back(log(profile[k].intensity));
	}

	// Fit a quadratic to (m/z, log I) points.  Only valid if I >
	// 0 and baseline of gaussian = 0.
	double a, b, c; polyfit3(a, b, c, mz, logI);

	// Use quadratic fit to get height, width, and position of gaussian peak.
	// DO NOT USE Height, erroneously estimates peak intensity 
	// double Height=exp( a - powk( c*( b/( 2*c ) ), 2) );
	// double Width=2.35703/(sqrt(2)*sqrt(-c));
	double Position=-b/(2*c);
	if(isfin_portable(Position) && fabs(Position) <= maxShift ) {
	  p.intensity = maxI; p.mz = mzCent + Position;
	}
	else { p.intensity = maxI; p.mz = mzCent; /*failed_fit++;*/ }
      }
      centroid.push_back(p);
    }
    //if(failed_fit > 0) cout << "failed_fit=" << failed_fit << endl;
  }

  
};

#endif
