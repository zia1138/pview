#ifndef __FFT_HPP__
#define __FFT_HPP__

#include <complex>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

namespace numeric {
  using namespace std;

  // Returns true of integer is a lower of 2. 
  inline bool isPow2(unsigned int v) { return v && !(v & (v - 1)); }

  // Returns the floor form of binary logarithm for a 32 bit integer
  inline int floorLog2(unsigned int v) {  unsigned r = 0;  while (v >>= 1) r++; return r; }

  inline int nearestPow2(int v) { if(isPow2(v)) return 1 << floorLog2(v); else return 1 << (floorLog2(v) + 1); }

  typedef complex<double> complex_t; // Use doubles for now.

  // Adapted from:
  // http://en.literateprograms.org/Cooley-Tukey_FFT_algorithm_%28C%29
  complex_t* DFT(complex_t* x, int N /* must be a power of 2 */, bool inverse) {
    complex_t* X = new complex_t[N];

    if (N == 1) { X[0] = x[0]; return X; }

    complex_t *e = new complex_t[N/2], *d = new complex_t[N/2];
    for(int k = 0; k < N/2; k++) {
      e[k] = x[2*k];
      d[k] = x[2*k + 1];
    }

    complex_t *E = DFT(e, N/2, inverse), *D = DFT(d, N/2, inverse);

    delete[] e; delete[] d;

    for(int k = 0; k < N/2; k++) {
      // Multiply entries of D by the twiddle factors e^(-2*pi*i/N * k) 
      if(inverse) D[k] = polar(1.0,  2.0*M_PI*k/N) * D[k];
      else        D[k] = polar(1.0, -2.0*M_PI*k/N) * D[k];
    }

    for(int k = 0; k < N/2; k++) {
      X[k]       = E[k] + D[k]; 
      X[k + N/2] = E[k] - D[k];
    }

    delete[] D; delete[] E;
    return X;
  }

}

#endif // __FFT_HPP__
