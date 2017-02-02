#ifndef HELIOPROP_CONSTANTS_H
#define HELIOPROP_CONSTANTS_H

#ifdef _OPENMP
#include <omp.h>
//#define OMP_NUM_THREADS 8
#endif

#include <cmath>

static const double kms2UAd = (24.0 * 60.0 * 60.0) / 149597870.691; // conversion factor from km/s to UA/d
static const double SunRadius = 0.005; // Radius of the Sun [ UA ] . When a particle falls in it is lost
static const double MinEn = 0.001;  // Minimal energy to be considered [ GeV ]
static const double c_light = 2.99792458e8 / 1e3 * kms2UAd;
static const double DegToRad = M_PI / 180.0;
static const double PiOver2 = M_PI / 2.0;
static const double TwoPi = 2.0 * M_PI;

inline double sgn(const double& x) {
	return (x > 0) - (x < 0);
}

inline double Heaviside(const double& x) {
	if (x > 0)
		return -1;
	if (x < 0)
		return 1;
	return 1e-30;
}

inline double ASin(const double& x) {
	if (x < -1.)
		return -M_PI / 2.;
	if (x >  1.)
		return  M_PI / 2.;
	return asin(x);
}

#endif

