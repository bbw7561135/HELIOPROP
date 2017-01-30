#ifndef HELIOPROP_CONSTANTS_H
#define HELIOPROP_CONSTANTS_H

#ifdef _OPENMP
#include <omp.h>
//#define OMP_NUM_THREADS 8
#endif

#include <cmath>
#include "TMath.h"

typedef enum {
	NUCLEUS, LEPTON
} PARTICLETYPE;

const double kms2UAd = (24.0 * 60.0 * 60.0) / 149597870.691; // conversion factor from km/s to UA/d
const double SunRadius = 0.005; // Radius of the Sun [ UA ] . When a particle falls in it is lost
const double MinEn = 0.001;  // Minimal energy to be considered [ GeV ]
const double Clight = TMath::C() / 1e3 * kms2UAd;

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

#endif

