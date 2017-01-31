#include "bfield.h"

void TBfield::set(double value) {
	set(value, value, value, value, value, value, value);
}

void TBfield::set(double r_, double theta_, double phi_, double qf_, double beta_, double momentum_, double charge_) {
	r = r_;
	theta = theta_;
	phi = phi_;
	qf = qf_;
	betavelocity = beta_;
	charge = charge_;
	momentum = momentum_;
	Vsw = Vsw_min; // *(1.0+pow(cos(theta),2)); //Implements Fichtner et al. A&A 308:248 (1996), which is better than Bobik et al. 2012 because it is differentiable
	Omega_Vsw = Omega / Vsw;
	gamma = Omega_Vsw * r * sin(theta);
	gamma_squared = pow(gamma, 2);
	psi = atan(gamma);
	if (theta < 30.0 * DegToRad() || theta > 150.0 * DegToRad())
		Kperp_factor = 10.0 * Kperp_factor_constant;
	else
		Kperp_factor = Kperp_factor_constant;
}
