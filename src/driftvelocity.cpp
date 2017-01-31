#include "bfield.h"

void TBfield::getDriftVelocity(std::vector<double>& vd) {
	double signfactors = sgn(charge) * Ac;
	double velocity = betavelocity * c_light;
	double LarmorRadius = (3.3 / 149597870691.0) * momentum / B_total() / fabs(charge);
	double commonfactor = 2.0 / 3.0 * velocity * (LarmorRadius / r) / pow(1 + gamma_squared, 3.0 / 2.0) * signfactors * Heaviside(theta - thetaprime());
	vd[0] = -gamma / tan(theta) * commonfactor;
	vd[1] = -gamma * vd[0];
	vd[2] = gamma * (2.0 + gamma_squared) * commonfactor;

	double LrL = closest_distance() / LarmorRadius;

	if (LrL <= 2) {
		double vns = (0.457 - 0.412 * LrL + 0.0915 * pow(LrL, 2)) * signfactors
				* velocity;
		double betaval = beta();
		vd[0] += (vns * (cos(betaval) * sin(psi)));
		vd[1] += (vns * (cos(betaval) * cos(psi)));
		vd[2] += (vns * (sin(betaval)));
	}
	double momfactor = 10.0 * pow(momentum / fabs(charge), 2)
	/ (1.0 + 10.0 * pow(momentum / fabs(charge), 2));
	vd[0] *= momfactor;
	vd[1] *= momfactor;
	vd[2] *= momfactor;
}

double TBfield::beta() const {
	if (alpha == 0) {
		return 0.0;
	}
	double b = atan(Omega_Vsw * r / sin(psi)
			* sqrt(pow(sin_alpha, 2) - pow(cos(thetaprime()), 2))
			/ sin(thetaprime()));
	return fabs(b) * sgn(cos(phi - phi_0 + Omega_Vsw * r));
}

double TBfield::thetaprime() const {
	return PiOver2() + ASin(sin_alpha * sin(phi - phi_0 + Omega_Vsw * r));
}

double TBfield::thetaprime(const double& r_, const double& phi_) const {
	return PiOver2() + ASin(sin_alpha * sin(phi_ - phi_0 + Omega_Vsw * r_));
}

double TBfield::thetaprime_phi() const {
	return sin_alpha * cos(phi - phi_0 + Omega_Vsw * r)
	/ sqrt(1.0 - pow(sin_alpha * sin(phi - phi_0 + Omega_Vsw * r), 2));
}  // to be revised

double TBfield::thetaprime_r() const {
	return Omega_Vsw * thetaprime_phi();
}  // to be revised

double TBfield::thetaprime_phi(const double& r_, const double& phi_) const {
	return sin_alpha * cos(phi_ - phi_0 + Omega_Vsw * r_)
	/ sqrt(1.0 - pow(sin_alpha * sin(phi_ - phi_0 + Omega_Vsw * r_), 2));
}  // to be revised

double TBfield::thetaprime_r(const double& r_, const double& phi_) const {
	return Omega_Vsw * thetaprime_phi(r_, phi_);
}  // to be revised

double TBfield::distance_from_HCS(const double* xx) {
	double thetap = thetaprime(xx[0], xx[1]);
	double rp = xx[0];
	double phip = xx[1];
	double xmxp2 = pow(r * sin(theta) * cos(phi) - rp * sin(thetap) * cos(phip), 2);
	double ymyp2 = pow(r * sin(theta) * sin(phi) - rp * sin(thetap) * sin(phip), 2);
	double zmzp2 = pow(r * cos(theta) - rp * cos(thetap), 2);
	return xmxp2 + ymyp2 + zmzp2;
}

double TBfield::deriv_distance_from_HCS(const double* xx, unsigned int up) {
	double thetap = thetaprime(xx[0], xx[1]);
	double rp = xx[0];
	double phip = xx[1];
	double xmxp = (r * sin(theta) * cos(phi) - rp * sin(thetap) * cos(phip));
	double ymyp = (r * sin(theta) * sin(phi) - rp * sin(thetap) * sin(phip));
	double zmzp = (r * cos(theta) - rp * cos(thetap));

	if (up == 0) {
		// Derivative in r direction
		return 2.0 * xmxp
				* (-sin(thetap) * cos(phip)
						- rp * cos(phip) * cos(thetap)
						* thetaprime_r(xx[0], xx[1]))
						+ 2.0 * ymyp
						* (-sin(thetap) * sin(phip)
								- rp * sin(phip) * cos(thetap)
								* thetaprime_r(xx[0], xx[1]))
								+ 2.0 * zmzp
								* (-cos(thetap)
										+ rp * sin(thetap) * thetaprime_r(xx[0], xx[1]));
	} else {
		// Derivative in phi direction
		return 2.0 * xmxp
				* (rp * sin(phip) * sin(thetap)
						- rp * cos(phip) * cos(thetap)
						* thetaprime_phi(xx[0], xx[1]))
						+ 2.0 * ymyp
						* (-rp * sin(thetap) * cos(phip)
								- rp * sin(phip) * cos(thetap)
								* thetaprime_phi(xx[0], xx[1]))
								+ 2.0 * zmzp * (rp * sin(thetap) * thetaprime_phi(xx[0], xx[1]));
	}

	return -1;
}

double TBfield::closest_distance() {

	if (alpha == 0) {
		return fabs(r * cos(theta));
	}

	/*ROOT::Math::GSLMinimizer min(ROOT::Math::kVectorBFGS2);

	min.SetMaxFunctionCalls(10000);
	min.SetMaxIterations(10000);
	min.SetTolerance(0.1);
	min.SetPrintLevel(0);

	ROOT::Math::GradFunctor f(this, &Bfield::distance_from_HCS,
			&Bfield::deriv_distance_from_HCS, 2);
	double step[2] = { 0.01, 0.01 };
	double variable[2] = { r + 1.0, phi };
	min.SetFunction(f);

	// Set the free variables to be minimized!
	min.SetLowerLimitedVariable(0, "rprime", variable[0], step[0], 0);
	min.SetLimitedVariable(1, "phiprime", variable[1], step[1], 0,
			TMath::TwoPi());

#ifdef _OPENMP
#pragma omp critical
#endif
	{
		min.Minimize();
	}

	const double *xs = min.X();

	return std::sqrt(distance_from_HCS(xs));*/

	return 0;
}



