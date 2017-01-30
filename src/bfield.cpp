#include "bfield.h"

// Drift velocity
void Bfield::GetVd(std::vector<double>& vd) {

	double signfactors = sgn(charge) * Ac;

	double velocity = betavelocity * Clight;

	double LarmorRadius = rL();
	double commonfactor = (2.0 / 3.0) * velocity * (LarmorRadius / r)
			/ pow(1 + gamma2, 3.0 / 2.0) * signfactors
			* Heaviside(theta - thetaprime());
	vd[0] = -gamma / tan(theta) * commonfactor;
	vd[1] = -gamma * vd[0];
	vd[2] = gamma * (2.0 + gamma2) * commonfactor;

	double LrL = closest_distance() / LarmorRadius;

	if (LrL <= 2) {
		double vns = (0.457 - 0.412 * LrL + 0.0915 * pow(LrL, 2)) * signfactors
				* velocity;
		double betaval = beta();
		vd[0] += (vns * (cos(betaval) * sin(psi)));
		vd[1] += (vns * (cos(betaval) * cos(psi)));
		vd[2] += (vns * (sin(betaval)));

	}

	double momfactor = 10.0 * pow(mom / fabs(charge), 2)
			/ (1.0 + 10.0 * pow(mom / fabs(charge), 2));
	vd[0] *= momfactor;
	vd[1] *= momfactor;
	vd[2] *= momfactor;

	return;
}

double Bfield::distance_from_HCS(const double* xx) {

	double thetap = thetaprime(xx[0], xx[1]);
	double rp = xx[0];
	double phip = xx[1];
	double xmxp2 = pow(r * sin(theta) * cos(phi) - rp * sin(thetap) * cos(phip),
			2);
	double ymyp2 = pow(r * sin(theta) * sin(phi) - rp * sin(thetap) * sin(phip),
			2);
	double zmzp2 = pow(r * cos(theta) - rp * cos(thetap), 2);

	return xmxp2 + ymyp2 + zmzp2;

}

double Bfield::deriv_distance_from_HCS(const double* xx, unsigned int up) {
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

double Bfield::closest_distance() {

	if (alpha == 0) {
		return fabs(r * cos(theta));
	}

	ROOT::Math::GSLMinimizer min(ROOT::Math::kVectorBFGS2);

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

	return sqrt(distance_from_HCS(xs));
}

void Bfield::GetKTensors(std::vector<double>& K, std::vector<double>& Kderiv) {

	double kappa_par = Kpar();
	double cpsi = cos(psi);
	double spsi = sin(psi);

	K[0] = kappa_par * (pow(cpsi, 2) + Kperp_factor * pow(spsi, 2)); // Krr
	K[1] = -kappa_par * (1.0 - Kperp_factor) * spsi * cpsi;             // Krphi
	K[2] = Kperp_factor * kappa_par;                              // Kthetatheta
	K[3] = kappa_par * (Kperp_factor * pow(cpsi, 2) + pow(spsi, 2));  // Kphiphi

	Kderiv[0] = 2.0 * r * K[0] + r * r * dKrrdr(); // d/dr (r^2*Krr)
	Kderiv[1] = dKrphidphi(); // d/dphi Krphi
	Kderiv[2] = cos(theta) * K[2] + sin(theta) * dKthetathetadtheta(); // d/dtheta (sin(theta)*Kthetatheta)
	Kderiv[3] = dKphiphidphi(); // d/dphi Kphiphi
	Kderiv[4] = K[1] + r * dKrphidr(); // d/dr (r*Krphi)
}

