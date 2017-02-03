#include "bfield.h"
#include "gsl_minimizer.h"

void Bfield::getDriftVelocity(std::vector<double>& vd) {
	double signfactors = sgn(charge) * Ac;
	double velocity = betavelocity * c_light;
	double LarmorRadius = (3.3 / 149597870691.0) * momentum / B_total()
																													/ fabs(charge);
	double commonfactor = 2.0 / 3.0 * velocity * (LarmorRadius / r)
																													/ pow(1 + gamma_squared, 3.0 / 2.0) * signfactors
																													* Heaviside(theta - thetaprime());
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

double Bfield::beta() const {
	if (alpha == 0) {
		return 0.0;
	}
	double b = atan(
			Omega_Vsw * r / sin(psi)
			* sqrt(pow(sin_alpha, 2) - pow(cos(thetaprime()), 2))
			/ sin(thetaprime()));
	return fabs(b) * sgn(cos(phi - phi_0 + Omega_Vsw * r));
}

double Bfield::thetaprime() const {
	return PiOver2 + ASin(sin_alpha * sin(phi - phi_0 + Omega_Vsw * r));
}

double Bfield::thetaprime(const double& r_, const double& phi_) const {
	return PiOver2 + ASin(sin_alpha * sin(phi_ - phi_0 + Omega_Vsw * r_));
}

double Bfield::thetaprime_phi() const {
	return sin_alpha * cos(phi - phi_0 + Omega_Vsw * r)
	/ sqrt(1.0 - pow(sin_alpha * sin(phi - phi_0 + Omega_Vsw * r), 2));
}  // to be revised

double Bfield::thetaprime_r() const {
	return Omega_Vsw * thetaprime_phi();
}  // to be revised

double Bfield::thetaprime_phi(const double& r_, const double& phi_) const {
	return sin_alpha * cos(phi_ - phi_0 + Omega_Vsw * r_)
	/ sqrt(1.0 - pow(sin_alpha * sin(phi_ - phi_0 + Omega_Vsw * r_), 2));
}  // to be revised

double Bfield::thetaprime_r(const double& r_, const double& phi_) const {
	return Omega_Vsw * thetaprime_phi(r_, phi_);
}  // to be revised

double Bfield::distance_from_HCS(const double& rp, const double& phip) {
	double thetap = thetaprime(rp, phip);
	double xmxp2 = pow(r * sin(theta) * cos(phi) - rp * sin(thetap) * cos(phip),
			2);
	double ymyp2 = pow(r * sin(theta) * sin(phi) - rp * sin(thetap) * sin(phip),
			2);
	double zmzp2 = pow(r * cos(theta) - rp * cos(thetap), 2);
	return xmxp2 + ymyp2 + zmzp2;
}

double Bfield::deriv_distance_from_HCS(const double& rp, const double& phip, unsigned int up) {
	double thetap = thetaprime(rp, phip);
	double xmxp = (r * sin(theta) * cos(phi) - rp * sin(thetap) * cos(phip));
	double ymyp = (r * sin(theta) * sin(phi) - rp * sin(thetap) * sin(phip));
	double zmzp = (r * cos(theta) - rp * cos(thetap));

	double sin_thp = sin(thetap);
	double cos_thp = cos(thetap);
	double sin_php = sin(phip);
	double cos_php = cos(phip);

	if (up == 0) {
		// Derivative in r direction
		double thetap_r = thetaprime_r(rp, phip);
		return 2.0 * xmxp * (-sin_thp * cos_php - rp * cos_php * cos_thp * thetap_r)
				+ 2.0 * ymyp * (-sin_thp * sin_php - rp * sin_php * cos_thp  * thetap_r)
				+ 2.0 * zmzp * (-cos_thp + rp * sin_thp * thetap_r);
	} else {
		// Derivative in phi direction
		double thetap_phi = thetaprime_phi(rp, phip);
		return 2.0 * xmxp * (rp * sin_php * sin_thp  - rp * cos_php * cos_thp * thetap_phi)
				+ 2.0 * ymyp * (-rp * sin_thp * cos_php - rp * sin_php * cos_thp * thetap_phi)
				+ 2.0 * zmzp * (rp * sin_thp * thetap_phi);
	}
	return -1;
}

double Bfield::closest_distance() {

	if (alpha == 0) {
		return fabs(r * cos(theta));
	}

	GSLMinimizer f;
	std::cout << f.minimize(this) << "\n";

	return 0;
}
