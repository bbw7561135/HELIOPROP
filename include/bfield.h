#ifndef HELIOPROP_BFIELD_H
#define HELIOPROP_BFIELD_H

#include <vector>

//#include "TMath.h"
//#include "Math/Minimizer.h"
//#include "Math/GSLMinimizer.h"
//#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "constants.h"
#include "TEnvironment.h"

class Bfield {

public:
	Bfield(double r0_,
			double B0_,
			double Ac_,
			double Omega_,
			double Vsw_,
			double alpha_,
			double phi0_,
			double lambda0_,
			double kpfact_,
			double deltaKpar_,
			double b_,
			double c_) :
				r0(r0_),
				Be(B0_),
				Ac(Ac_),
				Omega(Omega_),
				Vswmin(Vsw_),
				Vsw(Vsw_),
				alpha(alpha_),
				phi0(phi0_),
				lambda0(lambda0_),
				Kperp_factor0(kpfact_),
				Kperp_factor(kpfact_),
				deltaKpar(deltaKpar_),
				b(b_),
				c(c_) {
		sinalpha = sin(alpha);
		Vswmax = 760.0 * kms2UAd;
		OmegaVsw = Omega / Vsw;
		B0 = Be / sqrt(1.0 + pow(OmegaVsw * r0, 2));
		r = 0;
		theta = 0;
		phi = 0;
		gamma = OmegaVsw * r * sin(theta);
		gamma_squared = pow(gamma, 2);
		psi = GetPsi();
		qf = 0;
		charge = 0;
		mom = 0;
		betavelocity = 0;
	}

	~Bfield() {
	}

	inline double dGammadtheta() const {
		return gamma / tan(theta);
	} // gamma*cos(theta)*(3.0-pow(cos(theta),2))/(1.0+pow(cos(theta),2));

	inline double GetPsi() const {
		return atan(gamma);
	}

	inline double dPsidr() const {
		return 1.0 / (1.0 + gamma_squared) * gamma / r;
	}
	inline double dPsidtheta() const {
		return 1.0 / (1.0 + gamma_squared) * dGammadtheta();
	}

	inline double GetVsw() const {
		return Vsw;
	}
	// B components
	inline double B_r() {
		return Ac * B0 * pow(r0 / r, 2) * Heaviside(theta - thetaprime());
	}
	inline double B_phi() {
		return -B_r() * gamma;
	}
	inline double B_total() {
		return fabs(Heaviside(theta - thetaprime())) * B0 * pow(r0 / r, 2)
		* sqrt(1.0 + gamma_squared); //sqrt(pow(Br(r,theta,phi),2) + pow(Bphi(r,theta,phi),2));
	}
	inline double dBtotal_dr() {
		return -B_total() / r * (2.0 + gamma_squared) / (1.0 + gamma_squared);
	}
	inline double dBtotal_dphi() {
		return 0;
	}
	inline double dBtotal_dtheta() {
		return -B_total() * gamma / (1.0 + gamma_squared) * dGammadtheta();
	}

	// Drift velocity
	void GetVd(std::vector<double>& vd);
	double distance_from_HCS(const double* xx);
	double deriv_distance_from_HCS(const double* xx, unsigned int up);
	double closest_distance();

	inline double rL() {
		return (3.3 / 149597870691.0) * mom / B_total() / fabs(charge);
	}
	inline double thetaprime() {
		return PiOver2() + ASin(sinalpha * sin(phi - phi0 + OmegaVsw * r));
	}
	inline double thetaprime(const double& r_, const double& phi_) {
		return PiOver2() + ASin(sinalpha * sin(phi_ - phi0 + OmegaVsw * r_));
	}
	inline double thetaprime_phi() {
		return sinalpha * cos(phi - phi0 + OmegaVsw * r)
		/ sqrt(1.0 - pow(sinalpha * sin(phi - phi0 + OmegaVsw * r), 2));
	}  // to be revised
	inline double thetaprime_r() {
		return OmegaVsw * thetaprime_phi();
	}  // to be revised
	inline double thetaprime_phi(const double& r_, const double& phi_) {
		return sinalpha * cos(phi_ - phi0 + OmegaVsw * r_)
		/ sqrt(1.0 - pow(sinalpha * sin(phi_ - phi0 + OmegaVsw * r_), 2));
	}  // to be revised
	inline double thetaprime_r(const double& r_, const double& phi_) {
		return OmegaVsw * thetaprime_phi(r_, phi_);
	}  // to be revised
	inline double beta() {
		if (alpha == 0) {
			return 0.0;
		}
		double b = atan(OmegaVsw * r / sin(psi)
				* sqrt(pow(sinalpha, 2) - pow(cos(thetaprime()), 2))
				/ sin(thetaprime()));
		return fabs(b) * sgn(cos(phi - phi0 + OmegaVsw * r));
	}
	inline void SetPhi0(const double& phi0_) {
		phi0 = phi0_;
	}

	void GetKTensors(std::vector<double>& K, std::vector<double>& Kderiv);

	inline double Kpar() {

		//double lambda_par = lambda0*mom/(Btot()/B0);//(1.0+r/r0);
		double lambda_par = lambda0 / (B_total() / B0);  //(1.0+r/r0);
		//double lambda_par = lambda0;//*(1.0+r/r0);

		//if (mom >= 0.1*fabs(charge))
		double rig = mom / fabs(charge);

		//if (rig > 1.)
		lambda_par *= pow(rig, deltaKpar);

		// <!-- as in Potgieter paper: 1302.1284 -->
		if (b != 0 && c != 0)
			lambda_par *= pow(
					(pow(rig, c) + pow(4.0 / fabs(charge), c))
					/ (1.0 + pow(4.0 / fabs(charge), c)),
					(b - deltaKpar) / c);

		return lambda_par * betavelocity * Clight / 3.0;
	}

	inline double dKpardr() {
		//return 0;
		//return (mom >= 1) ? lambda0*mom/r0 : lambda0/r0;
		//return lambda0*mom/r0;
		return -Kpar() / B_total() * dBtotal_dr();
	}

	inline double dKpardphi() {
		//return 0;
		return -Kpar() / B_total() * dBtotal_dphi();

	}

	inline double dKpardtheta() {
		//return 0;
		return -Kpar() / B_total() * dBtotal_dtheta();
	}

	inline double dKrrdr() {
		return dKpardr() * (pow(cos(psi), 2) + Kperp_factor * pow(sin(psi), 2))
				+ Kpar() * (-1.0 + Kperp_factor) * sin(2.0 * psi) * dPsidr();
	}

	inline double dKrphidphi() {
		return -0.5 * sin(2.0 * psi) * (1.0 - Kperp_factor) * dKpardphi();
	}

	inline double dKrphidr() {
		return -0.5 * sin(2.0 * psi) * (1.0 - Kperp_factor) * dKpardr()
				- cos(2.0 * psi) * dPsidr() * (1.0 - Kperp_factor) * Kpar();
	}

	inline double dKthetathetadtheta() {
		return Kperp_factor * dKpardtheta();
	}  // to be revised

	inline double dKphiphidphi() {
		return dKpardphi()
				* (pow(sin(psi), 2) + Kperp_factor * pow(cos(psi), 2));
	}

	inline void Set(const double& r_, const double& theta_, const double& phi_,
			const double& qf_, double beta_, double mom_,
			const double& charge_) {
		r = r_;
		theta = theta_;
		phi = phi_;
		Vsw = Vswmin; //*(1.0+pow(cos(theta),2)); // Implements Fichtner et al. A&A 308:248 (1996), which is better than Bobik et al. 2012 because it is differentiable
		OmegaVsw = Omega / Vsw;
		gamma = OmegaVsw * r * sin(theta);
		gamma_squared = pow(gamma, 2);
		psi = GetPsi();
		if (theta < 30 * DegToRad() || theta > 150.0 * DegToRad())
			Kperp_factor = 10.0 * Kperp_factor0;
		else
			Kperp_factor = Kperp_factor0;
		qf = qf_;
		betavelocity = beta_;
		charge = charge_;
		mom = mom_;
		return;
	}

protected:
	double r0;
	double alpha;
	double phi0;
	double lambda0;
	double Kperp_factor;
	double Kperp_factor0;
	double Be;
	double B0;
	double Ac;
	double Omega;
	double Vsw;
	double Vswmin;
	double Vswmax;
	double r;
	double theta;
	double phi;
	double psi;
	double OmegaVsw;
	double sinalpha;
	double deltaKpar;
	double gamma;
	double gamma_squared;
	double qf;
	double betavelocity;
	double charge;
	double mom;
	double b;
	double c;
};

#endif

