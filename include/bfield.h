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

class TBfield {

public:
	TBfield() {
		r0 = 0;
		Be = 0;
		Ac = 0;
		Omega = 0;
		Vsw = 0;
		alpha = 0;
		phi_0 = 0;
		lambda_0 = 0;
		Kperp_factor_constant = 0;
		set(0);
	}

	TBfield(double r0_,
			double Be_,
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
				Be(Be_),
				Ac(Ac_),
				Omega(Omega_),
				Vsw(Vsw_),
				alpha(alpha_),
				phi_0(phi0_),
				lambda_0(lambda0_),
				Kperp_factor_constant(kpfact_),
				delta_Kpar(deltaKpar_),
				b(b_),
				c(c_) {
		sin_alpha = sin(alpha);
		Vsw_max = 760.0 * kms2UAd;
		Vsw_min = Vsw;
		B0 = Be / sqrt(1.0 + pow(Omega_Vsw * r0, 2));
		reference_rigity = 4.0; // GV
		set(0);
	}

	~TBfield() {
	}

	void set(double value) {
		set(value, value, value, value, value, value, value);
	}

	void set(double r_, double theta_, double phi_, double qf_, double beta_, double momentum_, double charge_) {
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

	inline double GetVsw() const {
		return Vsw;
	}

	// B components
	double dGamma_dtheta() const;
	double B_r() const;
	double B_phi() const;
	double B_total() const;
	double dBtotal_dr() const;
	double dBtotal_dphi() const;
	double dBtotal_dtheta() const;

	// Drift velocity
	void getDriftVelocity(std::vector<double>& vd);
	double thetaprime() const;
	double thetaprime(const double& r_, const double& phi_) const;
	double thetaprime_phi() const;
	double thetaprime_r() const;
	double thetaprime_phi(const double& r_, const double& phi_) const;
	double thetaprime_r(const double& r_, const double& phi_) const;
	double distance_from_HCS(const double* xx);
	double deriv_distance_from_HCS(const double* xx, unsigned int up);
	double closest_distance();

	/*inline double rL() {
		return (3.3 / 149597870691.0) * momentum / B_total() / fabs(charge);
	}*/

	inline double beta() {
		if (alpha == 0) {
			return 0.0;
		}
		double b = atan(Omega_Vsw * r / sin(psi)
				* sqrt(pow(sin_alpha, 2) - pow(cos(thetaprime()), 2))
				/ sin(thetaprime()));
		return fabs(b) * sgn(cos(phi - phi_0 + Omega_Vsw * r));
	}

	inline void SetPhi0(const double& phi0_) {
		phi_0 = phi0_;
	}

	// Diffusion Tensor
	double dPsi_dr() const;
	double dPsi_dtheta() const;
	double Kpar() const;
	double dKpar_dr() const;
	double dKpar_dphi() const;
	double dKpar_dtheta() const;
	double dKrr_dr() const;
	double dKrphi_dphi() const;
	double dKrphi_dr() const;
	double dKthetatheta_dtheta() const;
	double dKphiphi_dphi() const;
	void getKTensors(std::vector<double>& K, std::vector<double>& Kderiv);

protected:
	double r0;
	double alpha;
	double phi_0;
	double lambda_0;
	double Kperp_factor;
	double Kperp_factor_constant;
	double Be;
	double B0;
	double Ac;
	double Omega;
	double Vsw;
	double Vsw_min;
	double Vsw_max;
	double r;
	double theta;
	double phi;
	double psi;
	double Omega_Vsw;
	double sin_alpha;
	double delta_Kpar;
	double gamma;
	double gamma_squared;
	double qf;
	double betavelocity;
	double charge;
	double momentum;
	double b;
	double c;
	double reference_rigity;
};

#endif
