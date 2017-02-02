#ifndef HELIOPROP_BFIELD_H
#define HELIOPROP_BFIELD_H

#include <vector>

#include "constants.h"
#include "input.h"
#include "TEnvironment.h"

class Bfield {

public:
	Bfield() {
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

	Bfield(double Be_,
			double Ac_,
			double Omega_,
			double Vsw_,
			double alpha_,
			double phi0_,
			double lambda0_,
			double kpfact_,
			double deltaKpar_) :
				Be(Be_),
				Ac(Ac_),
				Omega(Omega_),
				Vsw(Vsw_),
				alpha(alpha_),
				phi_0(phi0_),
				lambda_0(lambda0_),
				Kperp_factor_constant(kpfact_),
				delta_parallel_low(deltaKpar_) {
		r0 = 1;
		sin_alpha = sin(alpha);
		Vsw_max = 760.0 * kms2UAd;
		Vsw_min = Vsw;
		B0 = Be / sqrt(1.0 + pow(Omega_Vsw * r0, 2));
		reference_rigity = 4.0; // GV
		set(0);
	}

	// r0, Be, Ac, Omega, Vsw, alpha, phi0, lambda0
	/*	double bp[12] = { 1,  // Earth position [UA]
		input->MagField,  ///sqrt(2.0), // Reference magnetic field [T]
		input->polarity, // Magnetic field polarity
		TwoPi() / 27.0, // Solar differential rotation rate [d^-1]
		400.0 * kms2UAd, // Velocity of solar wind (constant, radial) [ km/s --> UA/d ]
		input->tiltangle * DegToRad(), // Tilt angle of current sheet [deg --> rad]
		0, // Reference phi position of current sheet (basically not used, but useful to make the current sheet corotate with Sun)
		input->lambda_par, // parallel mean-free-path [ UA ]
		input->Kperp_factor, input->delta, input->b, input->c };*/

	Bfield(const Input& input) {
		r0 = 1.0;
		Be = input.MagField;
		Ac = input.polarity;
		Omega = TwoPi / 27.0;
		Vsw = 400.0 * kms2UAd;
		alpha = input.tiltangle * DegToRad;
		phi_0 = 0.0;
		lambda_0 = input.lambda_par;
		Kperp_factor_constant = input.Kperp_factor;
		delta_parallel_low = input.delta_low;
		delta_parallel_high = input.delta_hi;
		//
		r0 = 1;
		sin_alpha = sin(alpha);
		Vsw_max = 760.0 * kms2UAd;
		Vsw_min = Vsw;
		B0 = Be / sqrt(1.0 + pow(Omega_Vsw * r0, 2));
		reference_rigity = 4.0; // GV
		set(0);
	}

	~Bfield() {
	}

	void set(double value);
	void set(double r_, double theta_, double phi_, double qf_, double beta_, double momentum_, double charge_);

	inline double getVsw() const {
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
	double beta() const;
	double distance_from_HCS(const double* xx);
	double deriv_distance_from_HCS(const double* xx, unsigned int up);
	double closest_distance();

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
	double delta_parallel_low;
	double delta_parallel_high;
	double gamma;
	double gamma_squared;
	double qf;
	double betavelocity;
	double charge;
	double momentum;
	double reference_rigity;
};

#endif

