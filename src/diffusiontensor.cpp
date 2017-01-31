#include "bfield.h"

double TBfield::dPsi_dr() const {
	return 1.0 / (1.0 + gamma_squared) * gamma / r;
}

double TBfield::dPsi_dtheta() const {
	return 1.0 / (1.0 + gamma_squared) * dGamma_dtheta();
}

double TBfield::Kpar() const {

	//double lambda_par = lambda0*mom/(Btot()/B0);//(1.0+r/r0);
	double lambda_par = lambda_0 / (B_total() / B0);  //(1.0+r/r0);

	//if (mom >= 0.1*fabs(charge))
	double rigidity = momentum / fabs(charge);

	//if (rig > 1.)
	lambda_par *= pow(rigidity, delta_Kpar);

	// see Potgieter: 1302.1284
	if (b != 0 && c != 0) {
		double Pk_to_c = pow(reference_rigity / fabs(charge), c);
		lambda_par *= pow((pow(rigidity, c) + Pk_to_c) / (1.0 + Pk_to_c), (b - delta_Kpar) / c);
	}
	return lambda_par * betavelocity * c_light / 3.0;
}

double TBfield::dKpar_dr() const {
	return -Kpar() / B_total() * dBtotal_dr();
}

double TBfield::dKpar_dphi() const {
	return -Kpar() / B_total() * dBtotal_dphi();
}

double TBfield::dKpar_dtheta() const {
	return -Kpar() / B_total() * dBtotal_dtheta();
}

double TBfield::dKrr_dr() const {
	return dKpar_dr() * (pow(cos(psi), 2) + Kperp_factor * pow(sin(psi), 2))
			+ Kpar() * (-1.0 + Kperp_factor) * sin(2.0 * psi) * dPsi_dr();
}

double TBfield::dKrphi_dphi() const {
	return -0.5 * sin(2.0 * psi) * (1.0 - Kperp_factor) * dKpar_dphi();
}

double TBfield::dKrphi_dr() const {
	return -0.5 * sin(2.0 * psi) * (1.0 - Kperp_factor) * dKpar_dr()
			- cos(2.0 * psi) * dPsi_dr() * (1.0 - Kperp_factor) * Kpar();
}

double TBfield::dKthetatheta_dtheta() const {
	return Kperp_factor * dKpar_dtheta();
}  // to be revised

double TBfield::dKphiphi_dphi() const {
	return dKpar_dphi() * (pow(sin(psi), 2) + Kperp_factor * pow(cos(psi), 2));
}

void TBfield::getKTensors(std::vector<double>& K, std::vector<double>& Kderiv) {

	double kappa_par = Kpar();
	double cos_psi = cos(psi);
	double sin_psi = sin(psi);

	double Krr = kappa_par * (pow(cos_psi, 2) + Kperp_factor * pow(sin_psi, 2));
	double Krphi = -kappa_par * (1.0 - Kperp_factor) * sin_psi * cos_psi;
	double Kthetatheta = Kperp_factor * kappa_par;
	double Kphiphi = kappa_par * (Kperp_factor * pow(cos_psi, 2) + pow(sin_psi, 2));

	K[0] = Krr;
	K[1] = Krphi;
	K[2] = Kthetatheta;
	K[3] = Kphiphi;

	Kderiv[0] = 2.0 * r * Krr + r * r * dKrr_dr(); // d/dr (r^2*Krr)
	Kderiv[1] = dKrphi_dphi(); // d/dphi Krphi
	Kderiv[2] = cos(theta) * Kthetatheta + sin(theta) * dKthetatheta_dtheta(); // d/dtheta (sin(theta)*Kthetatheta)
	Kderiv[3] = dKphiphi_dphi(); // d/dphi Kphiphi
	Kderiv[4] = Krphi + r * dKrphi_dr(); // d/dr (r*Krphi)
}


