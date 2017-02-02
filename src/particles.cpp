#include "particles.h"

TPseudoParticle::TPseudoParticle(
		const PARTICLETYPE& hd,
		const double& A_,
		const double& rf_,
		const double& phif_,
		const double& thetaf_,
		const double& qf_,
		const double& tf_,
		const double& af_,
		const double& alphaf_,
		const double& charge_,
		double* bparam) :
		Anumber(A_),
		rf(rf_),
		phif(phif_),
		thetaf(thetaf_),
		qf(qf_),
		tf(tf_),
		af(af_),
		alphaf(alphaf_),
		iE(qf_),
		charge(charge_)
{
	if (hd == NUCLEUS) M = 0.938;
	else if (hd == LEPTON) M = 0.000511;
	bf = Bfield();//bparam[0],bparam[1],bparam[2],bparam[3],bparam[4], bparam[5], bparam[6], bparam[7], bparam[8], bparam[9], bparam[10], bparam[11]);
	dW = std::vector<double>(4,0);
	vd = std::vector<double>(3,0);
	Ktensor = std::vector<double>(4,0);
	Kderiv = std::vector<double>(5,0);
}

void TPseudoParticle::Evolve(const TEnvironment* env, const TRandomNumberGenerator& rn) {

	const double Rmax = env->GetRmax();
	const double dt0 = env->Getdt();

	while (rf < Rmax && rf > SunRadius) {
		bf.set(rf, thetaf, phif, qf, beta_velocity(qf), momentum(qf), charge);
		bf.getDriftVelocity(vd);
		bf.getKTensors(Ktensor, Kderiv);
		double Vsw = bf.getVsw();

		double dt = dt0*pow(rf,2)/Ktensor[0];
		double sqrtdt = sqrt(dt);

		tf += dt;

		for (int j = 0; j < 2; j++) {
			double R1 = std::sqrt(-2.0 * log(rn.getRandomNumber()));
			double R2 = TwoPi() * rn.getRandomNumber();
			dW[2*j] = R1 * cos(R2) * sqrtdt;
			dW[2*j+1] = R1 * sin(R2) * sqrtdt;
		}

		double Ar = pow(rf,-2)*Kderiv[0] + 1.0/(rf*sin(thetaf))*Kderiv[1] - Vsw - vd[0];

		double Aphi = pow(rf*sin(thetaf),-2)*Kderiv[3] + 1.0/(pow(rf,2)*sin(thetaf))*Kderiv[4] - vd[1]/(rf*sin(thetaf));

		double Atheta = 1.0/(pow(rf,2)*sin(thetaf))*Kderiv[2] - vd[2]/rf;

		double Gamma = (qf+2.0*M)/(qf+M);
		double AE = pow(rf,-2)/3.0*2.0*(rf*Vsw)*Gamma*qf;

		double secterm = sqrt(2.0/Ktensor[3])*Ktensor[1];
		double Wt = sqrt(2.0*Ktensor[2])/rf;
		double Wphi = sqrt(2.0*Ktensor[3])/(rf*sin(thetaf));

		rf += Ar*dt + sqrt(2.0*Ktensor[0]-pow(secterm,2))*dW[0] + secterm*dW[1];

		thetaf += (Atheta*dt + Wt*dW[2]);

		phif += (Aphi*dt + Wphi*dW[1]);

		qf += AE*dt;

		// Consistency conditions
		thetaf = fmod(thetaf, TwoPi());

		if (thetaf < 0) {
			thetaf *= -1.0;
			phif += M_PI;
		}

		if (thetaf > M_PI) {
			thetaf = TwoPi() - thetaf;
			phif += M_PI;
		}

		if (rf < 0) {
			phif += M_PI;
			thetaf = M_PI - thetaf;
			rf *= -1.0;
		}

		phif = fmod(phif, TwoPi());
		if (phif < 0)
			phif += TwoPi();
	}

	af += boundarycondition(qf,tf,rf);

	return ;

}


