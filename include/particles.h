#ifndef HELIOPROP_PARTICLES_H
#define HELIOPROP_PARTICLES_H

#include <cmath>
#include <vector>

//#include "TNtupleD.h"
//#include "TRandom3.h"

#include "constants.h"
#include "bfield.h"
#include "random.h"

class TPseudoParticle {
public:
	TPseudoParticle(const PARTICLETYPE& hd, const double& A_, const double& rf_,
			const double& phif_, const double& thetaf_, const double& qf_,
			const double& tf_, const double& af_, const double& alphaf_,
			const double& charge_, double* bparam);

	~TPseudoParticle() {
	}

	void Evolve(const TEnvironment* env, const TRandomNumberGenerator& rn);

	inline double momentum(const double& qf) {
		return Anumber * sqrt(pow(qf + M, 2) - pow(M, 2));
	}
	inline double beta_velocity(const double& q) {
		return sqrt(1.0 - 1.0 / pow(1.0 + q / M, 2));
	}

	inline double GetAlpha() const {
		return alphaf;
	}
	inline double GetA() const {
		return af;
	}

	inline double boundarycondition(const double& qf, const double& tf,
			const double& rf) {
		return (rf < SunRadius || qf < MinEn) ? 0 : 1; ///momentum(qf)*pow(pow(M,2)+pow(momentum(qf),2),-1.8);//12.14*beta_velocity(qf)*pow(qf+0.5*M,-2.6); // Boundary condition. Put only 0 or 1 to allow reweighting // 12.14*beta_velocity(qf)*pow(qf+0.5*M,-2.6);
	}

	//inline void Print(TNtupleD* nt, const double& pid, const double& q0) {
	//	nt->Fill(pid, q0, rf, thetaf, phif, qf, af, alphaf, tf, iE);
	//	return;
	//}

protected:
	double Anumber;
	double M;
	double rf;
	double phif;
	double thetaf;
	double qf;
	double tf;
	double af;
	double alphaf;
	double charge;
	double iE;
	std::vector<double> dW;
	std::vector<double> vd;
	std::vector<double> Ktensor;
	std::vector<double> Kderiv;
	Bfield bf;
};

#endif
