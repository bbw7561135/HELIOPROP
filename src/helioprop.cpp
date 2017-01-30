#include <iostream>
#include <fstream>

//#include "TFile.h"
//#include "TNtupleD.h"
//#include "TStopwatch.h"
//#include "TMath.h"
//#include "TRandom3.h"

#include "constants.h"
#include "particles.h"
#include "Bfield.h"
#include "input.h"

//#include "config.h" // deprecated
//#include "RConfigure.h"

int main(int argc, char** argv) {

	if (argc != 2) {
		std::cerr << "Usage: ./HelioMod <xml input file>" << "\n";
		exit(-1);
	}

	//TStopwatch ts;
	//ts.Start();

	Input* inp = new Input();
	inp->LoadFile(argv[1]);
	inp->Print();

	//TFile* outfile = new TFile(inp->outputfilename.c_str(), "RECREATE");

	//TNtupleD* nt = new TNtupleD("nt", "Positions", "pid:q0:rf:thetaf:phif:qf:af:alphaf:tf:iE");
	//TNtupleD* header = new TNtupleD("header", "Init data", "dt:Nparticles:NE:logEmin:logEmax");
	//TRandom3* rn = new TRandom3(inp->seed);
	TRandomNumberGenerator rn(inp->seed);

	std::cout << inp->seed << "\n";

	// r0, Be, Ac, Omega, Vsw, alpha, phi0, lambda0
	double bp[12] = { 1,  // Earth position [UA]
			inp->MagField,  ///sqrt(2.0), // Reference magnetic field [T]
			inp->polarity, // Magnetic field polarity
			TwoPi() / 27.0, // Solar differential rotation rate [d^-1]
			400.0 * kms2UAd, // Velocity of solar wind (constant, radial) [ km/s --> UA/d ]
			inp->tiltangle * DegToRad(), // Tilt angle of current sheet [deg --> rad]
			0, // Reference phi position of current sheet (basically not used, but useful to make the current sheet corotate with Sun)
			inp->lambda_par, // parallel mean-free-path [ UA ]
			inp->Kperp_factor, inp->delta, inp->b, inp->c };

	const double tl = 0;
	TEnvironment* env = new TEnvironment(tl,
			inp->Rmax /* Rmax = Radius of the HelioPause [ UA ] */,
			inp->dt /* Integration backward time step [ d ] */);

	const int Nparticles = inp->Nparticles;
	const double charge = double(inp->Znumber);

	double logEmin;
	double logEmax;
	int NE;

	std::vector<double> enPam;

	if (inp->Pamelamode == false) {
		if (inp->SingleEnergy) {
			NE = 1;
			logEmin = log10(inp->Emin);
			logEmax = log10(inp->Emin);
			enPam.push_back(inp->Emin);
		} else {
			logEmin = log10(inp->Emin);
			logEmax = log10(inp->Emax);
			NE = inp->NE;
			for (int iE = 0; iE < NE; iE++)
				enPam.push_back(
						pow(10,
								logEmin
								+ double(iE) / double(NE - 1)
								* (logEmax - logEmin)));
		}
	} else {
		std::ifstream Pamenergy("energy_Pamela.dat", std::ios::in);
		double Pamq0;
		while (Pamenergy >> Pamq0)
			enPam.push_back(Pamq0);
		Pamenergy.close();

		logEmin = log10(enPam.front());
		logEmax = log10(enPam.back());
		NE = enPam.size();
	}
	int particleID = 0;

	//header->Fill(env->Getdt(), double(Nparticles), double(NE), logEmin, logEmax);
	//header->Write();

#ifdef _OPENMP
#pragma omp parallel num_threads(OMP_NUM_THREADS)
#endif
	{
		for (int iE = 0; iE < NE; iE++) {

#pragma omp critical
			std::cout << "+++ iE = " << iE << " out of " << NE << "\n";
			double q0 = enPam[iE]; // Kinetic energy [GeV]

#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
			for (particleID = 0; particleID < Nparticles; particleID++) {

				//cout << "+++ particle number = " << particleID << " out of " << Nparticles << endl;

				double rf = 1.0; // Starting position (Earth)
				double phif = 0.0;
				double thetaf = PiOver2();
				double af = 0; // Final amplitude
				double alphaf = 1; // Final weight

				TPseudoParticle ps(inp->hd, double(inp->Anumber), rf, phif,
						thetaf, q0, tl, af, alphaf, charge, bp);
				ps.Evolve(env, rn);
#ifdef _OPENMP
#pragma omp critical
#endif
				//ps.Print(nt, double(iE), q0);

			}
		}
	}

	//nt->Write();
	//outfile->Close();
	delete env;
	delete inp;

	//ts.Stop();

	//std::cout << "Ended in " << ts.RealTime() << " s. Actual CPU time "
	//		<< ts.CpuTime() << "\n";

	return 0;
}

