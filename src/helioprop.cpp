#include <iostream>
#include <fstream>

#include "constants.h"
#include "particles.h"
#include "Bfield.h"
#include "energyaxis.h"
#include "input.h"
#include "Timer.h"

int main(int argc, char** argv) {
	Timer tmr;
	Input input;
	if (argc == 2) {
		input.load_file(argv[1]);
		TRandomNumberGenerator rn(input.seed);
		Bfield bfield(input);
		EnergyAxis energy(input);
		Environment env(0, input.Rmax, input.dt);

		std::cout << "Starting computation ---" << "\n";
#ifdef _OPENMP
#pragma omp parallel num_threads(OMP_NUM_THREADS)
#endif
		for (size_t iE = 0; iE < energy.get_size(); iE++) {
#ifdef _OPENMP
#pragma omp critical
#endif
			std::cout << "energy index " << iE << " out of " << energy.get_size() << "\n";
			double q0 = energy.at(iE); // Kinetic energy [GeV]
#ifdef _OPENMP
#pragma omp for schedule(dynamic) nowait
#endif
			for (size_t id = 0; id < input.particle_number; ++id) {
				double rf = 1.0; // Starting position (Earth)
				double phif = 0.0;
				double thetaf = PiOver2;
				double af = 0; // Final amplitude
				double alphaf = 1; // Final weight
				//TPseudoParticle ps(input->particle_type, double(input->Anumber), rf, phif,
				//		thetaf, q0, tl, af, alphaf, charge, bp);
				//	ps.Evolve(env, rn);
#ifdef _OPENMP
#pragma omp critical
#endif
				//ps.Print(nt, double(iE), q0);

			}
		}
	}
	else {
		std::cerr << "Usage: ./HELIOPROP <xml input file>" << "\n";
	}
	std::cout << "Helioprop ends in " << tmr.elapsed() << " seconds.\n";


	//TFile* outfile = new TFile(inp->outputfilename.c_str(), "RECREATE");
	//TNtupleD* nt = new TNtupleD("nt", "Positions", "pid:q0:rf:thetaf:phif:qf:af:alphaf:tf:iE");
	//TNtupleD* header = new TNtupleD("header", "Init data", "dt:Nparticles:NE:logEmin:logEmax");

#ifdef HELIOPROP

	const int Nparticles = input->particle_number;
	const double charge = double(input->Znumber);

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

				TPseudoParticle ps(input->particle_type, double(input->Anumber), rf, phif,
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
#endif
	return 0;
}

