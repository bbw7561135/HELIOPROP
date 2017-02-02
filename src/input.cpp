#include "input.h"

void Input::set_default_values() {
	output_filename = "output/default.helioprop";
	Anumber = 1;
	Znumber = 1;
	is_lepton = false;
	particle_number = 1000;
	kenergyn_min = 0.1;
	kenergyn_max = 1e3;
	kenergyn_size = 40;
	delta_low = 0.34;
	delta_hi = 0.34;
	reference_rigidity = 4;
	Kperp_factor = 1;
	lambda_par = 1;
	MagField = 0;
	polarity = 0;
	tiltangle = 20;
	Rmax = 100;
	dt = 1e-3;
}

void Input::set_output_filename() {
	int found = xml_filename.rfind('/');
	if (found == std::string::npos)
		found = -1;
	output_filename = "output/";
	output_filename += std::string(xml_filename, found + 1,
			xml_filename.size() - 4 - (found + 1));
	output_filename += ".helioprop";
}

void Input::load_file(const std::string& filename) {
	xml_filename = filename;
	set_output_filename();
	try {
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file( xml_filename.c_str() );
		if (!result) throw PROBINPUT;
		pugi::xml_node root = doc.child("Helioprop");
		if (!root) throw NOHELIOPROP;
		load_input(root);
		print();
	} catch (int e) {
		if (e == PROBINPUT)
			std::cerr << "Error reading XML card: " << xml_filename << "\n";
		if (e == NOHELIOPROP)
			std::cerr << "Error reading XML card: Root element Helioprop not found" << "\n";
	}
}

void Input::load_input(pugi::xml_node root) {
	seed = childIntValue(root, "RandomGeneratorSeed", seed);
	Rmax = childDoubleValue(root, "HeliosphereRadius", Rmax);
	dt = childDoubleValue(root, "TimeStep", dt);
	kenergyn_min = childDoubleValue(root, "KEnergyNucMin", kenergyn_min);
	kenergyn_max = childDoubleValue(root, "KEnergyNucMax", kenergyn_max);
	kenergyn_size = childIntValue(root, "KEnergyNucSize", kenergyn_size);
	particle_number = childIntValue(root, "ParticleNumber", 1000);
	Znumber = childIntValue(root, "ParticleCharge", 1);
	Anumber = childIntValue(root, "ParticleAtomicNumber", 1);
	is_lepton = (Anumber == 0);
	lambda_par = childDoubleValue(root, "ParallelMfp", lambda_par);
	MagField = childDoubleValue(root, "BfieldAtSun", MagField);
	delta_low = childDoubleValue(root, "Delta", delta_low);
	delta_low = childDoubleValue(root, "DeltaLow", delta_low);
	delta_hi = childDoubleValue(root, "Delta", delta_hi);
	delta_hi = childDoubleValue(root, "DeltaHigh", delta_hi);
	reference_rigidity = childDoubleValue(root, "ReferenceRigidity", reference_rigidity);
	Kperp_factor = childDoubleValue(root, "PerpendicularDiffusion", Kperp_factor);
	polarity = childDoubleValue(root, "Polarity", polarity);
	tiltangle = childDoubleValue(root, "TiltAngle", tiltangle);
}

void Input::print() {
	std::cout << "Setting list" << "\n";
	std::cout << "- Output file : " << output_filename << "\n";
	std::cout << "- RNG seed : " << seed << "\n";
	std::cout << "- time step = " << dt << "\n";
	std::cout << "- heliosphere radius = " << Rmax << "\n";
	std::cout << "- KEnergy size = " << kenergyn_size << "\n";
	std::cout << "- KEnergy max  = " << kenergyn_max << "\n";
	std::cout << "- KEnergy min  = " << kenergyn_min << "\n";
	std::cout << "- particle number = " << particle_number << "\n";
	std::cout << "- charge = " << Znumber << "\n";
	std::cout << "- atomic number = " << Anumber << "\n";
	if (is_lepton)
		std::cout << "- type = lepton\n";
	else
		std::cout << "- type = nucleus\n";
	std::cout << "- parallel mean free path = " << lambda_par << "\n";
	std::cout << "- B field at Sun = " << MagField << "\n";
	std::cout << "- delta = " << delta_low << " - " << delta_hi << "\n";
	std::cout << "- ref rigidity = " << reference_rigidity << "\n";
	std::cout << "- perpendicular diffusion = " << Kperp_factor << "\n";
	std::cout << "- polarity = " << polarity << "\n";
	std::cout << "- tilt angle = " << tiltangle << "\n";
	std::cout << "\n";
}
