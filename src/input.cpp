#include "input.h"

int Input::QueryIntAttribute(std::string obj, TiXmlElement* el) const {
	TiXmlElement* el1 = el->FirstChildElement(obj.c_str());
	int result = -1;
	if (el1)
		el1->QueryIntAttribute("value", &result);
	return result;
}

double Input::QueryDoubleAttribute(std::string obj, TiXmlElement* el) const {
	TiXmlElement* el1 = el->FirstChildElement(obj.c_str());
	double result = -1.0;
	if (el1)
		el1->QueryDoubleAttribute("value", &result);
	return result;
}

int Input::load_file(const std::string inputfilename) {

	TiXmlDocument* doc = new TiXmlDocument(inputfilename.c_str());

	if (doc->LoadFile(inputfilename.c_str(), TIXML_DEFAULT_ENCODING)) {

		TiXmlElement* el = doc->FirstChildElement("HelioProp");
		if (el) {
			int found = inputfilename.rfind('/');
			if (found == std::string::npos)
				found = -1;
			output_filename = "output/";
			output_filename += std::string(inputfilename, found + 1,
					inputfilename.size() - 4 - (found + 1));
			output_filename += ".root";

			Pamelamode = false;
			TiXmlElement* el1 = el->FirstChildElement("Pamela");
			if (el1)
				Pamelamode = true;

			SingleEnergy = false;
			el1 = el->FirstChildElement("SingleEnergy");
			if (el1)
				SingleEnergy = true;

			seed = QueryIntAttribute("Seed", el);

			Rmax = QueryDoubleAttribute("Rmax", el);
			dt = QueryDoubleAttribute("Deltat", el);

			kenergy_min = QueryDoubleAttribute("Emin", el);
			kenergy_max = QueryDoubleAttribute("Emax", el);
			kenergy_size = QueryIntAttribute("NE", el);

			particle_number = QueryIntAttribute("Nparticles", el);
			Znumber = QueryIntAttribute("Charge", el);

			el1 = el->FirstChildElement("Mass");
			if (el1) {
				el1->QueryIntAttribute("value", &Anumber);
				std::string gmod = "Nucleus";
				//el1->QueryStringAttribute("type", &gmod);
				if (gmod == "Nucleus")
					particle_type = NUCLEUS;
				else if (gmod == "Lepton")
					particle_type = LEPTON;
			}

			lambda_par = QueryDoubleAttribute("LambdaPar", el);

			MagField = QueryDoubleAttribute("B0", el);
			delta = QueryDoubleAttribute("Delta", el);

			rig_break = false;
			b = 0.;
			c = 0.;
			el1 = el->FirstChildElement("Break");
			if (el1) {
				rig_break = true;
				b = QueryDoubleAttribute("b", el1);
				c = QueryDoubleAttribute("c", el1);
			}

			Kperp_factor = QueryDoubleAttribute("Kperp", el);

			polarity = QueryDoubleAttribute("Polarity", el);
			tiltangle = QueryDoubleAttribute("TiltAngle", el);

		} else {
			std::cerr << "Problems with file format. No HelioProp tag?" << "\n";
			exit (NOHELIOPROP);
		}
	} else {
		std::cerr << "Problems with the input file " << inputfilename << "\n";
		exit (PROBINPUT);
	}

	delete doc;

	return 0;
}

void Input::print() {

	std::cout << "Settings are " << "\n";
	std::cout << "Output file " << output_filename << "\n";
	std::cout << "NE = " << kenergy_size << "\n";
	std::cout << "Emax = " << kenergy_max << "\n";
	std::cout << "Emin = " << kenergy_min << "\n";
	std::cout << "Rmax = " << Rmax << "\n";
	std::cout << "Delta t = " << dt << "\n";
	std::cout << "How many particles? " << particle_number << "\n";
	std::cout << "Particle type = ";
	if (particle_type == NUCLEUS)
		std::cout << "nucleus" << "\n";
	else if (particle_type == LEPTON)
		std::cout << "lepton" << "\n";
	std::cout << "Charge = " << Znumber << "\n";
	std::cout << "Mass = " << Anumber << "\n";
	std::cout << "Parallel mean free path = " << lambda_par << "\n";
	std::cout << "Delta = " << delta << "\n";
	std::cout << "Perpendicular diffusion = " << Kperp_factor << "\n";
	std::cout << "Polarity = " << polarity << "\n";
	std::cout << "Tilt Angle = " << tiltangle << "\n";
}
