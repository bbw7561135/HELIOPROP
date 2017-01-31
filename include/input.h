#ifndef HELIOPROP_INPUT_H
#define HELIOPROP_INPUT_H

#include <iostream>
#include <string>

#include "constants.h"
#include "errorcode.h"
#include "tinystr.h"
#include "tinyxml.h"

class Input {

private:
	int QueryIntAttribute(std::string, TiXmlElement*) const;
	double QueryDoubleAttribute(std::string, TiXmlElement*) const;

public:
	Input(const std::string filename) {
		load_file(filename);
		print();
	}
	~Input() {
	}

	int load_file(const std::string);
	void print();

	std::string output_filename;

	int Anumber;
	int Znumber;
	unsigned int seed;

	double Kperp_factor;
	double lambda_par;
	double delta;
	double MagField;

	double polarity;
	double tiltangle;

	double kenergy_min;
	double kenergy_max;
	size_t kenergy_size;

	size_t particle_number;
	PARTICLETYPE particle_type;

	double Rmax;
	double dt;

	bool Pamelamode;
	bool SingleEnergy;

	bool rig_break;
	double b;
	double c;
};

#endif
