#ifndef HELIOPROP_INPUT_H
#define HELIOPROP_INPUT_H

#include <iostream>
#include <string>

#include "constants.h"
#include "errorcode.h"
#include "tinystr.h"
#include "tinyxml.h"
#include "pugi_wrapper.h"

class Input {

private:
	int QueryIntAttribute(std::string, TiXmlElement*) const;
	double QueryDoubleAttribute(std::string, TiXmlElement*) const;

public:
	Input() {
		set_default_values();
	}
	~Input() {
	}

	void set_default_values();
	void load_file(const std::string& filename);
	void load_input(pugi::xml_node root);
	void set_output_filename();
	void print();

	std::string xml_filename;
	std::string output_filename;

	int seed;

	double Kperp_factor;
	double lambda_par;
	double delta;
	double MagField;
	double polarity;
	double tiltangle;
	double Rmax;
	double dt;

	double kenergyn_min;
	double kenergyn_max;
	size_t kenergyn_size;

	size_t particle_number;
	PARTICLETYPE particle_type;
	unsigned int Anumber;
	int Znumber;

	bool rig_break;
	double b;
	double c;
};

#endif
