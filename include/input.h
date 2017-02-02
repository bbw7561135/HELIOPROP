#ifndef HELIOPROP_INPUT_H
#define HELIOPROP_INPUT_H

#include <iostream>
#include <string>

#include "constants.h"
#include "errorcode.h"
#include "pugi_wrapper.h"

class Input {
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

	double Kperp_factor;
	double lambda_par;
	double delta_low;
	double delta_hi;
	double reference_rigidity;
	double MagField;
	double polarity;
	double tiltangle;
	double Rmax;
	double dt;
	double kenergyn_min;
	double kenergyn_max;
	size_t kenergyn_size;
	size_t particle_number;
	bool is_lepton;
	unsigned int Anumber;
	int Znumber;
	int seed;
};

#endif
