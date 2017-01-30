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
	Input() {
	}
	~Input() {
	}

	int LoadFile(const std::string);
	void Print();

	std::string outputfilename;

	int Anumber;
	int Znumber;
	unsigned int seed;

	double Kperp_factor;
	double lambda_par;
	double delta;
	double MagField;

	double polarity;
	double tiltangle;

	double Emin;
	double Emax;
	int NE;
	int Nparticles;
	PARTICLETYPE hd;

	double Rmax;
	double dt;

	bool Pamelamode;
	bool SingleEnergy;

	bool rig_break;
	double b;
	double c;
};

#endif
