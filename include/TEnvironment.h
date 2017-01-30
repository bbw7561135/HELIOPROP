#ifndef HELIOPROP_TENVIRONMENT_H_
#define HELIOPROP_TENVIRONMENT_H_

class TEnvironment {
public:
	TEnvironment() :
		Tend(0), Rmax(0), dt(0) {
	}
	TEnvironment(double Tend_, double Rmax_, double dt_) :
		Tend(Tend_), Rmax(Rmax_), dt(dt_) {
	}
	~TEnvironment() {
	}
	inline double GetTend() const {
		return Tend;
	}
	inline double GetRmax() const {
		return Rmax;
	}
	inline double Getdt() const {
		return dt;
	}
protected:
	double Tend;
	double dt;
	double Rmax;
};

#endif /* HELIOPROP_TENVIRONMENT_H_ */
