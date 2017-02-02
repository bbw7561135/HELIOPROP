#ifndef HELIOPROP_TENVIRONMENT_H_
#define HELIOPROP_TENVIRONMENT_H_

class Environment {
public:
	Environment() :
		t_end(0), r_max(0), dt(0) {
	}
	Environment(const double& tend_, const double& rmax_, const double& dt_) :
		t_end(tend_), r_max(rmax_), dt(dt_) {
	}
	~Environment() {
	}
	inline double get_tend() const {
		return t_end;
	}
	inline double get_rmax() const {
		return r_max;
	}
	inline double get_dt() const {
		return dt;
	}
protected:
	double t_end;
	double dt;
	double r_max;
};

#endif /* HELIOPROP_TENVIRONMENT_H_ */
