#ifndef HELIOPROP_UTILITIES_H_
#define HELIOPROP_UTILITIES_H_

inline double Heaviside(const double& x) {
	if (x > 0)
		return -1;
	if (x < 0)
		return 1;
	return 1e-30;
}

inline double sgn(const double& x) {
	return (x > 0) - (x < 0);
}

inline double ASin(const double& x) {
	if (x < -1.)
		return -M_PI / 2.;
	if (x >  1.)
		return  M_PI / 2.;
	return asin(x);
}

#endif /* INCLUDE_UTILITIES_H_ */
