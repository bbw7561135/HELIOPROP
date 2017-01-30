#ifndef HELIOPROP_RANDOM_H_
#define HELIOPROP_RANDOM_H_

#include <gsl/gsl_rng.h>

class TRandomNumberGenerator {
public:
	TRandomNumberGenerator(unsigned long int seed) {
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
		gsl_rng_set (r, seed);

	}
	~TRandomNumberGenerator() {
		gsl_rng_free (r);
	}
	double getRandomNumber() const {
		return gsl_rng_uniform (r);
	}
private:
	const gsl_rng_type * T;
	gsl_rng * r;
};

#endif /* INCLUDE_RANDOM_H_ */
