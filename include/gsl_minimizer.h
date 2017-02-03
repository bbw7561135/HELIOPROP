#ifndef HELIOPROP_GSLMINIMIZER_H_
#define HELIOPROP_GSLMINIMIZER_H_

#include "gsl/gsl_multimin.h"

double gsl_f(const gsl_vector *v, void *params) {
	double rp, phip;
	Bfield *bfield = (Bfield *) params;

	rp = gsl_vector_get(v, 0);
	phip = gsl_vector_get(v, 1);

	return bfield->distance_from_HCS(rp, phip);
}

void gsl_df(const gsl_vector *v, void *params, gsl_vector *df) {
	double rp, phip;
	Bfield *bfield = (Bfield *) params;

	rp = gsl_vector_get(v, 0);
	phip = gsl_vector_get(v, 1);

	gsl_vector_set(df, 0, bfield->deriv_distance_from_HCS(rp, phip, 0)); // dr
	gsl_vector_set(df, 1, bfield->deriv_distance_from_HCS(rp, phip, 1)); // dphi
}

void gsl_fdf (const gsl_vector *x, void *params,
		double *f, gsl_vector *df)
{
	*f = gsl_f(x, params);
	gsl_df(x, params, df);
}

class GSLMinimizer {
public:
	GSLMinimizer() {
		iter = 0;
		max_iter = 1000;
		status = 0;
		n_params = 2;
		step_size = 0.01;
		tol = 1e-2;
		absolute_tolerance = 1e-3;
		T = gsl_multimin_fdfminimizer_vector_bfgs;
		s = gsl_multimin_fdfminimizer_alloc(T, 2);
	}
	~GSLMinimizer() {
		gsl_multimin_fdfminimizer_free(s);
	}
	double minimize(Bfield* ptr) {
		gsl_vector *x;
		gsl_multimin_function_fdf func;

		func.f = &gsl_f;
		func.df = &gsl_df;
		func.fdf = &gsl_fdf;
		func.n = n_params;
		func.params = ptr;

		/* Starting point, x = (5,7) */
		x = gsl_vector_alloc(2);
		gsl_vector_set(x, 0, 5.0);
		gsl_vector_set(x, 1, 7.0);

		gsl_multimin_fdfminimizer_set(s, &func, x, step_size, tol);

		do {
			iter++;
			status = gsl_multimin_fdfminimizer_iterate(s);
			if (status)
				break;
			status = gsl_multimin_test_gradient(s->gradient,
					absolute_tolerance);
		} while (status == GSL_CONTINUE && iter < max_iter);

		if (status == GSL_SUCCESS)
			std::cout << "Minimum found!" << "\n";

		std::cout << iter << "\t" << gsl_vector_get(s->x, 0) << "\t"
				<< gsl_vector_get(s->x, 1) << "\n";
		std::cout << s->f << "\n";

		gsl_vector_free (x);

		return s->f;
	}
protected:
	double step_size;
	double tol;
	double absolute_tolerance;
	size_t n_params;
	size_t iter;
	size_t max_iter;
	int status;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
};

#ifdef BOH
/*GSLIntegrator::~GSLIntegrator() {
 gsl_integration_workspace_free (w);
 }*/

/*double GSLIntegrator::los(MapBuilder* ptr, const size_t& ipix, const double& maxDistance) {
 params.mapBuilder = ptr;
 params.ipix = ipix;
 gsl_function F;
 F.params = &params;
 F.function = &gslClassWrapper;
 double minDistance = 0.;
 int status = gsl_integration_qag (&F, minDistance, maxDistance, epsabs, epsrel, integration_rule, workspace_size, w, &result, &error);
 return result;
 }*/


/*ROOT::Math::GSLMinimizer min(ROOT::Math::kVectorBFGS2);

 min.SetMaxFunctionCalls(10000);
 min.SetPrintLevel(0);

 ROOT::Math::GradFunctor f(this, &Bfield::distance_from_HCS,
 &Bfield::deriv_distance_from_HCS, 2);
 double step[2] = { 0.01, 0.01 };
 double variable[2] = { r + 1.0, phi };
 min.SetFunction(f);

 // Set the free variables to be minimized!
 min.SetLowerLimitedVariable(0, "rprime", variable[0], step[0], 0);
 min.SetLimitedVariable(1, "phiprime", variable[1], step[1], 0,
 TMath::TwoPi());

 #ifdef _OPENMP
 #pragma omp critical
 #endif
 {
 min.Minimize();
 }

 const double *xs = min.X();

 return std::sqrt(distance_from_HCS(xs));*/

//min.SetMaxIterations(10000);
	size_t max_iter = 10000;
	//min.SetTolerance(0.1);
	double tol = 1e-2;
	double step_size = 0.01;
	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	/* Position of the minimum (1,2), scale factors
	 10,20, height 30. */
	double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };

	gsl_vector *x;
	gsl_multimin_function_fdf F;

	F.n = 2;
	F.f = my_f;
	F.df = my_df;
	F.fdf = my_fdf;
	F.params = par;

	/* Starting point, x = (r + 1.0, phi) */
	x = gsl_vector_alloc(2);
	gsl_vector_set(x, 0, r + 1.0);
	gsl_vector_set(x, 1, phi);

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc(T, 2);

	gsl_multimin_fdfminimizer_set(s, &F, x, step_size, tol);

	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 1e-3);

		if (status == GSL_SUCCESS)
			printf("Minimum found at:\n");

		std::cout << iter << "\t" << gsl_vector_get(s->x, 0) << "\t"
				<< gsl_vector_get(s->x, 1) << "\t" << s->f << "\n";
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);
#endif

#endif /* INCLUDE_GSLMINIMIZER_H_ */
