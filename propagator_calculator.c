#include <stdio.h>
#include "action.c"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#define NP 20

void integral_vegas(size_t, double*, double*, double*, double*, double*);
void defineX_0(double, double, double*);

void integral_vegas(size_t dim, double x_0s[], double lower[], double upper[], double vals[], double error[])
{
	FILE *fp = fopen("1D_x^4_potential_propagator.txt", "w");
	
	double res, err, x0;
	double A = pow((0.5 * MASS * M_1_PI)/LATTICE_SPACING, 0.5*N);
	double standard_qm_val = 0.;
	size_t calls = 1.e6;
	size_t warmCalls = 1.e5; //calls for warm up
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus2);
	unsigned long seed = 1UL;
	gsl_rng_set (r, seed);
	for(int i = 0; i < NP; i++)
	{
		x0 = x_0s[i];
		standard_qm_val = exp(-x0 * x0 - TIME * 0.5)/M_SQRTPI;
		gsl_monte_function G = {&integrand, dim, &x0};
		gsl_monte_vegas_state *sv = gsl_monte_vegas_alloc(dim); //allocates space for vegas algorithm structure
		gsl_monte_vegas_init (sv);
		gsl_monte_vegas_integrate (&G, lower, upper, dim, warmCalls, r, sv, &res, &err); //initial warm up
		do 
		{
			gsl_monte_vegas_integrate (&G, lower, upper, dim, calls, r, sv, &res, &err);
			fflush(stdout);
		}
		while(fabs (gsl_monte_vegas_chisq (sv) - 1.0) > 0.2); //loop until vegas chi^2 is close to 1.
		
		vals[i] = res;
		error[i] = err;
		/* FOR HARMONIC OSCILLATOR POTENTIAL
		printf("%f %f %f %f\n", x0, A*res, standard_qm_val, err);
		fprintf(fp, "%f %f %f %f\n", x0, A*res, standard_qm_val, err);
		*/
		//FOR x^4 POTENTIAL
		printf("%f %f %f\n", x0, A*res, err);
		fprintf(fp, "%f %f %f\n", x0, A*res, err);
		
		gsl_monte_vegas_free(sv); //frees memory of VEGAS algorithm structure
	}
		
	
	gsl_rng_free(r); //frees space of random number generator
	fclose(fp);
}
	

void defineX_0(double lower, double upper, double x_0[])
{
	double x0 = 0.; //mimics parameter r in integrand
	double xStep = (upper-lower)/NP;
	for(int i = 0; i < NP; i++)
	{
			x0 = lower + xStep*i; //increases dist by diststep
			x_0[i] = x0;
	}
}


int main(void)
{	
	
	//5 ~= infinity since the integrand is somewhat of a gaussian
	//bounds are for paths x_1 -> x_(N-1)
	size_t dim = N - 1;
	double lower_bounds[] = {-5., -5., -5., -5., -5., -5., -5};
	double upper_bounds[] = {5., 5., 5., 5., 5., 5., 5};
	double x_0s[NP];
	double xmin = 0., xmax = 2.1;
	double values[NP], errors[NP];
	defineX_0(xmin, xmax, x_0s);
	integral_vegas(dim, x_0s, lower_bounds, upper_bounds, values, errors);
	
	return 0;
}
