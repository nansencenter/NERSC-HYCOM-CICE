#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* --------------------------------------------------------------------------
 * Sample usage from Fortran:
 * integer :: seed, n
 * integer(kind=8) state ! or any array >= 8 bytes
 * real(kind=8) :: a, b, x(100)
 *
 * n = 100
 * seed = 1234
 * call gsl_drandinitialize(seed, state)
 * gsl_drandgamma(n, a, b, state, x)
 *
 * NOTE: use environment GSL_RNG_TYPE to select the random number generator,
 *       Mersenne-Twister is the default.
 * --------------------------------------------------------------------------
 */

#if defined(_AIX)
#define GSL_DRANDINITIALIZE gsl_drandinitialize
#define GSL_DRANDGAMMA gsl_drandgamma
#elif defined(__linux) || defined(DOUBLE_UNDERSCORE)
#define GSL_DRANDINITIALIZE gsl_drandinitialize_
#define GSL_DRANDGAMMA gsl_drandgamma_
#else
#define GSL_DRANDINITIALIZE gsl_drandinitialize_
#define GSL_DRANDGAMMA gsl_drandgamma_
#endif


void
GSL_DRANDINITIALIZE(unsigned *seed, void **state)
{
	const gsl_rng_type *T;
	gsl_rng *r;

	
	/*
	 * Choose a generator and a seed as defined by the environment variables
	 * GSL_RNG_TYPE and GSL_RNG_SEED. Defaults to mt19937 and 0.
	 */
	T = gsl_rng_env_setup();

	/* Allocate resources and instantaniate the chosen generator */
	r = gsl_rng_alloc (T);

	/* Seed the newly created generator (overriding GSL_RNG_SEED) */
	gsl_rng_set (r, (unsigned long)*seed);

	*state = r;
}


void
GSL_DRANDGAMMA(unsigned *n, double *a, double *b, void **state, double *x)
{
	unsigned i;
	gsl_rng *r = (gsl_rng *)*state;

	/* Generate N variates from the Gamma distribution */
	for (i = 0; i < *n; i++)
		x[i] = gsl_ran_gamma (r, *a, *b);
}
