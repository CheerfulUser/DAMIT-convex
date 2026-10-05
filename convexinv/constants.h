#define POINTS_MAX         3000             /* max number of data points in one lc. */
/* raised from 2000: TESS visits are continuous for days, so a single session can run to
   several thousand points. (3550) Link has 2,061 in one visit and convexinv refused the
   whole object. Only ~0.1% of the catalogue is affected, but those are long, dense,
   high-quality lightcurves -- exactly the ones worth modelling.
   Kept modest: these are STATIC arrays, and 8000/60000 segfaulted on stack overflow.
   Sampled max single visit is 1,955, so 3000 has headroom without risking that. */
#define MAX_N_OBS         20000             /* max number of data points */
/* raised from 10000: 0.2% of sampled objects exceed it across all visits (max seen 12,899). */
#define MAX_LC            10000             /* max number of lightcurves (raised from 100; small static arrays) */
#define MAX_LINE_LENGTH    1000             /* max length of line in the input file */
#define MAX_N_FAC          1000             /* max number of facets */
#define MAX_N_ITER         1000             /* maximum number of iterations */
#define MAX_N_PAR           300             /* maximum number of parameters */
#define MAX_LM               15             /* maximum degree and order of sph. harm. */
#define N_PHOT_PAR            3             /* maximum number of parameters in scattering  law */
#define EPSILON               0             /* precision parameter */
#define TINY                  1e-8          /* precision parameter for mu, mu0*/
#define ALAMDA_START	      0.001
#define ALAMDA_COEFF         10	
#define A_ELL_INIT	      1.05	    /* initial ellipsoid semiaxis 'a' */	
#define B_ELL_INIT	      1		    /* initial ellipsoid semiaxis 'b' */	
#define C_ELL_INIT	      0.95	    /* initial ellipsoid semiaxis 'c' */	

#define PI                    3.14159265358979323846

#define DEG2RAD      (PI / 180)
#define RAD2DEG      (180 / PI)
