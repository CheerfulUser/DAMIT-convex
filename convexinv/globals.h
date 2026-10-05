#include <stdio.h>

#include "constants.h"

extern int Lmax, Mmax, Niter, Lastcall,
           Ncoef, Numfac, Lcurves, Nphpar,
           Lpoints[MAX_LC+1], Inrel[MAX_LC+1],
	   Deallocate, Nyorp;   /* Nyorp = 1: YORP term is parameter Ncoef+Nphpar+6 */
    
extern double Yorp_now,   /* current YORP d(omega)/dt (rad/day^2) used by matrix() */
              Ochisq, Chisq, Alamda, Alamda_incr, Alamda_start, Phi_0, Scale,
              Area[MAX_N_FAC+1], Darea[MAX_N_FAC+1], Sclnw[MAX_LC+1], 
	      *Yout, *Weight, *Resid,   /* Resid: weighted-fit residual / sigma per point (robust mode) */ /* heap arrays, one per data point (+3 convexity rows); Weight = chi^2 weight */
              Fc[MAX_N_FAC+1][MAX_LM+1], Fs[MAX_N_FAC+1][MAX_LM+1], 
	      Tc[MAX_N_FAC+1][MAX_LM+1], Ts[MAX_N_FAC+1][MAX_LM+1], 
	      Dsph[MAX_N_FAC+1][MAX_N_PAR+1], Dg[MAX_N_FAC+1][MAX_N_PAR+1],   
              Nor[MAX_N_FAC+1][4], Blmat[4][4],
              Pleg[MAX_N_FAC+1][MAX_LM+1][MAX_LM+1],
              Dblm[3][4][4];
    
	   

