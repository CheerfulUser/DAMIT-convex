/* this computes brightness and its derivatives w.r.t. parameters */

/*
Copyright (C) 2006  Mikko Kaasalainen, Josef Durech

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
#include "declarations.h"
#include "constants.h"

/* srow != NULL (batched mrqcof): instead of the shape derivatives dyda[1..Ncoef], write this point's
   facet terms srow[f-1] = dA_f s_f (0 for facets not both lit and visible); the caller forms all
   points' shape derivatives as one matrix product srow x Dg and scales each row by Scale */
double bright_s(double ee[], double ee0[], double t, double cg[],
                double dyda[], int ncoef, double *srow)
{
   int ncoef0, i, j, k, nc, icl, ils, iy,
       incl[MAX_N_FAC+1];
 
   double cos_alpha, br, cl, cls, alpha, sum, dnom, wj, *dg, a1[4], a2[4], s_cl, s_ls,
          e[4], e0[4],
          php[N_PHOT_PAR+1], dphp[N_PHOT_PAR+1],
	  mu[MAX_N_FAC+1], mu0[MAX_N_FAC+1+1], s[MAX_N_FAC+1], dbr[MAX_N_FAC+1],
          dsmu[MAX_N_FAC+1], dsmu0[MAX_N_FAC+1],
	  de[4][4], de0[4][4], tmat[4][4],
	  dtm[4][4][4];
   
   /* parameter positions from the number of shape coeffs., not from ncoef (= all parameters),
      so an optional YORP term can follow Lommel-Seeliger without moving anything */
   ncoef0 = Ncoef + 3;                 /* omega */
   icl = Ncoef + Nphpar + 4;           /* log Lambert */
   ils = Ncoef + Nphpar + 5;           /* Lommel-Seeliger */
   iy = Ncoef + Nphpar + 6;            /* YORP d(omega)/dt, rad/day^2, when Nyorp = 1 */
   cl = exp(cg[icl]); /* Lambert */
   cls = cg[ils];       /* Lommel-Seeliger */
   Yorp_now = Nyorp ? cg[iy] : 0;
   cos_alpha = dot_product(ee, ee0);
   alpha = acos(cos_alpha);
   for (i = 1; i <= Nphpar; i++)
      php[i] = cg[ncoef0+i];

   phasec(dphp,alpha,php); /* computes also Scale */

   matrix(cg[ncoef0],t,tmat,dtm);

   br = 0;
   /* Directions (and ders.) in the rotating system */
   for (i = 1; i <= 3; i++)
   {
      e[i] = 0;
      e0[i] = 0;
      for (j = 1; j <= 3; j++)
      {
         e[i] = e[i] + tmat[i][j] * ee[j];
         e0[i] = e0[i] + tmat[i][j] * ee0[j];
         de[i][j] = 0;
         de0[i][j] = 0;
         for (k = 1; k <= 3; k++)
	 {
            de[i][j] = de[i][j] + dtm[j][i][k] * ee[k];
            de0[i][j] = de0[i][j] + dtm[j][i][k] * ee0[k];
         }
      }
   } 

   for (j = 1; j <= 3; j++)
      a1[j] = a2[j] = 0;
   s_cl = s_ls = 0;
   /*Integrated brightness (phase coeff. used later) */
   for (i = 1; i <= Numfac; i++)
   {
      incl[i] = 0;
      mu[i] = e[1] * Nor[i][1] + e[2] * Nor[i][2] + e[3] * Nor[i][3];
      mu0[i] = e0[1] * Nor[i][1] + e0[2] * Nor[i][2] + e0[3] * Nor[i][3];
      if((mu[i] > TINY) && (mu0[i] > TINY)) 
      {
         double q0, q1, am;
         incl[i] = 1;
         dnom = mu[i] + mu0[i];
         s[i] = mu[i] * mu0[i] * (cl + cls / dnom);
         br = br + Area[i] * s[i];
         q0 = mu0[i] / dnom;
         q1 = mu[i] / dnom;
         dsmu[i] = cls * q0 * q0 + cl * mu0[i];
         dsmu0[i] = cls * q1 * q1 + cl * mu[i];
         dbr[i] = Darea[i] * s[i];
         /* rotation and scattering-law sums in the same pass (were two more facet loops) */
         for (j = 1; j <= 3; j++)
         {
            a1[j] += Area[i] * dsmu[i] * Nor[i][j];
            a2[j] += Area[i] * dsmu0[i] * Nor[i][j];
         }
         am = mu[i] * mu0[i] * Area[i];
         s_cl += am;
         s_ls += Area[i] * mu[i] * mu0[i] / dnom;
      }
    }

   /* Derivatives of brightness w.r.t. g-coeffs */
   /* facet-outer, coefficient-inner: Dg[j][.] is contiguous, so the inner loop vectorises and
      the visibility test runs once per facet rather than once per facet per coefficient
      (92% of the run time was here; same sums, different order) */
   nc = ncoef0 - 3;
   if (srow != NULL)
   {
      for (j = 1; j <= Numfac; j++)
         srow[j - 1] = incl[j] ? dbr[j] : 0;
   }
   else
   {
   for (i = 1; i <= nc; i++)
      dyda[i] = 0;
   for (j = 1; j <= Numfac; j++)
      if (incl[j] == 1)
      {
         dg = Dg[j];
         wj = dbr[j];
         for (i = 1; i <= nc; i++)
            dyda[i] += wj * dg[i];
      }
   for (i = 1; i <= nc; i++)
      dyda[i] = Scale * dyda[i];
   }
/*   printf("%f \n", dyda[1]);   */
   /* Ders. of brightness w.r.t. rotation parameters */
   /* sum_i Area_i (dsmu_i N_i.de_k + dsmu0_i N_i.de0_k) = sum_j de_jk a1_j + de0_jk a2_j with
      a1 = sum_i Area_i dsmu_i N_i, a2 = sum_i Area_i dsmu0_i N_i: one facet pass, not three */
   for (k = 1; k <= 3; k++)
   {
      sum = 0;
      for (j = 1; j <= 3; j++)
         sum += de[j][k] * a1[j] + de0[j][k] * a2[j];
      dyda[ncoef0-3+k] = Scale * sum;
   }
   /* phi = omega t + Phi_0 + Yorp t^2 / 2 enters only through phi, so d/dYorp = (t / 2) d/domega */
   if (Nyorp)
      dyda[iy] = 0.5 * t * dyda[ncoef0];
   
   /* Ders. of br. w.r.t. phase function params. */
   for(i = 1; i <= Nphpar; i++)
      dyda[ncoef0+i] = br * dphp[i];

   /* Ders. of br. w.r.t. cl, cls */
   dyda[icl] = Scale * s_cl * cl;
   dyda[ils] = Scale * s_ls;

   /* Scaled brightness */
   br = br * Scale;
   
   return(br);
}

double bright(double ee[], double ee0[], double t, double cg[],
              double dyda[], int ncoef)
{
   return bright_s(ee, ee0, t, cg, dyda, ncoef, NULL);
}
