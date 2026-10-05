/* this is a modified veresion of the routine from 
   Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "declarations.h"
#include "constants.h"
/* standard CBLAS dgemm, declared here so no BLAS header is needed (GNU gcc cannot parse Apple's
   Accelerate headers); the same symbol is provided by Accelerate (macOS) and OpenBLAS */
enum { CblasRowMajor = 101 };
enum { CblasNoTrans = 111 };
void cblas_dgemm(int order, int transa, int transb, int m, int n, int k, double alpha,
                 const double *a, int lda, const double *b, int ldb, double beta, double *c, int ldc);
enum { CblasTrans = 112, CblasLower = 122 };
void cblas_dsyrk(int order, int uplo, int trans, int n, int k, double alpha,
                 const double *a, int lda, double beta, double *c, int ldc);
void cblas_dgemv(int order, int trans, int m, int n, double alpha, const double *a, int lda,
                 const double *x, int incx, double beta, double *y, int incy);

double mrqcof(double **x1, double **x2, double x3[], double y[], 
              double sig[], double a[], int ia[], int ma, 
	      double **alpha, double beta[], double (*funcs)())
{
   int mfit,i,j,k,l, np, np1, np2, jp, ic;

   double xx1[4], xx2[4],dy,sig2i, ymod,
	  dave[MAX_N_PAR+1], 
	  coef, ave = 0, trial_chisq;
   /* per-lightcurve scratch, on the heap and sized to the longest lightcurve: as stack arrays of
      POINTS_MAX x MAX_N_PAR they capped lightcurve length and overflowed the stack when raised */
   static double *ytemp = NULL, **dytemp = NULL;
   static int cap_pts = 0, cap_par = 0;
   int maxp = 0;
   /* batched shape derivatives: S (points x facets) holds each point's dA_f s_f, D = S Dg */
   static double *S = NULL, *D = NULL, *sc = NULL;
   static int cap_s = 0;
   int nfit_c, np0, pt;
   /* normal equations per lightcurve: Jw (points x fitted), rows sqrt(w / sig^2) dy/da over the
      fitted parameters only; alpha += Jw^T Jw (dsyrk), beta += Jw^T rw (dgemv) */
   static double *Jw = NULL, *rw = NULL, *Cm = NULL, *bv = NULL;
   static int cap_j = 0, cap_m = 0;
   static int *fidx = NULL;
   double swt;

   for (i = 1; i <= Lcurves; i++)
      if (Lpoints[i] > maxp) maxp = Lpoints[i];
   if ((maxp > cap_pts) || (ma > cap_par))
   {
      if (dytemp != NULL)
      {
         for (i = 0; i <= cap_pts; i++) free(dytemp[i]);
         free(dytemp); free(ytemp);
      }
      free(S); free(D); free(sc); cap_s = 0;
      cap_pts = maxp > cap_pts ? maxp : cap_pts;
      cap_par = ma > cap_par ? ma : cap_par;
      ytemp = (double *) malloc((cap_pts + 1) * sizeof(double));
      dytemp = (double **) malloc((cap_pts + 1) * sizeof(double *));
      for (i = 0; i <= cap_pts; i++)
         dytemp[i] = (double *) malloc((cap_par + 1) * sizeof(double));
   }

   if ((S == NULL) || (cap_s < cap_pts * Numfac))
   {
      free(S); free(D); free(sc);
      cap_s = cap_pts * Numfac;
      S = (double *) malloc((size_t) cap_s * sizeof(double));
      D = (double *) malloc((size_t) cap_pts * Ncoef * sizeof(double));
      sc = (double *) malloc((size_t) (cap_pts + 1) * sizeof(double));
   }
   nfit_c = Ncoef;

   /* N.B. curv and blmatrix called outside bright 
      because output same for all points */
   curv(a);

   blmatrix(a[Ncoef+1],a[Ncoef+2]);

   mfit=0;
   for (j = 1; j <= ma; j++)
      if (ia[j]) mfit++;
   if ((Jw == NULL) || (cap_j < cap_pts * mfit) || (cap_m < mfit))
   {
      free(Jw); free(rw); free(Cm); free(bv); free(fidx);
      cap_j = (cap_pts > 3 ? cap_pts : 3) * mfit;
      cap_m = mfit;
      Jw = (double *) malloc((size_t) cap_j * sizeof(double));
      rw = (double *) malloc((size_t) (cap_pts + 4) * sizeof(double));
      Cm = (double *) malloc((size_t) mfit * mfit * sizeof(double));
      bv = (double *) malloc((size_t) mfit * sizeof(double));
      fidx = (int *) malloc((size_t) (mfit + 1) * sizeof(int));
   }
   j = 0;
   for (l = 1; l <= ma; l++)
      if (ia[l]) fidx[j++] = l;
   for(j = 1; j <= mfit; j++)
   {
      for (k = 1; k <= j; k++)
         alpha[j][k]=0;
      beta[j]=0;
   }
   trial_chisq = 0;
   np = 0;
   np1 = 0;
   np2 = 0;

   for (i = 1; i <= Lcurves; i++)
   {
      if (Inrel[i] == 1) /* is the LC relative? */
      {
         ave = 0;
         for (l = 1; l <= ma; l++)
            dave[l]=0;
      }
      np0 = np;
      if (i < Lcurves)
      {
         /* every point's brightness and non-shape derivatives, with its facet row in S ... */
         for (jp = 1; jp <= Lpoints[i]; jp++)
         {
            pt = np0 + jp;
            for (ic = 1; ic <= 3; ic++) /* position vectors */
            {
               xx1[ic] = x1[pt][ic];
               xx2[ic] = x2[pt][ic];
            }
            ytemp[jp] = bright_s(xx1, xx2, x3[pt], a, dytemp[jp], ma, S + (size_t) (jp - 1) * Numfac);
            sc[jp] = Scale;
         }
         /* ... then all the shape derivatives at once: D = S Dg, row-scaled by each point's Scale */
         cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Lpoints[i], nfit_c, Numfac,
                     1.0, S, Numfac, &Dg[1][1], MAX_N_PAR + 1, 0.0, D, nfit_c);
         for (jp = 1; jp <= Lpoints[i]; jp++)
            for (l = 1; l <= nfit_c; l++)
               dytemp[jp][l] = sc[jp] * D[(size_t) (jp - 1) * nfit_c + l - 1];
      }
      for (jp = 1; jp <= Lpoints[i]; jp++)
      {
         np++;
         if (i == Lcurves)
	    ytemp[jp] = conv(jp,dytemp[jp],ma);
         ymod = ytemp[jp];
	    
         if (Inrel[i] == 1)
            ave = ave + ymod;
     
         if (Inrel[i] == 1)
            for (l = 1; l <= ma; l++)
               dave[l] = dave[l] + dytemp[jp][l];
         /* save lightcurves */
	 
         if (Lastcall == 1) 
	    Yout[np] = ymod;
      } /* jp, lpoints */

   if (Lastcall != 1)
   {
      for (jp = 1; jp <= Lpoints[i]; jp++)
      {
         np1++;
         if (Inrel[i] == 1) 
         {
            coef = sig[np1] * Lpoints[i] / ave;
            for (l = 1; l <= ma; l++)
               dytemp[jp][l] = coef * (dytemp[jp][l] - ytemp[jp] * dave[l] / ave);
            ytemp[jp] = coef * ytemp[jp];
            /* Set the size scale coeff. deriv. explicitly zero for relative lcurves */
            dytemp[jp][1] = 0;
         }
      }

      for (jp = 1; jp <= Lpoints[i]; jp++)
      {
         ymod = ytemp[jp];
         np2++;
         sig2i = Weight[np2] / (sig[np2] * sig[np2]);
         dy = y[np2] - ymod;
         Resid[np2] = dy / sig[np2];
         swt = sqrt(sig2i);
         for (j = 0; j < mfit; j++)
            Jw[(size_t) (jp - 1) * mfit + j] = swt * dytemp[jp][fidx[j]];
         rw[jp - 1] = swt * dy;
         trial_chisq = trial_chisq + dy * dy * sig2i;
      } /* jp */
      cblas_dsyrk(CblasRowMajor, CblasLower, CblasTrans, mfit, Lpoints[i], 1.0, Jw, mfit, 0.0, Cm, mfit);
      cblas_dgemv(CblasRowMajor, CblasTrans, Lpoints[i], mfit, 1.0, Jw, mfit, rw, 1, 0.0, bv, 1);
      for (j = 0; j < mfit; j++)
      {
         for (k = 0; k <= j; k++)
            alpha[j + 1][k + 1] += Cm[(size_t) j * mfit + k];
         beta[j + 1] += bv[j];
      }
     } /* Lastcall != 1 */
         
     if ((Lastcall == 1) && (Inrel[i] == 1))
        Sclnw[i] = Scale * Lpoints[i] * sig[np]/ave;

   } /* i,  lcurves */

   for (j = 2; j <= mfit; j++)
      for (k = 1; k <= j-1; k++)
         alpha[k][j] = alpha[j][k];

   return trial_chisq;
   
}

