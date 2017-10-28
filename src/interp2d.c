#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "prototype.h"

/***************************************************************
 * 2d spline interpolation via GSL
 ***************************************************************/
void interp2d(double *x1a, double *x2a, double **ya, int m, int n, double x1, double x2, double *y){
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();

    gsl_spline *splinerow = gsl_spline_alloc (gsl_interp_cspline, n);

    gsl_spline *splinecol = gsl_spline_alloc (gsl_interp_cspline, m);

                   	
    int j;
    double *ymtmp;
 
    ymtmp= (double*) (malloc(sizeof(double) * m));
    
    for (j=0;j<m;j++) {
        gsl_spline_init (splinerow, x2a, ya[j], n);
        ymtmp[j] = gsl_spline_eval (splinerow, x2, acc);
    }
  
    gsl_spline_init (splinecol, x1a, ymtmp, m);
 
    *y = gsl_spline_eval (splinecol, x1, acc);

    free(ymtmp);
    gsl_spline_free (splinerow);
    gsl_spline_free (splinecol);
    gsl_interp_accel_free (acc);

}
