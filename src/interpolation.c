/***************************************************************
 * Uses the pre-computed tables and sets interpolated 
 * File interpolation.c returns from the generated tables the required 
values of the integral.
 ***************************************************************/

/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/***************************************************************
 * GSL header files
 ***************************************************************/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "constants.h"
#include "table_settings.h"
#include "prototype.h"
#include "allvars.h"

/***************************************************************
 * Cached interpolation objects.
 *
 * interpolation() runs once per grid-point per timestep and, via its 14
 * interpolation_2D() calls, used to allocate/free ~40 gsl_spline +
 * gsl_interp_accel objects on every call. Since all splines are size
 * INTERPOINTS (and the 2D helper's row/col splines are sized to the (m,n)
 * passed in), we build them once, re-init the data on each call, and free
 * only at shutdown (interpolation_free, from memory_free_all).
 *
 * Reusing a gsl_interp_accel across different data is safe and gives
 * bit-identical results: the accelerator is only a search hint and always
 * returns the correct bracketing index.
 ***************************************************************/
static gsl_interp_accel *ip_acc    = NULL;
static gsl_spline       *ip_spline = NULL;   /* size INTERPOINTS */

static gsl_interp_accel *ip2d_acc   = NULL;
static gsl_spline       *ip2d_row   = NULL;  static int ip2d_n     = 0;
static gsl_spline       *ip2d_col   = NULL;  static int ip2d_m     = 0;
static double           *ip2d_ymtmp = NULL;  static int ip2d_mcap  = 0;

void interpolation_free(void){
    if(ip_spline){   gsl_spline_free(ip_spline);        ip_spline = NULL; }
    if(ip_acc){      gsl_interp_accel_free(ip_acc);     ip_acc    = NULL; }
    if(ip2d_row){    gsl_spline_free(ip2d_row);         ip2d_row  = NULL; ip2d_n = 0; }
    if(ip2d_col){    gsl_spline_free(ip2d_col);         ip2d_col  = NULL; ip2d_m = 0; }
    if(ip2d_acc){    gsl_interp_accel_free(ip2d_acc);   ip2d_acc  = NULL; }
    if(ip2d_ymtmp){  free(ip2d_ymtmp);                  ip2d_ymtmp= NULL; ip2d_mcap = 0; }
}

void interpolation(double nHX1, double nHeX2, double nHX13, double radius){


    int i,j;
    double y1, y2, y3;


    if (!ip_spline){
        ip_acc    = gsl_interp_accel_alloc();
        ip_spline = gsl_spline_alloc(gsl_interp_cspline, INTERPOINTS);
    }
    gsl_interp_accel *acc    = ip_acc;
    gsl_spline       *spline = ip_spline;

    
    double radNorm = 4 * M_PI * radius * radius * KPC * KPC; // normalization due to 4*pi*r^2 [r in cm]         


    double nHX1_val, nHeX2_val, nHX13_val;     /* table value wanted for */ 

    /* Picking up integral values to be used in the ODE */
    int back_step = (int) ((INTERPOINTS - 1) / 2);
  
    nHX1_val = ((int) (nHX1 * (1 / TABLERES))) / (1 / TABLERES) - back_step * TABLERES;    
    
    for (i = 0; i < INTERPOINTS; i++)
        nHx1a[i] = xa[i] = nHX1_val + TABLERES * i;
  
    nHeX2_val = ((int) (nHeX2 * (1 / TABLERES))) / (1 / TABLERES) - back_step * TABLERES;
    
    for (i = 0; i < INTERPOINTS; i++)
        nHex2a[i] = nHeX2_val + TABLERES * i;

    nHX13_val = ((int) (nHX13 * (1 / TABLERES))) / (1 / TABLERES) - back_step * TABLERES;
    
    for (i = 0; i < INTERPOINTS; i++)
        nHx13a[i] = nHX13_val + TABLERES * i;
  
    /***************************************************************
     * The above three statments calculates the array which is 
     * centred on the value for which you need the table value; 
     * In this case you have 2*TABLERES because Interpoints = 5. 
     * And thus you need 2 points on either side of the value 
     * for which you need the answer 
     ***************************************************************/
  

    /***************************************************************
     * First Set - Ionization
     ***************************************************************/
    
    /* First Intergral */
    for (i = 0; i < INTERPOINTS; i++)
        ya[i] = gsl_vector_get( e1h1_p1, (long int) ((xa[i] - LOWLIM) / TABLERES) + 2);
  
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP2a[i][j] = gsl_matrix_get(  matrix_e1h1_p2, 
                                          (long int) ((nHx1a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM)/ TABLERES) + 2);
  
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get( matrix_e1h1_p3, 
                                         (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                         (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
    gsl_spline_init (spline, xa, ya, INTERPOINTS);
    y1 = gsl_spline_eval (spline, nHX1, acc);
  
    interpolation_2D(nHx1a, nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1, nHeX2, &y2);  
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
  
  
    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
    
    fuku_e1h1[iGrid] =  (y1 + y2 + y3) / (radNorm);
      
     
    /* Second Integral */
     for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP2a[i][j] = gsl_matrix_get(  matrix_ehe1_p2,
                                          (long int) ((nHx1a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM)/ TABLERES) + 2);
  
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  matrix_ehe1_p3, 
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
    interpolation_2D(nHx1a, nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1, nHeX2, &y2);
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
  
    y2 = fabs (y2);
    y3 = fabs (y3);

    fuku_ehe1[iGrid] = (y2 + y3) / (radNorm);
  
    /* Third Integral */
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  matrix_ehe2_p3, 
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
  
    y3 = fabs (y3);

    fuku_ehe2[iGrid] = y3 / (radNorm);

  
    /***************************************************************
     * Second set - Temperatures
     ***************************************************************/

    /* First Intergral */
    for (i = 0; i < INTERPOINTS; i++)
        ya[i] = gsl_vector_get(  temp_e1h1_p1, (long int) ((xa[i] - LOWLIM) / TABLERES) + 2);

    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP2a[i][j] = gsl_matrix_get(  temp_matrix_e1h1_p2, 
                                          (long int) ((nHx1a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  temp_matrix_e1h1_p3,
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
  
    gsl_spline_init (spline, xa, ya, INTERPOINTS);
    y1 = gsl_spline_eval (spline, nHX1, acc);
    
    interpolation_2D(nHx1a,  nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1,  nHeX2, &y2);
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
    
    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);

    integral_H1[iGrid] = (y1 + y2 + y3) / ( radNorm);
  
    /* Second Integral */
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP2a[i][j] = gsl_matrix_get(  temp_matrix_ehe1_p2,
                                          (long int) ((nHx1a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);  
  
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  temp_matrix_ehe1_p3,
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);

    interpolation_2D(nHx1a,  nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1,  nHeX2, &y2);
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);  
  
    y2 = fabs (y2);
    y3 = fabs (y3);

    integral_He1[iGrid] = (y2 + y3) / ( radNorm );
  
    
    /* Third Integral */
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  temp_matrix_ehe2_p3,
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);  
  
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);  
  
    y3 = fabs (y3);

    integral_He2[iGrid] = (y3) / ( radNorm );
  

    /***************************************************************
     * Third Set - Compton
     ***************************************************************/  
    
    /* first integral */  
    for (i = 0; i < INTERPOINTS; i++)
        ya[i] = gsl_vector_get( comp1_p1, (long int) ((xa[i] - LOWLIM) / TABLERES) + 2 );
  
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP2a[i][j] = gsl_matrix_get(  comp1_matrix_p2,
                                          (long int) ((nHx1a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
 
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get( comp1_matrix_p3,
                                         (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                         (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
    gsl_spline_init (spline, xa, ya, INTERPOINTS);
    y1 = gsl_spline_eval (spline, nHX1, acc);
  
    interpolation_2D(nHx1a, nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1, nHeX2, &y2);
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);

    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
  
    comp_integ1[iGrid] = (y1 + y2 + y3) / ( radNorm );

    /* second integral */
    for (i = 0; i < INTERPOINTS; i++)
        ya[i] = gsl_vector_get( comp2_p1, (long int) ((xa[i] - LOWLIM) / TABLERES) + 2);
 
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP2a[i][j] = gsl_matrix_get( comp2_matrix_p2,
                                         (long int) ((nHx1a[i] - LOWLIM) / TABLERES) + 2,
                                         (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
 
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get( comp2_matrix_p3,
                                         (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                         (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
 
    gsl_spline_init (spline, xa, ya, INTERPOINTS);
    y1 = gsl_spline_eval (spline, nHX1, acc);
    
    interpolation_2D(nHx1a,  nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1,  nHeX2, &y2);
    interpolation_2D(nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);

    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
    

    comp_integ2[iGrid] = (y1 + y2 + y3) / ( radNorm);

    /* cached objects (ip_spline / ip_acc) are freed at shutdown, not here */
}

/***************************************************************
 * 2d spline interpolation via GSL
 ***************************************************************/
void interpolation_2D(double *x1a, double *x2a, double **ya, int m, int n, double x1, double x2, double *y){

    /* lazily (re)allocate the cached splines; sizes are constant in practice
     * (INTERPOINTS), so this reallocates at most once. */
    if (!ip2d_acc) ip2d_acc = gsl_interp_accel_alloc();
    if (ip2d_n != n){ if(ip2d_row) gsl_spline_free(ip2d_row); ip2d_row = gsl_spline_alloc(gsl_interp_cspline, n); ip2d_n = n; }
    if (ip2d_m != m){ if(ip2d_col) gsl_spline_free(ip2d_col); ip2d_col = gsl_spline_alloc(gsl_interp_cspline, m); ip2d_m = m; }

    gsl_interp_accel *acc = ip2d_acc;
    gsl_spline *splinerow = ip2d_row;
    gsl_spline *splinecol = ip2d_col;

    int j;
    double *ymtmp;

    if (ip2d_mcap < m){ free(ip2d_ymtmp); ip2d_ymtmp = (double*) malloc(sizeof(double) * m); ip2d_mcap = m; }
    ymtmp = ip2d_ymtmp;

    for (j=0;j<m;j++) {
        gsl_spline_init (splinerow, x2a, ya[j], n);
        ymtmp[j] = gsl_spline_eval (splinerow, x2, acc);
    }
  
    gsl_spline_init (splinecol, x1a, ymtmp, m);
 
    *y = gsl_spline_eval (splinecol, x1, acc);

    /* cached objects are reused across calls and freed at shutdown */
}

