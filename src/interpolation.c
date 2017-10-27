/*! File interpolation.c returns from the generated tables the required 
values of the integral.*/


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

void interpolation(double nHX1, double nHeX2, double nHX13, double radius){

    
    int i,j;    
    double y1, y2, y3;
    
    
    gsl_interp_accel *acc    = gsl_interp_accel_alloc ();
    gsl_spline       *spline = gsl_spline_alloc (gsl_interp_cspline, INTERPOINTS);    

    
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
  
    interp2d (nHx1a, nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1, nHeX2, &y2);  
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
  
  
    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
    
    //fuku_e1h1[iGrid] = Ag * (y1 + y2 + y3) / (radius * radius * kpc2Mpc);
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
  
    interp2d (nHx1a, nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1, nHeX2, &y2);
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
  
    y2 = fabs (y2);
    y3 = fabs (y3);
//     fuku_ehe1[iGrid] = Ag * (y2 + y3) / (radius * radius * kpc2Mpc);
    fuku_ehe1[iGrid] = (y2 + y3) / (radNorm);
  
    /* Third Integral */
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  matrix_ehe2_p3, 
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);
  
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
  
    y3 = fabs (y3);
//     fuku_ehe2[iGrid] = Ag * (y3) / (radius * radius * kpc2Mpc);
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
    
    interp2d (nHx1a,  nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1,  nHeX2, &y2);
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);
    
    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
  
//   integral_H1[iGrid] = Ag * (y1 + y2 + y3) / (radius * radius * kpc2Mpc);
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

    interp2d (nHx1a,  nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1,  nHeX2, &y2);
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);  
  
    y2 = fabs (y2);
    y3 = fabs (y3);

//     integral_He1[iGrid] = Ag * (y2 + y3) / (radius * radius * kpc2Mpc);
    integral_He1[iGrid] = (y2 + y3) / ( radNorm );
  
    
    /* Third Integral */
    for (i = 0; i < INTERPOINTS; i++)
        for (j = 0; j < INTERPOINTS; j++)
            yP3a[i][j] = gsl_matrix_get(  temp_matrix_ehe2_p3,
                                          (long int) ((nHx13a[i] - LOWLIM) / TABLERES) + 2,
                                          (long int) ((nHex2a[j] - LOWLIM) / TABLERES) + 2);  
  
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);  
  
    y3 = fabs (y3);
//   integral_He2[iGrid] = Ag * (y3) / (radius * radius * kpc2Mpc);
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
  
    interp2d (nHx1a, nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1, nHeX2, &y2);
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);

    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
  
//  comp_integ1[iGrid] = Ag * (y1 + y2 + y3) / (radius * radius * kpc2Mpc);
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
    
    interp2d (nHx1a,  nHex2a, yP2a, INTERPOINTS, INTERPOINTS, nHX1,  nHeX2, &y2);
    interp2d (nHx13a, nHex2a, yP3a, INTERPOINTS, INTERPOINTS, nHX13, nHeX2, &y3);

    y1 = fabs (y1);
    y2 = fabs (y2);
    y3 = fabs (y3);
    
//     comp_integ2[iGrid] = Ag * (y1 + y2 + y3) / (radius * radius * kpc2Mpc);
    comp_integ2[iGrid] = (y1 + y2 + y3) / ( radNorm);

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

}

