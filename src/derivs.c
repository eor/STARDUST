/***************************************************************
 * Sets up the system of differential equations as described in 
 * Fukugita & Kawasaki (Ref [2])
 ***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***************************************************************
 * GSL header files
 ***************************************************************/
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

/***************************************************************
 * Constants, function declarations, etc
 ***************************************************************/
#include "constants.h"

#include "allvars.h"
#include "prototype.h"

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))


/***************************************************************
 * compute derivatives
 ***************************************************************/
int derivs (double t, const double y[], double dydt[], void * p){
  
    /* extract the parameters */
    struct ode_params * params = (struct ode_params *)p;

    double fe1h1        = (params->fe1h1);   /*    fuku_e1h1[iGrid]     */
    double fehe1        = (params->fehe1);   /*    fuku_ehe1[iGrid]     */
    double fehe2        = (params->fehe2);   /*    fuku_ehe2[iGrid]     */
    double integralH1   = (params->intH1);   /*    integral_H1[iGrid]   */
    double integralHe1  = (params->intHe1);   /*    integral_He1[iGrid]  */
    double integralHe2  = (params->intHe2);   /*    integral_He2[iGrid]  */
    
    
    double y0 = y[0]; 
    double y1 = y[1];
    double y2 = y[2];
    double y3 = y[3];
    

    
    // sanity checks. number densities and T_e should not be negative
    y0 = GSL_MAX_DBL( y0, 0. ); // n_H_II
    y1 = GSL_MAX_DBL( y1, 0. ); // n_He_II 
    y2 = GSL_MAX_DBL( y2, 0. ); // n_He_III 
    y3 = GSL_MAX_DBL( y3, 0. ); // T_e 
        
    

    double n_Hn, n_Hen;

    double a2h2, R1c;

    double bhe1, bhe2, b1h1;
    double ahe2, ahe3;
    double zhe2;

    double zeta_H1, zeta_He1, zeta_He2;
    double eta_H2, eta_He2, eta_He3;
    double psi_H1, psi_He1, psi_He2;
    double thetaff, w_He2, mucon = 1.24;
    
    double T_CMB0 = myConfig.cosmoTCMB;
    
    // TODO change names , do so also in jac.c
    

    
    /* eq(B6) in Ref [2] */
    a2h2 = 2.6e-13 * pow(y3 / 1.e4, -0.8);                                                    // alpha2HII
    
    /* eq(B1) in Ref [2] */
    b1h1 = 5.85e-11 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-1.578e5 / y3);    // beta1HI
    
    /* eq(B3) in Ref [2] */
    bhe1 = 2.38e-11 * pow(y3, .5) * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-2.853e5 / y3);     // betaHeI    
    
    /* eq(B4) in Ref [2] */
    bhe2 = 5.68e-12 * pow(y3, .5) * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-6.315e5 / y3);     // betaHeII
    
    /* eq(B7) in Ref [2] */
    ahe2 = 1.50e-10 * pow(y3, -.6353);
    
    /* eq(B11) in Ref [2] */
    ahe3 = 3.36e-10 * pow(y3, -.5) * pow(y3 / 1e3, -.2) * pow(1 + pow(y3 / 4e6, .7), -1.);
    
    /* eq(B10) in Ref [2] */
    zhe2 = 1.9e-3 * pow(y3, -1.5) * exp(-4.7e5 / y3) * (1 + 0.3 * exp(-9.4e4 / y3));
    
    /* eq(B28) in Ref [2] */
    thetaff = 1.42e-27 * 1.1 * pow(y3, .5);
    
    /* eq(B17) in Ref [2] */
    zeta_H1 = 1.27e-21 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-1.58e5 / y3);    

    /* eq(B18) in Ref [2] */
    zeta_He1 = 9.38e-22 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-2.85e5 / y3);
    
    /* eq(B20) in Ref [2] */
    zeta_He2 = 4.95e-22 * pow(y3, .5) * pow(1 + pow(y3 / 1.e5, .5), -1.) * exp(-6.31e5 / y3);
    
    /* eq(B21) in Ref [2] */
    eta_H2 = 6.5e-27 * pow(y3, .5) * pow(y3 / 1e3, -.2) * pow(1 + pow(y3 / 1.e6, .7), -1.);
    
    /* eq(B22) in Ref [2] */
    eta_He2 = 1.55e-26 * pow(y3, .3647);
    
    /* eq(B24) in Ref [2] */
    eta_He3 = 3.48e-26 * pow(y3, .5) * pow(y3 / 1e3, -.2) * pow(1 + pow(y3 / 4.e6, .7), -1.);
    
    /* eq(B23) in Ref [2] */
    w_He2 = 1.24e-13 * pow(y3, -1.5) * exp(-4.7e5 / y3) * (1 + .3 * exp(-9.4e4 / y3));
    
    /* eq(B25) in Ref [2] */
    psi_H1 = 7.5e-19 * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-1.18e5 / y3);
    
    /* eq(B26) in Ref [2] */
    psi_He1 = 9.1e-27 * pow(y3, -.1687) * pow(1 + pow(y3 / 1e5, .5),-1.) 
                * exp(-1.31e4 / y3) * (1.0 * y0 + 1.0 * y1 + 2.0 * y2)
                * y1 / (n_He0 * OverDensity * pow((1 + z), 3.0) - y1);
    
    /* eq(B27) in Ref [2] */
    psi_He2 =  5.54e-17 * pow(y3, -.397) * pow(1 + pow(y3 / 1e5, .5), -1.) * exp(-4.73e5 / y3);
    


    n_Hn  = n_H0  * OverDensity * pow3(1 + z) - y0;
    n_Hen = n_He0 * OverDensity * pow3(1 + z) - y1 - y2;
    
     
    
 
    /* Explicitly regularize */

    
    n_Hn  = GSL_MAX_DBL( n_Hn,  1.e-55 );   
    n_Hen = GSL_MAX_DBL( n_Hen, 1.e-55 );
    n_Hn  = GSL_MIN_DBL( n_Hn,  n_H0  * OverDensity * pow3(1 + z) ); 
    n_Hen = GSL_MIN_DBL( n_Hen, n_He0 * OverDensity * pow3(1 + z) );
    
    #ifdef STROEMGRENTEST
    n_Hen = 0.0;
    psi_He1 = 0.0; //?
    psi_He2 = 0.0; //?  check!
    integralHe1 = 0.0;
    integralHe2 = 0.0;
    fehe1 = 0.0;
    fehe2 = 0.0;
    #endif    



    R1c = b1h1 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) + fe1h1;

    /* eq(26) in Ref [2] */
    dydt[0] = R1c * n_Hn - a2h2 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y0;
  
    /* eq(29) in Ref [2] */
    dydt[1] =   n_Hen * fehe1 
                + bhe1 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * n_Hen 
                - bhe2 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y1 
                - ahe2 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y1 
                + ahe3 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y2 
                - zhe2 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y1;
    

    dydt[2] =   y1 * fehe2 
                + bhe2 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y1 
                - ahe3 * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * y2;
    
    /* eq(36) in Ref [2] */
    dydt[3] = ( (n_Hn * integralH1 + n_Hen * integralHe1 + y1 * integralHe2)
                - ((1.0 * y0 + 1.0 * y1 + 2.0 * y2) * (zeta_H1 * n_Hn + zeta_He1 * n_Hen + zeta_He2 * y1)) 
                - ((1.0 * y0 + 1.0 * y1 + 2.0 * y2) * (eta_H2 * y0 + eta_He2 * y1 + eta_He3 * y2)) 
                - (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * w_He2 * y2 
                - ( (1.0 * y0 + 1.0 * y1 + 2.0 * y2) * (psi_H1 * n_Hn + psi_He1 * n_Hen + psi_He2 * y1) ) 
                - (y3 - T_CMB0 * (1 + z)) * lamcons * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) 
                - thetaff * (y0 + y1 + 4 * y2) * (1.0 * y0 + 1.0 * y1 + 2.0 * y2) 
                - 7.5 * hubble (z) * (k_BeV * y3 * nB0 * pow((1 + z), 3.0) / mucon)
              ) * mucon / (k_BeV * nB0 * pow3(1 + z) * 1.5);

  
              
    #ifdef STROEMGRENTEST
    dydt[1] = dydt[2] = 0.0;
    #endif             
              
              
    return GSL_SUCCESS;
}
