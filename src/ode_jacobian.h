#include <stdio.h>
#include <math.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "constants.h"
#include "allvars.h"
#include "prototype.h"

/***************************************************************
 * Etc
 ***************************************************************/
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

/*************************************************************** 
 * 
 * We compute the jacobian needed for the Rosenbrock ODE solver.
 *
 * Terms of the Jacobian are computed analytically using MATHEMATICA
 * and then coverted to ANSI C format.
 *
 * We need the Jacobian for the set of coupled ODEs found in Fukugita
 * & Kawasaki (1994), eq-26 to 36.
 *
 * Therefore, the terms that need to be computed are,
 *  J[1][1] = d/dnH2(d/dt nH2)   J[1][2] = d/dnHe2(d/dt nH2)  ...
 *  J[2][1] = d/dnH2(d/dt nHe2)  J[2][2] = d/dnHe2(d/dt nHe2) ...
 *   .
 *   .
 *  J[4][1] = d/dnH2(d/dt T)..............J[4][4]=d/dT(d/dt T)
 *
 *
 **************************************************************/
  
struct ode_system_jacobi
{
    
#ifdef EULER
    void operator()( const vector_type &y  , matrix_type &J , const double & /* t */  ) const
#else
     void operator()( const vector_type &y  , matrix_type &J , const double & /* t */ , vector_type &dfdt ) 
#endif    
    {
        
        
        
        
        double fe1h1        = (params.fe1h1);   /*    fuku_e1h1[iGrid]     */        
        double fehe1        = (params.fehe1);   /*    fuku_ehe1[iGrid]     */        
        double fehe2        = (params.fehe2);   /*    fuku_ehe2[iGrid]     */        
        double integralH1   = (params.intH1);   /*    integral_H1[iGrid]   */        
        double integralHe1  = (params.intHe1);   /*    integral_He1[iGrid]  */        
        double integralHe2  = (params.intHe2);   /*    integral_He2[iGrid]  */
        double localOD      = (params.localOD); /*     over_densities[iGrid] */
                
        /* For all terms in equation-(26) of Fukugita & Kawasaki (1994) */
        double dT_a2h2, dT_b1h1;

        /* For all terms in equation-(29 & 30) of Fukugita & Kawasaki (1994) */
        double dT_bhe1, dT_bhe2, dT_ahe2, dT_ahe3, dT_zhe2;

        /* For all terms in equation-(36) of Fukugita & Kawasaki (1994) */
        double dT_zetah1, dT_zetahe1, dT_zetahe2;
        double dT_etah2, dT_etahe2, dT_etahe3;
        double dT_whe2;
        double dT_psih1, dT_psihe1, dT_psihe2;
        double dT_thetaff;

        double n_e,n_Hn, n_Hen;

        double a2h2;

        double bhe1, bhe2, b1h1;
        double ahe2, ahe3;
        double zhe2;

        double zeta_H1, zeta_He1, zeta_He2;
        double eta_H2, eta_He2, eta_He3;
        double psi_H1, psi_He1, psi_He2;
        double thetaff, w_He2, mucon = 1.24;

        //---eq(B6) Ref(2)
        a2h2    = 2.6e-13  * pow (y[3] / 1.e4, -0.8);

        //---eq(B1) Ref(2)
        b1h1    = 5.85e-11 * pow (y[3], .5) * pow (1 + pow (y[3] / 1.e5, .5), -1.) * exp (-1.578e5 / y[3]);

        //---eq(B3) Ref(2)
        bhe1    = 2.38e-11 * pow (y[3], .5) * pow (1 + pow (y[3] / 1e5, .5), -1.) * exp (-2.853e5 / y[3]);

        //---eq(B4) Ref(2)
        bhe2    = 5.68e-12 * pow (y[3], .5) * pow (1 + pow (y[3] / 1e5, .5), -1.) * exp (-6.315e5 / y[3]);

        //---eq(B7) Ref(2)
        ahe2    = 1.50e-10 * pow (y[3], -.6353);

        //---eq(B11) Ref(2)
        ahe3    = 3.36e-10 * pow (y[3], -.5) * pow (y[3] / 1e3, -.2) * pow (1 + pow (y[3] / 4e6, .7), -1.);

        //---eq(B10) Ref(2)
        zhe2    = 1.9e-3   * pow (y[3], -1.5) * exp (-4.7e5 / y[3]) * (1 + 0.3 * exp (-9.4e4 / y[3]));

        //---eq(B28) Ref(2)
        thetaff = 1.42e-27 * 1.1 * pow (y[3], .5);

        //---eq(B17) Ref(2)
        zeta_H1 = 1.27e-21 * pow (y[3], .5) * pow (1 + pow (y[3] / 1.e5, .5), -1.) * exp (-1.58e5 / y[3]);
        
        //---eq(B18) Ref(2)
        zeta_He1= 9.38e-22 * pow (y[3], .5) * pow (1 + pow (y[3] / 1.e5, .5), -1.) * exp (-2.85e5 / y[3]);

        //---eq(B20) Ref(2)
        zeta_He2= 4.95e-22 * pow (y[3], .5) * pow (1 + pow (y[3] / 1.e5, .5), -1.) * exp (-6.31e5 / y[3]);

        //---eq(B21) Ref(2)
        eta_H2  = 6.5e-27 * pow (y[3], .5) * pow (y[3] / 1e3, -.2) * pow (1 + pow (y[3] / 1.e6, .7), -1.);

        //---eq(B22) Ref(2)   
        eta_He2 = 1.55e-26 * pow (y[3], .3647);

        //---eq(B24) Ref(2)
        eta_He3 = 3.48e-26 * pow (y[3], .5) * pow (y[3] / 1e3, -.2) * pow (1 + pow (y[3] / 4.e6, .7), -1.);

        //---eq(B23) Ref(2)
        w_He2   = 1.24e-13 * pow (y[3], -1.5) * exp (-4.7e5 / y[3]) * (1 + .3 * exp (-9.4e4 / y[3]));

        //---eq(B25) Ref(2)
        psi_H1  = 7.5e-19 * pow (1 + pow (y[3] / 1e5, .5), -1.) * exp (-1.18e5 / y[3]);

        //---eq(B26) Ref(2)
        psi_He1 = 9.1e-27 * pow (y[3], -.1687) * pow (1 + pow (y[3] / 1e5, .5), -1.) * exp (-1.31e4 / y[3]) * (1.0 * y[0] + 1.0 * y[1] + 2.0 * y[2])
                    * y[1] / (n_He0 * localOD * pow ((1 + z), 3.0) - y[1]);

        //---eq(B27) Ref(2)
        psi_He2 = 5.54e-17 * pow (y[3], -.397) * pow (1 + pow (y[3] / 1e5, .5), -1.) * exp (-4.73e5 / y[3]);


        // densities
        n_Hn    = n_H0  * localOD * pow3( 1. + z ) - y[0];
        n_Hen   = n_He0 * localOD * pow3( 1. + z ) - y[1] - y[2];
        n_e     = y[0] + y[1] + 2*y[2];

#ifdef STROEMGRENTEST
        n_Hen = 0.0;
        psi_He1 = 0.0;
        psi_He2 = 0.0;
#endif    
    
    
    /* Here is gets messy. 
     * We are using pre-computed derivatives (found with MATHEMATICA) 
     * ported to ANSI C.
     * 
     * TODO: code refactoring ^^
     * 
     */

    dT_a2h2 = -3.29657e-10 / pow (y[3], 1.8);


    dT_b1h1 = -9.24966e-14 / ( exp (157800. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5) , 2) )
             +9.2313e-6   / (exp (157800. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) 
             * pow(y[3], 1.5)) + 2.925e-11 / ( exp(157800. / y[3]) * (1 + 0.0031622 * pow (y[3], 0.5)) * pow (y[3], 0.5) );



  dT_bhe1 = -3.7631e-14 / (exp (285300. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2))
           +6.79014e-6 / (exp (285300. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3], 1.5)) + 1.19e-11 / (exp (285300. / y[3]) * (1 +
                                                                   0.003162 *
                                                                   pow (y[3], 0.5))  * pow (y[3], 0.5));



  dT_bhe2 =
    -8.9808e-15 / (exp (631500. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    3.58692e-6 / (exp (631500. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3],
                       1.5)) +
    2.8399999999999996e-12 / (exp (631500. / y[3]) *
                              (1 + 0.003162 * pow (y[3], 0.5)) * pow (y[3], 0.5));



  dT_ahe2 = -9.5295e-11 / pow (y[3], 1.6353);

  dT_ahe3 =
    -9.3634e-10 / ((1 + 0.000024 * pow (y[3], 0.7)) * pow (y[3], 1.7)) -
    2.2386e-14 / (pow (1 + 0.00002390881249475094 * pow (y[3], 0.7), 2) *
                  pow (y[3], 1.));




  dT_zhe2 = 53.58 / (exp (564000. / y[3]) * pow (y[3], 3.5)) +
    (893. * (1 + 0.3 / exp (94000. / y[3]))) / (exp (470000. / y[3]) *
                                             pow (y[3],
                                                  3.5)) - (0.00285 * (1 +
                                                                      0.3 /
                                                                      exp
                                                                      (94000.
                                                                       /
                                                                       y[3]))) /
    (exp (470000. / y[3]) * pow (y[3], 2.5));



  dT_zetah1 =
    -2.0080e-24 / (exp (158000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    2.0066e-16 / (exp (158000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3],
                       1.5)) + 6.35e-22 / (exp (158000. / y[3]) * (1 +
                                                                0.003162 *
                                                                pow (y[3],
                                                                     0.5)) *
                                           pow (y[3], 0.5));




  dT_zetahe1 =
    -1.4831e-24 / (exp (285000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    2.6733e-16 / (exp (285000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3],
                       1.5)) + 4.69e-22 / (exp (285000. / y[3]) * (1 +
                                                                0.003162 *
                                                                pow (y[3],
                                                                     0.5)) *
                                           pow (y[3], 0.5));



  dT_zetahe2 =
    -7.8266e-25 / (exp (631000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    3.12345e-16 / (exp (631000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                   pow (y[3],
                        1.5)) + 2.475e-22 / (exp (631000. / y[3]) * (1 +
                                                                  0.003162 *
                                                                  pow (y[3],
                                                                       0.5)) *
                                             pow (y[3], 0.5));



  dT_etah2 =
    7.7631e-27 / ((1 + 0.0000631 * pow (y[3], 0.7)) * pow (y[3], 0.7)) -
    1.1429e-30 / (pow (1 + 0.0000631 * pow (y[3], 0.7), 2) *
                  pow (y[3], 5.5511e-17));


  dT_etahe2 = 5.65285e-27 / pow (y[3], 0.6353);


  dT_etahe3 =
    4.156e-26 / ((1 + 0.000024 * pow (y[3], 0.7)) * pow (y[3], 0.7)) -
    2.31865e-30 / (pow (1 + 0.000024 * pow (y[3], 0.7), 2) * pow (y[3], 5.551e-17));



  dT_whe2 = 3.4968e-9 / (exp (564000. / y[3]) * pow (y[3], 3.5)) +
    (5.828e-8 * (1 + 0.3 / exp (94000. / y[3]))) /
    (exp (470000. / y[3]) * pow (y[3], 3.5)) -
    (1.86e-13 * (1 + 0.3 / exp (94000. / y[3]))) / (exp (470000. / y[3]) *
                                                 pow (y[3], 2.5));



  dT_psih1 =
    8.85e-14 / (exp (118000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                pow (y[3],
                     2)) - 1.18585e-21 / (exp (118000. / y[3]) * pow (1 +
                                                                   0.003162 *
                                                                   pow (y[3],
                                                                        0.5),
                                                                   2) *
                                          pow (y[3], 0.5));



  dT_psihe1 = 1.1921e-22 / (exp (13100. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                           pow (y[3],
                                2.1687)) -
    1.53517e-27 / (exp (13100. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                   pow (y[3],
                        1.1687)) - 1.4388e-29 / (exp (13100. / y[3]) * pow (1 +
                                                                         0.003162
                                                                         *
                                                                         pow
                                                                         (y[3],
                                                                          0.5),
                                                                         2) *
                                                 pow (y[3], 0.6687));




  dT_psihe2 =
    2.62042e-11 / (exp (473000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                   pow (y[3],
                        2.397)) - 2.12e-17 / (exp (473000. / y[3]) * (1 +
                                                                   0.003162 *
                                                                   pow (y[3],
                                                                        0.5))
                                              * pow (y[3],
                                                     1.397)) -
    8.75951e-20 / (exp (473000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2) *
                   pow (y[3], 0.897));



  dT_thetaff = 7.1e-28 / pow (y[3], 0.5);


  /* End of computing derivatives of complicated functions */


/* Computing terms corresponding to the FIRST row of the Jacobian
     i.e., J[1][1] =d/dnH2(d/dT_ nH2); J[1][2]=d/dnHe2(d/dt nH2)
           J[1][3] =d/dnHe3(d/dt nH2); J[1][4]=d/dT(d/dt nH2)
*/

    double dnH2_nH2eq, dnHe2_nH2eq, dnHe3_nH2eq, dT_nH2eq;

    dnH2_nH2eq      = b1h1 * (n_Hn - n_e) - fe1h1 - a2h2 * (n_e + y[0]);

    dnHe2_nH2eq     = n_Hn * b1h1 - a2h2 * y[0];

    dnHe3_nH2eq     = 2 * dnHe2_nH2eq;

    dT_nH2eq        = dT_b1h1 * n_e * n_Hn - dT_a2h2 * n_e * y[0];


  /* Computing terms corresponding to the SECOND row of the Jacobian
     i.e., J[2][1] =d/dnH2(d/dt nHe2); J[2][2]=d/dnHe2(d/dt nHe2)
     J[2][3] =d/dnHe3(d/dt nHe2); J[2][4]=d/dT(d/dt nHe2)
   */

    double dnH2_nHe2eq, dnHe2_nHe2eq, dnHe3_nHe2eq, dT_nHe2eq;

    dnH2_nHe2eq     =  bhe1 * n_Hen - y[1] * (bhe2 + ahe2 + zeta_He2) + ahe3 * y[2];

    dnHe2_nHe2eq    = - fehe1 + bhe1 * (n_Hen - n_e) 
                      - (y[1] + n_e) * (bhe2 + ahe2 + zeta_He2) + ahe3 * y[2];

    dnHe3_nHe2eq    = - fehe1 + bhe1 * (2 * n_Hen - n_e) 
                      - 2 * y[1] *(bhe2 + ahe2 + zeta_He2) 
                      + 2 * ahe3 * y[2] + ahe3 * n_e;

    dT_nHe2eq = n_e * (dT_bhe1 * n_Hen - dT_bhe2 * y[1] - dT_ahe2 * y[1]
                        + dT_ahe3 * y[2] - dT_zhe2 * y[1]);


  /* Computing terms corresponding to the THIRD row of the Jacobian
     i.e., J[3][1] =d/dnH2(d/dt nHe3); J[3][2]=d/dnHe2(d/dt nHe3)
     J[3][3] =d/dnHe3(d/dt nHe3); J[3][4]=d/dT(d/dt nHe3)
   */

    double dnH2_nHe3eq, dnHe2_nHe3eq, dnHe3_nHe3eq, dT_nHe3eq;

    dnH2_nHe3eq  = bhe2 * y[1] - ahe3 * y[2];

    dnHe2_nHe3eq = fehe2 + bhe2 * (y[1] + n_e) - ahe3 * y[2];

    dnHe3_nHe3eq = 2 * dnH2_nHe3eq - ahe3 * n_e;

    dT_nHe3eq    = n_e * (dT_bhe2 * y[1] - ahe3 * y[2]);


  /* Computing terms corresponding to the FOURTH row of the Jacobian
     i.e., J[4][1] =d/dnH2(d/dt T); J[4][2]=d/dnHe2(d/dt T)
     J[4][3] =d/dnHe3(d/dt T); J[4][4]=d/dT(d/dt T)
   */

    double dnH2_Teq, dnHe2_Teq, dnHe3_Teq, dT_Teq;
    double T_CMB0 = myConfig.cosmoTCMB;

    dnH2_Teq =  - integralH1
                - (zeta_H1 * n_Hn + zeta_He1 * n_Hen + zeta_He2 * y[1]) 
                + n_e * zeta_H1
                - (eta_H2 * y[0] + eta_He2 * y[1] + eta_He3 * y[2]) 
                - eta_H2 * n_e
                - w_He2 * y[2]
                - (psi_H1 * n_Hn + psi_He1 * n_Hen + psi_He2 * y[1]) 
                + n_e * psi_H1
                - lamcons * (y[3] - T_CMB0 * (1 + z))
                - thetaff * (n_e + y[0] + y[1] + 4 * y[2]);


    dnHe2_Teq = - integralHe1
                - (zeta_H1 * n_Hn + zeta_He1 * n_Hen + zeta_He2 * y[1]) + n_e * zeta_He1
                - (eta_H2 * y[0] + eta_He2 * y[1] + eta_He3 * y[2]) 
                - eta_He2 * n_e
                - w_He2 * y[2]
                - (psi_H1 * n_Hn + psi_He1 * n_Hen + psi_He2 * y[1]) 
                + n_e * psi_He1
                - lamcons * (y[3] - T_CMB0 * (1 + z))
                - thetaff * (n_e + y[0] + y[1] + 4 * y[2]);


    dnHe3_Teq = - integralHe2
                - 2 * (zeta_H1 * n_Hn + zeta_He1 * n_Hen + zeta_He2 * y[1]) 
                + n_e * zeta_He1 
                - 2 * (eta_H2 * y[0] + eta_He2 * y[1] + eta_He3 * y[2]) 
                - eta_He3 * n_e 
                - 2 * w_He2 * y[2] 
                - w_He2 * n_e 
                - (psi_H1 * n_Hn + psi_He1 * n_Hen + psi_He2 * y[1]) 
                + n_e * psi_He1 
                - 2 * lamcons * (y[3] - T_CMB0 * (1 + z)) 
                - thetaff * (4 * n_e + 2 * (y[0] + y[1] + 4 * y[2]));



    dT_Teq =    - (n_e * (dT_zetah1 * n_Hn + dT_zetahe1 * n_Hen + dT_zetahe2 * y[1]))
                - (n_e * (dT_etah2 * y[0] + dT_etahe2 * y[1] + dT_etahe3 * y[2]))
                - n_e * w_He2 * y[2] 
                - (n_e * (dT_psih1 * n_Hn + dT_psihe1 * n_Hen + dT_psihe2 * y[1]))
                - lamcons * n_e - dT_thetaff * (y[0] + y[1] + 4 * y[2]) * n_e
                - 7.5 * hubble (z) * (k_BeV * nB0 * pow ((1 + z), 3.0) / mucon);

                
                

    double T_denom;

    T_denom      = k_BeV * nB0 * pow3(1 + z) * 1.5;

    dnH2_Teq    /= T_denom;
    dnHe2_Teq   /= T_denom;
    dnHe3_Teq   /= T_denom;
    dT_Teq      /= T_denom;
    
            
    // And finally, here is the jacobian   
    J( 0 , 0 ) = MYR * dnH2_nH2eq;
    J( 0 , 1 ) = MYR * dnHe2_nH2eq;
    J( 0 , 2 ) = MYR * dnHe3_nH2eq;
    J( 0 , 3 ) = MYR * dT_nH2eq;

    J( 1 , 0 ) = MYR * dnH2_nHe2eq;
    J( 1 , 1 ) = MYR * dnHe2_nHe2eq;
    J( 1 , 2 ) = MYR * dnHe3_nHe2eq;
    J( 1 , 3 ) = MYR * dT_nHe2eq;

    J( 2 , 0 ) = MYR * dnH2_nHe3eq;
    J( 2 , 1 ) = MYR * dnHe2_nHe3eq;
    J( 2 , 2 ) = MYR * dnHe3_nHe3eq;
    J( 2 , 3 ) = MYR * dT_nHe3eq;

    J( 3 , 0 ) = MYR * dnH2_Teq;         
    J( 3 , 1 ) = MYR * dnHe2_Teq;        
    J( 3 , 2 ) = MYR * dnHe3_Teq;        
    J( 3 , 3 ) = MYR * dT_Teq;           

    // There is no explicit time-dependence        
#ifndef EULER
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    dfdt[2] = 0.0;
    dfdt[3] = 0.0;             
#endif    
    
    
    
#ifdef STROEMGRENTEST
        J( 0 , 1 ) = J( 0 , 2 ) = J( 0 , 3 ) = 0.0;
        J( 1 , 0 ) = J( 1 , 1 ) = J( 1 , 2 ) = J( 1 , 3 ) = 0.0;
        J( 2 , 0 ) = J( 2 , 1 ) = J( 2 , 2 ) = J( 2 , 3 ) = 0.0;
        J( 3 , 1 ) = J( 3 , 2 ) = 0.0;
#endif    

    
    
        
        
    // Testing   
//     if (iGrid==1 || iGrid==800){
//         printf("--------------------------------------\n");
//         printf("jacobian: dfdy[0][0] = %le at iGrid = %d \n",dnH2_nH2eq, iGrid);
//         printf("jacobian: dfdy[0][1] = %le \n",dnHe2_nH2eq);
//         printf("jacobian: dfdy[0][2] = %le \n",dnHe3_nH2eq);
//         printf("jacobian: dfdy[0][3] = %le \n",dT_nH2eq);
//                 
//         printf("jacobian: dfdy[1][0] = %le \n",dnH2_nHe2eq);
//         printf("jacobian: dfdy[1][1] = %le \n",dnHe2_nHe2eq);
//         printf("jacobian: dfdy[1][2] = %le \n",dnHe3_nHe2eq);
//         printf("jacobian: dfdy[1][3] = %le \n",dT_nHe2eq);
//                 
//         printf("jacobian: dfdy[2][0] = %le \n",dnH2_nHe3eq);
//         printf("jacobian: dfdy[2][1] = %le \n",dnHe2_nHe3eq);
//         printf("jacobian: dfdy[2][2] = %le \n",dnHe3_nHe3eq);
//         printf("jacobian: dfdy[2][3] = %le \n",dT_nHe3eq);
//                 
//         printf("jacobian: dfdy[3][0] = %le \n",dnH2_Teq);
//         printf("jacobian: dfdy[3][1] = %le \n",dnHe2_Teq);
//         printf("jacobian: dfdy[3][2] = %le \n",dnHe3_Teq);
//         printf("jacobian: dfdy[3][3] = %le \n",dT_Teq);
//    }
        
 
    } // operator
}; // struct
