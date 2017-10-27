#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

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

  /* We compute the jacobian needed to use  GSL - ODE solvers.

     Terms of the Jacobian are computed analytically using MATHEMATICA
     and then coverted to ANSI C format.

     We need the Jacobian is for the set of coupled ODEs found in Fukugita
     & Kawasaki (1994), eq-26 to 36.

     Therefore, the terms that need to be computed are,
     J[1][1] = d/dnH2(d/dt nH2); J[1][2]=d/dnHe2(d/dt nH2)...............
     J[2][1] = d/dnH2(d/dt nHe2);J[2][2] = d/dnHe2(d/dt nHe2).............
     .
     .
     J[4][1] = d/dnH2(d/dt T)..............J[4][4]=d/dT(d/dt T)


     NOTE: The matrix is set according to GSL conventions. But you can
     change this.  int (* jac) (double t, const double y[], double * dfdy,
     double dfdt[], void * params); This function should store the vector
     of derivative elements df_i(t,y,params)/dt in the array dfdt and the
     Jacobian matrix J_{ij} in the array dfdy, regarded as a row-ordered
     matrix J(i,j) = dfdy[i * dimension + j] where dimension is the
     dimension of the system. The function should return GSL_SUCCESS if the
     calculation was completed successfully. Any other return value
     indicates an error.

     Some of the simpler solver algorithms do not make use of the
     Jacobian matrix, so it is not always strictly necessary to provide
     it (the jacobian element of the struct can be replaced by a null
     pointer for those algorithms). However, it is useful to provide
     the Jacobian to allow the solver algorithms to be
     interchanged%Gâ€”%@the best algorithms make use of the
     Jacobian.*** REFER GSL MANUAL, 3rd Edition Mark Galassi et al.,***

   */

int jac (double t, const double y[], double *dfdy, double dfdt[], void *p){


    /* Extract parameters required for calculation below */
  
    struct ode_params * params  = (struct ode_params *)p;

    double fuku_e1h1    = (params->fe1h1);    /*    fuku_e1h1[iGrid]     */
    double fuku_ehe1    = (params->fehe1);    /*    fuku_ehe1[iGrid]     */
    double fuku_ehe2    = (params->fehe2);    /*    fuku_ehe2[iGrid]     */
    double integral_H1  = (params->intH1);    /*    integral_H1[iGrid]   */
    double integral_He1 = (params->intHe1);   /*    integral_He1[iGrid]  */
    double integral_He2 = (params->intHe2);   /*    integral_He2[iGrid]  */


    /* For all terms in equation-(26) of Fukugita & Kawasaki (1994) */
    double dTa2h2, dTb1h1;

    /* For all terms in equation-(29 & 30) of Fukugita & Kawasaki (1994) */
    double dTbhe1, dTbhe2, dTahe2, dTahe3, dTzhe2;

    /* For all terms in equation-(36) of Fukugita & Kawasaki (1994) */
    double dTzetah1, dTzetahe1, dTzetahe2;
    double dTetah2, dTetahe2, dTetahe3;
    double dTwhe2;
    double dTpsih1, dTpsihe1, dTpsihe2;
    double dTthetaff;


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
                * y[1] / (n_He0 * OverDensity * pow ((1 + z), 3.0) - y[1]);

    //---eq(B27) Ref(2)
    psi_He2 = 5.54e-17 * pow (y[3], -.397) * pow (1 + pow (y[3] / 1e5, .5), -1.) * exp (-4.73e5 / y[3]);


    // densities
    n_Hn    = n_H0  * OverDensity * pow3( 1. + z ) - y[0];
    n_Hen   = n_He0 * OverDensity * pow3( 1. + z ) - y[1] - y[2];
    n_e     = y[0] + y[1] + 2*y[2];

   #ifdef STROEMGRENTEST
    n_Hen = 0.0;
    psi_He1 = 0.0;
    psi_He2 = 0.0;
    #endif    
    
    
//     if (n_Hn<0) printf("jac.c: n_Hn=%e\n", n_Hn );
//     if (n_e<0)  printf("jac.c: n_e =%e\n", n_e );
//     if (y[0]<0 || y[1]<0 || y[2]<0 || y[3]<0)printf("jac.c: y0=%e\ty1=%e\ty2=%e\ty3=%e\n", y[0], y[1], y[2], y[3] );
//     
    
    
    
    /* Here is gets messy. 
     * We are using pre-computed derivatives (found with MATHEMATICA) 
     * ported to ANSI C.
     * 
     * TODO: code refactorization ^^
     * 
     */

    dTa2h2 = -3.29657e-10 / pow (y[3], 1.8);


    dTb1h1 = -9.24966e-14 / ( exp (157800. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5) , 2) )    
             +9.2313e-6   / (exp (157800. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) 
             * pow(y[3], 1.5)) + 2.925e-11 / ( exp(157800. / y[3]) * (1 + 0.0031622 * pow (y[3], 0.5)) * pow (y[3], 0.5) );



  dTbhe1 =
    -3.7631e-14 / (exp (285300. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    6.79014e-6 / (exp (285300. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3], 1.5)) + 1.19e-11 / (exp (285300. / y[3]) * (1 +
                                                                   0.003162 *
                                                                   pow (y[3],
                                                                        0.5))
                                              * pow (y[3], 0.5));



  dTbhe2 =
    -8.9808e-15 / (exp (631500. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    3.58692e-6 / (exp (631500. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3],
                       1.5)) +
    2.8399999999999996e-12 / (exp (631500. / y[3]) *
                              (1 + 0.003162 * pow (y[3], 0.5)) * pow (y[3], 0.5));



  dTahe2 = -9.5295e-11 / pow (y[3], 1.6353);

  dTahe3 =
    -9.3634e-10 / ((1 + 0.000024 * pow (y[3], 0.7)) * pow (y[3], 1.7)) -
    2.2386e-14 / (pow (1 + 0.00002390881249475094 * pow (y[3], 0.7), 2) *
                  pow (y[3], 1.));




  dTzhe2 = 53.58 / (exp (564000. / y[3]) * pow (y[3], 3.5)) +
    (893. * (1 + 0.3 / exp (94000. / y[3]))) / (exp (470000. / y[3]) *
                                             pow (y[3],
                                                  3.5)) - (0.00285 * (1 +
                                                                      0.3 /
                                                                      exp
                                                                      (94000.
                                                                       /
                                                                       y[3]))) /
    (exp (470000. / y[3]) * pow (y[3], 2.5));



  dTzetah1 =
    -2.0080e-24 / (exp (158000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    2.0066e-16 / (exp (158000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3],
                       1.5)) + 6.35e-22 / (exp (158000. / y[3]) * (1 +
                                                                0.003162 *
                                                                pow (y[3],
                                                                     0.5)) *
                                           pow (y[3], 0.5));




  dTzetahe1 =
    -1.4831e-24 / (exp (285000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    2.6733e-16 / (exp (285000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                  pow (y[3],
                       1.5)) + 4.69e-22 / (exp (285000. / y[3]) * (1 +
                                                                0.003162 *
                                                                pow (y[3],
                                                                     0.5)) *
                                           pow (y[3], 0.5));



  dTzetahe2 =
    -7.8266e-25 / (exp (631000. / y[3]) * pow (1 + 0.003162 * pow (y[3], 0.5), 2)) +
    3.12345e-16 / (exp (631000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                   pow (y[3],
                        1.5)) + 2.475e-22 / (exp (631000. / y[3]) * (1 +
                                                                  0.003162 *
                                                                  pow (y[3],
                                                                       0.5)) *
                                             pow (y[3], 0.5));



  dTetah2 =
    7.7631e-27 / ((1 + 0.0000631 * pow (y[3], 0.7)) * pow (y[3], 0.7)) -
    1.1429e-30 / (pow (1 + 0.0000631 * pow (y[3], 0.7), 2) *
                  pow (y[3], 5.5511e-17));


  dTetahe2 = 5.65285e-27 / pow (y[3], 0.6353);


  dTetahe3 =
    4.156e-26 / ((1 + 0.000024 * pow (y[3], 0.7)) * pow (y[3], 0.7)) -
    2.31865e-30 / (pow (1 + 0.000024 * pow (y[3], 0.7), 2) * pow (y[3], 5.551e-17));



  dTwhe2 = 3.4968e-9 / (exp (564000. / y[3]) * pow (y[3], 3.5)) +
    (5.828e-8 * (1 + 0.3 / exp (94000. / y[3]))) /
    (exp (470000. / y[3]) * pow (y[3], 3.5)) -
    (1.86e-13 * (1 + 0.3 / exp (94000. / y[3]))) / (exp (470000. / y[3]) *
                                                 pow (y[3], 2.5));



  dTpsih1 =
    8.85e-14 / (exp (118000. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
                pow (y[3],
                     2)) - 1.18585e-21 / (exp (118000. / y[3]) * pow (1 +
                                                                   0.003162 *
                                                                   pow (y[3],
                                                                        0.5),
                                                                   2) *
                                          pow (y[3], 0.5));



  dTpsihe1 = 1.1921e-22 / (exp (13100. / y[3]) * (1 + 0.003162 * pow (y[3], 0.5)) *
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




  dTpsihe2 =
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



  dTthetaff = 7.1e-28 / pow (y[3], 0.5);


  /* End of computing derivatives of complicated functions */


/* Computing terms correspoding to the FIRST row of the Jacobian 
     i.e., J[1][1] =d/dnH2(d/dt nH2); J[1][2]=d/dnHe2(d/dt nH2)
           J[1][3] =d/dnHe3(d/dt nH2); J[1][4]=d/dT(d/dt nH2)
*/

    double dnH2_nH2eq, dnHe2_nH2eq, dnHe3_nH2eq, dT_nH2eq;

    dnH2_nH2eq      = b1h1 * (n_Hn - n_e) - fuku_e1h1 - a2h2 * (n_e + y[0]);

    dnHe2_nH2eq     = n_Hn * b1h1 - a2h2 * y[0];

    dnHe3_nH2eq     = 2 * dnHe2_nH2eq;

    dT_nH2eq        = dTb1h1 * n_e * n_Hn - dTa2h2 * n_e * y[0];


  /* Computing terms correspoding to the SECOND row of the Jacobian 
     i.e., J[2][1] =d/dnH2(d/dt nHe2); J[2][2]=d/dnHe2(d/dt nHe2)
     J[2][3] =d/dnHe3(d/dt nHe2); J[2][4]=d/dT(d/dt nHe2)
   */

    double dnH2_nHe2eq, dnHe2_nHe2eq, dnHe3_nHe2eq, dT_nHe2eq;

    dnH2_nHe2eq     =  bhe1 * n_Hen - y[1] * (bhe2 + ahe2 + zeta_He2) + ahe3 * y[2];

    dnHe2_nHe2eq    = - fuku_ehe1 + bhe1 * (n_Hen - n_e) 
                      - (y[1] + n_e) * (bhe2 + ahe2 + zeta_He2) + ahe3 * y[2];

    dnHe3_nHe2eq    = - fuku_ehe1 + bhe1 * (2 * n_Hen - n_e) 
                      - 2 * y[1] *(bhe2 + ahe2 + zeta_He2) 
                      + 2 * ahe3 * y[2] + ahe3 * n_e;

    dT_nHe2eq = n_e * (dTbhe1 * n_Hen - dTbhe2 * y[1] - dTahe2 * y[1]
                        + dTahe3 * y[2] - dTzhe2 * y[1]);


  /* Computing terms correspoding to the THIRD row of the Jacobian 
     i.e., J[3][1] =d/dnH2(d/dt nHe3); J[3][2]=d/dnHe2(d/dt nHe3)
     J[3][3] =d/dnHe3(d/dt nHe3); J[3][4]=d/dT(d/dt nHe3)
   */

    double dnH2_nHe3eq, dnHe2_nHe3eq, dnHe3_nHe3eq, dT_nHe3eq;

    dnH2_nHe3eq  = bhe2 * y[1] - ahe3 * y[2];

    dnHe2_nHe3eq = fuku_ehe2 + bhe2 * (y[1] + n_e) - ahe3 * y[2];

    dnHe3_nHe3eq = 2 * dnH2_nHe3eq - ahe3 * n_e;

    dT_nHe3eq    = n_e * (dTbhe2 * y[1] - ahe3 * y[2]);


  /* Computing terms correspoding to the FOURTH row of the Jacobian 
     i.e., J[4][1] =d/dnH2(d/dt T); J[4][2]=d/dnHe2(d/dt T)
     J[4][3] =d/dnHe3(d/dt T); J[4][4]=d/dT(d/dt T)
   */

    double dnH2_Teq, dnHe2_Teq, dnHe3_Teq, dT_Teq;
    double T_CMB0 = myConfig.cosmoTCMB;

    dnH2_Teq =  - integral_H1
                - (zeta_H1 * n_Hn + zeta_He1 * n_Hen + zeta_He2 * y[1]) 
                + n_e * zeta_H1
                - (eta_H2 * y[0] + eta_He2 * y[1] + eta_He3 * y[2]) 
                - eta_H2 * n_e
                - w_He2 * y[2]
                - (psi_H1 * n_Hn + psi_He1 * n_Hen + psi_He2 * y[1]) 
                + n_e * psi_H1
                - lamcons * (y[3] - T_CMB0 * (1 + z))
                - thetaff * (n_e + y[0] + y[1] + 4 * y[2]);


    dnHe2_Teq = - integral_He1
                - (zeta_H1 * n_Hn + zeta_He1 * n_Hen + zeta_He2 * y[1]) + n_e * zeta_He1
                - (eta_H2 * y[0] + eta_He2 * y[1] + eta_He3 * y[2]) 
                - eta_He2 * n_e
                - w_He2 * y[2]
                - (psi_H1 * n_Hn + psi_He1 * n_Hen + psi_He2 * y[1]) 
                + n_e * psi_He1
                - lamcons * (y[3] - T_CMB0 * (1 + z))
                - thetaff * (n_e + y[0] + y[1] + 4 * y[2]);


    dnHe3_Teq = - integral_He2
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



    dT_Teq =    - (n_e * (dTzetah1 * n_Hn + dTzetahe1 * n_Hen + dTzetahe2 * y[1])) 
                - (n_e * (dTetah2 * y[0] + dTetahe2 * y[1] + dTetahe3 * y[2])) 
                - n_e * w_He2 * y[2] 
                - (n_e * (dTpsih1 * n_Hn + dTpsihe1 * n_Hen + dTpsihe2 * y[1])) 
                - lamcons * n_e - dTthetaff * (y[0] + y[1] + 4 * y[2]) * n_e 
                - 7.5 * hubble (z) * (k_BeV * nB0 * pow ((1 + z), 3.0) / mucon);

  /*

     This is basically the Jacobian. We will use GSL matrix_set 
     to place it in the correct format.

     dfdy[0][0] = dnH2_nH2eq;
     dfdy[0][1] = dnHe2_nH2eq;
     dfdy[0][2] = dnHe3_nH2eq;
     dfdy[0][3] = dT_nH2eq;

     dfdy[1][1] = dnH2_nHe2eq;
     dfdy[1][1] = dnHe2_nHe2eq;
     dfdy[1][2] = dnHe3_nHe2eq;
     dfdy[1][3] = dT_nHe2eq;

     dfdy[2][1] = dnH2_nHe3eq;
     dfdy[2][1] = dnHe2_nHe3eq;
     dfdy[2][2] = dnHe3_nHe3eq;
     dfdy[2][3] = dT_nHe3eq;

     dfdy[3][1] = dnH2_Teq / (mucon / (k_BeV * nB0 * pow ((1 + z), 3.0) * 1.5));
     dfdy[3][1] = dnHe2_Teq / (mucon / (k_BeV * nB0 * pow ((1 + z), 3.0) * 1.5));
     dfdy[3][2] = dnHe3_Teq / (mucon / (k_BeV * nB0 * pow ((1 + z), 3.0) * 1.5));
     dfdy[3][3] = dT_Teq / (mucon / (k_BeV * nB0 * pow ((1 + z), 3.0) * 1.5));


   */

    double T_denom;

    T_denom      = k_BeV * nB0 * pow3(1 + z) * 1.5;

    dnH2_Teq    /= T_denom;
    dnHe2_Teq   /= T_denom;
    dnHe3_Teq   /= T_denom;
    dT_Teq      /= T_denom;



  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
  gsl_matrix *m = &dfdy_mat.matrix;


  gsl_matrix_set (m, 0, 0, dnH2_nH2eq);
  gsl_matrix_set (m, 0, 1, dnHe2_nH2eq);
  gsl_matrix_set (m, 0, 2, dnHe3_nH2eq);
  gsl_matrix_set (m, 0, 3, dT_nH2eq);

  gsl_matrix_set (m, 1, 0, dnH2_nHe2eq);
  gsl_matrix_set (m, 1, 1, dnHe2_nHe2eq);
  gsl_matrix_set (m, 1, 2, dnHe3_nHe2eq);
  gsl_matrix_set (m, 1, 3, dT_nHe2eq);

  gsl_matrix_set (m, 2, 0, dnH2_nHe3eq);
  gsl_matrix_set (m, 2, 1, dnHe2_nHe3eq);
  gsl_matrix_set (m, 2, 2, dnHe3_nHe3eq);
  gsl_matrix_set (m, 2, 3, dT_nHe3eq);

  gsl_matrix_set (m, 3, 0, dnH2_Teq);
  gsl_matrix_set (m, 3, 1, dnHe2_Teq);
  gsl_matrix_set (m, 3, 2, dnHe3_Teq);
  gsl_matrix_set (m, 3, 3, dT_Teq);

  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[3] = 0.0;                /* There is no explicit time-dependence */


  /* Testing */
  
//   printf("jac.c: dfdy[0][0] = %le \n",dnH2_nH2eq);
//   printf("jac.c: dfdy[0][1] = %le \n",dnHe2_nH2eq);
//   printf("jac.c: dfdy[0][2] = %le \n",dnHe3_nH2eq);
//   printf("jac.c: dfdy[0][3] = %le \n",dT_nH2eq);
//           
//   printf("jac.c: dfdy[1][0] = %le \n",dnH2_nHe2eq);
//   printf("jac.c: dfdy[1][1] = %le \n",dnHe2_nHe2eq);
//   printf("jac.c: dfdy[1][2] = %le \n",dnHe3_nHe2eq);
//   printf("jac.c: dfdy[1][3] = %le \n",dT_nHe2eq);
//           
//   printf("jac.c: dfdy[2][0] = %le \n",dnH2_nHe3eq);
//   printf("jac.c: dfdy[2][1] = %le \n",dnHe2_nHe3eq);
//   printf("jac.c: dfdy[2][2] = %le \n",dnHe3_nHe3eq);
//   printf("jac.c: dfdy[2][3] = %le \n",dT_nHe3eq);
//           
//   printf("jac.c: dfdy[3][0] = %le \n",dnH2_Teq);
//   printf("jac.c: dfdy[3][1] = %le \n",dnHe2_Teq);
//   printf("jac.c: dfdy[3][2] = %le \n",dnHe3_Teq);
//   printf("jac.c: dfdy[3][3] = %le \n",dT_Teq);
//   
  //exit(0);
  


  return GSL_SUCCESS;


}
