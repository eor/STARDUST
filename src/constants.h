/***************************************************************
 * Physical constants
 ***************************************************************/

// TODO: rename all to comply with naming scheme!

#ifndef CONSTANTS_H

#define CONSTANTS_H


/* vel. of light  [cm/sec] */
#define LIGHTVEL 2.9979e10   

/* two state decay rate from 2S for H-atoms */
#define LRATE 8.227                 

/* Electron mass  [grams] */
#define MASSELECTRON 9.10938188e-28 

/* Electron mass  [eV] */
#define MASSELECTRONeV  511e6  

/* Proton mass [grams] */
#define MASSPROTON  1.672e-24

/* Boltzmann Constant [erg/kelvin] */
#define k_B 1.3807e-16          

/* Boltzmann Constant [eV/kelvin] */
#define k_BeV 8.617e-05         

/* eV --> ergs */
#define eVCGS 1.60217646e-12       

/* Threshold for hydrogen ionization 13.6 [eV] */
#define HIONIZEeV  13.6

/* Threshold for helium one e- ionization 24.6 [eV] */
#define He1IONIZEeV  24.6

/* Threshold for helium full ionization 54.4 [eV] */
#define He2IONIZEeV  54.4

/* v=Lyman-alpha frequency (Actually h $\nu$ /k) */
#define LYALPHAhvk 118408.90           

/* Refer (Fukugita & Kawasaki, 1994) eq.26 */
#define K_CONS  7.1541e-17       

/* Mpc-->cms */
#define MPC  3.086e24            

/* Kpc-->cms */
#define KPC  3.086e21  

/* kilometres to centimetres conversion */
#define KM2CM  1.0e5

/* Million year-->seconds */
#define MYR  3.153e13            

/* Planck's Constant [ergs sec] */
#define PLANCKCONSTCGS  6.265e-27            

/* Planck's Constant [eV sec] */
#define PLANCKCONSTeV  4.14e-15

/* thompson cross section [cm^2] */
#define THOMSON_ELEC_CROSS  6.6524e-25

/* Einstein A Coeffiecient [sec^-1] */
#define EINSTEIN_A  2.85e-15     

/* Mass of SUN [grams] */
#define MSUNCGS  1.99e33 

/* Gravitational constant [cm^3 grams^-1 s^-2] */
#define GRAVCONST  6.67e-8               

/* Integral Blackbody [erg cm^-2 K-4 sec^-1] */
#define STEFANBOLTZCONST  5.6705e-5       

/* Solar Luminosity [ergs sec^-1] */
#define LSUN 3.89e33            

/* Electron Charge [CGS units] */
#define ELEC_CHARGE  4.803206814e-10       

/* Lyman-alpha frequency at line-centre*/
#define LYALPHAFREQ  2.6095e15          

/* Hydrogen density (total = protons + H-atoms) at z = 0 [ number cm^-3 ]*/
#define  n_H0  1.9e-7 

/* Helium density total at z = 0 [ number cm^-3 ] */
#ifdef STROEMGRENTEST
    #define n_He0  0.0
    #define STROEMGRENTEMP 1.0e4
#else
    #define n_He0  1.5e-8
#endif    

/* All at z=0 [ number cm^-3 ] */
#define  nB0  2.5e-7

/* Redshift at which kinetic temperature of the IGM equals background CMB */
#define  zTkinEQTCMB 250


#endif
