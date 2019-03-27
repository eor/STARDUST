/***************************************************************
 * About:
 * Declarations of all global variables.
 * Instances of these variables are provided in allvars.c
 * 
 ***************************************************************/

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/***************************************************************
 * STARDUST version
 ***************************************************************/
#define SD_VERSION 0.1234

/***************************************************************
 * config struct
 ***************************************************************/
extern struct global_config{       
    
    /* paths */
    char    pathOutDir[256];
    char    pathSED[256];
    char    pathDensity[256];
    char    pathID[256];    
    
    /* simulation */
    double  sourceELow;
    double  sourceEHigh;     
    double  sourceLifetime;  

    double  haloMass;
    
    double  redshiftLow;        // to be removed at some point
    double  redshiftHigh;    
    double  redshiftStride;      
    
    /* cosmology */
    double  cosmoOmegaM;        //  total matter density 
    double  cosmoOmegaL;        //  Dark energy / cosmological constant
    double  cosmoOmegaB;        //  Baryon density
    double  cosmoH0;            //  Hubble parameter [km/sec/Mpc] at z = 0
    double  cosmoH100;          //  cosmoH0/100 at  z = 0
    double  cosmoSigma8;        //  power spectrum normalization (overdensity in 8kpc spherical vol)
    double  cosmoTauThom;       //  Optical depth of Thomson scattering
    double  cosmoTCMB;          //  CMB temperature at z = 0    
    
    /*general setting */
    int     settingsDebug;   
  
    double  settingsRMax;    
    double  settingsRStart;  
    
    double  settingsDeltaR;  
    double  settingsDeltaT;  

    double  settingsWriteT;
    
} myConfig;


/***************************************************************
 * SED related
 ***************************************************************/
extern double *Lambda, *Energy;
extern int type; 



/***************************************************************
 * all things logging & file writing
 ***************************************************************/
extern const char  *debugLogFile;
extern const char  *mainLogFile;

extern FILE *fpDebugLog;
extern FILE *fpMainLog;

extern int DEBUG;    

/* Stores the output file name */
extern char *outFile; 


/***************************************************************
 * physics
 ***************************************************************/
extern double drKpc;

extern float OverDensity;

/* Just a constant evaluted in rt.c */
extern double lamcons;


/***************************************************************
 * computing grids
 ***************************************************************/

/* array to store the background overdensity along the LOS */
extern double *over_densities;

/* An array that stores the free-electron densities at every grid */
extern double *ne;

/* Arrays that stores the neutral hydrogen and helium densities at every grid point */
extern double *n_H1;
extern double *n_He1;

/* iGrid is the radial grid index for the our computing grids */
extern int iGrid;

/* Number of grid points of our computing grids*/
extern int numGridPoints;     


/*  Arrays that store the densities of the different species at all 
    grid points. Needs to be initialized with the right BCs. */
extern double  *n_H2, *n_He2, *n_He3, *T_e;

/* tmp storage for the fractions of the different species */
extern double *x_HI, *x_HII,  *x_HeI, *x_HeII, *x_HeIII;  

/* The spin and brightness temperatures at grid locations */
extern double *T_spin, *T_brig;



/***************************************************************
 * variables for solving the ODE 
 ***************************************************************/

/* Constants that stores that 4*PI*R^2 factor */
extern double Ag;

/* the current physical radius at which the computations are performed */
extern double radius;


/* Current redshift */
extern double z;

/* Value of the redshift at the previous timestep */
extern double zPrevious;

/* column densities */
extern double nHX1, nHeX2, nHX13;    

extern struct ode_params{
    
    double fe1h1;    /*    fuku_e1h1[iGrid]     */
    double fehe1;    /*    fuku_ehe1[iGrid]     */
    double fehe2;    /*    fuku_ehe2[iGrid]     */
    double intH1;    /*    integral_H1[iGrid]   */
    double intHe1;   /*    integral_He1[iGrid]  */
    double intHe2;   /*    integral_He2[iGrid]  */
    double localOD;  /*   over_densities[iGrid] */

} params;
/* params used at a given radius includes the following:
 *
 * - Ionization rates as given in Ref[2] for solution of the ODE
 * - Heating rates used to compute the temperature evolution
 * - Compton heating / cooling in the temp-evol ODE added on top of the
 *   equations found in Ref[2] (CURRENTLY EXCLUDED)
 */


/* Structure used by Findz function to locate a particular redshift */
extern struct Findz_params{
    
    double cummtime; /* time since the source has turned ON */
    double zsrc;   /* redshift at which the source turns ON */
    
} findzparams;

/***************************************************************
 * everything related to tabulated integrals & interpolation
 ***************************************************************/

/* Heating rates used to compute the temperature evolution */
extern double *integral_H1, *integral_He1, *integral_He2;
/* The values of the above are taken from pre-computed tables generated
   using table_temp.c */

extern double *comp_integ1, *comp_integ2;
/* Compton heating / cooling in the temp-evol ODE added on top of the
   equations found in Ref[2] */

/* Ionization rates as given in Ref[2] for solution of the ODE */
extern double *fuku_e1h1, *fuku_ehe1, *fuku_ehe2;
/* The values of the above are taken from pre-computed tables generated
   using table_ion.c */

/* variables for interpolation */
extern double *xa, *nHx1a, *nHex2a, *nHx13a, *ya, **yP2a, **yP3a;

/* tabulated values of integral */
extern gsl_vector *e1h1_p1,*e1h1_p2,*e1h1_p3;                          
extern gsl_vector *ehe1_p2,*ehe1_p3,*ehe2_p3;                  
                                                        
extern gsl_vector *temp_e1h1_p1,*temp_e1h1_p2,*temp_e1h1_p3;   
extern gsl_vector *temp_ehe1_p2,*temp_ehe1_p3,*temp_ehe2_p3;   

extern gsl_vector *comp1_p1,*comp1_p2,*comp1_p3;
extern gsl_vector *comp2_p1,*comp2_p2,*comp2_p3;        

/* Converting to matrix form */
extern gsl_matrix *matrix_e1h1_p2,*matrix_e1h1_p3;
extern gsl_matrix *matrix_ehe1_p2,*matrix_ehe1_p3;
extern gsl_matrix *matrix_ehe2_p3;

extern gsl_matrix *temp_matrix_e1h1_p2,*temp_matrix_e1h1_p3;   
extern gsl_matrix *temp_matrix_ehe1_p2,*temp_matrix_ehe1_p3;   
extern gsl_matrix *temp_matrix_ehe2_p3;

extern gsl_matrix *comp1_matrix_p2,*comp1_matrix_p3;           
extern gsl_matrix *comp2_matrix_p2,*comp2_matrix_p3; 


#endif /* allvars.h */
