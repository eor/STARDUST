/***************************************************************
 * About:
 * provides instances of all global variables.
 * Declarations can be found in allvars.h
 * 
 ***************************************************************/

#include <stdio.h>
#include <gsl/gsl_rng.h>

#include "constants.h"
#include "allvars.h"

//TODO: Lots of variables should be renamed here

/***************************************************************
 * config
 ***************************************************************/
struct global_config myConfig;


/***************************************************************
 * SED related
 ***************************************************************/
double  *Lambda, *Energy;  // TODO: change this!
int     type = 2;


/***************************************************************
 * file writing and logging related
 ***************************************************************/
const char  *debugLogFile = "log_debug";
const char  *mainLogFile  = "log_main";

FILE *fpDebugLog = 0;
FILE *fpMainLog  = 0;
int DEBUG = 0;    

 /* Stores the output file name */
char *outFile; 
  

/***************************************************************
 * physics
 ***************************************************************/
double drKpc = KPC;

/* Square of the clumping factor */
float OverDensity = 1.;

/* Just a constant evaluated in rt.c */
double lamcons;


/***************************************************************
 * computing grids
 ***************************************************************/

/* array to store the background overdensity along the LOS */
double *over_densities;

/* An array that stores the free-electron densities at every grid */
double *ne;

/* Arrays that store the neutral hydrogen and helium densities at every grid */
double *n_H1; 
double *n_He1;

/* iGrid is the index to a particular grid */
int iGrid;

/* Number of grid points in Space  */
int numGridPoints;  

/*  Arrays that store the densities of the different species at all 
    grid points. Needs to be initialized with the right BCs. */
double *n_H2, *n_He2, *n_He3;

/* tmp storage for the fractions of the different species */
double *x_HI, *x_HII,  *x_HeI, *x_HeII, *x_HeIII;  

/* The spin and brightness temperatures at grid locations */
double *T_e, *T_spin, *T_brig;


/***************************************************************
 * variables for solving the ODE 
 ***************************************************************/

/* Constants that stores that 4*PI*R^2 factor that is computed below */
double Ag;

/* Radius at which the computations are performed */
double radius;

/* Current redshift */
double z;

/* Value of the redshift at the previous time step */
double zPrevious;

/* column densities */
double nHX1, nHeX2, nHX13;    



struct ode_params params;
/* params used at a given radius includes the following:
 *
 * - Ionization rates as given in Ref[2] for solution of the ODE
 * - Heating rates used to compute the temperature evolution
 * - Compton heating / cooling in the temp-evol ODE added on top of the
 *   equations found in Ref[2] (CURRENTLY EXCLUDED)
 */


/* Parameters passed to the Findz function in function.c */
struct Findz_params findzparams;


/***************************************************************
 * everything related to tabulated integrals & interpolation
 ***************************************************************/

/* Heating rates used to compute the temperature evolution */
double *integral_H1, *integral_He1, *integral_He2;
/* The values of the above are taken from pre-computed tables generated
   using table_temp.c */

double *comp_integ1, *comp_integ2;
/* Compton heating / cooling in the temp-evol ODE added on top of the
   equations found in Ref[2] */

/* Ionization rates as given in Ref[2] for solution of the ODE */
double *fuku_e1h1, *fuku_ehe1, *fuku_ehe2;
/* The values of the above are taken from pre-computed tables generated
   using table_ion.c */

/* variables for interpolation */
double *xa, *nHx1a, *nHex2a, *nHx13a, *ya, **yP2a, **yP3a;

/* tabulated values of integral */
gsl_vector *e1h1_p1,*e1h1_p2,*e1h1_p3;                          
gsl_vector *ehe1_p2,*ehe1_p3,*ehe2_p3;                  
                                                        
gsl_vector *temp_e1h1_p1,*temp_e1h1_p2,*temp_e1h1_p3;   
gsl_vector *temp_ehe1_p2,*temp_ehe1_p3,*temp_ehe2_p3;   

gsl_vector *comp1_p1,*comp1_p2,*comp1_p3;
gsl_vector *comp2_p1,*comp2_p2,*comp2_p3;        

/* Converting to matrix form */
gsl_matrix *matrix_e1h1_p2,*matrix_e1h1_p3;
gsl_matrix *matrix_ehe1_p2,*matrix_ehe1_p3;
gsl_matrix *matrix_ehe2_p3;

gsl_matrix *temp_matrix_e1h1_p2,*temp_matrix_e1h1_p3;   
gsl_matrix *temp_matrix_ehe1_p2,*temp_matrix_ehe1_p3;   
gsl_matrix *temp_matrix_ehe2_p3;

gsl_matrix *comp1_matrix_p2,*comp1_matrix_p3;	        
gsl_matrix *comp2_matrix_p2,*comp2_matrix_p3; 


/***************************************************************
 * graveyard
 ***************************************************************/

// int r_index;
/* r_index stores the maximum number of grids you have i.e.,
   r_index = (r_max - r_start)/dr ; were dr=resolution of the run 
   And iGrid is the index to a particular grid */

//extern double tot_time,z_ON;
/* Used in functions.h. Basically is the time to which all grids have 
   time evolved, i.e., tot_time = n times dt; n belongs to integers */