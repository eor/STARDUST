/***************************************************************
 * 
 * this file contains all function prototypes of the code
 * 
 ***************************************************************/

#ifndef ALLVARS_H
#include "allvars.h"
#endif


/***************************************************************
 * rt.c - main run
 ***************************************************************/

void rt_main_run();
void rt_initialize_grids(double);
int rt_check_is_last_time_step(double, double);

/***************************************************************
 * config.c - config related
 ***************************************************************/
void config_load_from_file(char *fileName);

/***************************************************************
 * SED-related functions
 ***************************************************************/
void sed_read_file(char *SEDFileName);

double sed_compute_norm();

// Calculates the normalization of the SED 
// double Norm (double );       
  
// return interpolation of SED (in log) for given E
double sed_get_log_value(double , void *); 

// return interpolation of SED for given E
double sed_get_value(double );


double sed_estimate_stroemgren_radius( double redshift, double time );

/***************************************************************
 * logging
 ***************************************************************/
void log_start();
void log_close();
void log_time(int type);

/***************************************************************
 * functions.c - misc utils
 ***************************************************************/


 
/* Calculates H(z)  */
double hubble (double );

/* Returns the redshift given time*/
double Findz (double , void *);
         
/* Given redshift computes time */
double time_from_redshift(double , void *);


double cross_sec_e1h1(double );
double cross_sec_ehe1(double );
double cross_sec_ehe2(double );
  

/***************************************************************
 * derivs.c and jac.c functions
 ***************************************************************/
         
/* Sets up the system of differential equations */      
int derivs (double t, const double y[], double dydt[], void *params);
         
/* Computes the jacobian */
int jac (double t, const double y[], double *dfdy, double dfdt[], void *p);


/***************************************************************
 * table_* functions
 ***************************************************************/

/* Integrals for ionizaton */
int table_ion (void );
         
/* Integrals for temperature */
int table_temp (void );
         
/* Integrals for compton term */
int table_compton (void );
         
/* table reading routine */
void table_read_all(void );
         
/**Prototypes for function to make the tables of the integrals**/

double tau_e1h1_p1_log(double ,void * params);  /* p1 => part-(1), */
double tau_e1h1_p2_log(double ,void * params);  /* p2 => part-(2) etc., */
double tau_e1h1_p3_log(double ,void * params);
double tau_e1h1_p4_log(double ,void * params);

double tau_ehe1_p2_log(double ,void * params);
double tau_ehe1_p3_log(double ,void * params);
double tau_ehe1_p4_log(double ,void * params);

double tau_ehe2_p3_log(double ,void * params);
double tau_ehe2_p4_log(double ,void * params);
/*Integrals for ionizaton above */

double temp_e1h1_p1_log(double ,void * params);  /* p1 => part-(1), */
double temp_e1h1_p2_log(double ,void * params);  /* p2 => part-(2) etc., */
double temp_e1h1_p3_log(double ,void * params);
double temp_e1h1_p4_log(double ,void * params);

double temp_ehe1_p2_log(double ,void * params);
double temp_ehe1_p3_log(double ,void * params);
double temp_ehe1_p4_log(double ,void * params);

double temp_ehe2_p3_log(double,void * params);
double temp_ehe2_p4_log(double,void * params);
/* Integrals for temperature above */

double compton_heating1_p1(double ,void * params);
double compton_heating1_p2(double ,void * params);
double compton_heating1_p3(double ,void * params);
double compton_heating1_p4(double ,void * params);

double compton_heating2_p1(double ,void * params);
double compton_heating2_p2(double ,void * params);
double compton_heating2_p3(double ,void * params);
double compton_heating2_p4(double ,void * params);        
/* Integrals for compton term above */

/***************************************************************
 * functions in utils.c 
 ***************************************************************/

/* contruct paths */
char* utils_concat_path(char *dir, char *file);
char* utils_concat_path_noID(char *dir, char *file);

/* fallback help */
void utils_print_help();    

/* standard greeting */
void utils_print_hello();   

/***************************************************************
 * functions in interpolation.c
 ***************************************************************/

/* interpolates the table of integrals*/
void interpolation(double, double, double ,double );
         

/* The 2-D interpolation routine */
void interpolation_2D(double*, double*, double** , int , int , double , double , double *);

/***************************************************************
 * etc
 ***************************************************************/

/* Allocates all memory required */
void memory_allocate_all(void );

/* Frees all allocated memory */
void memory_free_all(void );

