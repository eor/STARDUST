/***************************************************************
 * libraries
 ***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/***************************************************************
 * GSL header files
 ***************************************************************/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

/***************************************************************
 * SD headers
 ***************************************************************/
#include "constants.h"
#include "allvars.h"
#include "prototype.h"
#include "table_settings.h"
#include "log.h"

/***************************************************************
 * Etc
 ***************************************************************/
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))


/***************************************************************
 * given a redshift the function below computes the time
 ***************************************************************/
double time_from_redshift(double z, void * pars){

    double omegaM   = myConfig.cosmoOmegaM;
    double omegaL   = myConfig.cosmoOmegaL;

    return 1. / ( (1. + z) * sqrt( pow2(1 + z) * ( 1. + omegaM * z ) - z * (2. + z) * omegaL) );

}

/***************************************************************
 * return the redshift for a particular time
 ***************************************************************/
double Findz (double rshift, void * params){


    struct Findz_params *p = (struct Findz_params *) params;

    double CummTime = p->cummtime;
    CummTime *= MYR;
    double z_of_src = p->zsrc;
 
    double Hubb;
    Hubb = myConfig.cosmoH0 * KM2CM / MPC;     
 
    if (rshift == z_of_src)
        return (CummTime - 0.);
    else{
        
        gsl_integration_workspace *wbal = gsl_integration_workspace_alloc (1000);

        double result, error;

        gsl_function Func;
        Func.function = &time_from_redshift;
        Func.params = 0;

        gsl_integration_qags (&Func, rshift, z_of_src, 0, 1e-7, 1000, wbal, &result, &error);
        
        //TODO: relErr absErr, different use rk8 instead of qags

        gsl_integration_workspace_free (wbal);

        return (CummTime - (result / Hubb));

    }// else statement not needed, fk
}

/***************************************************************
 * calculate the Hubble constant at a given redshift
 ***************************************************************/
double hubble (double z){

    double Hubble;    
    double omegaM   = myConfig.cosmoOmegaM;
    double hubble0  = myConfig.cosmoH0;
    
    Hubble = hubble0 * KM2CM / MPC;     
    return Hubble * pow( omegaM * pow3( 1.+z ) + (1.-omegaM) , .5 );
  
}



/***************************************************************
 * functions to calculate various cross-sections
 ***************************************************************/
double cross_sec_e1h1(double E){
    /* ?? HI (n=1 --> free) */ 
    double sig_0,F,E0,x,y,y0,y1,yw,yaa,p,Mb;
    
    Mb      = 1.e-18;
    sig_0   = 5.475e4;
    E0      = 4.298e-1; 
    y0      = 0.;
    y1      = 0.;
    yw      = 0.;
    yaa     = 3.288e1;
    p       = 2.963;
    x       = E/E0 - y0;
    y       = sqrt(x*x+y1*y1); 
    F       = ((x-1)*(x-1)+yw*yw)*pow(y,(0.5*p-5.5))*pow(1+sqrt(y/yaa),-p);
    
    return sig_0*F*Mb;
}


double cross_sec_ehe1(double E){
    
    double sig_0,F,E0,x,y,y0,y1,yw,yaa,p,Mb;

    Mb      = 1.e-18;
    sig_0   = 9.492e2;
    E0      = 1.361e1; 
    y0      = 4.434e-1;
    y1      = 2.136;
    yw      = 2.039; 
    yaa     = 1.469;
    p       = 3.188;
    x       = E/E0 - y0;
    y       = sqrt(x*x+y1*y1); 
    F       = ((x-1)*(x-1)+yw*yw)*pow(y,(0.5*p-5.5))*pow(1+sqrt(y/yaa),-p);
    
    return sig_0*F*Mb;
}

double cross_sec_ehe2(double E){
    
    double sig_0,F,E0,x,y,y0,y1,yw,yaa,p,Mb;
    
    Mb      = 1.e-18;
    sig_0   = 1.369e4;
    E0      = 1.72; 
    y0      = 0.;
    y1      = 0.;
    yw      = 0.;
    yaa     = 3.288e1;
    p       = 2.963;
    x       = E/E0 - y0;
    y       = sqrt(x*x+y1*y1); 
    F       = ((x-1)*(x-1)+yw*yw)*pow(y,(0.5*p-5.5))*pow(1+sqrt(y/yaa),-p);
    
    return sig_0 * F * Mb;
}


