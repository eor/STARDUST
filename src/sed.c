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
 * global variables needed only in the scope of this file
 ***************************************************************/
int SEDLineCount    = 0;
double SEDNorm      = 0;



/***************************************************************
 * for a given energy return (interpolated) SED value
 ***************************************************************/
double sed_get_value(double E){

    
#ifdef STROEMGRENTEST 
//        if (E>=13.6 && E<14.1)
//        if (E>=13.6 && E<13.8)

        if (E>=13.6 && E<14.6)
            return SEDNorm * 1.0;
        else 
            return 1e-100; 
#else  

  /* We need to interpolate between the Lambda from the array to
     finds its corresponding Energy */
    double interpolSED; 
//     double dy, dy2;

    gsl_interp_accel *lookupCache  = gsl_interp_accel_alloc();
//    gsl_spline       *splineObject = gsl_spline_alloc (gsl_interp_cspline, SEDLineCount);
     gsl_spline       *splineObject = gsl_spline_alloc (gsl_interp_linear, SEDLineCount);
    
    gsl_spline_init (splineObject, Lambda, Energy, SEDLineCount);

    interpolSED = gsl_spline_eval (splineObject, E, lookupCache);
//     dy  = gsl_spline_eval_deriv(splineObject, E, lookupCache);
//     dy2 = gsl_spline_eval_deriv2(splineObject, E, lookupCache);

    gsl_spline_free (splineObject);
    gsl_interp_accel_free (lookupCache);    


    return interpolSED;
#endif
    
}

/***************************************************************
 * SED log interpolation
 ***************************************************************/
double sed_get_log_value(double E, void *pars){

   
    /***************************************************************
     * For a given energy this function returns an interpolated 
     * SED value in log. 
     ***************************************************************/
 
    E = exp(E);
    double interpolValue;

    gsl_interp_accel *lookupCache    = gsl_interp_accel_alloc();
    gsl_spline       *spline = gsl_spline_alloc(gsl_interp_cspline, SEDLineCount);

    gsl_spline_init(spline, Lambda, Energy, SEDLineCount);

    interpolValue = gsl_spline_eval(spline, E, lookupCache);

    gsl_spline_free (spline);
    gsl_interp_accel_free (lookupCache);
  
    return interpolValue; 

}


/***************************************************************
 * compute normalization for a given mass
 ***************************************************************/
// double  Norm (double mass){
// 
//     /* 
//      * mass here is just the parameter which is used to determine normalization 
//      *
//      *  fk: deprecated, will be deleted.
//      */
// 
//     log_debug("Entering Norm()");
// 
// 
//     gsl_integration_workspace *wnorm = gsl_integration_workspace_alloc (1000);
// 
//     double result, error;
//     double alpha =0.;
//     gsl_function Func;
//     Func.function = &sed_get_log_value;
//     Func.params = &alpha;
//     
//     log_debug(" we are runing gsl_integration_qags (sed_get_log_value)");
//     
//     gsl_integration_qags( &Func, log(myConfig.sourceELow), log(myConfig.sourceEHigh), 0, 1e-5, 1000,
//                             wnorm, &result, &error );
// 
//     gsl_integration_workspace_free (wnorm);
// 
//     return (mass * SEDNorm) / (result * 4 * M_PI * pow (MPC, 2.0));
//     
// 
// 
// }


/***************************************************************
 * read spectrum of the radiating source
 ***************************************************************/
void sed_read_file(char *SEDFileName){
    
    /***************************************************************
     * This function reads the provided SED file.
     * It assumes that the first line is a header, i.e.
     * it will be ignored.
     ***************************************************************/
    
    if(DEBUG)log_debug("Entering sed_read_file(...)");

    FILE *fp;
    char buf[256];
    double lam, egy;
    int i;
    
    /* open file */    
    printf(" Trying to read SED from file \'%s\' \n",SEDFileName);
    fp = fopen(SEDFileName,"r");
 
    if(!fp){ 
        printf("ERROR: Could not read SED file %s \n", SEDFileName);
        exit(1); 
    }

    /* Reading Header */
    if( !fgets(buf, sizeof buf, fp) ){
        printf(" Error. Could not header from file '%s'. Exiting.\n", SEDFileName);
        exit(1);        
    }    
    
    if(DEBUG) log_debug("File header: %s", buf);       
    
    /* Get line SEDLineCount */
    while(!feof(fp)){
        if( !fscanf(fp,"%le %le \n",&lam,&egy)){
            printf(" Error. Could not read from file '%s'. Exiting.\n", SEDFileName);
            exit(1);      
            
        }
        SEDLineCount++;
    }    
    
    if(DEBUG)log_debug("SEDLineCount = %d", SEDLineCount);        
     
    
    /* Go back to the beginning */    
    rewind(fp);    
    
    /* skip header again */
    if( !fgets(buf, sizeof buf, fp) ){
        printf(" Error. Could not read header file '%s'. Exiting.\n", SEDFileName);
        exit(1);        
    }
    
    /* allocate memory */ 
    double *LambdaT, *EnergyT;
    LambdaT = (double*) malloc(SEDLineCount * sizeof(double));      // these two will be removed
    EnergyT = (double*) malloc(SEDLineCount * sizeof(double));

    Lambda = (double*) malloc(SEDLineCount * sizeof(double));       // global, TODO: rename those two
    Energy = (double*) malloc(SEDLineCount * sizeof(double));   
        
    
    // legacy stuff: the following will be thrown away:
    if(type!=2){
        for(i=0;i<SEDLineCount;i++){
            
            if( !(fscanf(fp,"%le %le \n", &LambdaT[i],&EnergyT[i]) ) ){
                printf(" Error. Could not read from file '%s'. Exiting.\n", SEDFileName);
                exit(1);                                
            }
            
            
            /* Converting Lambda (in Angstroms) to energy*/
            /* and converting energy from erg/sec to eV/sec */
            Lambda[SEDLineCount-1-i] = PLANCKCONSTeV * LIGHTVEL/ (LambdaT[i] * 1.0e-8) ; /* [eV] */
            
            Energy[SEDLineCount-1-i] = pow(10.0,EnergyT[i])/eVCGS;
            /* Because the table is in LOG units and we reverse the order because */
            /* the interpolation scheme likes things to be monotonically increasing */
        }
    } // legacy read function. TODO remove 
 
 
    /* read data */
    if(type==2){        
        for(i=0; i<SEDLineCount; i++){            
             if( !( fscanf( fp,"%le %le \n", &Lambda[i],&Energy[i] ) ) ){
                    printf(" Error. Could not read from file '%s'. Exiting.\n", SEDFileName);
                    log_error("Could not read from file '%s'. Exiting.", SEDFileName);
                    log_close();
                    exit(1);
            } 
        }      
    }
 
    if(DEBUG){        
        log_debug("Printing SED:");
        log_debug("# Photon E [eV]\tTotal Intensity [eV/s/eV]");           
        for(i=0;i<SEDLineCount;i++){ log_debug(" %e\t%e",Lambda[i], Energy[i]); }  
        log_debug("End of SED file");        
    }
    
    
    /* check energy range of the SED. part 1. 
     * If the user-provided range [ELow, EHigh] is out of bounds, throw an error. */
    if(myConfig.sourceELow< Lambda[0] || myConfig.sourceEHigh> Lambda[SEDLineCount-1]){
        printf("ERROR: Provided energy range [sourceELow, sourceEHigh] is out of bounds. Exiting.\n");
        log_error("Provided energy range [sourceELow, sourceEHigh] is out of bounds.");
        log_close();
        exit(1);
    }    

    /* Check energy range of the SED. part 2. 
     * If the provided E range  is smaller than the range in SED file, the 
     * unnecessary entries are removed.  
     * 
     * We use an el-cheapo method to resize the dynamic SED arrays. 
     * Not elegant but it seems to do the job:
     * 
     * 1. find indexes of energies that come closest to the provided limits     * 
     */      
    
    double  absDeltaELow, absDeltaEHigh;
    double  tmpDeltaELow    = 0.0;
    double  tmpDeltaEHigh   = 0.0;    
    int     indexELow       = 0;
    int     indexEHigh      = 0;
    
    absDeltaELow  = fabs(myConfig.sourceELow  - Lambda[0]);        
    absDeltaEHigh = fabs(myConfig.sourceEHigh - Lambda[0]);
        
    for(i=0; i<SEDLineCount; i++){
        
        tmpDeltaELow  = fabs(myConfig.sourceELow - Lambda[i]);        
        tmpDeltaEHigh = fabs(myConfig.sourceEHigh - Lambda[i]);          
        
        if (tmpDeltaELow<=absDeltaELow && Lambda[i]<13.6){
            absDeltaELow = tmpDeltaELow;
            indexELow = i;   
        }
        
        if(tmpDeltaEHigh<=absDeltaEHigh){
            absDeltaEHigh = tmpDeltaEHigh;
            indexEHigh = i;
        }
        
    }
    printf(" Energy range closest to user input is [%e, %e]\n", Lambda[indexELow], Lambda[indexEHigh]);
    
    /* 2. count elements  */
    int newSEDLineCount = indexEHigh - indexELow+1;
    if(newSEDLineCount != SEDLineCount){printf(" Resizing SED arrays accordingly\n");}
    
    /* 3. allocate memory for tmp E and lambda arrays here */
    double *tmpLambda = (double*) malloc(newSEDLineCount * sizeof(double)); 
    double *tmpEnergy = (double*) malloc(newSEDLineCount * sizeof(double)); 
    
    /* 4. copy content */
    for(i=0; i<newSEDLineCount; i++){        
        tmpLambda[i] = Lambda[i+indexELow];
        tmpEnergy[i] = Energy[i+indexELow];        
    }
    
    /* 5. resize original arrays and copy content back */
    Lambda = (double*)realloc(Lambda, newSEDLineCount * sizeof(double));
    Energy = (double*)realloc(Energy, newSEDLineCount * sizeof(double));
    
    for(i=0; i<newSEDLineCount; i++){  
        Lambda[i] = tmpLambda[i];
        Energy[i] = tmpEnergy[i];
    }
    
    SEDLineCount = newSEDLineCount;
    
 
//     for(i=0;i<SEDLineCount;i++){
//         Energy[i] = Energy[i]/6.2415e11 ;
//     }
//     
    /* compute normalization  */
    SEDNorm = sed_compute_norm();   // 
    printf(" Integral over SED (normalization) : %le eV/sec\n", SEDNorm);

    
    /* clean up */
    free(tmpLambda);
    free(tmpEnergy);  
    free(LambdaT);
    free(EnergyT);
    fclose(fp);
}


/***************************************************************
 * compute SED normalization 
 ***************************************************************/
double sed_compute_norm(){
    
    /***************************************************************
     * This function computes the normalization of the SED
     * in the (approximated) user-provided limits.
     * It uses the GSL QAGS integrator. If the desired relative
     * tolerance cannot be achieved, it will be increased and then
     * integration will be re-run.
     ***************************************************************/
    
    /* logging */
    if(DEBUG){
        log_debug("Entering sed_compute_norm() to compute SEDNorm");
        log_debug("Integrating provided SED in the limits of (%le,%le) [eV] ", myConfig.sourceELow, myConfig.sourceEHigh);
    }       
    
#ifdef STROEMGRENTEST
    /* here we by-pass the integration for the Strömgen sphere test (delta-function SED). */    
    double maxE = 0.0;
    int i;
    for(i=0;i<SEDLineCount;i++){
        if (Energy[i]>maxE)
            maxE = Energy[i];                        
    } 
  
    return maxE;
#else   
    
    /* save error handler to old_handler, then switch it off */
    gsl_error_handler_t * oldErrorHandler = gsl_set_error_handler_off();    
    
    /* set up integration workspace, variables and function */
    gsl_integration_workspace *wnorm = gsl_integration_workspace_alloc(1000);

    double result, error;
    double alpha      =  0.;                            // alpha is only here for legacy reasons, TODO: get rid of it!
    double absError   = 1e-4;
    double relError   = 1e-6;                           // relative error
    double lowerLimit = log(Lambda[0]);
    double upperLimit = log(Lambda[SEDLineCount-1]);
    
    gsl_function  myFunction;
    myFunction.function     = &sed_get_log_value;
    myFunction.params       = &alpha;
   
    int status = 1;     
    
    while(status){

        /* if successful, 0 is returned. If the relError is too small, status will be 13 */
        status     = gsl_integration_qags(&myFunction, lowerLimit , upperLimit, absError, relError, 1000, wnorm, &result, &error);        
        relError  *= 5;

        if(status){            
            if(DEBUG){
                    log_debug("gsl_integration_qags status = %d", status);                    
                    log_debug("SEDNorm integration: Increased relError to %e",relError);                                          
            }
        }
    }

    gsl_set_error_handler(oldErrorHandler); // reset error handler (might be unnecessary).
    gsl_integration_workspace_free (wnorm);

    if(DEBUG){
        log_debug("Integrated SED = %e", result);
    }
    
    return result; 
  
    /* references for:
     * the integrator: https://www.gnu.org/software/gsl/manual/html_node/QAGS-adaptive-integration-with-singularities.html
     * error handling: https://lists.gnu.org/archive/html/help-gsl/2006-04/msg00055.html   
     */ 
#endif
}


/***************************************************************
 * compute an estimated radius of the Stroemgren's sphere
 ***************************************************************/
double sed_estimate_stroemgren_radius( double redshift, double time ){
    
        // find density
      double density = n_H0  * OverDensity * pow3(1 + redshift);
      
      // find number of emitted ionizing photons 
      double ionPhotonN = (SEDNorm/13.6) * (time)*MYR;  
                    
      // find ionized volume
      double volume = ionPhotonN / density; 

      // volume [cm^3] --> radius [kpc]
      double radius = pow(0.75 * volume / 3.14, 1./3.) / KPC;
      
      return radius;  
    
    
}


